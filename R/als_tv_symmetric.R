####
# symmetric gaussian rw-smoothed time-varying als via l-bfgs.
# under B_t = A_t the data term ||Y - A X A^T||^2 is quartic in A, so the
# row-wise bcd machinery used in the directed case is not an exact
# conditional minimiser. instead optimise the whole trajectory
# {A_t}_{t=2..T} via l-bfgs, with sign canonicalisation by frobenius
# inner-product dp (operating on A only since B is identical).
####

#' L-BFGS solver for symmetric Gaussian TV-ALS (Stage 2)
#'
#' Optimizes the trajectory `{A_t}_{t=2..T}` directly under the quartic data
#' fit + RW(1) penalty + level ridge. The data Hessian is approximated by
#' L-BFGS; the RW penalty Hessian (block-tridiagonal) and level ridge are
#' added to the gradient analytically.
#'
#' @param Y 4D array [n, n, p, Tt]
#' @param lambda RW penalty
#' @param mu level penalty
#' @param max_iter L-BFGS iteration cap
#' @param tol_obj relative objective convergence tolerance
#' @param init list with A_init [n, n, Tt]; M_init [n, n, p]
#' @param verbose logical
#' @return list with A, M, objective trajectory, convergence diagnostics
#' @keywords internal
#' @noRd
.dbn_als_tv_symmetric_lbfgs <- function(Y, lambda, mu,
                                          max_iter = 200L, tol_obj = 1e-5,
                                          init = NULL,
                                          family = "gaussian",
                                          Y_obs = NULL,
                                          z_outer_iter = 5L,
                                          verbose = FALSE) {
	dims <- dim(Y); n <- dims[1]; p <- dims[3]; Tt <- dims[4]
	stopifnot(dims[1] == dims[2])  # must be square for symmetric

	# build response mask
	Omega <- array(1, dim = c(n, n, Tt))
	for (t in 1:Tt) diag(Omega[, , t]) <- 0
	for (r in 1:p) for (t in 1:Tt) {
		na_mask <- is.na(Y[, , r, t])
		Omega[, , t][na_mask] <- 0
	}
	Y_filled <- Y; Y_filled[is.na(Y_filled)] <- 0

	# initialize from static fit if not given
	if (is.null(init)) {
		stat_fit <- dbn_als(Y, family = "gaussian", symmetric = TRUE,
			bootstrap = 0, verbose = FALSE)
		# for symmetric static, fit$A and fit$B should be equal
		A0 <- stat_fit$A[[1]][, , 1]
		# force symmetric
		A0 <- 0.5 * (A0 + t(A0))
		A_init <- array(rep(c(A0), Tt), dim = c(n, n, Tt))
		M_init <- array(NA, dim = c(n, n, p))
		for (r in 1:p) M_init[, , r] <- stat_fit$M[, , r, 1]
		init <- list(A_init = A_init, M_init = M_init)
	}

	# working state: Phi_t = Z_t - M (for Gaussian, Z = Y)
	# we freeze M at init (matching directed Stage 1A behavior)
	A <- init$A_init
	M <- init$M_init

	# average centered state over p (assume p=1 for now)
	Phi <- array(0, dim = c(n, n, Tt))
	for (t in 1:Tt) {
		Phi[, , t] <- Y_filled[, , 1, t] - M[, , 1]
		Phi[, , t][Omega[, , t] == 0] <- 0
	}
	# boundary multiplier
	S <- Tt - 1L

	# symmetry parameterization: parameterize A by its upper triangle
	# (including diagonal) as a length n*(n+1)/2 vector per time slice,
	# and rebuild the full A via A = U + U^T - diag(diag(U)). this keeps
	# L-BFGS on the symmetric manifold automatically.
	idx_upper <- which(upper.tri(matrix(0, n, n), diag = TRUE))
	npar_per_t <- length(idx_upper)  # n*(n+1)/2

	pack_A <- function(A_arr) {
		# symmetrize each slice (defensive), then extract upper-tri entries
		v <- numeric(npar_per_t * S)
		for (s in 1:S) {
			A_sym <- 0.5 * (A_arr[, , s + 1L] + t(A_arr[, , s + 1L]))
			v[((s - 1L) * npar_per_t + 1L):(s * npar_per_t)] <- A_sym[idx_upper]
		}
		v
	}
	unpack_A <- function(v) {
		A_arr <- array(0, dim = c(n, n, Tt))
		A_arr[, , 1] <- diag(n)
		for (s in 1:S) {
			vs <- v[((s - 1L) * npar_per_t + 1L):(s * npar_per_t)]
			A_s <- matrix(0, n, n)
			A_s[idx_upper] <- vs
			# mirror to full symmetric matrix
			A_s <- A_s + t(A_s) - diag(diag(A_s))
			A_arr[, , s + 1L] <- A_s
		}
		A_arr
	}

	# objective for a candidate A trajectory (A is guaranteed symmetric)
	obj_fn <- function(v) {
		A_arr <- unpack_A(v)
		data_fit <- 0
		for (s in 1:S) {
			A_t <- A_arr[, , s + 1L]
			pred <- A_t %*% Phi[, , s] %*% t(A_t)
			resid <- Omega[, , s + 1L] * (Phi[, , s + 1L] - pred)
			data_fit <- data_fit + 0.5 * sum(resid^2)
		}
		rw <- 0
		if (S >= 2L) for (s in 2:S) {
			d <- A_arr[, , s + 1L] - A_arr[, , s]
			rw <- rw + 0.5 * lambda * sum(d^2)
		}
		lev <- 0
		for (s in 1:S) lev <- lev + 0.5 * mu * sum(A_arr[, , s + 1L]^2)
		data_fit + rw + lev
	}

	# gradient w.r.t. the symmetric parameterization.
	# strategy: compute the unconstrained gradient g_full w.r.t. the n*n
	# matrix A (as before), then transform to the symmetric param. For a
	# scalar function f(A) restricted to symmetric A:
	#   d f(sym(U)) / d U = g_full + t(g_full) - diag(diag(g_full + t(g_full))) / 2
	# more carefully: A = U + U^T - diag(U), so dA/dU_{ij} for i<j is
	# (E_{ij} + E_{ji}); dA/dU_{ii} = E_{ii}. So:
	#   df/dU_{ij} (i<j) = g_full_{ij} + g_full_{ji}    (off-diagonal upper)
	#   df/dU_{ii}        = g_full_{ii}                 (diagonal)
	# we compute g_full as before and then collapse onto the upper triangle.
	grad_fn <- function(v) {
		A_arr <- unpack_A(v)
		g_full <- array(0, dim = c(n, n, Tt))
		for (s in 1:S) {
			A_t <- A_arr[, , s + 1L]
			X_t <- Phi[, , s]
			Y_t <- Phi[, , s + 1L]
			Omg <- Omega[, , s + 1L]
			R_t <- Omg * (A_t %*% X_t %*% t(A_t) - Y_t)
			gA <- R_t %*% A_t %*% t(X_t) + t(R_t) %*% A_t %*% X_t
			gA <- gA + mu * A_t
			if (s >= 2L) gA <- gA + lambda * (A_t - A_arr[, , s])
			if (s <= S - 1L) gA <- gA - lambda * (A_arr[, , s + 2L] - A_t)
			g_full[, , s + 1L] <- gA
		}
		# collapse to upper-triangular parameter gradient
		v_grad <- numeric(npar_per_t * S)
		for (s in 1:S) {
			G <- g_full[, , s + 1L]
			# symmetric param gradient: off-diag upper = G[i,j] + G[j,i];
			# diagonal = G[i,i]
			G_sym <- G + t(G)
			diag(G_sym) <- diag(G)
			v_grad[((s - 1L) * npar_per_t + 1L):(s * npar_per_t)] <- G_sym[idx_upper]
		}
		v_grad
	}

	v0 <- pack_A(A)
	obj_init <- obj_fn(v0)
	if (verbose) cli::cli_inform("Symmetric L-BFGS init: obj = {sprintf('%.4f', obj_init)}")

	# for binary/ordinal: wrap L-BFGS inside an outer Z-update loop (ECM).
	# for gaussian: single L-BFGS run.
	if (family == "gaussian" || is.null(Y_obs)) {
		opt <- optim(v0, fn = obj_fn, gr = grad_fn, method = "L-BFGS-B",
			control = list(maxit = max_iter, factr = 1 / tol_obj))
		A_final <- unpack_A(opt$par)
	} else {
		# ECM loop: alternate L-BFGS on A trajectory and Z update
		A_curr <- A
		v_curr <- v0
		for (k in seq_len(z_outer_iter)) {
			opt <- optim(v_curr, fn = obj_fn, gr = grad_fn, method = "L-BFGS-B",
				control = list(maxit = ceiling(max_iter / z_outer_iter), factr = 1 / tol_obj))
			A_curr <- unpack_A(opt$par)
			v_curr <- opt$par
			# Z update from current theta predictions
			for (t in 1:Tt) {
				if (t == 1L) {
					theta_t <- M[, , 1] + Phi[, , 1]
				} else {
					A_t <- A_curr[, , t]
					theta_t <- M[, , 1] + A_t %*% Phi[, , t - 1L] %*% t(A_t)
				}
				Y_filled[, , 1, t] <- .expected_Z(Y_obs[, , 1, t], theta_t, family)
				Phi[, , t] <- Y_filled[, , 1, t] - M[, , 1]
				Phi[, , t][Omega[, , t] == 0] <- 0
			}
		}
		A_final <- A_curr
	}

	# sign canonicalization: walk t=3..T, choose s_t in {-1,+1} to minimize
	# ||s_t A_t - A_{t-1}||_F^2 (anchor s_2 = +1)
	for (t in 3:Tt) {
		ip <- sum(A_final[, , t] * A_final[, , t - 1L])
		if (ip < 0) A_final[, , t] <- -A_final[, , t]
	}

	# compute final objective and stationarity
	obj_final <- obj_fn(pack_A(A_final))

	list(
		A = A_final,
		M = M,
		Phi = Phi,
		Omega = Omega,
		objective = c(obj_init, obj_final),
		converged = opt$convergence == 0L,
		n_iter = if (!is.null(opt$counts)) opt$counts[1] else NA,
		lambda = lambda, mu = mu,
		opt_message = opt$message
	)
}
