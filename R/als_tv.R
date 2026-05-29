####
# rw(1)-smoothed time-varying als for the directed gaussian case.
# the outer loop is: closed-form m-update, block-tridiagonal a- and
# b-trajectory updates (thomas solve in c++), z imputation, scale-gauge
# normalisation, sign-continuity dp, convergence checks.
####

#' Closed-form M-update conditional on operators and Z
#'
#' For each (i, j, r), the weighted mean of (Z_{t,ij,r} - A_t Phi_{t-1} B_t^T)
#' across observed transitions, plus prior shrinkage toward M_0 with variance g2.
#'
#' @param Z 4D latent array [n_row, n_col, p, Tt]
#' @param A_arr 3D cube [n_row, n_row, Tt]
#' @param B_arr 3D cube [n_col, n_col, Tt]
#' @param Phi 4D centered-state array [n_row, n_col, p, Tt] (current iterate)
#' @param Omega 3D mask array [n_row, n_col, Tt] (per-time mask; shared across p)
#' @param M_0 prior mean for M [n_row, n_col, p]
#' @param g2 prior variance for M
#' @return updated M [n_row, n_col, p]
#' @keywords internal
#' @noRd
.tv_als_update_M <- function(Z, A_arr, B_arr, Phi, Omega, M_0, g2) {
	dims <- dim(Z)
	n_row <- dims[1]; n_col <- dims[2]; p <- dims[3]; Tt <- dims[4]
	M_new <- array(0, dim = c(n_row, n_col, p))

	# weight per cell: count of observed transitions + 1/g2
	# response at t=1 contributes (Z_1 - M) with weight = mask_1
	# response at t>=2 contributes (Z_t - A_t (Z_{t-1} - M) B_t^T - M)... but
	# this couples M into A_t M B_t^T which is intractable in closed form.
	# pragmatic simplification: use the "Theta_t = M + Phi_t" decomposition
	# where Phi_t is the current centered-state iterate. Then E[Z_t - Phi_t] = M
	# given the current Phi. This is the "M | Phi" conditional and is exactly
	# what the MCMC's m-step uses (see R/models.R update_mu_dynamic).

	# numerator: sum over t of mask_{t,ij} * (Z_{t,ij,r} - Phi_{t,ij,r})
	# denominator: sum over t of mask_{t,ij} + 1/g2
	for (r in 1:p) {
		num <- matrix(0, n_row, n_col)
		den <- matrix(0, n_row, n_col)
		for (t in 1:Tt) {
			mask_t <- Omega[, , t]
			Z_t <- Z[, , r, t]
			Phi_t <- Phi[, , r, t]
			diff_t <- Z_t - Phi_t
			diff_t[is.na(diff_t)] <- 0
			num <- num + mask_t * diff_t
			den <- den + mask_t
		}
		num <- num + M_0[, , r] / g2
		den <- den + 1 / g2
		M_new[, , r] <- num / den
	}
	M_new
}

#' Construct the lag-state Phi^state (handles NA on the structural diagonal)
#'
#' Per the plan section3, the lag state cannot contain NA values (they would
#' propagate through M_t M_t^T). Default rule: zero out the structural
#' diagonal of square networks.
#'
#' @param Phi 4D centered-state array
#' @param Omega 3D mask (1 = observed, 0 = missing or structural diagonal)
#' @return 4D state array with no NAs
#' @keywords internal
#' @noRd
.tv_als_make_state <- function(Phi, Omega) {
	state <- Phi
	# zero out where mask = 0 (structural NA in state means "no contribution")
	dims <- dim(state)
	p <- dims[3]
	Tt <- dims[4]
	for (r in 1:p) for (t in 1:Tt) {
		state_rt <- state[, , r, t]
		state_rt[Omega[, , t] == 0] <- 0
		state_rt[is.na(state_rt)] <- 0
		state[, , r, t] <- state_rt
	}
	state
}

#' Penalty-minimizing scale gauge for the directed asymmetric case
#'
#' Solves min_u P(u) where A_t(u) = exp(u_t) A_t, B_t(u) = exp(-u_t) B_t
#' and P(u) is the RW + level penalty. Cheap Newton on the (Tt-1)-dim problem.
#'
#' @param A_arr cube [n_row, n_row, Tt]
#' @param B_arr cube [n_col, n_col, Tt]
#' @param lambda_A,lambda_B,mu_A,mu_B penalties
#' @param max_iter Newton iters
#' @return list with rescaled A_arr, B_arr
#' @keywords internal
#' @noRd
.tv_als_gauge_penalty_min <- function(A_arr, B_arr, lambda_A, lambda_B,
                                       mu_A, mu_B, max_iter = 5L) {
	Tt <- dim(A_arr)[3]
	S <- Tt - 1L
	if (S < 1L) return(list(A_arr = A_arr, B_arr = B_arr))

	# precompute scalars used in the penalty.
	# let aA_s = ||A_{s+1}||_F^2, bB_s = ||B_{s+1}||_F^2 (s = 0..S-1)
	# let inner_A_s = <A_{s+1}, A_s>_F, inner_B_s = <B_{s+1}, B_s>_F  (s = 1..S-1)
	aA <- numeric(S); bB <- numeric(S)
	inner_A <- numeric(S); inner_B <- numeric(S)  # inner_A[s] = <A_{s+2}, A_{s+1}>
	for (s in seq_len(S)) {
		aA[s] <- sum(A_arr[, , s + 1L]^2)
		bB[s] <- sum(B_arr[, , s + 1L]^2)
	}
	if (S >= 2L) for (s in seq_len(S - 1L)) {
		inner_A[s] <- sum(A_arr[, , s + 1L] * A_arr[, , s + 2L])
		inner_B[s] <- sum(B_arr[, , s + 1L] * B_arr[, , s + 2L])
	}

	# objective (additive constant dropped):
	# P(u) = sum_s [0.5 mu_A exp(2 u_s) aA_s + 0.5 mu_B exp(-2 u_s) bB_s]
	#      + sum_{s=1..S-1} [0.5 lambda_A (exp(2u_s) aA_s + exp(2 u_{s-1}) aA_{s-1}
	#                                       - 2 exp(u_s + u_{s-1}) inner_A_{s-1})
	#                       + 0.5 lambda_B (exp(-2u_s) bB_s + exp(-2 u_{s-1}) bB_{s-1}
	#                                       - 2 exp(-u_s - u_{s-1}) inner_B_{s-1})]
	# use 1-based indexing: u has length S (one u per t=2..T).

	# gradient and Hessian (tridiagonal): Newton with a numerical Hessian
	# is adequate at this scale (S is small).
	u <- rep(0, S)
	objective_u <- function(u_vec) {
		# level
		level_A <- sum(0.5 * mu_A * exp(2 * u_vec) * aA)
		level_B <- sum(0.5 * mu_B * exp(-2 * u_vec) * bB)
		# RW
		rw_A <- 0; rw_B <- 0
		if (S >= 2L) {
			for (s in seq_len(S - 1L)) {
				# diff between A_{s+2} (u_{s+1}) and A_{s+1} (u_s) -- but 1-based
				u_curr <- u_vec[s + 1L]; u_prev <- u_vec[s]
				rw_A <- rw_A + 0.5 * lambda_A * (exp(2 * u_curr) * aA[s + 1L] +
					exp(2 * u_prev) * aA[s] - 2 * exp(u_curr + u_prev) * inner_A[s])
				rw_B <- rw_B + 0.5 * lambda_B * (exp(-2 * u_curr) * bB[s + 1L] +
					exp(-2 * u_prev) * bB[s] - 2 * exp(-u_curr - u_prev) * inner_B[s])
			}
		}
		level_A + level_B + rw_A + rw_B
	}

	# newton via L-BFGS-B for robustness on this small problem. If the
	# gauge optimization itself fails we keep u = 0 (which is always a
	# feasible no-op rebalance) -- the fit is then identical to "no gauge".
	# surface the error class so a *systematic* failure is visible.
	if (max_iter > 0L && S >= 1L) {
		opt <- tryCatch(
			optim(u, objective_u, method = "L-BFGS-B",
				control = list(maxit = max_iter, factr = 1e10)),
			error = function(e) {
				cli::cli_warn(c(
					"Penalty-minimizing gauge step failed; keeping the no-rescale fit.",
					"x" = "{conditionMessage(e)}",
					"i" = "Try {.code gauge = \"equal_norm_conditional\"} or {.code gauge = \"none\"} if this happens repeatedly."
				))
				NULL
			}
		)
		if (!is.null(opt) && is.finite(opt$value)) u <- opt$par
	}

	# apply the gauge to t = 2..T (penalty-minimizing optimization runs over
	# these S = T-1 indices).
	for (s in seq_len(S)) {
		A_arr[, , s + 1L] <- exp(u[s]) * A_arr[, , s + 1L]
		B_arr[, , s + 1L] <- exp(-u[s]) * B_arr[, , s + 1L]
	}
	# rebalance t = 1 to equal Frobenius norm. A_1 / B_1 are not used in any
	# transition prediction (t = 1 has no lag), so rescaling them does not
	# change the data fit. But un-rebalanced t = 1 slices left over from the
	# static warm-start cause `dbn_operator(fit)` to show wildly different
	# ||A_1|| vs ||B_1|| in dbn_operator(fit).
	nA1 <- sqrt(sum(A_arr[, , 1]^2))
	nB1 <- sqrt(sum(B_arr[, , 1]^2))
	if (nA1 > 1e-10 && nB1 > 1e-10) {
		r1 <- sqrt(nB1 / nA1)
		A_arr[, , 1] <- A_arr[, , 1] * r1
		B_arr[, , 1] <- B_arr[, , 1] / r1
	}
	list(A_arr = A_arr, B_arr = B_arr, u = u)
}

#' Cheap equal-norm rebalance (conditional-accept fallback)
#'
#' Sets ||A_t||_F = ||B_t||_F per t. Always preserves data fit; may
#' increase the RW penalty. The caller is responsible for accepting only
#' if the full objective does not increase.
#' @keywords internal
#' @noRd
.tv_als_gauge_equal_norm <- function(A_arr, B_arr) {
	Tt <- dim(A_arr)[3]
	for (t in 2:Tt) {
		nA <- sqrt(sum(A_arr[, , t]^2))
		nB <- sqrt(sum(B_arr[, , t]^2))
		if (nA > 1e-10 && nB > 1e-10) {
			r <- sqrt(nB / nA)
			A_arr[, , t] <- A_arr[, , t] * r
			B_arr[, , t] <- B_arr[, , t] / r
		}
	}
	list(A_arr = A_arr, B_arr = B_arr)
}

#' Sign continuity via 2-state Viterbi DP
#'
#' Solves
#'   min_{s_2..s_T in {-1,+1}} sum_{t=3..T} lambda_A * ||s_t A_t - s_{t-1} A_{t-1}||^2
#'                                          + lambda_B * ||s_t B_t - s_{t-1} B_{t-1}||^2
#' anchored at s_2 = +1.
#' @keywords internal
#' @noRd
.tv_als_sign_dp <- function(A_arr, B_arr, lambda_A, lambda_B) {
	Tt <- dim(A_arr)[3]
	S <- Tt - 1L  # number of operator time points t=2..T
	if (S < 2L) return(list(A_arr = A_arr, B_arr = B_arr, signs = c(1)))

	# DP: V_s[k] = min cost to reach time s with sign k (k=1 for -1, k=2 for +1)
	# anchor: V_0 fixed at sign +1
	# link cost c_s(s_prev, s_curr) = lambda_A * ||s_curr A_curr - s_prev A_prev||^2
	#                                 + lambda_B * ||s_curr B_curr - s_prev B_prev||^2
	# note: ||sA - s'A'||^2 = ||A||^2 + ||A'||^2 - 2 s s' <A, A'>
	#   so the dependence is only on s * s' (-1 or +1).
	# so min_{s_curr} depends on whether s_curr * s_prev = +1 or -1.

	sign_vals <- c(-1, +1)
	V <- matrix(Inf, nrow = S, ncol = 2)
	backptr <- matrix(0L, nrow = S, ncol = 2)
	# anchor: at s=1 (i.e. t=2), force sign +1
	V[1, 2] <- 0  # +1
	V[1, 1] <- Inf  # -1 disallowed at anchor

	for (s in 2:S) {
		A_curr <- A_arr[, , s + 1L]
		A_prev <- A_arr[, , s]
		B_curr <- B_arr[, , s + 1L]
		B_prev <- B_arr[, , s]
		# precompute inner products (sign-independent factors)
		na2_curr <- sum(A_curr^2)
		na2_prev <- sum(A_prev^2)
		nb2_curr <- sum(B_curr^2)
		nb2_prev <- sum(B_prev^2)
		ip_A <- sum(A_curr * A_prev)
		ip_B <- sum(B_curr * B_prev)
		for (k_curr in 1:2) for (k_prev in 1:2) {
			s_curr <- sign_vals[k_curr]
			s_prev <- sign_vals[k_prev]
			link_A <- na2_curr + na2_prev - 2 * s_curr * s_prev * ip_A
			link_B <- nb2_curr + nb2_prev - 2 * s_curr * s_prev * ip_B
			cost <- V[s - 1, k_prev] + lambda_A * link_A + lambda_B * link_B
			if (cost < V[s, k_curr]) {
				V[s, k_curr] <- cost
				backptr[s, k_curr] <- k_prev
			}
		}
	}

	# backtrack
	signs <- integer(S)
	signs[S] <- which.min(V[S, ])
	if (S > 1L) for (s in (S - 1L):1L) signs[s] <- backptr[s + 1L, signs[s + 1L]]
	sign_seq <- sign_vals[signs]
	# apply
	for (s in seq_len(S)) {
		A_arr[, , s + 1L] <- sign_seq[s] * A_arr[, , s + 1L]
		B_arr[, , s + 1L] <- sign_seq[s] * B_arr[, , s + 1L]
	}
	list(A_arr = A_arr, B_arr = B_arr, signs = sign_seq)
}

#' Stage 1A: directed Gaussian RW-smoothed TV-ALS at fixed lambda
#'
#' @param Y 4D array [n_row, n_col, p, Tt]; the working response (for Gaussian, Y = Z)
#' @param family character, must be "gaussian" in Stage 1A
#' @param symmetric logical, must be FALSE in Stage 1A
#' @param lambda RW penalty (scalar)
#' @param mu level penalty (scalar)
#' @param max_iter outer-loop iteration cap
#' @param tol_obj relative objective decrease tolerance
#' @param tol_par per-period parameter movement tolerance
#' @param gauge one of "penalty_min", "equal_norm_conditional", "none"
#' @param init list with A_init, B_init, M_init; if NULL, warm-start from static
#' @param prior_M_mean,prior_M_var hyperparameters for M (default zero/1.0)
#' @param verbose integer/logical
#' @return list with A, B, M, Phi, objective trajectory, convergence diagnostics
#' @keywords internal
#' @noRd
.dbn_als_tv_fixed_lambda <- function(Y, family, symmetric,
                                       lambda, mu,
                                       max_iter = 200L, warn_iter = 100L,
                                       tol_obj = 1e-5, tol_par = 1e-3,
                                       gauge = c("penalty_min", "equal_norm_conditional", "none"),
                                       init = NULL,
                                       prior_M_mean = NULL, prior_M_var = 1.0,
                                       check_monotonicity = TRUE,
                                       Y_obs = NULL,
                                       z_update = c("mean", "sample"),
                                       verbose = FALSE) {
	gauge <- match.arg(gauge)
	z_update <- match.arg(z_update)
	stopifnot(!symmetric)
	# for binary/ordinal, Y is the LATENT Z working response and Y_obs is
	# the original observed Y. The latent Z is re-imputed via .expected_Z
	# at the end of each outer iteration.

	dims <- dim(Y)
	n_row <- dims[1]; n_col <- dims[2]; p <- dims[3]; Tt <- dims[4]

	# TV-ALS needs T >= 2 (at least one transition); T = 1 makes the
	# `for (t in 2:Tt)` loop reverse to `c(2, 1)` and crash.
	if (Tt < 2L) {
		cli::cli_abort(c(
			"TV-ALS requires at least 2 time points (got T = {Tt}).",
			"i" = "Use {.code lambda = \"static\"} for a single-period fit."
		))
	}

	# build response masks: Omega[i,j,t] = 1 where observed; 0 else.
	# for square networks, structural diagonal is 0.
	Omega <- array(1, dim = c(n_row, n_col, Tt))
	if (n_row == n_col) for (t in 1:Tt) diag(Omega[, , t]) <- 0
	# also propagate NA from the input
	for (r in 1:p) for (t in 1:Tt) {
		na_mask <- is.na(Y[, , r, t])
		Omega[, , t][na_mask] <- 0
	}
	# average mask over relations (we use a single mask across p)
	# if the user's mask varies across relations, we conservatively zero
	# where ANY relation is missing. (This is a limitation but matches the
	# rest of the package's per-time-but-not-per-relation mask convention.)

	# build a "filled" working response with NAs replaced by 0 (they'll be
	# multiplied by the mask anyway in the data-fit term)
	Y_filled <- Y
	Y_filled[is.na(Y_filled)] <- 0

	# initialize from static ALS warm-start if no init provided. Pass the
	# original Y (with NAs intact) to the static fit so its internal NA
	# handling can fire correctly; the TV-internal Y_filled is only used
	# for the per-iteration solves below.
	if (is.null(init)) {
		init <- .tv_als_warm_start_static(Y, family, symmetric, Omega, n_row, n_col, p, Tt)
	}

	# validate init
	A_arr <- init$A_init
	B_arr <- init$B_init
	M     <- init$M_init
	stopifnot(identical(dim(A_arr), c(n_row, n_row, Tt)))
	stopifnot(identical(dim(B_arr), c(n_col, n_col, Tt)))
	stopifnot(identical(dim(M), c(n_row, n_col, p)))

	# prior on M
	if (is.null(prior_M_mean)) prior_M_mean <- array(0, dim = c(n_row, n_col, p))
	g2 <- prior_M_var

	# state: Phi[,,r,t] = Y[,,r,t] - M[,,r] for Gaussian
	make_Phi <- function(M) {
		Phi <- Y_filled
		for (r in 1:p) for (t in 1:Tt) Phi[, , r, t] <- Y_filled[, , r, t] - M[, , r]
		Phi
	}

	# for the TV solver we need a 3D Phi cube (n_row, n_col, Tt) by averaging
	# over relations (or running per-relation -- for Stage 1A we assume p=1
	# for simplicity and the multi-p extension goes through unchanged).
	make_Phi3d <- function(Phi) {
		if (p == 1L) return(Phi[, , 1, ])
		apply(Phi, c(1, 2, 4), mean)
	}

	Phi <- make_Phi(M)
	Phi3d <- make_Phi3d(Phi)
	# state version: replace NAs / masked entries with 0
	Phi3d_state <- Phi3d
	for (t in 1:Tt) {
		st <- Phi3d_state[, , t]
		st[Omega[, , t] == 0] <- 0
		st[is.na(st)] <- 0
		Phi3d_state[, , t] <- st
	}

	# initial objective
	obj_prev <- compute_tv_als_objective_cpp(Phi3d_state, A_arr, B_arr, Omega,
		lambda, lambda, mu, mu)

	obj_history <- numeric(0)
	obj_history <- c(obj_history, obj_prev)
	converged <- FALSE
	if (verbose) cli::cli_inform("Iter 0: obj = {round(obj_prev, 6)}")

	for (iter in seq_len(max_iter)) {
		A_prev <- A_arr
		B_prev <- B_arr

		# (a) M update -- DEFERRED: for Gaussian observed data, the M
		# conditional given Z and Phi is degenerate because Phi = Z - M by
		# construction. The package convention (matching the static SIR and
		# the MCMC's update_mu_dynamic) is to fix M from preprocessing.
		# stage 6 may revisit a bilinear-residual M update if needed.
		# M <- .tv_als_update_M(Y_filled, A_arr, B_arr, Phi, Omega, prior_M_mean, g2)
		# phi already constructed from initial M; do not recompute here.

		# (c) A trajectory
		A_arr <- solve_A_trajectory_cpp(Phi3d_state, B_arr, Omega, lambda, mu)

		# (d) B trajectory
		B_arr <- solve_B_trajectory_cpp(Phi3d_state, A_arr, Omega, lambda, mu)

		# (e) Z imputation for binary/ordinal: deterministic conditional mean
		if (family != "gaussian" && !is.null(Y_obs) && z_update == "mean") {
			# compute current predicted mean Theta_t = M + A_t Phi_{t-1} B_t^T
			# note: Phi3d_state[, , t] is the centered state; the predicted Z
			# at time t (conditional on operators) is M + (pred at t).
			# we re-impute Z at each cell; the resulting Y_filled is updated.
			for (r in 1:p) {
				for (t in 1:Tt) {
					if (t == 1L) {
						theta_t <- M[, , r] + Phi3d_state[, , 1]
					} else {
						pred_t <- A_arr[, , t] %*% Phi3d_state[, , t - 1L] %*% t(B_arr[, , t])
						theta_t <- M[, , r] + pred_t
					}
					Y_filled[, , r, t] <- .expected_Z(Y_obs[, , r, t], theta_t, family)
				}
			}
			# update Phi3d_state from new Y_filled
			Phi <- make_Phi(M)
			Phi3d <- make_Phi3d(Phi)
			Phi3d_state <- Phi3d
			for (t in 1:Tt) {
				st <- Phi3d_state[, , t]
				st[Omega[, , t] == 0] <- 0
				st[is.na(st)] <- 0
				Phi3d_state[, , t] <- st
			}
		}

		# (f) gauge
		if (gauge == "penalty_min") {
			gauged <- .tv_als_gauge_penalty_min(A_arr, B_arr, lambda, lambda, mu, mu)
			A_arr <- gauged$A_arr; B_arr <- gauged$B_arr
		} else if (gauge == "equal_norm_conditional") {
			gauged <- .tv_als_gauge_equal_norm(A_arr, B_arr)
			# conditional accept: only if objective doesn't increase
			obj_try <- compute_tv_als_objective_cpp(Phi3d_state, gauged$A_arr,
				gauged$B_arr, Omega, lambda, lambda, mu, mu)
			if (obj_try <= obj_prev + 1e-10) {
				A_arr <- gauged$A_arr; B_arr <- gauged$B_arr
			}
		}

		# (g) sign continuity
		signed <- .tv_als_sign_dp(A_arr, B_arr, lambda, lambda)
		A_arr <- signed$A_arr; B_arr <- signed$B_arr

		# (h) convergence checks
		obj_curr <- compute_tv_als_objective_cpp(Phi3d_state, A_arr, B_arr,
			Omega, lambda, lambda, mu, mu)
		obj_history <- c(obj_history, obj_curr)

		# defensive monotonicity check (gaussian case only; with Z imputation
		# the objective on the re-imputed Z is different from the previous Z
		# objective, so monotonicity is not guaranteed there)
		if (check_monotonicity && family == "gaussian") {
			if (obj_curr > obj_prev * (1 + 1e-8) + 1e-10) {
				cli::cli_warn(c(
					"TV-ALS objective increased between iterations (iter {iter}: {obj_prev} -> {obj_curr}).",
					"i" = "This may indicate a bug in the gauge step or mask handling."
				))
			}
		}

		rel_obj <- abs(obj_curr - obj_prev) / (1 + abs(obj_prev))
		# per-period param movement (scale-aware)
		max_delta <- 0
		for (t in 2:Tt) {
			da <- sqrt(sum((A_arr[, , t] - A_prev[, , t])^2)) /
				(1 + sqrt(sum(A_prev[, , t]^2)))
			db <- sqrt(sum((B_arr[, , t] - B_prev[, , t])^2)) /
				(1 + sqrt(sum(B_prev[, , t]^2)))
			max_delta <- max(max_delta, da + db)
		}

		if (verbose && (iter %% 10 == 0 || iter == 1)) {
			cli::cli_inform("Iter {iter}: obj = {sprintf('%.6f', obj_curr)}, rel = {sprintf('%.3e', rel_obj)}, delta = {sprintf('%.3e', max_delta)}")
		}

		if (rel_obj < tol_obj && max_delta < tol_par) {
			converged <- TRUE
			if (verbose) cli::cli_inform("Converged at iter {iter}.")
			break
		}

		obj_prev <- obj_curr
	}

	if (!converged) {
		# unconditional warning: convergence failure is a fit-quality signal
		# the user should always see, irrespective of `verbose`.
		cli::cli_warn(c(
			"TV-ALS did not converge in {max_iter} iterations.",
			"i" = "Increase {.arg tv_max_iter} (currently {max_iter}) or loosen {.arg tv_tol_obj}.",
			"i" = "Inspect {.code fit$meta$tv} for the objective trajectory."
		))
	}

	# final stationarity residual
	stat_res <- .tv_als_stationarity_residual(Phi3d_state, A_arr, B_arr, Omega,
		lambda, lambda, mu, mu)

	list(
		A = A_arr, B = B_arr, M = M,
		Phi = Phi3d_state,
		Omega = Omega,
		objective = obj_history,
		converged = converged,
		n_iter = iter,
		stationarity_residual = stat_res,
		lambda = lambda, mu = mu
	)
}

#' Static-ALS warm start, broadcast across all t
#' @keywords internal
#' @noRd
.tv_als_warm_start_static <- function(Y_filled, family, symmetric, Omega,
                                       n_row, n_col, p, Tt) {
	# use the existing static ALS solver (dbn_als with no bootstrap, no lambda)
	# A more careful version would call als_core_fit directly to avoid the
	# preprocessing overhead -- skip for now.
	# replace masked entries with 0 (they have zero weight in static ALS too).
	Y_for_static <- Y_filled

	# we need to give static ALS a Y array. For Stage 1A (gaussian, p=1) we
	# can just call dbn_als with bootstrap = 0.
	static_fit <- dbn_als(Y_for_static, family = family, symmetric = symmetric,
		bootstrap = 0, verbose = FALSE)
	# extract static A, B (in static fit, these are [n, n, Tt] cubes but with
	# every slice the same). Just take the first.
	A0 <- static_fit$A[[1]][, , 1]
	B0 <- static_fit$B[[1]][, , 1]
	M0 <- static_fit$M[, , , 1, drop = FALSE]
	dim(M0) <- c(n_row, n_col, p)

	A_arr <- array(rep(c(A0), Tt), dim = c(n_row, n_row, Tt))
	B_arr <- array(rep(c(B0), Tt), dim = c(n_col, n_col, Tt))
	# slice 1 (t=1) is unused but we tile it with the static est anyway
	list(A_init = A_arr, B_init = B_arr, M_init = M0)
}

#' Compute stationarity residual norm
#' @keywords internal
#' @noRd
.tv_als_stationarity_residual <- function(Phi, A_arr, B_arr, Omega,
                                          lambda_A, lambda_B, mu_A, mu_B) {
	Tt <- dim(Phi)[3]
	S <- Tt - 1
	if (S < 1) return(0)
	# gradient norm: sum_t ||grad_A F||^2 + ||grad_B F||^2
	grad_norm2 <- 0
	for (s in 1:S) {
		A_t <- A_arr[, , s + 1]
		B_t <- B_arr[, , s + 1]
		Phi_prev <- Phi[, , s]
		Phi_t <- Phi[, , s + 1]
		resid <- Omega[, , s + 1] * (Phi_t - A_t %*% Phi_prev %*% t(B_t))
		# gradient wrt A_t. the data-fit term f = 0.5 * ||Phi - A X B^T||^2
		# gives df/dA = -(resid %*% B) %*% t(X) with X = Phi_prev. the
		# matrix orientation here handles n_row != n_col (bipartite).
		grad_A <- -(resid %*% B_t) %*% t(Phi_prev) + mu_A * A_t
		if (s >= 2) grad_A <- grad_A + lambda_A * (A_t - A_arr[, , s])
		if (s <= S - 1) grad_A <- grad_A + lambda_A * (A_t - A_arr[, , s + 2])
		# gradient wrt B_t: df/dB = -t(resid) %*% A %*% X
		grad_B <- -t(resid) %*% A_t %*% Phi_prev + mu_B * B_t
		if (s >= 2) grad_B <- grad_B + lambda_B * (B_t - B_arr[, , s])
		if (s <= S - 1) grad_B <- grad_B + lambda_B * (B_t - B_arr[, , s + 2])
		grad_norm2 <- grad_norm2 + sum(grad_A^2) + sum(grad_B^2)
	}
	param_norm2 <- 0
	for (s in 1:S) {
		param_norm2 <- param_norm2 + sum(A_arr[, , s + 1]^2) + sum(B_arr[, , s + 1]^2)
	}
	sqrt(grad_norm2) / (1 + sqrt(param_norm2))
}
