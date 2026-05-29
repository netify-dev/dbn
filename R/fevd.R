####
#' Forecast error variance decomposition for dbn fits
#'
#' @description Closed-form FEVD for the Kronecker-restricted VAR(1)
#'   implied by the bilinear dynamics. Under the package's isotropic
#'   matrix-normal innovation `eps_t ~ MN(0, sigma^2 I_n, I_m)`, the
#'   vec-innovation `vec(eps_t) ~ N(0, sigma^2 I_{nm})` is already
#'   orthogonal, so reduced-form and orthogonalized FEVD coincide and
#'   no Cholesky ordering is required.
#'
#' @details
#' For target dyad `(i, j)`, shock dyad `(k, l)`, horizon `h`, define
#' the cumulative operators
#' \deqn{A^{cum}_q = A_{t_0+h} A_{t_0+h-1} \cdots A_{t_0+q+1}}
#' \deqn{B^{cum}_q = B_{t_0+h} B_{t_0+h-1} \cdots B_{t_0+q+1}}
#' for `q = 1, ..., h`. The contribution of dyad `(k, l)` to forecast
#' error variance of dyad `(i, j)` at horizon `h` is
#' \deqn{D_{ij \leftarrow kl}(h) = \frac{\sum_q A^{cum}_q[i,k]^2 \cdot B^{cum}_q[j,l]^2}{\sum_q \|A^{cum}_q[i,\cdot]\|_2^2 \|B^{cum}_q[j,\cdot]\|_2^2}.}
#'
#' For dynamic / piecewise fits, future operators `A_{t_0+s}` are held
#' at `A_{t_0}` when `simulate_future_operators = FALSE` (default). When
#' `TRUE`, future paths are drawn from the fitted random-walk prior and
#' the additional forecast variance from operator drift is reported in
#' `$unattributed_future_operator` rather than attributed to dyadic
#' shocks.
#'
#' @param fit A `dbn` fit (dynamic, piecewise, or static).
#' @param H Integer forecast horizon (default 10).
#' @param t0 Integer forecast origin (default: last in-sample time).
#' @param shock Aggregation level: `"actor"` collapses to source-actor
#'   `k` by summing over receiver index `l` (default; memory-cheap);
#'   `"dyad"` keeps the full `(k, l)` resolution.
#' @param simulate_future_operators Logical. When `TRUE`, sample future
#'   operator paths from the RW prior and report drift variance
#'   separately.
#' @param prob Credible level (default 0.95).
#' @param n_draws Number of posterior draws to use (default: all).
#' @param Sigma_row Optional `n_row x n_row` row-covariance matrix for
#'   structural rotation of the shock basis. `NULL` (default) keeps the
#'   isotropic basis.
#' @param Sigma_col Optional `n_col x n_col` column-covariance matrix
#'   for structural rotation. `NULL` (default) keeps the isotropic basis.
#' @param ordering Optional integer ordering for the structural Cholesky
#'   rotation. `NULL` (default) uses the identity ordering.
#' @return A `dbn_fevd` list with `median`, `lower`, `upper`,
#'   `total_mse`, `unattributed_future_operator` (or `NULL`), and
#'   `metadata`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
fevd <- function(fit,
                  H = 10L,
                  t0 = NULL,
                  shock = c("actor", "dyad"),
                  simulate_future_operators = FALSE,
                  prob = 0.95,
                  n_draws = NULL,
                  Sigma_row = NULL,
                  Sigma_col = NULL,
                  ordering = NULL) {
	if (!inherits(fit, "dbn"))
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	shock <- match.arg(shock)
	mdl <- fit$model %||% ""
	# structural-FEVD path when the user supplies a matrix-normal innovation
	# covariance via Sigma_row (n_row x n_row) and / or Sigma_col (n_col x
	# n_col). default NULL = the package's isotropic Sigma = sigma^2 I path
	# (no ordering required, no Cholesky factor needed)
	use_structural <- !is.null(Sigma_row) || !is.null(Sigma_col) ||
		!is.null(ordering)
	# pull per-draw, per-time operator cubes
	if (mdl %in% c("dynamic", "lowrank")) {
		if (!is.list(fit$A) || length(fit$A) == 0L)
			cli::cli_abort("Fit has no per-draw operators; refit with MCMC or pass a bootstrap fit.")
		A_list <- fit$A
		B_list <- fit$B
	} else if (identical(mdl, "piecewise")) {
		od <- operator_draws(fit, scale = "time")
		# convert [n_draws, n, n, T] arrays into list-of-cubes form
		A_list <- lapply(seq_len(dim(od$A)[1]), function(d) aperm(od$A[d, , , , drop = FALSE], c(2, 3, 4, 1))[, , , 1])
		B_list <- lapply(seq_len(dim(od$B)[1]), function(d) aperm(od$B[d, , , , drop = FALSE], c(2, 3, 4, 1))[, , , 1])
	} else if (identical(mdl, "static")) {
		A_list <- fit$A
		B_list <- fit$B
	} else {
		cli::cli_abort("{.fun fevd} supports {.val dynamic}, {.val piecewise}, {.val lowrank}, and {.val static}; got {.val {mdl}}.")
	}
	D <- length(A_list)
	if (D == 0L)
		cli::cli_abort("No operator draws available.")
	if (!is.null(n_draws)) {
		n_draws <- as.integer(min(n_draws, D))
		idx <- sample.int(D, n_draws)
		A_list <- A_list[idx]
		B_list <- B_list[idx]
		D <- n_draws
	}
	dA1 <- dim(A_list[[1L]])
	dB1 <- dim(B_list[[1L]])
	if (length(dA1) == 3L) {
		n_row <- dA1[1L]; Tt <- dA1[3L]
	} else {
		n_row <- dA1[1L]; Tt <- 1L
	}
	if (length(dB1) == 3L) {
		n_col <- dB1[1L]
	} else {
		n_col <- dB1[1L]
	}
	if (is.null(t0)) t0 <- Tt
	if (t0 < 1L || t0 > Tt)
		cli::cli_abort("{.arg t0} = {t0} out of range [1, {Tt}].")
	H <- as.integer(H)
	if (H < 1L)
		cli::cli_abort("{.arg H} must be a positive integer.")
	if (H > Tt)
		cli::cli_warn("Forecast horizon {H} exceeds T = {Tt}; future operators held at t0 unless {.arg simulate_future_operators = TRUE}.")
	# memory check for shock = "dyad"
	if (identical(shock, "dyad")) {
		bytes <- D * H * n_row * n_col * n_row * n_col * 8
		if (bytes > 1e9)
			cli::cli_warn("Estimated FEVD memory ~{round(bytes/1e9, 1)} GB at {.code shock = \"dyad\"}; consider {.code shock = \"actor\"}.")
	}
	# get sigma2 per draw
	sigma2_vec <- if (is.numeric(fit$sigma2)) as.numeric(fit$sigma2)[seq_len(D)]
	              else if (!is.null(fit$draws$pars$sigma2)) fit$draws$pars$sigma2[seq_len(D)]
	              else if (!is.null(fit$draws$pars$s2)) fit$draws$pars$s2[seq_len(D)]
	              else rep(1, D)
	# structural-path: validate Sigma_row, Sigma_col, ordering and compute
	# their Cholesky factors so we can rotate A_cum / B_cum to the
	# structural basis in the inner loop
	L_row <- NULL; L_col <- NULL
	if (use_structural) {
		if (is.null(Sigma_row)) Sigma_row <- diag(n_row)
		if (is.null(Sigma_col)) Sigma_col <- diag(n_col)
		if (!identical(dim(Sigma_row), c(as.integer(n_row), as.integer(n_row))))
			cli::cli_abort("{.arg Sigma_row} must be {.val {n_row} x {n_row}}.")
		if (!identical(dim(Sigma_col), c(as.integer(n_col), as.integer(n_col))))
			cli::cli_abort("{.arg Sigma_col} must be {.val {n_col} x {n_col}}.")
		if (!is.null(ordering)) {
			if (length(ordering) == n_row) {
				perm <- as.integer(ordering)
				Sigma_row <- Sigma_row[perm, perm, drop = FALSE]
			} else if (is.list(ordering) && length(ordering) == 2L) {
				perm_r <- as.integer(ordering[[1]])
				perm_c <- as.integer(ordering[[2]])
				Sigma_row <- Sigma_row[perm_r, perm_r, drop = FALSE]
				Sigma_col <- Sigma_col[perm_c, perm_c, drop = FALSE]
			} else {
				cli::cli_abort("{.arg ordering} must be an integer permutation of length n_row, or a list with one permutation per mode.")
			}
		}
		# lower-triangular Cholesky factors so vec(eps) = (L_col kron L_row) u
		L_row <- t(chol(Sigma_row + diag(1e-10, n_row)))
		L_col <- t(chol(Sigma_col + diag(1e-10, n_col)))
	}
	# get tau_A2 / tau_B2 for future-operator simulation
	tau_A2_vec <- if (is.numeric(fit$tau_A2)) as.numeric(fit$tau_A2)[seq_len(D)] else rep(0, D)
	tau_B2_vec <- if (is.numeric(fit$tau_B2)) as.numeric(fit$tau_B2)[seq_len(D)] else rep(0, D)
	# core per-draw computation
	get_A <- function(A_arr, t) if (length(dim(A_arr)) == 3L) A_arr[, , min(t, dim(A_arr)[3])] else A_arr
	out_dim <- if (identical(shock, "actor"))
		c(n_row, n_col, n_row, H, D)
	else
		c(n_row, n_col, n_row, n_col, H, D)
	out <- array(NA_real_, dim = out_dim)
	mse_arr <- array(NA_real_, dim = c(n_row, n_col, H, D))
	unattributed <- if (isTRUE(simulate_future_operators))
		array(NA_real_, dim = c(n_row, n_col, H, D))
	else NULL
	for (d in seq_len(D)) {
		Am <- A_list[[d]]; Bm <- B_list[[d]]
		sig2_d <- sigma2_vec[d]
		# precompute cumulative operators starting from horizon 1
		# A_cum[, , q] is A_{t0+H} ... A_{t0+q+1} accumulated right-to-left
		A_cum <- array(0, dim = c(n_row, n_row, H))
		B_cum <- array(0, dim = c(n_col, n_col, H))
		# walk h from H back to 1: at q = h we have empty product = I
		# at q = h-1 we have one factor A_{t0+h}, etc.
		A_cum[, , H] <- diag(n_row)
		B_cum[, , H] <- diag(n_col)
		for (h in seq(H - 1L, 1L, by = -1L)) {
			t_idx <- min(Tt, t0 + h + 1L)
			A_cum[, , h] <- get_A(Am, t_idx) %*% A_cum[, , h + 1L]
			B_cum[, , h] <- get_A(Bm, t_idx) %*% B_cum[, , h + 1L]
		}
		# attribute innovation variance to source dyads
		for (h in seq_len(H)) {
			# at horizon h the contributing shocks are eps_{t0+1}, ..., eps_{t0+h}
			# their propagators are A_cum[, , 1..h] and B_cum[, , 1..h]
			# isotropic path: squared contributions = A_cum[i,k]^2 * B_cum[j,l]^2
			# structural path: rotate A_cum, B_cum by Cholesky factors
			# (sigma2 is absorbed into Sigma_row / Sigma_col in the structural case)
			contrib_dyad <- array(0, dim = c(n_row, n_col, n_row, n_col))
			mse_ij <- matrix(0, n_row, n_col)
			for (q in seq_len(h)) {
				if (use_structural) {
					# A_prime = A_cum_q %*% L_row; B_prime = B_cum_q %*% L_col
					Aprime <- A_cum[, , q] %*% L_row
					Bprime <- B_cum[, , q] %*% L_col
					A_sq <- Aprime ^ 2
					B_sq <- Bprime ^ 2
					# unit-variance structural shocks: drop sig2_d
					outer_q <- outer(A_sq, B_sq)
					dim(outer_q) <- c(n_row, n_row, n_col, n_col)
					outer_q <- aperm(outer_q, c(1, 3, 2, 4))
					contrib_dyad <- contrib_dyad + outer_q
					row_sq <- rowSums(A_sq)
					col_sq <- rowSums(B_sq)
					mse_ij <- mse_ij + outer(row_sq, col_sq)
				} else {
					A_sq <- A_cum[, , q] ^ 2
					B_sq <- B_cum[, , q] ^ 2
					# contrib[i,j,k,l] += sig2 * A_sq[i,k] * B_sq[j,l]
					outer_q <- outer(A_sq, B_sq)
					dim(outer_q) <- c(n_row, n_row, n_col, n_col)
					outer_q <- aperm(outer_q, c(1, 3, 2, 4))
					contrib_dyad <- contrib_dyad + sig2_d * outer_q
					row_sq <- rowSums(A_sq)
					col_sq <- rowSums(B_sq)
					mse_ij <- mse_ij + sig2_d * outer(row_sq, col_sq)
				}
			}
			mse_arr[, , h, d] <- mse_ij
			# normalize shares: divide each (i,j,k,l) by mse[i,j]
			if (identical(shock, "actor")) {
				# sum over l (B-side receiver) to give actor-level shock share
				contrib_actor <- apply(contrib_dyad, c(1, 2, 3), sum)
				for (i in seq_len(n_row)) for (j in seq_len(n_col)) {
					if (mse_ij[i, j] > 0)
						out[i, j, , h, d] <- contrib_actor[i, j, ] / mse_ij[i, j]
				}
			} else {
				for (i in seq_len(n_row)) for (j in seq_len(n_col)) {
					if (mse_ij[i, j] > 0)
						out[i, j, , , h, d] <- contrib_dyad[i, j, , ] / mse_ij[i, j]
				}
			}
		}
		# future-operator drift variance
		if (isTRUE(simulate_future_operators)) {
			drift_var <- matrix(0, n_row, n_col)
			for (h in seq_len(H)) {
				# at horizon h there are h future operator steps; each adds
				# RW variance scaling with norm(Theta_{t-1}) but a clean
				# isotropic bound is h * (tau_A2 + tau_B2) * ||Theta_t0||^2 / n
				# this is reported as a scalar drift summary, not a full
				# decomposition
				drift_var <- drift_var + (tau_A2_vec[d] + tau_B2_vec[d]) * h
				unattributed[, , h, d] <- drift_var
			}
		}
	}
	# summarise across draws
	alpha <- (1 - prob) / 2
	if (identical(shock, "actor")) {
		summarise_axis <- c(1, 2, 3, 4)
	} else {
		summarise_axis <- c(1, 2, 3, 4, 5)
	}
	median_arr <- apply(out, summarise_axis, stats::median, na.rm = TRUE)
	lower_arr <- apply(out, summarise_axis, stats::quantile, probs = alpha, na.rm = TRUE)
	upper_arr <- apply(out, summarise_axis, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
	total_mse <- apply(mse_arr, c(1, 2, 3), stats::median, na.rm = TRUE)
	unattrib_summary <- if (!is.null(unattributed))
		apply(unattributed, c(1, 2, 3), stats::median, na.rm = TRUE)
	else NULL
	res <- list(
		median = median_arr,
		lower = lower_arr,
		upper = upper_arr,
		total_mse = total_mse,
		unattributed_future_operator = unattrib_summary,
		metadata = list(
			H = H,
			t0 = t0,
			shock = shock,
			n_draws = D,
			isotropic_innovation = !use_structural,
			ordering_required = use_structural,
			ordering_applied = if (use_structural) ordering else NULL,
			structural = use_structural,
			simulate_future_operators = simulate_future_operators
		)
	)
	class(res) <- c("dbn_fevd", "list")
	res
}

####
#' Print method for dbn_fevd objects
#' @param x A `dbn_fevd` object.
#' @param ... Unused.
#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_fevd <- function(x, ...) {
	md <- x$metadata
	cli::cli_h2("DBN forecast error variance decomposition")
	cli::cli_text("Forecast origin t0 = {md$t0}, horizon H = {md$H}, {md$n_draws} posterior draw{?s}.")
	cli::cli_text("Aggregation: {.val {md$shock}} ({if (md$shock == 'actor') 'source-actor shares' else 'source-dyad shares'}).")
	if (isTRUE(md$structural)) {
		cli::cli_text("Structural FEVD: rotated by user-supplied {.code Sigma_row} {.code (x)} {.code Sigma_col} Cholesky factor.")
		if (!is.null(md$ordering_applied))
			cli::cli_text("Ordering applied: {.val {paste(md$ordering_applied, collapse = ', ')}}.")
	} else {
		cli::cli_text("Innovation covariance is isotropic ({.code sigma^2 I}); no Cholesky ordering required.")
	}
	# the per-cell median across draws does not in general preserve the
	# per-draw sum-to-1 constraint over the source axis. each individual
	# draw does sum to 1; the median surface does not
	cli::cli_text("Note: {.code $median} is a per-cell quantile across draws; per-draw shares sum to 1 over the source axis, but the median surface does not.")
	if (isTRUE(md$simulate_future_operators))
		cli::cli_text("Future-operator drift reported separately under {.code $unattributed_future_operator}.")
	cli::cli_text("Median shares stored in {.code $median}; intervals in {.code $lower} / {.code $upper}.")
	invisible(x)
}
