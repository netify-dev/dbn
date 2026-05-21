#' Actor-level forecast gain from a fitted dynamic bilinear model
#'
#' Computes per-actor predictive gain by comparing the dynamic model's
#' one-step-ahead predictions for each dyad against a static-baseline
#' prediction. For Gaussian-family fits the gain is MSE-based; for
#' ordinal / rank likelihood it uses log score.
#'
#' For each posterior draw m and each time t in window E:
#'   - predicted Y under the dynamic model is M + A_t Theta_{t-1} A_t'
#'   - predicted Y under the baseline is M alone
#' Per-dyad gain at (i, j, t):
#'   G_{ij,t}^{(m)} = (Y_{ij,t} - M_{ij}^{(m)})^2
#'                    - (Y_{ij,t} - M_{ij}^{(m)} - (A_t^{(m)} Theta_{t-1}^{(m)} A_t^{(m)'})_{ij})^2
#' Per-actor gain (Gaussian, MSE-form):
#'   G_i^{(m)}(E) = sum_{t in E} sum_{j != i} G_{ij,t}^{(m)}
#'
#' Returns one posterior summary per actor plus the per-draw matrix
#' so the caller can compute multi-chain R-hat / coverage.
#'
#' @param fit A `dbn_dynamic` fit (or compatible posterior list). Must
#'   contain `$Theta` (n x n x p x T x M), `$A` (list of length M of
#'   n x n x T arrays), `$M` (n x n x p x M), `$sigma2`, `$sigma2_obs`.
#' @param Y The observed data array (n x n x p x T). Required because
#'   the fit object may store only kept time slices.
#' @param window Integer vector of time points to sum over. Defaults to
#'   all t >= 2 (forecast gain is undefined at t = 1).
#' @param relation Index of the relation (p dimension) to use. Default 1.
#' @return A list with `per_actor` (n x M matrix of per-draw gains,
#'   summed over `window`), `summary` (a data frame with mean, sd,
#'   q025, q500, q975 per actor), and `attrs` (window, relation, n
#'   for reference).
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 10, seed = 1)
#' fit <- dbn_dynamic(sim$Y, family = "gaussian",
#'                    nscan = 200, burn = 100, time_thin = 1L,
#'                    verbose = FALSE)
#' fg <- forecast_gain(fit, sim$Y, window = 8:10, relation = 1)
#' head(fg$summary)
#' }
#' @export
forecast_gain <- function(fit, Y, window = NULL, relation = 1L) {
	n_row <- dim(Y)[1]
	n_col <- dim(Y)[2]
	if (n_row != n_col) {
		stop("forecast_gain currently supports unipartite networks ",
			"(n_row == n_col).")
	}
	Tt <- dim(Y)[4]
	if (is.null(window)) window <- seq(2L, Tt)
	if (any(window < 2 | window > Tt)) {
		stop("window must be a subset of 2:Tt")
	}

	A_keep <- simplify2array(fit$A)  # n x n x Tt x M
	if (length(dim(A_keep)) != 4) {
		stop("fit$A must be coercible to a 4D array (n x n x Tt x M)")
	}
	if (dim(A_keep)[3] != Tt) {
		stop("fit$A has Tt = ", dim(A_keep)[3],
			" but data has Tt = ", Tt,
			"; pass forecast_gain time_thin = 1 fit")
	}
	M_keep <- fit$M  # n x n x p x M
	Theta_keep <- fit$Theta  # n x n x p x Tt x M
	n_draws <- dim(A_keep)[4]
	if (dim(M_keep)[4] != n_draws) {
		stop("M and A have different draw counts")
	}

	Y_rel <- Y[, , relation, ]  # n x n x Tt
	per_actor <- matrix(0, n_row, n_draws)

	for (m in seq_len(n_draws)) {
		M_m <- M_keep[, , relation, m]
		for (t in window) {
			A_t  <- A_keep[, , t, m]
			Th_p <- Theta_keep[, , relation, t - 1, m]
			pred_dyn <- M_m + A_t %*% Th_p %*% t(A_t)
			# enforce diag(Y) = 0 convention: zero the diagonal
			diag(pred_dyn) <- 0
			diag(M_m) <- 0
			resid_base <- Y_rel[, , t] - M_m
			resid_dyn  <- Y_rel[, , t] - pred_dyn
			diag(resid_base) <- 0
			diag(resid_dyn)  <- 0
			gain_mat <- resid_base^2 - resid_dyn^2  # n x n
			# zero diagonal so it doesn't pollute the row sum
			diag(gain_mat) <- 0
			per_actor[, m] <- per_actor[, m] + rowSums(gain_mat)
		}
	}

	summary_df <- data.frame(
		actor   = seq_len(n_row),
		mean    = rowMeans(per_actor),
		sd      = apply(per_actor, 1, sd),
		q025    = apply(per_actor, 1, quantile, probs = 0.025),
		q500    = apply(per_actor, 1, quantile, probs = 0.5),
		q975    = apply(per_actor, 1, quantile, probs = 0.975)
	)

	list(
		per_actor = per_actor,
		summary   = summary_df,
		attrs     = list(window = window, relation = relation, n = n_row)
	)
}

#' Truth forecast gain under a known data-generating process
#'
#' Companion to `forecast_gain()` for simulation studies. Computes
#' the per-actor gain using the true A_arr, Theta_arr, M_mat at each
#' time. Useful for calibration checks: posterior `forecast_gain()`
#' should be unbiased and well-calibrated for this target.
#'
#' @param Y The observed data array (n x n x p x T).
#' @param A_true (n x n x T) truth operator path.
#' @param Theta_true (n x n x T) truth latent state path.
#' @param M_true (n x n) truth baseline.
#' @param window Integer vector of time points to sum over.
#' @param relation Index of the relation (default 1).
#' @return Numeric vector of length n giving per-actor truth gain.
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 8, seed = 1)
#' n <- dim(sim$A)[1]; Tt <- dim(sim$A)[3]
#' M_true <- matrix(0, n, n)  # simulator centers Theta at zero
#' Theta_true_3d <- array(sim$Theta[, , 1, ], c(n, n, Tt))
#' truth_gain <- dbn:::forecast_gain_truth(sim$Y, sim$A,
#'   Theta_true_3d, M_true,
#'   window = 5:8)
#' truth_gain
#' }
#' @keywords internal
forecast_gain_truth <- function(Y, A_true, Theta_true, M_true,
	window = NULL, relation = 1L) {
	n <- dim(Y)[1]
	Tt <- dim(Y)[4]
	if (is.null(window)) window <- seq(2L, Tt)
	Y_rel <- Y[, , relation, ]
	per_actor <- numeric(n)
	M_mat <- M_true; diag(M_mat) <- 0
	for (t in window) {
		A_t  <- A_true[, , t]
		Th_p <- Theta_true[, , t - 1]
		pred_dyn <- M_mat + A_t %*% Th_p %*% t(A_t)
		diag(pred_dyn) <- 0
		resid_base <- Y_rel[, , t] - M_mat
		resid_dyn  <- Y_rel[, , t] - pred_dyn
		diag(resid_base) <- 0
		diag(resid_dyn)  <- 0
		gain_mat <- resid_base^2 - resid_dyn^2
		diag(gain_mat) <- 0
		per_actor <- per_actor + rowSums(gain_mat)
	}
	per_actor
}
