#' Actor-level forecast gain from a fitted dynamic bilinear model
#'
#' Computes per-actor predictive gain by comparing the dynamic model's
#' one-step-ahead predictions for each dyad against a static-baseline
#' prediction. For Gaussian-family fits the gain is MSE-based; for
#' ordinal / rank likelihood it uses log score.
#'
#' For each posterior draw m and each time t in window E:
#'   - predicted Theta under the dynamic model is
#'     `M + A_t (Theta_{t-1} - M) B_t'` (the one-step state prediction
#'     under the bilinear AR; the FFBS sampler centers the state on M);
#'   - predicted Theta under the baseline is M alone.
#' Per-dyad gain at (i, j, t):
#'   G_{ij,t}^{(m)} = (Y_{ij,t} - M_{ij}^{(m)})^2
#'                    - (Y_{ij,t} - pred_dyn_{ij,t}^{(m)})^2
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
#' @param ... Reserved. Passing `H = ...` here triggers a directive
#'   pointer to `predict(fit, H = ...)` instead, since `forecast_gain()`
#'   is an in-sample predictive-gain measure, not a multi-step forecast.
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
#' @author Tosin Salau and Shahryar Minhas
#' @export
forecast_gain <- function(fit, Y = NULL, window = NULL, relation = 1L, ...) {
	# catch the common confusion: `forecast_gain` is IN-SAMPLE predictive
	# gain over a chosen window, not multi-step forecasting. Users who reach
	# for an `H = ...` argument want `predict(fit, H = ...)`.
	extra <- list(...)
	if ("H" %in% names(extra) || "horizon" %in% names(extra)) {
		cli::cli_abort(c(
			"{.fun forecast_gain} does not take an {.arg H} (horizon) argument.",
			"i" = "It compares one-step in-sample predictions against a static baseline over {.arg window}.",
			"i" = "For multi-step forecasts use {.code predict(fit, H = N)} instead."
		))
	}
	if (length(extra) > 0L) {
		cli::cli_warn(c(
			"Ignoring {length(extra)} unrecognized argument{?s} to {.fun forecast_gain}: {.arg {names(extra)}}.",
			"i" = "Allowed arguments: {.arg fit}, {.arg Y}, {.arg window}, {.arg relation}."
		))
	}
	# refuse on non-gaussian families: prediction is on the latent gaussian
	# scale while Y is on the observed scale ({0,1} or ordinal category),
	# so squared-difference loss against raw Y is meaningless.
	fam_name <- if (is.list(fit$family)) fit$family$name else fit$family
	if (!is.null(fam_name) && !identical(fam_name, "gaussian")) {
		cli::cli_abort(c(
			"{.fun forecast_gain} currently supports only {.code family = \"gaussian\"} fits.",
			"x" = "This fit has {.code family = \"{fam_name}\"}.",
			"i" = "The squared-error gain compares predictions on the latent scale to observations on the {fam_name} scale, which is not a meaningful target.",
			"i" = "For binary / ordinal in-sample predictive checks, use {.fun posterior_predict_dbn} and {.fun plot_ppc_density} / {.fun plot_ppc_ecdf}."
		))
	}
	# default Y to the data stored on the fit
	if (is.null(Y)) {
		if (is.null(fit$Y)) {
			cli::cli_abort(c(
				"{.arg Y} is required and {.code fit$Y} was not stored on this fit.",
				"i" = "Pass the original data array explicitly."
			))
		}
		Y <- fit$Y
	}
	n_row <- dim(Y)[1]
	n_col <- dim(Y)[2]
	if (n_row != n_col) {
		cli::cli_abort(c(
			"{.fun forecast_gain} currently supports unipartite (square) networks.",
			"x" = "Got dimensions {.val {n_row}} x {.val {n_col}}.",
			"i" = "Compute per-actor in/out gain manually via {.fun predict} for bipartite cases."
		))
	}
	Tt <- dim(Y)[4]
	if (is.null(window)) window <- seq(2L, Tt)
	if (any(window < 2 | window > Tt)) {
		stop("window must be a subset of 2:Tt")
	}

	A_keep <- simplify2array(fit$A)  # n x n x Tt x M
	B_keep <- simplify2array(fit$B)  # n x n x Tt x M
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
			B_t  <- B_keep[, , t, m]
			Th_p <- Theta_keep[, , relation, t - 1, m]
			# bilinear-AR one-step prediction:
			# E[Theta_t | Theta_{t-1}, M, A, B] = M + A_t (Theta_{t-1} - M) B_t'
			pred_dyn <- M_m + A_t %*% (Th_p - M_m) %*% t(B_t)
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

	# na.rm = TRUE on summaries: NA observations in Y propagate to
	# per_actor rows for actors absent in any window time-point.
	summary_df <- data.frame(
		actor   = seq_len(n_row),
		mean    = rowMeans(per_actor, na.rm = TRUE),
		sd      = apply(per_actor, 1, sd, na.rm = TRUE),
		q025    = apply(per_actor, 1, quantile, probs = 0.025, na.rm = TRUE),
		q500    = apply(per_actor, 1, quantile, probs = 0.5, na.rm = TRUE),
		q975    = apply(per_actor, 1, quantile, probs = 0.975, na.rm = TRUE)
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
#' truth_gain <- forecast_gain_truth(sim$Y, sim$A,
#'   Theta_true_3d, M_true,
#'   window = 5:8)
#' truth_gain
#' }
#' @export
#' @keywords internal
forecast_gain_truth <- function(Y, A_true, Theta_true, M_true,
	window = NULL, relation = 1L) {
	# validate required inputs upfront with cli messages; without these,
	# NULL inputs trip a cryptic base-R seq() error several lines down.
	for (nm in c("Y", "A_true", "Theta_true", "M_true")) {
		val <- get(nm)
		if (is.null(val)) {
			cli::cli_abort(c(
				"{.arg {nm}} must not be NULL.",
				"i" = "{.fn forecast_gain_truth} needs all four ground-truth arrays: {.var Y}, {.var A_true}, {.var Theta_true}, {.var M_true}."
			))
		}
	}
	if (is.null(dim(Y)) || length(dim(Y)) != 4L) {
		cli::cli_abort(c(
			"{.arg Y} must be a 4D array [n_row, n_col, p, time].",
			"x" = "Got {.cls {class(Y)[1]}} with dim {.val {if (is.null(dim(Y))) length(Y) else paste(dim(Y), collapse = 'x')}}."
		))
	}
	n <- dim(Y)[1]
	Tt <- dim(Y)[4]
	if (is.null(window)) window <- seq(2L, Tt)
	Y_rel <- Y[, , relation, ]
	per_actor <- numeric(n)
	M_mat <- M_true; diag(M_mat) <- 0
	for (t in window) {
		A_t  <- A_true[, , t]
		Th_p <- Theta_true[, , t - 1]
		# match the corrected forecast_gain() formula: deviation-form prediction
		pred_dyn <- M_mat + A_t %*% (Th_p - M_mat) %*% t(A_t)
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
