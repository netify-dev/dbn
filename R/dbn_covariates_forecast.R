####
# covariate-aware H-step forecast
#
# the dynamic-covariate model:
#   theta_t = L_t + A_t (theta_{t-1} - L_{t-1}) B_t' + eps
# for forecasts we need L_t at horizons T+1, ..., T+H. the user supplies
# these via a new_covariates argument: a dbn_covariates object describing
# the design at the same K terms but at length-H time stamps.
####

####
#' Covariate-aware forecast
#'
#' @description H-step-ahead forecasts that thread the exogenous covariate
#'   path through the bilinear lag operand. Each posterior draw of the
#'   sampler is propagated forward H steps using random-walk innovations
#'   on A_t / B_t and Gaussian state noise on the centered residual.
#'
#' @param fit A \code{dbn_covariates_fit} produced by
#'   \code{dbn(..., covariates = ...)}.
#' @param H Forecast horizon (positive integer).
#' @param new_covariates A \code{dbn_covariates} object whose terms match
#'   the fit's terms; the design at each of the new H time stamps is
#'   used to build \code{L_{T+1}, ..., L_{T+H}}. Time-invariant terms
#'   may be reused from the fit's \code{$covariates} slot.
#' @param draws Number of posterior draws to use (default 100, sampled
#'   with replacement from the stored chain).
#' @param summary One of \code{"none"} (default; return full draws array),
#'   \code{"mean"}, \code{"ci"} (return list with mean / lower / upper).
#' @param level Coverage for the \code{summary = "ci"} bands (default 0.95).
#' @param seed Optional RNG seed.
#' @return Either a 5D array \code{[n_row, n_col, p, H, draws]}, a 4D
#'   posterior-mean array, or a list of three 4D arrays (\code{summary =
#'   "ci"}).
#' @author Tosin Salau and Shahryar Minhas
#' @export
forecast_with_covariates <- function(fit, H, new_covariates,
                                      draws = 100L,
                                      summary = c("none", "mean", "ci"),
                                      level = 0.95,
                                      seed = NULL) {
	if (!inherits(fit, "dbn_covariates_fit"))
		cli::cli_abort(c(
			"{.arg fit} must be a {.cls dbn_covariates_fit}.",
			"i" = "Fit with {.code dbn(..., covariates = ...)}."
		))
	if (!inherits(new_covariates, "dbn_covariates"))
		cli::cli_abort("{.arg new_covariates} must be a {.cls dbn_covariates} object.")
	if (length(H) != 1L || !is.numeric(H) || !is.finite(H) || H < 1 || H != round(H))
		cli::cli_abort("{.arg H} must be a single positive integer.")
	summary <- match.arg(summary)
	if (!is.null(seed)) set.seed(seed)

	n_row <- fit$dims$n_row
	n_col <- fit$dims$n_col
	p     <- fit$dims$p
	Tt    <- fit$dims$Tt
	H     <- as.integer(H)

	covars_new <- .dbn_validate_covariates(new_covariates, n_row, n_col, p, H,
		symmetric = isTRUE(fit$symmetric))
	# the term schema must match the fit's covariate object
	if (!identical(sort(covars_new$terms$name), sort(fit$covariates$terms$name)))
		cli::cli_abort(c(
			"{.arg new_covariates} terms must match the fit's terms.",
			"i" = "Fit terms: {.val {fit$covariates$terms$name}}",
			"i" = "New terms: {.val {covars_new$terms$name}}"
		))
	# permute new term order to match fit's beta column order
	term_order <- match(fit$covariates$terms$name, covars_new$terms$name)
	# build a permuted index for the design materialiser by reordering entries
	covars_new$entries <- covars_new$entries[term_order]
	covars_new$terms   <- covars_new$terms[term_order, , drop = FALSE]

	# pull L at last observed time from the fit (needed as L_{T} for the
	# first step's lag operand). it lives in fit$L if stored
	have_L <- !is.null(fit$L) && is.array(fit$L) && length(dim(fit$L)) == 4L
	n_saved <- length(fit$A)
	idx <- sample(n_saved, draws, replace = TRUE)
	have_theta <- !is.null(fit$Theta) && is.array(fit$Theta) &&
		length(dim(fit$Theta)) == 5L
	if (!have_theta)
		cli::cli_abort(c(
			"Fit does not store {.code fit$Theta}; cannot forecast.",
			"i" = "Refit with default settings ({.code store_theta = TRUE})."
		))
	t_last_theta <- dim(fit$Theta)[4]

	tau_A2_avail <- !all(is.na(fit$tau_A2))
	tau_B2_avail <- !all(is.na(fit$tau_B2))
	propagate_operator <- tau_A2_avail && tau_B2_avail

	Theta_pred <- array(0, c(n_row, n_col, p, H, draws))

	for (s in seq_len(draws)) {
		dr <- idx[s]
		beta_dr <- as.numeric(fit$beta[dr, ])
		# safe sigma2 index: bootstrap-expanded fits have length-1 sigma2
		sigma2  <- if (length(fit$sigma2) >= dr) fit$sigma2[dr] else fit$sigma2[1]
		if (!is.finite(sigma2)) sigma2 <- 1
		# operator state at T
		t_last_A <- dim(fit$A[[dr]])[3]
		t_last_B <- dim(fit$B[[dr]])[3]
		A_curr <- fit$A[[dr]][, , t_last_A]
		B_curr <- fit$B[[dr]][, , t_last_B]
		if (propagate_operator) {
			tauA <- sqrt(max(fit$tau_A2[dr], 0))
			tauB <- sqrt(max(fit$tau_B2[dr], 0))
		}
		# theta at T (last observed time)
		Theta_T <- fit$Theta[, , , t_last_theta, dr, drop = FALSE]
		dim(Theta_T) <- c(n_row, n_col, p)
		# L at T (for the first step's lag operand)
		L_T <- if (have_L) fit$L[, , t_last_theta, dr] else
			.dbn_build_linear_predictor_t(fit$covariates, beta_dr,
				dim(fit$L)[3])
		L_prev <- L_T
		Theta_prev <- Theta_T

		for (h in seq_len(H)) {
			# build L_{T+h} from new design using this draw's beta
			L_h <- .dbn_build_linear_predictor_t(covars_new, beta_dr, h)
			# propagate operators via RW innovation
			if (propagate_operator) {
				if (tauA > 0)
					A_curr <- A_curr + matrix(stats::rnorm(n_row * n_row, 0, tauA),
						n_row, n_row)
				if (tauB > 0)
					B_curr <- B_curr + matrix(stats::rnorm(n_col * n_col, 0, tauB),
						n_col, n_col)
			}
			Theta_new <- array(0, c(n_row, n_col, p))
			for (rel in seq_len(p)) {
				dev_prev <- Theta_prev[, , rel] - L_prev
				Theta_new[, , rel] <- L_h + A_curr %*% dev_prev %*% t(B_curr) +
					sqrt(sigma2) *
					matrix(stats::rnorm(n_row * n_col), n_row, n_col)
			}
			Theta_pred[, , , h, s] <- Theta_new
			Theta_prev <- Theta_new
			L_prev <- L_h
		}
	}

	# map latent forecast to observation scale per family
	fam_name <- if (is.list(fit$family)) fit$family$name else fit$family
	if (identical(fam_name, "binary")) {
		Theta_pred[] <- stats::pnorm(Theta_pred)
	}

	if (identical(summary, "mean")) {
		return(apply(Theta_pred, 1:4, mean))
	}
	if (identical(summary, "ci")) {
		alpha <- (1 - level) / 2
		mean_arr  <- apply(Theta_pred, 1:4, mean)
		lower_arr <- apply(Theta_pred, 1:4, stats::quantile, probs = alpha)
		upper_arr <- apply(Theta_pred, 1:4, stats::quantile, probs = 1 - alpha)
		return(list(mean = mean_arr, lower = lower_arr, upper = upper_arr,
			level = level))
	}
	Theta_pred
}
