####
# s3 method dispatch for dbn objects
####

#' S3 Methods for DBN Objects
#'
#' @description These methods dispatch to model-specific implementations
#' @name dbn-methods
#' @details
#' The following S3 methods are available for dbn objects:
#' \itemize{
#'   \item \code{print}: Display a concise summary of the dbn model
#'   \item \code{plot}: Create model-specific diagnostic plots
#'   \item \code{summary}: Generate detailed model summaries
#'   \item \code{predict}: Generate predictions or simulations
#' }
NULL

#' Plot Diagnostics for a Fitted DBN Model
#'
#' @description Creates model-specific diagnostic plots for a fitted DBN
#'   object. The type of plot depends on the model variant:
#'   \itemize{
#'     \item **Static**: trace plots, posterior histograms, and a network
#'       summary of the estimated B matrix.
#'     \item **Dynamic**: trace plots for variance parameters and, if
#'       available, time-varying A/B summaries.
#'     \item **Lowrank**: trace plots, estimated factor trajectories
#'       (\eqn{\alpha_t}), and the posterior mean node-loading matrix U.
#'     \item **HMM**: regime probabilities over time, the estimated
#'       transition matrix, and MCMC trace plots.
#'     \item **Piecewise**: trace plots for each regime block.
#'   }
#'
#' @param x A fitted `dbn` object returned by [dbn()].
#' @param ... Additional arguments passed to model-specific plot functions
#'   (e.g., `alpha` for edge significance in the static model).
#' @return A ggplot2 object (or arranged multi-panel plot) is printed and
#'   returned invisibly.
#' @seealso [dbn()], [plot_trace()], [check_convergence()], [summary.dbn()]
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' plot(fit)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method plot dbn
plot.dbn <- function(x, ...) {
	if (!inherits(x, "dbn")) cli::cli_abort("{.arg x} must be a {.cls dbn} object.")
	if (isTRUE(x$meta$sampler_used %in% c("als", "als_tv", "als_piecewise"))) {
		cli::cli_abort(c(
			"{.fun plot} draws MCMC trace diagnostics, which need a posterior chain.",
			"x" = "This is a point-estimate fit from {.fun dbn_als} -- there is no chain to plot.",
			"i" = "Fit with {.fun dbn} for trace plots, or inspect {.code fit$A} / {.code fit$Theta} or {.fun compute_irf} output directly.",
			"i" = "Use {.code dyad_path(fit, i, j)}, {.code net_snapshot(fit, t)}, or {.code plot_theta(fit, time)} for time-resolved views with bootstrap intervals."
		))
	}
	plot_fun <- switch(x$model,
		static = "plot_static",
		dynamic = "plot_dynamic",
		lowrank = "plot_lowrank",
		hmm = "plot_hmm",
		piecewise = "plot_piecewise",
		cli::cli_abort("Unknown model: {.val {x$model}}.")
	)
	do.call(plot_fun, list(x, ...))
}

#' Summarize a Fitted DBN Model
#'
#' @description Prints a structured summary of a fitted DBN model including
#'   data dimensions, posterior means and 95\% credible intervals for scalar
#'   parameters, and model-specific details (e.g., transition matrix for HMM,
#'   block structure for piecewise, DIC for Gaussian models).
#'
#' @param object A fitted `dbn` object returned by [dbn()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `object`. The summary is printed to the console.
#' @seealso [dbn()], [param_summary()], [check_convergence()], [plot.dbn()]
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' summary(fit)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method summary dbn
summary.dbn <- function(object, ...) {
	if (!inherits(object, "dbn")) cli::cli_abort("{.arg object} must be a {.cls dbn} object.")
	summary_fun <- switch(object$model,
		static = "summary_static",
		dynamic = "summary_dynamic",
		lowrank = "summary_lowrank",
		hmm = "summary_hmm",
		piecewise = "summary_piecewise",
		cli::cli_abort("Unknown model: {.val {object$model}}.")
	)
	# the model-specific summaries print their cli block and return the fit;
	# wrap in invisible() so the returned <dbn> object does not auto-print
	# (which would fire print.dbn and emit a second, redundant block).
	invisible(do.call(summary_fun, list(object, ...)))
}

#' Predict from a Fitted DBN Model
#'
#' @description Generates predictions from a fitted DBN model. With no
#'   forecasting arguments, returns posterior predictive replications of the
#'   observed data for model checking. With `H` (and optionally `draws` and
#'   `summary`), propagates the estimated dynamics forward to produce
#'   multi-step-ahead forecasts.
#'
#' @details
#' **Posterior predictive distribution (default):** Simulates new datasets
#' from the fitted model to compare against the observed data. Use
#' [plot_ppc_ecdf()] or [plot_ppc_density()] for visual checks.
#'
#' **Forecasting:** Pass `H` (horizon), `draws` (number of posterior draws
#' to use), and optionally `summary = "mean"` to average across draws. The
#' forecast propagates \eqn{A_t}, \eqn{B_t}, and process noise forward from
#' the final observed time point.
#'
#' @param object A fitted `dbn` object returned by [dbn()].
#' @param ... Model-specific arguments:
#'   \describe{
#'     \item{H}{Forecast horizon (number of steps ahead). Triggers
#'       forecasting mode.}
#'     \item{draws}{Number of posterior draws to use (for both forecasting
#'       and posterior predictive checks). Default depends on the entry
#'       point used downstream.}
#'     \item{summary}{If `"mean"`, return the pointwise mean across draws.
#'       If `NULL` (default), return the full array of simulated
#'       trajectories.}
#'     \item{seed}{Random seed for reproducibility.}
#'   }
#' @return For forecasting (`H` specified): an array of dimensions
#'   `[n_row, n_col, p, H]` (if `summary = "mean"`) or
#'   `[n_row, n_col, p, H, draws]`.
#'   For posterior predictive checks: a list of class `"dbn_ppd"` containing
#'   replicated datasets.
#' @seealso [dbn()], [posterior_predict_dbn()], [plot_ppc_ecdf()],
#'   [plot_ppc_density()]
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 6886)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'            nscan = 200, burn = 100, verbose = FALSE)
#'
#' # forecasting: 5 steps ahead
#' fc <- predict(fit, H = 5, draws = 50, summary = "mean")
#' dim(fc)
#'
#' # posterior predictive check
#' ppd <- predict(fit, draws = 10)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method predict dbn
predict.dbn <- function(object, ...) {
	if (!inherits(object, "dbn")) cli::cli_abort("{.arg object} must be a {.cls dbn} object.")
	dots <- list(...)
	# catch a non-NULL `newdata = ...` directly: predict.dbn cannot accept
	# new panel data because the model conditions on the fitted latent path
	# Theta_{1..T}. A user passing newdata expects out-of-sample
	# prediction; silently falling back to a PPD on the training data
	# is misleading. an explicit newdata = NULL is a no-op and passes
	# through.
	if ("newdata" %in% names(dots) && !is.null(dots$newdata)) {
		cli::cli_abort(c(
			"{.fun predict.dbn} does not support {.arg newdata}.",
			"i" = "For H-step-ahead forecasts (post-sample): {.code predict(fit, H = h, draws = D)}.",
			"i" = "For predictions on a new panel: refit with {.code dbn()} on the concatenated array; the latent path is conditional on the observed Theta history.",
			"i" = "For posterior-predictive replicates of the training data: {.code predict(fit, draws = D)} (no {.arg newdata} needed)."
		))
	}
	dots$newdata <- NULL
	# validate `summary=` against the user-facing vocabulary (none / mean /
	# ci / draws); reject unknowns directively. downstream-normalise "draws"
	# to "none" so the simulator sees its native names.
	if ("summary" %in% names(dots)) {
		s_user <- dots$summary
		valid_summary <- c("none", "mean", "ci", "draws")
		if (!is.character(s_user) || length(s_user) != 1L || is.na(s_user) ||
		    !(s_user %in% valid_summary)) {
			cli::cli_abort(c(
				"{.arg summary} must be one of {.val {valid_summary}}.",
				"x" = "Got {.val {s_user}}.",
				"i" = "{.val none} or {.val draws}: return all per-draw forecasts.",
				"i" = "{.val mean}: return posterior-mean forecasts.",
				"i" = "{.val ci}: return a list of {.var mean}/{.var lower}/{.var upper} arrays."
			))
		}
		# user-facing alias: "draws" reads naturally but the downstream
		# simulator accepts only "none" / "mean" / "ci".
		if (identical(s_user, "draws")) dots$summary <- "none"
	}
	if (object$model == "hmm") return(predict_hmm(object, ...))
	if (object$model == "lowrank") return(predict_lowrank(object, ...))
	# alias common synonyms for the forecast horizon onto `H`:
	# `h` (lowercase typo), `horizon`, `n.ahead` (vars / forecast convention).
	for (alias in c("h", "horizon", "horizons", "n.ahead")) {
		if (alias %in% names(dots) && !("H" %in% names(dots))) {
			dots$H <- dots[[alias]]
			dots[[alias]] <- NULL
		} else if (alias %in% names(dots) && "H" %in% names(dots)) {
			cli::cli_warn(c(
				"Both {.arg H} and {.arg {alias}} were specified; using {.arg H}.",
				"i" = "Pass only one."
			))
			dots[[alias]] <- NULL
		}
	}
	# always validate `draws =` when supplied; both the forecast path and the
	# ppd path consume it downstream and crash on 0 / NULL / NA.
	if ("draws" %in% names(dots)) {
		d_val <- dots$draws
		if (length(d_val) != 1L || !is.numeric(d_val) || !is.finite(d_val) ||
		    d_val < 1 || d_val != round(d_val)) {
			cli::cli_abort(c(
				"{.arg draws} must be a single positive integer.",
				"x" = "Got {.val {d_val}}."
			))
		}
	}
	is_forecast <- "H" %in% names(dots)
	if (is_forecast) {
		H_val <- dots$H
		if (length(H_val) != 1L || !is.numeric(H_val) || !is.finite(H_val) ||
		    H_val < 1 || H_val != round(H_val)) {
			cli::cli_abort(c(
				"{.arg H} (forecast horizon) must be a single positive integer.",
				"x" = "Got {.val {H_val}}."
			))
		}
	}
	if (!is_forecast && length(dots) > 0) {
		known_ppd <- c("draws", "seed", "summary")
		unknown <- setdiff(names(dots), known_ppd)
		if (length(unknown) > 0) {
			cli::cli_warn(c(
				"Unrecognized argument{?s}: {.arg {unknown}}.",
				"i" = "For forecasting use {.arg H}, {.arg draws}, or {.arg summary}.",
				"i" = "Falling back to posterior predictive distribution."
			))
		}
	}
	if (is_forecast) {
		if (object$model == "static") {
			cli::cli_abort(c(
				"H-step-ahead forecasting is not available for the {.val static} model.",
				"i" = "The static operator is time-invariant; use {.val dynamic} or {.val piecewise} for forecasts.",
				"i" = "For posterior predictive replication of the observed data, call {.code predict()} without {.arg H}."
			))
		}
		# bootstrap-aware path: expand A/B/M draws so the existing forecast
		# code sees R draws and produces honest bootstrap intervals.
		su <- object$meta$sampler_used; is_als <- length(su) == 1L && !is.na(su) && su %in% c("als", "als_tv")
		if (is_als && !is.null(object$bootstrap) && !isTRUE(object$meta$bootstrap_expanded)) {
			object <- .bootstrap_expand(object)
			# ALS+bootstrap forecasts hold (A_T, B_T) constant across the
			# horizon: variance comes from the (A, B, M) bootstrap spread,
			# not from operator-evolution RW innovation. bands are
			# conservative at h = 1 and increasingly anti-conservative as
			# h grows. fire once per session (idempotent on repeat calls)
			# to avoid spamming production pipelines that forecast in loops.
			if (!isTRUE(getOption("dbn.als_bands_disclaimer_fired", FALSE))) {
				cli::cli_warn(c(
					"ALS+bootstrap forecast bands hold the operator constant across the horizon.",
					"i" = "Forecast variance reflects (A, B, M) bootstrap uncertainty but NOT the random-walk innovation that would propagate {.code A_t -> A_t+1}.",
					"i" = "Bands are conservative at h = 1 and increasingly anti-conservative as h grows.",
					"i" = "For honest forecast intervals with operator-evolution variance, refit with {.code dbn(model = \"dynamic\")} (MCMC).",
					"i" = "This message fires only once per session; clear with {.code options(dbn.als_bands_disclaimer_fired = FALSE)}."
				))
				options(dbn.als_bands_disclaimer_fired = TRUE)
			}
		} else if (is_als || isFALSE(object$meta$uncertainty_available)) {
			cli::cli_warn(c(
				"Forecasting from a point-estimate fit ({.fun dbn_als}).",
				"i" = "The forecast carries no posterior uncertainty; every draw is identical.",
				"i" = "Refit with {.code bootstrap = N} or with {.fun dbn} for forecast intervals."
			))
		}
		predict_fun <- switch(object$model,
			dynamic = "simulate_dynamic",
			piecewise = "simulate_piecewise",
			cli::cli_abort("Forecasting is not supported for model {.val {object$model}}.")
		)
		# the downstream simulators accept summary in {"none", "mean"}. expose a
		# higher-level `summary = "ci"` shortcut here: pull draws, compute
		# {lower, mean, upper} per (i, j, p, h), and return a list of three
		# 4D arrays plus the metadata so users don't have to compute quantiles
		# manually.
		summary_user <- dots$summary
		want_ci <- is.character(summary_user) && length(summary_user) == 1L &&
			summary_user == "ci"
		# accept `level = N` (forecast-package convention, scalar percent in
		# (1, 100)) as a way to control CI width; converts to probs.
		# `summary_probs = c(lo, hi)` is the lower-level alternative.
		# refuse both supplied at once + refuse vector `level` + refuse a
		# probability-scale level (< 1) so a user passing 0.95 expecting 95%
		# gets a directive abort instead of a silent 0.95% interval.
		ci_probs <- c(0.025, 0.975)
		level_supplied <- "level" %in% names(dots)
		probs_supplied <- "summary_probs" %in% names(dots)
		if (level_supplied && probs_supplied) {
			cli::cli_abort(c(
				"Both {.arg level} and {.arg summary_probs} were supplied.",
				"i" = "Pass only one: {.arg level} is the forecast-package convention (percent), {.arg summary_probs} is the explicit c(lo, hi) form."
			))
		}
		if (level_supplied) {
			lvl <- dots$level
			if (!is.numeric(lvl) || length(lvl) != 1L || !is.finite(lvl)) {
				cli::cli_abort(c(
					"{.arg level} must be a single finite number on the percent scale, i.e. in (1, 100).",
					"x" = "Got {.val {lvl}}.",
					"i" = "For multiple bands use {.code summary_probs} with the desired c(lo, hi) twice, or pick one band."
				))
			}
			if (lvl < 1) {
				cli::cli_abort(c(
					"{.arg level} must be on the {.emph percent} scale (1..100), not the probability scale (0..1).",
					"x" = "Got {.val {lvl}}, which would request a {sprintf('%.2g%%', lvl)} interval.",
					"i" = "Pass {.code level = {lvl * 100}} for a {sprintf('%.0f%%', lvl * 100)} interval, or use {.code summary_probs = c({(1 - lvl) / 2}, {1 - (1 - lvl) / 2})} directly."
				))
			}
			if (lvl >= 100) {
				cli::cli_abort(c(
					"{.arg level} must be strictly less than 100.",
					"x" = "Got {.val {lvl}}."
				))
			}
			lvl_pick <- lvl / 100
			ci_probs <- c((1 - lvl_pick) / 2, 1 - (1 - lvl_pick) / 2)
			dots$level <- NULL
			# `level=` implies `summary = "ci"` unless the user said otherwise.
			if (is.null(summary_user)) {
				dots$summary <- "ci"
				summary_user <- "ci"
				want_ci <- TRUE
			}
		}
		if (probs_supplied) {
			sp <- dots$summary_probs
			if (!is.numeric(sp) || length(sp) != 2L || any(!is.finite(sp)) ||
			    any(sp <= 0) || any(sp >= 1) || sp[1] >= sp[2]) {
				cli::cli_abort(c(
					"{.arg summary_probs} must be a numeric c(lo, hi) with 0 < lo < hi < 1.",
					"x" = "Got {.val {sp}}."
				))
			}
			ci_probs <- sp
			dots$summary_probs <- NULL
		}
		if (want_ci) dots$summary <- "none"  # force per-draw output upstream
		fc <- do.call(predict_fun, c(list(object), dots))
		if (want_ci) {
			# fc shape: [n_row, n_col, p, H, draws]
			fc_arr <- fc
			# warn on a degenerate CI: <30 draws gives uninformative quantiles
			# (and <2 produces lower == upper).
			n_draws_fc <- if (length(dim(fc_arr)) == 5L) dim(fc_arr)[5] else NA_integer_
			if (is.finite(n_draws_fc) && n_draws_fc < 30L) {
				cli::cli_warn(c(
					"{.code summary = \"ci\"} requested with only {.val {n_draws_fc}} draws.",
					"i" = "Quantile-based intervals from fewer than ~30 draws are unstable; raise {.arg draws} (200+ is common).",
					"i" = "With {.val 1} draw the {.var lower}/{.var upper} reduce to the point estimate."
				))
			}
			fc_mean <- apply(fc_arr, 1:4, mean)
			fc_lo   <- apply(fc_arr, 1:4, stats::quantile, probs = ci_probs[1], na.rm = TRUE)
			fc_hi   <- apply(fc_arr, 1:4, stats::quantile, probs = ci_probs[2], na.rm = TRUE)
			# stats::quantile / mean strip dimnames; reapply via the same
			# decorator that the per-draw / summary='mean' path uses so users
			# can read off (sender, receiver, relation, horizon).
			fc_mean <- .dbn_decorate_forecast(fc_mean, object, H = as.integer(dots$H))
			fc_lo   <- .dbn_decorate_forecast(fc_lo,   object, H = as.integer(dots$H))
			fc_hi   <- .dbn_decorate_forecast(fc_hi,   object, H = as.integer(dots$H))
			fc <- list(mean = fc_mean, lower = fc_lo, upper = fc_hi)
			class(fc) <- c("dbn_forecast_ci", "list")
			attr(fc, "H") <- as.integer(dots$H)
			attr(fc, "model") <- object$model
			attr(fc, "family") <- if (is.list(object$family)) object$family$name else object$family
			attr(fc, "summary_kind") <- "ci"
			attr(fc, "probs") <- ci_probs
			return(fc)
		}
		# decorate the forecast array with sane dimnames so users don't have
		# to keep horizon / actor labels separately. only sets dimnames if the
		# simulator didn't already provide them.
		fc <- .dbn_decorate_forecast(fc, object, H = as.integer(dots$H))
		attr(fc, "H") <- as.integer(dots$H)
		attr(fc, "model") <- object$model
		attr(fc, "family") <- if (is.list(object$family)) object$family$name else object$family
		attr(fc, "summary_kind") <- dots$summary %||% "draws"
		class(fc) <- c("dbn_forecast", class(fc))
		return(fc)
	}
	ppd_args <- dots[intersect(names(dots), c("draws", "seed"))]
	# get_family may legitimately fail for old fits that didn't record family
	# metadata; NULL fallback selects the simulator-based predict path below.
	# the fall-through is intentional, but if get_family fails for a *new*
	# fit that should have family info, the error is a real bug -- surface it.
	fam <- tryCatch(get_family(object), error = function(e) {
		if (!is.null(object$family)) {
			cli::cli_warn(c(
				"{.fun get_family} failed on a fit that has {.code fit$family} set.",
				"x" = "{conditionMessage(e)}",
				"i" = "Falling back to model-specific simulator for {.fun predict}."
			))
		}
		NULL
	})
	if (!is.null(fam)) {
		return(do.call(posterior_predict_dbn, c(list(object), ppd_args)))
	}
	predict_fun <- switch(object$model,
		static = "simulate_static",
		dynamic = "simulate_dynamic",
		piecewise = "simulate_piecewise",
		cli::cli_abort("Unknown model: {.val {object$model}}.")
	)
	do.call(predict_fun, c(list(object), ppd_args))
}

#' Forecast method for DBN objects
#'
#' @description Wrapper around [predict.dbn()] that defaults to forecasting
#'   semantics: a numeric horizon is required, and the return is the
#'   H-step-ahead forecast array, not a posterior-predictive replicate.
#'   Provided for users who reach for the `forecast` package's vocabulary.
#'
#' @param object A `dbn` object.
#' @param h Forecast horizon (integer >= 1). Aliased to `H` internally.
#' @param ... Passed to [predict.dbn()] (e.g., `draws`, `seed`, `summary`).
#' @return The same forecast array `predict(object, H = h, ...)` would
#'   return.
#' @seealso [predict.dbn()]
#' @examples
#' \donttest{
#' if (requireNamespace("forecast", quietly = TRUE)) {
#'   sim <- simulate_dynamic_dbn(n = 6, time = 8, seed = 1)
#'   fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'              nscan = 200, burn = 100, verbose = FALSE)
#'   fc <- forecast::forecast(fit, h = 2, draws = 50, summary = "mean")
#'   dim(fc)
#' }
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
forecast.dbn <- function(object, h = 1L, ...) {
	dots <- list(...)
	# accept H= as an alias for h= (predict.dbn's canonical name). Doing
	# this here avoids the base-R partial-match collision where `H = 2`
	# in the caller's `...` would otherwise match both `h` (via prefix)
	# and the inner `H = as.integer(h)` keyword.
	if ("H" %in% names(dots)) {
		if (!missing(h)) {
			cli::cli_warn(c(
				"Both {.arg h} and {.arg H} were specified; using {.arg H}.",
				"i" = "Pass only one."
			))
		}
		h <- dots$H
		dots$H <- NULL
	}
	if (length(h) != 1L || !is.numeric(h) || !is.finite(h) || h < 1 || h != round(h)) {
		cli::cli_abort(c(
			"{.arg h} (forecast horizon) must be a single positive integer.",
			"x" = "Got {.val {h}}."
		))
	}
	do.call(predict.dbn, c(list(object = object, H = as.integer(h)), dots))
}

#' Plot method for DBN objects (internal wrapper)
#'
#' @description Internal wrapper that calls [plot.dbn()]. Use `plot(fit)`
#'   directly; this bare wrapper is retained for internal package use.
#' @param x A `dbn` object
#' @param ... Additional arguments passed to the S3 method
#' @return Invisibly returns the plot object
#' @keywords internal
#' @noRd
plot_dbn <- function(x, ...) plot.dbn(x, ...)

#' Summary method for DBN objects (internal wrapper)
#'
#' @description Internal wrapper that calls [summary.dbn()]. Use
#'   `summary(fit)` directly; this bare wrapper is retained for internal
#'   package use.
#' @param object A `dbn` object
#' @param ... Additional arguments passed to the S3 method
#' @return A summary object specific to the model type
#' @keywords internal
#' @noRd
summary_dbn <- function(object, ...) summary.dbn(object, ...)

####
# print method for dbn objects
####

#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method print dbn
print.dbn <- function(x, ...) {
	if (!inherits(x, "dbn")) cli::cli_abort("{.arg x} must be a {.cls dbn} object.")

	# header
	cat("Dynamic Bilinear Network (DBN) Model\n")
	cat(rep("-", 40), "\n", sep = "")

	# model type: prefer the sampler-specific label when the ALS pipeline ran.
	# the label names what actually ran (which operator parameterisation),
	# not the dispatch routing target. Avoids the beginner-confusing "STATIC"
	# vs `?dbn_als` "fast point-estimate of the dynamic model" mismatch by
	# explicitly naming the operator parameterisation.
	model_label <- toupper(x$model)
	su <- x$meta$sampler_used %||% ""
	if (identical(su, "als")) model_label <- "ALS point estimate (time-invariant operator A, B)"
	else if (identical(su, "als_tv")) model_label <- "ALS point estimate (time-varying operators A_t, B_t)"
	else if (identical(su, "als_piecewise")) model_label <- "ALS point estimate (block-constant operators A_k, B_k)"
	cat("Model Type: ", model_label, "\n", sep = "")

	# distribution family
	if (!is.null(x$family)) {
		fam_name <- if (is.list(x$family) && !is.null(x$family$name)) x$family$name
					else if (is.character(x$family)) x$family
					else "unknown"
		cat("Family: ", fam_name, "\n", sep = "")
	}

	# dimensions
	if (!is.null(x$dims)) {
		cat("\nData Dimensions:\n")
		if (isTRUE(x$dims$is_bipartite)) {
			cat("  Senders: ", x$dims$n_row, "\n", sep = "")
			cat("  Receivers: ", x$dims$n_col, "\n", sep = "")
		} else {
			cat("  Nodes: ", x$dims$n_row, "\n", sep = "")
		}
		cat("  Relations: ", x$dims$p, "\n", sep = "")
		if (!is.null(x$dims$Tt)) {
			cat("  Time points: ", x$dims$Tt, "\n", sep = "")
		}
	}

	# fitting settings: an ALS exploratory fit has no MCMC, so report it as a
	# point estimate rather than printing "Iterations: 0"
	su <- x$meta$sampler_used; is_als <- length(su) == 1L && !is.na(su) && su %in% c("als", "als_tv")
	if (is_als) {
		cat("\nEstimator: ALS exploratory fit (point estimate, no MCMC)\n")
		if (!is.null(x$settings$max_iter)) {
			cat("  Max ALS iterations: ", x$settings$max_iter, "\n", sep = "")
		}
		if (!is.null(x$meta$converged)) {
			cat("  Converged: ", x$meta$converged, "\n", sep = "")
		}
	} else if (!is.null(x$settings)) {
		cat("\nMCMC Settings:\n")
		cat("  Iterations: ", x$settings$nscan, "\n", sep = "")
		cat("  Burn-in: ", x$settings$burn, "\n", sep = "")
		cat("  Thinning: ", x$settings$odens, "\n", sep = "")
		# settings$draws is only populated for static/dynamic; HMM and lowrank
		# record the kept-draw count under meta$draws instead
		n_saved <- x$settings$draws %||% x$meta$draws
		cat("  Saved draws: ", n_saved %||% "(not recorded)", "\n", sep = "")
	}

	# model-specific details
	if (x$model == "hmm" && !is.null(x$R)) {
		cat("\nRegimes: ", x$R, "\n", sep = "")
	} else if (x$model == "lowrank" && !is.null(x$rank)) {
		cat("\nRank: ", x$rank, "\n", sep = "")
	} else if (x$model == "piecewise" && !is.null(x$blocks)) {
		cat("\nBlocks: ", x$blocks$K, "\n", sep = "")
		cat("Boundaries: ", paste(x$blocks$boundaries, collapse = ", "), "\n", sep = "")
	}

	# stored components
	cat("\nModel Components:\n")
	component_names <- names(x)
	exclude <- c("model", "family", "dims", "settings", "meta", "diagnostics")
	components <- setdiff(component_names, exclude)

	for (comp in components) {
		if (!is.null(x[[comp]])) {
			cat("  ", comp, ": ", sep = "")
			if (is.list(x[[comp]])) {
				cat("list of length ", length(x[[comp]]), "\n", sep = "")
			} else if (is.array(x[[comp]]) || is.matrix(x[[comp]])) {
				dims <- dim(x[[comp]])
				cat("[", paste(dims, collapse = " x "), "]\n", sep = "")
			} else if (is.vector(x[[comp]])) {
				cat("vector of length ", length(x[[comp]]), "\n", sep = "")
			} else {
				cat(class(x[[comp]])[1], "\n", sep = "")
			}
		}
	}

	# diagnostics
	if (!is.null(x$diagnostics) && !is.null(x$diagnostics$dic)) {
		cat("\nModel Diagnostics:\n")
		cat("  DIC: ", round(x$diagnostics$dic, 2), "\n", sep = "")
	}

	# print a one-line convergence headline using a light ESS heuristic
	# on the kept sigma2 / tau_A2 / tau_B2 traces, and point at
	# check_convergence() for full diagnostics
	if (!is_als && !identical(x$model, "static")) {
		# surface helper failures via cli_warn rather than NULL-swallowing
		tryCatch({
			ess_min <- .dbn_print_min_ess(x)
			if (is.finite(ess_min)) {
				lvl <- if (ess_min < 50) "low" else if (ess_min < 200) "modest" else "ok"
				lab <- switch(lvl,
					low    = "LOW (<50)",
					modest = "modest (50-200)",
					ok     = "ok (>=200)"
				)
				cat("\nConvergence:\n")
				cat("  min ESS across (sigma2, tau_A2, tau_B2): ", round(ess_min), " -- ", lab, "\n", sep = "")
				if (lvl != "ok") {
					cat("  Run check_convergence(fit) for full diagnostics.\n")
				}
			}
		}, error = function(e) {
			cli::cli_warn(c(
				"Print-method convergence headline failed on this fit.",
				"x" = "{conditionMessage(e)}",
				"i" = "The fit itself is fine; call {.fn check_convergence} directly."
			))
		})
	}

	invisible(x)
}

#' @keywords internal
.dbn_print_min_ess <- function(x) {
	pars <- c("sigma2", "tau_A2", "tau_B2")
	vals <- numeric(0)
	for (p in pars) {
		v <- x[[p]]
		if (is.numeric(v) && length(v) > 5L) {
			# crude ESS = N * (1 - rho_1) / (1 + rho_1), assuming geometric
			# autocorrelation decay (Geweke-style)
			v <- v[is.finite(v)]
			if (length(v) < 10L || stats::sd(v) == 0) next
			ac1 <- stats::acf(v, lag.max = 1L, plot = FALSE)$acf[2L, , 1L]
			if (!is.finite(ac1)) next
			ess <- length(v) * (1 - ac1) / (1 + ac1)
			if (is.finite(ess) && ess > 0) vals <- c(vals, ess)
		}
	}
	if (length(vals) == 0L) return(NA_real_)
	min(vals)
}
####
#' S3 Generic for Sampler Description
#'
#' @description Describe a fitted sampler (ALS, MCMC, etc.) and its uncertainty properties.
#'
#' @param object A fitted model object (dbn, dbn_boot, etc.)
#' @param verbose Logical. If TRUE, print full description.
#' @param ... Additional arguments for methods.
#'
#' @return Invisibly, information about the sampler.
#'
#' @keywords internal
#' @export
sampler_describe <- function(object, verbose = TRUE) {
	UseMethod("sampler_describe", object)
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
sampler_describe.dbn <- function(object, verbose = TRUE) {
	sampler <- object$meta$sampler_used %||% NA_character_
	# infer the sampler tier from the model name when meta$sampler_used is
	# missing (older fits sometimes don't have it stamped)
	if (is.na(sampler) || identical(sampler, "unknown")) {
		mdl <- object$model %||% ""
		sampler <- switch(mdl,
			static    = "mcmc_static",
			dynamic   = "mcmc_dynamic",
			piecewise = "mcmc_piecewise",
			hmm       = "mcmc_hmm",
			lowrank   = "mcmc_lowrank",
			"unknown"
		)
	}
	if (isTRUE(verbose)) {
		cli::cli_h3("Sampler Information")
		cli::cli_text("Type: {.val {sampler}}")
		# TV-ALS specifics
		if (identical(sampler, "als_tv") && !is.null(object$meta$tv)) {
			tv <- object$meta$tv
			cli::cli_text("RW-smoothed time-varying ALS")
			cli::cli_text("  lambda used: {.val {sprintf('%.4g', tv$lambda_used %||% NA_real_)}}")
			cli::cli_text("  mu (level ridge): {.val {sprintf('%.4g', tv$mu %||% NA_real_)}}")
			if (!is.null(tv$cv)) {
				cli::cli_text("  CV-selected lambda (1-SE rule); see meta$tv$cv for path")
			}
			cli::cli_text("  Converged: {.val {tv$tv_converged}} in {tv$tv_n_iter} iter")
		}
		unc <- object$meta$uncertainty_available
		if (length(unc) == 1L && !is.na(unc)) {
			cli::cli_text("Uncertainty available: {.val {unc}}")
		}
		if (!is.null(object$meta$approximation_note)) {
			cli::cli_text("Note: {object$meta$approximation_note}")
		}
		als <- object$meta$als_solve_summary
		if (!is.null(als) && als$n_solves > 0) {
			cli::cli_text(
				"ALS linear solves: {als$n_qr} QR, {als$n_ridge} ridge ",
				"({als$n_rank_deficient} rank-deficient)")
			wk <- als$worst_kappa
			kappa_txt <- if (is.finite(wk)) sprintf("%.2e", wk) else "Inf (rank-deficient)"
			cli::cli_text("Worst design condition number: {kappa_txt}")
			if (als$n_ridge > 0) {
				cli::cli_text("Median ridge penalty: {sprintf('%.2e', als$median_lambda)}")
				cli::cli_text(
					"{.strong Note:} ridge regularization was applied; ",
					"A/B uncertainty may be sensitive to initialization.")
			}
		}
	}
	invisible(list(
		sampler = sampler,
		uncertainty_available = object$meta$uncertainty_available,
		approximation_note = object$meta$approximation_note
	))
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
sampler_describe.dbn_boot <- function(object, verbose = TRUE) {
	cli::cli_h3("Bootstrap Sampler")
	cli::cli_text("Type: {.field {object$type}}")
	cli::cli_text("Valid replicates: {.val {object$n_valid}}/{.val {object$n_total}}")
	cli::cli_text("Standard errors available for A, B, M")
	invisible(object)
}

####
# the als exploration and bootstrap routines below adapt the iterative
# block-coordinate-descent estimator of the social influence regression
# model (minhas & hoff 2025, political analysis) to dynamic bilinear
# networks.
####
#' Fast Exploratory ALS Fitting for Dynamic Bilinear Networks
#'
#' @description
#' Fits a dynamic bilinear network model via alternating least squares (ALS),
#' providing a fast point estimate without MCMC. Useful for exploring data,
#' diagnosing model fit, and choosing initial values for full Bayesian inference.
#' No credible intervals are produced.
#'
#' The ALS estimator adapts the iterative block coordinate descent approach of
#' the Social Influence Regression (SIR) model of Minhas and Hoff (2025) to the
#' dynamic bilinear network setting.
#'
#' @param data Numeric array or file path. Same format as [dbn()]:
#'   3-dimensional `[actors, actors, time]` (single relation) or
#'   4-dimensional `[actors, actors, relations, time]` (multiple relations).
#'   Diagonal entries should be `NA`.
#' @param family Character string: `"gaussian"` (default), `"ordinal"`, or `"binary"`.
#'   For ordinal/binary, ALS is applied to latent continuous scores from [shared_preprocess()].
#' @param symmetric Logical. If TRUE, enforce B = A (symmetric/undirected networks).
#' @param max_iter Maximum ALS iterations (default: 100).
#' @param tol Convergence tolerance: relative change in the ALS objective
#'   (default: 1e-4).
#' @param ridge Asymmetric ridge penalty, as a fraction of the mean
#'   lagged-state energy (default: 3e-3). With both `A` and `B` free the lagged
#'   state is near low rank, so the two operators are not separately identified;
#'   the ridge pins them to a unique, bounded, balanced estimate. Larger values
#'   shrink the operators and trade fit for stability. Ignored when
#'   `symmetric = TRUE` (the symmetric operator is identified without it).
#' @param verbose Logical. If TRUE, print iteration progress (default: TRUE).
#' @param seed Random seed for reproducibility. Default \code{NULL} inherits
#'   the caller's RNG state; pass an integer for a reproducible fit.
#'
#' @return
#' A list of class `"dbn"` whose layout matches a full MCMC fit, but with:
#' - Single pseudo-draw (n_keep = 1)
#' - No MCMC variance parameters (`tau_A2`, `tau_B2`, and `g2` are all `NA`).
#'   `sigma2` is populated with an OLS-style residual-variance estimate from
#'   the converged ALS fit; it is not a posterior draw.
#' - `meta$sampler_used = "als"`
#' - `meta$uncertainty_available = FALSE`
#'
#' One semantic difference from an MCMC fit is important: the `Theta` element
#' of an ALS fit holds the one-step \emph{fitted prediction}
#' \eqn{A Z_{t-1} B^\top}, not a smoothed latent state. A full [dbn()] fit's
#' `Theta` is the model's latent state. The two are different quantities, so
#' do not use an ALS fit's `Theta` where a latent-state estimate is expected.
#' Relatedly, the ALS operator layout matches an MCMC fit:
#' \itemize{
#'   \item \code{fit$A} is a \strong{list of length 1}, with the single element
#'     a \code{[n_row, n_row, T]} cube. Access slice \code{t} with
#'     \code{fit$A[[1]][, , t]}, not \code{fit$A[, , t]}.
#'   \item \code{fit$B} has the same layout for the receiver operator.
#'   \item For \code{lambda = "static"} all time slices of \code{A[[1]]} are
#'     identical (one operator over time).
#'   \item For \code{lambda = numeric} or \code{lambda = "cv"} the slices vary
#'     across t (RW-smoothed time-varying operator).
#'   \item For \code{breaks = c(...)} the slices are piecewise-constant within
#'     each block.
#' }
#'
#' Compatible with downstream tools: [compute_irf()], [dbn_compute_snr()],
#' [dbn_coupling_rank_probs()], etc. On a point-estimate fit, tools that
#' report posterior probabilities (e.g. [dbn_coupling_rank_probs()]) return
#' degenerate 0/1 values rather than genuine probabilities.
#'
#' @section Relationship to MCMC:
#' `dbn_als(..., bootstrap = N)` is designed as a **fast frequentist
#' approximation** to the dynamic bilinear structure estimated by the full
#' Bayesian sampler. In Gaussian models with aligned penalties and a
#' resampling scheme that preserves the network and temporal dependence,
#' the bootstrap distribution approximates uncertainty in the same
#' *invariant* summaries the Bayesian model emphasises -- fitted values,
#' the per-time interaction summary \eqn{W_t = A_t B_t^\top}, actor
#' coupling, and the timing of structural shifts. It is **not** a
#' substitute for full posterior inference. Treat the intervals as
#' *bootstrap intervals*, not credible intervals; in our calibrations they
#' are systematically narrower than MCMC credible intervals in weakly
#' identified settings (small `n`, small `T`, large `r`, sparse data).
#'
#' To align ALS smoothing with the MCMC prior, the empirical-Bayes pilot
#' \eqn{\hat\lambda = \hat\sigma^2 / \hat\tau^2} corresponds approximately
#' to the MCMC random-walk innovation variance via
#' \deqn{\lambda \approx \sigma^2 / \tau_A^2.}
#' When `lambda = "eb"` (default for TV) or a numeric value, bump it
#' upward by a factor of 2-5 if the MCMC posterior of \eqn{\tau_A^2} is
#' much smaller than the pilot suggests; the ALS point estimate then
#' tracks the MCMC posterior mean more closely.
#'
#' For Bayesian credible intervals, run MCMC via [dbn()].
#' Use `dbn_als()` for fast exploration, model selection across many
#' candidate breakpoints, and as a warm start (`als_refit_warmstart()`).
#'
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}. The alternating least
#' squares estimator implemented here adapts the iterative block coordinate
#' descent estimation routine introduced in that work.
#'
#' @seealso [dbn()], [compute_irf()]
#'
#' @param bootstrap Integer. If \code{> 0}, runs that many bootstrap
#'   replicates after the point-estimate fit and attaches the result as a
#'   \code{$bootstrap} component of class \code{dbn_boot}, with bootstrap
#'   standard errors and 95\% intervals for \code{A}, \code{B}, and \code{M}.
#'   When the bootstrap is run, \code{meta$uncertainty_available} is set to
#'   \code{TRUE} and \code{meta$uncertainty_source} to \code{"bootstrap"}.
#'   Default \code{0} (no bootstrap, point estimate only). 200 or more is
#'   recommended when intervals are wanted; \code{[1, 30)} fits but warns.
#' @param bootstrap_type Character: \code{"block"} (resamples time indices
#'   with replacement) or \code{"parametric"} (simulates from the fitted
#'   model). Default \code{"block"}. Used only when \code{bootstrap > 0}.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 8, time = 10, seed = 1)
#' # sim$Z is the continuous series for family = "gaussian";
#' # sim$Y is the ordinal discretization for family = "ordinal"
#' fit_als <- dbn_als(sim$Z, family = "gaussian")
#' sampler_describe(fit_als)
#'
#' # also get bootstrap CIs in the same call
#' fit_als_ci <- dbn_als(sim$Z, family = "gaussian", bootstrap = 200)
#' fit_als_ci$bootstrap$ci_A_lo[1:3]
#' }
#'
#' @param lambda Smoothing penalty for time-varying ALS. Options:
#'   \describe{
#'     \item{\code{"static"} (default)}{Time-invariant SIR fit (Minhas & Hoff
#'       2025). One operator pair (A, B) for the entire panel; result tiled
#'       across time for layout compatibility.}
#'     \item{\code{"eb_corrected"}}{Single fit at the bias-corrected
#'       empirical-Bayes pilot
#'       \eqn{\hat\lambda = \hat\sigma^2 / \hat\tau^2_{\text{corr}}}.
#'       The raw windowed-static pilot estimate of \eqn{\tau^2} is biased
#'       upward by the sampling variance of the per-window static
#'       \eqn{\hat A}, \eqn{\hat B}. The corrected pilot subtracts that
#'       measurement-noise contribution, giving a smoothing strength that
#'       reflects true operator drift instead of within-window noise. Fast.
#'       \code{lambda_A} and \code{lambda_B} are computed separately and
#'       combined via precision-weighted average; the per-side values are
#'       recorded under \code{meta$tv$pilot}.}
#'     \item{\code{"cv_local"}}{Three-point local CV refinement around the
#'       bias-corrected pilot. Costs three warm-started fits, not a full
#'       grid. Recommended when the pilot is suspected of being off by a
#'       small constant factor.}
#'     \item{\code{"cv"}}{Rolling-origin cross-validation across a log-spaced
#'       \eqn{\lambda} grid centered on the EB pilot; selects via the 1-SE
#'       rule. The recommended principled default; slower (one fit per
#'       \eqn{\lambda} per held-out transition).}
#'     \item{\code{"path"}}{Fit the full \eqn{\lambda}-path via continuation
#'       (large-to-small \eqn{\lambda}, warm-starting each fit from the
#'       previous). Returns all fits in \code{meta$tv$path$fits}; the middle
#'       \eqn{\lambda} is reported as the headline fit.}
#'     \item{Numeric}{Single non-negative finite scalar: fit at that fixed
#'       \eqn{\lambda}. Larger \eqn{\lambda} -> smoother (more static-like)
#'       trajectory.}
#'   }
#'   Non-static options activate the RW(1)-smoothed time-varying ALS solver,
#'   producing genuinely time-varying \code{(A_t, B_t)} (slices 2..T of the
#'   returned cubes vary across time; slice 1 is the unestimated initial
#'   state).
#' @param mu Level ridge for TV-ALS (Frobenius penalty on \eqn{A_t, B_t}).
#'   If \code{NULL}, set automatically to \code{eta * h_bar} where
#'   \code{h_bar} is a data-curvature scale computed once at initialization.
#'   Required for identifiability of the asymmetric (directed) trajectory.
#' @param eta Default scaling of the level ridge when \code{mu = NULL}.
#' @param breaks Optional integer vector of breakpoint time indices for
#'   piecewise-constant ALS. Each breakpoint `b` separates block `[1..b]`
#'   from block `[b+1..T]`. When supplied, the function fits an independent
#'   static ALS on each block and assembles a piecewise-constant TV cube.
#'   Mutually exclusive with smoothing `lambda` (the latter is ignored if
#'   both are given). Returns a fit with `meta$sampler_used = "als_piecewise"`.
#' @param blocks Alias for `breaks` (matches the `dbn(piecewise, blocks=)`
#'   argument). Supply one or the other, not both.
#' @param parallel Logical (default `TRUE`). If `TRUE` and a non-sequential
#'   `future::plan` is active, the bootstrap replicates run in parallel via
#'   `future.apply::future_lapply`. Set to `FALSE` to force a serial loop.
#'   Falls back to sequential automatically under `devtools::load_all`.
#' @param gauge Scale gauge for TV-ALS (asymmetric only): \code{"penalty_min"}
#'   (default; minimizes the penalty in the gauge-ambiguous direction via a
#'   small Newton subproblem in the per-t scalar gauge), \code{"equal_norm_conditional"}
#'   (cheap rebalance, accepted only if it doesn't increase the objective),
#'   or \code{"none"}.
#' @param tv_max_iter Outer-loop iteration cap for TV-ALS (default 200).
#' @param tv_tol_obj Relative objective convergence tolerance (default 1e-5).
#' @param tv_tol_par Per-period parameter movement tolerance (default 1e-3).
#'
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_als <- function(data,
                    family = c("gaussian", "ordinal", "binary"),
                    symmetric = FALSE,
                    max_iter = 100,
                    tol = 1e-4,
                    ridge = 3e-3,
                    bootstrap = 0,
                    bootstrap_type = c("block", "parametric"),
                    lambda = "static",
                    mu = NULL,
                    eta = 1e-3,
                    gauge = c("penalty_min", "equal_norm_conditional", "none"),
                    breaks = NULL,
                    blocks = NULL,  # accepted as alias for `breaks` (matches dbn(piecewise, blocks=))
                    tv_max_iter = 200L,
                    tv_tol_obj = 1e-5,
                    tv_tol_par = 1e-3,
                    verbose = TRUE,
                    parallel = TRUE,
                    seed = NULL) {

	family <- match.arg(family)
	bootstrap_type <- match.arg(bootstrap_type)
	gauge <- match.arg(gauge)
	# cross-API alias: `dbn(model="piecewise", blocks=...)` uses `blocks=`,
	# while `dbn_als(breaks=...)` is the native name. Accept either; warn if
	# both are given.
	if (!is.null(blocks)) {
		if (!is.null(breaks)) {
			cli::cli_abort(c(
				"{.arg breaks} and {.arg blocks} are aliases; supply only one.",
				"i" = "{.arg breaks} is the {.fun dbn_als} native name; {.arg blocks} is the {.fun dbn} (piecewise) name."
			))
		}
		breaks <- blocks
	}
	# seed = NULL inherits the caller's RNG state so two back-to-back fits
	# give Monte-Carlo variability; pass an integer for reproducibility.
	.dbn_restore_seed <- .use_seed_locally(seed)
	on.exit(.dbn_restore_seed(), add = TRUE)
	if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
		cli::cli_abort("{.arg ridge} must be a single non-negative number.")
	}
	if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
	    max_iter != round(max_iter) || max_iter < 1) {
		cli::cli_abort("{.arg max_iter} must be a single positive whole number, got {.val {max_iter}}.")
	}
	if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
		cli::cli_abort("{.arg tol} must be a single positive number, got {.val {tol}}.")
	}
	if (!is.numeric(bootstrap) || length(bootstrap) != 1L || !is.finite(bootstrap) ||
	    bootstrap != round(bootstrap) || bootstrap < 0) {
		cli::cli_abort("{.arg bootstrap} must be a single non-negative whole number (0 = no bootstrap), got {.val {bootstrap}}.")
	}
	# parallel must be a scalar logical: NULL / NA / numeric would otherwise
	# silently slip through and break the future_lapply / lapply branch.
	if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
		cli::cli_abort(c(
			"{.arg parallel} must be a single non-NA logical ({.code TRUE} or {.code FALSE}).",
			"x" = "Got {.cls {class(parallel)[1]}} {.val {parallel}}.",
			"i" = "Pass {.code parallel = FALSE} to force a sequential bootstrap, or {.code parallel = TRUE} (the default) to route through {.fn future.apply::future_lapply}."
		))
	}
	bootstrap <- as.integer(bootstrap)
	if (bootstrap > 0L && bootstrap < 30L) {
		# the inner bootstrap helper would re-fire the same warning tagged
		# as `R`; suppress its duplicate via dbn.in_bootstrap.
		.dbn_prev_in_boot_outer <- getOption("dbn.in_bootstrap", FALSE)
		options(dbn.in_bootstrap = TRUE)
		on.exit(options(dbn.in_bootstrap = .dbn_prev_in_boot_outer), add = TRUE)
		cli::cli_warn(c(
			"{.arg bootstrap} = {bootstrap} is small for bootstrap inference.",
			"i" = "Standard errors and intervals from fewer than ~30 replicates are unstable; 200 or more is recommended."
		))
	}

	# load data from file path or use the array directly. A length-1
	# character is a path; a character array is malformed data, not a path.
	if (is.character(data) && length(data) == 1L) {
		cli::cli_inform("Loading data from: {.path {data}}")
		env <- new.env()
		load(data, envir = env)
		Y <- env$Y
		if (is.null(Y)) cli::cli_abort("Data file must contain object {.var Y}")
	} else {
		if (inherits(data, "formula")) {
			cli::cli_abort(c(
				"{.fun dbn_als} does not accept a formula interface ({.code dbn_als(Y ~ x)}).",
				"x" = "Got a {.cls formula} for {.arg data}.",
				"i" = "Pass the network data directly as a 3D ({.code [actors, actors, time]}) or 4D ({.code [actors, actors, relations, time]}) numeric array.",
				"i" = "Exogenous covariates are not currently supported."
			))
		}
		if (is.character(data)) {
			cli::cli_abort(c(
				"{.arg data} must be a numeric array or a single file path.",
				"x" = "Got a character array of length {length(data)}."
			))
		}
		if (inherits(data, c("igraph", "network"))) {
			cli::cli_abort(c(
				"{.fun dbn_als} does not accept {.cls {class(data)[1]}} objects directly.",
				"i" = "Coerce to a 3D adjacency array first (one slice per time period)."
			))
		}
		# accept and coerce a list of matrices (a common form from sliding-window
		# analyses or per-period observation containers). Each element must be
		# the same shape; we stack along the time dimension.
		if (is.list(data) && !is.data.frame(data) && !is.null(data) && length(data) > 0L &&
		    all(vapply(data, is.matrix, logical(1)))) {
			shapes <- vapply(data, function(m) paste(dim(m), collapse = "x"), character(1))
			if (length(unique(shapes)) != 1L) {
				cli::cli_abort(c(
					"{.arg data} is a list of matrices but they have inconsistent shapes.",
					"i" = "All elements must have the same dim; got {.val {unique(shapes)}}."
				))
			}
			d_one <- dim(data[[1]])
			Y <- array(unlist(data, use.names = FALSE),
				dim = c(d_one[1], d_one[2], 1, length(data)))
			cli::cli_inform("Coerced {.cls list} of {length(data)} matrices to a 4D array ({d_one[1]} x {d_one[2]} x 1 x {length(data)}).")
		} else if (is.list(data) && !is.data.frame(data)) {
			cli::cli_abort(c(
				"{.arg data} is a list but not a list-of-matrices.",
				"x" = "Found element types: {.val {unique(vapply(data, function(e) class(e)[1], character(1)))}}.",
				"i" = "Pass a 3D / 4D numeric array, or a list where every element is a matrix of the same shape."
			))
		} else {
			Y <- data
		}
	}

	# convert 3D to 4D if needed
	if (length(dim(Y)) == 3) {
		dim_orig <- dim(Y)
		Y <- array(Y, dim = c(dim_orig[1], dim_orig[2], 1, dim_orig[3]))
	} else if (length(dim(Y)) != 4) {
		cli::cli_abort(c(
			"{.arg data} must be a 3D or 4D numeric array.",
			"x" = "Got class {.cls {class(data)[1]}} with dim {.val {paste(dim(data), collapse = 'x')}}."
		))
	}
	warn_filled_diagonal(Y)

	# preprocess data
	pre <- shared_preprocess(Y, family = family)
	Z <- pre$Z  # latent continuous scores
	M <- pre$M  # baseline
	dims <- pre$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	p <- dims$p
	Tt <- dims$Tt

	if (symmetric && n_row != n_col) {
		cli::cli_abort("Symmetric model requires square network (n_row == n_col).")
	}

	# ALS needs at least one transition (T >= 2); the static `for (t in 2:Tt)`
	# loop reverses to `c(2, 1)` otherwise and crashes.
	if (Tt < 2L) {
		cli::cli_abort(c(
			"{.fun dbn_als} requires at least 2 time points (got T = {Tt}).",
			"i" = "ALS estimates the bilinear operator from t -> t+1 transitions; a single period supplies no transitions."
		))
	}

	# piecewise ALS branch: if `breaks` is given, fit independent static ALS
	# models per block and assemble into a piecewise-constant TV-cube. Exits
	# early without touching the lambda machinery.
	if (!is.null(breaks)) {
		if (is.character(lambda) && lambda != "static") {
			cli::cli_warn(c(
				"{.arg breaks} and {.code lambda = \"{lambda}\"} are both specified; the smoothing penalty is ignored when piecewise breaks are given.",
				"i" = "Drop {.arg breaks} or use {.code lambda = \"static\"} explicitly."
			))
		}
		return(.dbn_als_piecewise(Y, family = family, symmetric = symmetric,
			breaks = breaks, ridge = ridge, max_iter = max_iter,
			verbose = verbose, seed = seed,
			bootstrap = bootstrap, bootstrap_type = bootstrap_type,
			parallel = parallel))
	}

	if (verbose) {
		cli::cli_h3("ALS Exploration")
		cli::cli_text("Network: {n_row} x {n_col}, {p} relation{?s}, {Tt} time point{?s}")
		cli::cli_text("Tolerance: {tol}, Max iterations: {max_iter}")
	}

	# initialize A and B: identity or small random
	A_est <- diag(n_row) + 0.01 * matrix(rnorm(n_row^2), n_row, n_row)
	if (symmetric) {
		A_est <- 0.5 * (A_est + t(A_est))  # symmetrize
		B_est <- A_est
	} else {
		B_est <- diag(n_col) + 0.01 * matrix(rnorm(n_col^2), n_col, n_col)
	}

	# extract latent Z: from 4D [n_row, n_col, p, Tt] to list of T matrices
	# Z[, , 1, t] is the first relation at time t; use all relations averaged
	Z_list <- list()
	Z_prev_list <- list()
	for (t in 1:Tt) {
		# average latent scores across relations (alternative: Z[, , 1, t] for first relation only)
		Z_t_all <- Z[, , , t]  # n_row x n_col x p
		Z_list[[t]] <- apply(Z_t_all, c(1, 2), mean)  # Average to get n_row x n_col matrix
		if (t > 1) {
			Z_prev_list[[t]] <- Z_list[[t-1]]
		}
	}

	# fit the bilinear operator by ALS. dbn_als() and the bootstrap
	# refit share als_core_fit(): monotone line-search descent for the
	# symmetric operator, penalized extrapolated alternating least squares
	# for the asymmetric operator. Both are hardened against the (A, B) ->
	# (cA, c^-1 B) scale invariance that otherwise makes the operators
	# diverge.
	fit <- als_core_fit(Z_list, Z_prev_list, Tt, n_row, n_col, symmetric,
	                    A_init = A_est, B_init = B_est,
	                    max_iter = max_iter, tol = tol,
	                    ridge_frac = ridge, verbose = isTRUE(verbose))
	A_est       <- fit$A
	B_est       <- fit$B
	converged   <- fit$converged
	solve_log   <- fit$solve_log
	n_iter_run  <- fit$iterations %||% max_iter

	if (!converged && verbose) {
		cli::cli_warn("ALS did not converge after {max_iter} iterations.")
	}

	# summarize ALS linear-solve conditioning (from the final iteration)
	als_solve_summary <- summarize_als_solves(solve_log)
	if (als_solve_summary$n_ridge > 0 || als_solve_summary$n_rank_deficient > 0) {
		A_zero <- all(abs(A_est) < sqrt(.Machine$double.eps))
		B_zero <- all(abs(B_est) < sqrt(.Machine$double.eps))
		if (A_zero || B_zero) {
			which_zero <- if (A_zero && B_zero) "A and B" else if (A_zero) "A" else "B"
			cli::cli_warn(c(
				"ALS produced a near-zero {which_zero} matrix.",
				"i" = "The design is rank-deficient; the influence structure is not identified from these data."
			))
		} else if (als_solve_summary$n_solves > 0) {
			frac <- als_solve_summary$n_rank_deficient / als_solve_summary$n_solves
			if (frac > 0.05) {
				cli::cli_warn(c(
					"{als_solve_summary$n_rank_deficient} of {als_solve_summary$n_solves} ALS updates were rank-deficient; ridge regularization was applied.",
					"i" = "Bootstrap and credible uncertainty for A and B may be understated for these data."
				))
			}
		}
	}

	# estimate Theta_t from current A, B, Z
	Theta_est <- array(NA, dim = c(n_row, n_col, p, Tt))
	for (t in 1:Tt) {
		Theta_t <- A_est %*% Z_list[[t]] %*% t(B_est)
		for (r in 1:p) {
			Theta_est[, , r, t] <- Theta_t
		}
	}

	# estimate M from residuals
	M_est <- M  # keep as-is from preprocessing

	# estimate sigma2 from residuals
	resid_sq <- 0
	n_resid <- 0
	for (t in 2:Tt) {
		pred_t <- A_est %*% Z_prev_list[[t]] %*% t(B_est)
		resid_t <- Z_list[[t]] - pred_t
		resid_sq <- resid_sq + sum(resid_t^2)
		n_resid <- n_resid + length(resid_t)
	}
	sigma2_est <- resid_sq / max(n_resid, 1)
	if (sigma2_est < 1e-6) sigma2_est <- 1.0

	# build output object (single pseudo-draw, n_keep = 1)
	n_keep <- 1L
	time_keep <- seq_len(Tt)

	# ============================================================
	# RW-smoothed time-varying ALS branch (Stage 1A).
	# when `lambda` is a numeric scalar, use the static (A_est, B_est)
	# computed above as the warm-start and run the TV solver to produce
	# genuinely time-varying (A_t, B_t).
	#
	# lambda options: "static" (time-invariant SIR), numeric (fixed lambda),
	# error out with a clear "not yet implemented" message.
	# ============================================================
	tv_meta <- NULL
	pilot <- NULL
	use_tv <- FALSE
	lambda_numeric <- NA_real_

	# reject malformed lambda values up front (-1, NA, c(1,2), Inf would
	# otherwise fall through to a static fit silently)
	if (!is.null(lambda)) {
		if (is.character(lambda)) {
			allowed_lambda <- c("static", "eb_corrected", "cv_local", "cv", "path")
			if (length(lambda) != 1L || !(lambda %in% allowed_lambda)) {
				cli::cli_abort(c(
					"{.code lambda = {.val {lambda}}} is not recognized.",
					"i" = "Allowed: {.val static} (default), {.val eb_corrected}, {.val cv_local}, {.val cv}, {.val path}, or a single non-negative finite number."
				))
			}
		} else if (is.numeric(lambda)) {
			if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
				cli::cli_abort(c(
					"{.arg lambda} must be a single non-negative finite number, or one of the string choices.",
					"x" = "Got {.val {lambda}}.",
					"i" = "Use {.val static} for the time-invariant fit; positive numerics give RW-smoothed time-varying ALS."
				))
			}
		} else {
			cli::cli_abort(c(
				"{.arg lambda} must be character or numeric.",
				"x" = "Got class {.cls {class(lambda)}}."
			))
		}
	}

	# resolve mu before dispatching CV/path/eb; mu = NULL in the CV path
	# silently crashes every fold and returns the largest grid value.
	if (is.character(lambda) && lambda != "static") {
		if (is.null(mu)) {
			h_bar <- mean(sapply(1:Tt, function(t) {
				M_t <- Z_list[[t]] %*% t(B_est)
				sum(M_t^2) / n_row
			}))
			mu <- eta * h_bar
		}
	}

	if (is.character(lambda)) {
		if (lambda == "static") {
			use_tv <- FALSE
		} else {
			# compute pilot, build grid, then dispatch
			pilot <- .tv_als_pilot_lambda(Y, family, symmetric, n_row, n_col, p, Tt)
			lambda_pilot <- pilot$lambda_pilot
			lambda_grid <- .tv_als_make_lambda_grid(lambda_pilot)

			if (lambda == "eb_corrected") {
				# single fit at the bias-corrected pilot
				use_tv <- TRUE
				lambda_numeric <- lambda_pilot
				if (verbose) cli::cli_inform("lambda = \"eb_corrected\" -> pilot = {sprintf('%.4g', lambda_pilot)} (raw = {sprintf('%.4g', pilot$lambda_raw)}, correction factor = {sprintf('%.2f', pilot$correction_factor)})")
			} else if (lambda == "cv_local") {
				# three-point local CV around the corrected pilot
				use_tv <- TRUE
				lambda_candidates <- lambda_pilot * c(0.5, 1, 2)
				lambda_numeric <- .tv_als_select_lambda_local(
					Y, family, symmetric, lambda_candidates,
					mu, eta, gauge, tv_max_iter, tv_tol_obj, tv_tol_par,
					verbose = isTRUE(verbose)
				)
				if (verbose) cli::cli_inform("lambda = \"cv_local\" -> selected {sprintf('%.4g', lambda_numeric)} from {paste(sprintf('%.3g', lambda_candidates), collapse = ', ')}")
			} else if (lambda == "cv") {
				# CV-select from the grid, then fit at lambda_1se.
				# auto-extend the grid if CV hits the boundary: doubles the range
				# in the indicated direction and re-runs (max 2 extensions).
				# if CV picks the boundary of the lambda grid, extend the grid
				# in that direction and re-run (up to 2 extensions)
				if (verbose) cli::cli_inform("lambda = \"cv\" -> rolling-origin CV across {length(lambda_grid)} values")
				# warn upfront when CV is combined with bootstrap: the
				# combined runtime is roughly (CV reps) x (bootstrap reps),
				# which can hang for many minutes with no feedback if the
				# user did not anticipate the multiplicative cost
				if (isTRUE(bootstrap > 0L)) {
					cli::cli_inform(c(
						"i" = "Combining {.code lambda = \"cv\"} with {.code bootstrap = {bootstrap}}: total cost is approximately {length(lambda_grid)} (CV) {.field x} {bootstrap} (bootstrap) fits.",
						"i" = "Each fit is the full ALS run; expect minutes to an hour. Silence this with {.code options(dbn.suppress_cv_progress = TRUE)} or pre-select lambda first."
					))
				}
				cv <- .tv_als_select_lambda_cv(Y, family, symmetric, lambda_grid,
					mu, eta, gauge, tv_max_iter, tv_tol_obj, tv_tol_par,
					verbose = isTRUE(verbose))
				ext_count <- 0L
				while (isTRUE(cv$boundary_hit) && ext_count < 2L) {
					ext_count <- ext_count + 1L
					if (cv$idx_1se == 1L) {
						# smoothest end -- extend upward (larger lambda)
						extend_factor <- 10
						lambda_top <- lambda_grid[1] * extend_factor
						lambda_grid <- exp(seq(log(lambda_top), log(lambda_grid[length(lambda_grid)]),
							length.out = length(lambda_grid) + 4L))
						if (isTRUE(verbose)) cli::cli_alert_info("CV hit smoothest boundary; auto-extending grid upward to lambda_max = {.val {signif(lambda_top, 4)}} ({ext_count}/2).")
					} else {
						# roughest end -- extend downward (smaller lambda)
						extend_factor <- 10
						lambda_bot <- lambda_grid[length(lambda_grid)] / extend_factor
						lambda_grid <- exp(seq(log(lambda_grid[1]), log(lambda_bot),
							length.out = length(lambda_grid) + 4L))
						if (isTRUE(verbose)) cli::cli_alert_info("CV hit roughest boundary; auto-extending grid downward to lambda_min = {.val {signif(lambda_bot, 4)}} ({ext_count}/2).")
					}
					cv <- .tv_als_select_lambda_cv(Y, family, symmetric, lambda_grid,
						mu, eta, gauge, tv_max_iter, tv_tol_obj, tv_tol_par,
						verbose = isTRUE(verbose))
				}
				# use cli_alert_danger AFTER all per-fold noise so it survives
				# the warning blizzard. Distinct from the in-CV cli_warn which
				# can get muffled.
				if (isTRUE(cv$boundary_hit) && isTRUE(verbose)) {
					cli::cli_alert_danger(c(
						"CV still at grid boundary after {ext_count} auto-extensions: lambda_1se = {.val {signif(cv$lambda_1se, 4)}} at idx {cv$idx_1se} of {length(cv$lambda_grid)}.",
						"i" = "The CV loss is monotone over the grid range; the chosen lambda may not be the true optimum.",
						"i" = "Inspect {.code fit$meta$tv$cv$losses} -- if monotone, the data may have minimal time-varying signal."
					))
				}
				use_tv <- TRUE
				lambda_numeric <- cv$lambda_1se
				if (verbose) cli::cli_inform("CV selected lambda = {sprintf('%.4g', lambda_numeric)} (min: {sprintf('%.4g', cv$lambda_min)})")
				tv_meta <- list(cv = cv, cv_n_extensions = ext_count)
			} else if (lambda == "path") {
				# fit full path and return
				if (verbose) cli::cli_inform("lambda = \"path\" -> fitting {length(lambda_grid)} values via continuation")
				path <- .tv_als_lambda_path(Y, family, symmetric, lambda_grid,
					mu, eta, gauge, tv_max_iter, tv_tol_obj, tv_tol_par,
					verbose = isTRUE(verbose))
				# for "path", return the lambda-path object on top of a chosen point
				# default chosen point: middle lambda
				idx_mid <- ceiling(length(lambda_grid) / 2)
				tv_meta <- list(path = path, path_idx = idx_mid)
				use_tv <- TRUE
				lambda_numeric <- lambda_grid[idx_mid]
			}
		}
	} else if (is.numeric(lambda)) {
		use_tv <- TRUE
		lambda_numeric <- lambda
	}
	if (use_tv) {
		# choose mu via h_bar from the warm-start design (only if not already
		# resolved above for the CV/path/eb pre-dispatch)
		if (is.null(mu)) {
			# h_bar = median over t of (1/n) tr(M_t M_t^T), using the static
			# warm-start (one M_t per t; in static mode all M_t equal M_global,
			# so compute it once from the data)
			h_bar <- mean(sapply(1:Tt, function(t) {
				M_t <- Z_list[[t]] %*% t(B_est)
				sum(M_t^2) / n_row
			}))
			mu <- eta * h_bar
		}

		# validate mu and eta; negative values would silently drop into a
		# chol() failure deeper in the solver
		if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) || mu < 0) {
			cli::cli_abort("{.arg mu} must be a single non-negative finite number.")
		}
		if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta < 0) {
			cli::cli_abort("{.arg eta} must be a single non-negative finite number.")
		}

		# warm-start: broadcast static A, B across t
		A_init <- array(rep(c(A_est), Tt), dim = c(n_row, n_row, Tt))
		B_init <- array(rep(c(B_est), Tt), dim = c(n_col, n_col, Tt))
		M_init <- array(NA, dim = c(n_row, n_col, p))
		# M_est is [n_row, n_col, p]; assign per-relation slice into M_init[,,r]
		for (r in 1:p) M_init[, , r] <- M_est[, , r]

		# stage 4: binary/ordinal supported via deterministic Z re-imputation
		# in the outer loop (.expected_Z). The TV solver receives the latent
		# Z as the "working response" and Y_obs as the actual observed Y.
		Y_for_tv <- Y
		Y_obs_for_tv <- NULL
		if (family != "gaussian") {
			# for binary: Z_init from the truncated normal mean given theta_0 = 0
			# for ordinal: Z_init from the empirical-quantile rule
			Y_obs_for_tv <- Y
			Z_init_arr <- Y
			if (family == "binary") {
				# binary: rough init -- Z = sign(Y - 0.5)
				Z_init_arr <- ifelse(is.na(Y), NA, ifelse(Y > 0.5, 1, -1))
			} else if (family == "ordinal") {
				# ordinal: convert categories to standard normal quantiles
				nonna <- !is.na(Y)
				Y_int <- as.integer(round(Y))
				cats <- sort(unique(Y_int[nonna]))
				K <- length(cats)
				if (K >= 2) {
					rank_centers <- qnorm((seq_along(cats) - 0.5) / K)
					Z_init_arr[nonna] <- rank_centers[match(Y_int[nonna], cats)]
				}
			}
			Y_for_tv <- Z_init_arr
		}

		# if a lambda-path was already computed (lambda = "path"), pull the
		# selected fit out of the path; otherwise fit fresh at lambda_numeric.
		if (!is.null(tv_meta$path)) {
			tv_fit <- tv_meta$path$fits[[tv_meta$path_idx]]
		} else if (symmetric) {
			# stage 2 (Gaussian) / Stage 5 (binary/ordinal): symmetric solver
			sym_fit <- .dbn_als_tv_symmetric_lbfgs(
				Y_for_tv, lambda = lambda_numeric, mu = mu,
				max_iter = tv_max_iter, tol_obj = tv_tol_obj,
				init = list(A_init = A_init, M_init = M_init),
				family = family,
				Y_obs = Y_obs_for_tv,
				verbose = isTRUE(verbose)
			)
			# pack into the same shape as directed TV fit
			tv_fit <- list(
				A = sym_fit$A,
				B = sym_fit$A,  # B = A in symmetric case
				M = sym_fit$M,
				Phi = sym_fit$Phi,
				Omega = sym_fit$Omega,
				objective = sym_fit$objective,
				converged = sym_fit$converged,
				n_iter = sym_fit$n_iter,
				stationarity_residual = NA_real_,  # not computed for L-BFGS
				lambda = sym_fit$lambda,
				mu = sym_fit$mu
			)
		} else {
			tv_fit <- .dbn_als_tv_fixed_lambda(
				Y_for_tv, family = family, symmetric = symmetric,
				lambda = lambda_numeric, mu = mu,
				max_iter = tv_max_iter, tol_obj = tv_tol_obj, tol_par = tv_tol_par,
				gauge = gauge,
				init = list(A_init = A_init, B_init = B_init, M_init = M_init),
				Y_obs = Y_obs_for_tv,
				z_update = "mean",
				verbose = isTRUE(verbose)
			)
		}

		# replace static estimates with TV estimates for the output
		A_est <- tv_fit$A
		B_est <- tv_fit$B
		# tv_fit$M is [n_row, n_col, p]
		M_est <- tv_fit$M
		# theta from TV: Theta_t = M + Phi^state_t * A_t Phi^state_{t-1} B_t^T ... actually
		# we re-derive Theta_t = M + A_t Phi^state_{t-1} B_t^T (one-step prediction)
		# (matches the original SIR static formulation, just per-t now)
		Theta_est <- array(NA, dim = c(n_row, n_col, p, Tt))
		Phi_state <- tv_fit$Phi
		for (r in 1:p) {
			Theta_est[, , r, 1] <- M_est[, , r] + Phi_state[, , 1]
			for (t in 2:Tt) {
				pred_t <- A_est[, , t] %*% Phi_state[, , t - 1] %*% t(B_est[, , t])
				Theta_est[, , r, t] <- M_est[, , r] + pred_t
			}
		}
		# residual variance under TV fit
		resid_sq2 <- 0; n_resid2 <- 0
		for (t in 2:Tt) {
			pred_t <- A_est[, , t] %*% Phi_state[, , t - 1] %*% t(B_est[, , t])
			resid_t <- Phi_state[, , t] - pred_t
			resid_sq2 <- resid_sq2 + sum(resid_t^2)
			n_resid2 <- n_resid2 + length(resid_t)
		}
		sigma2_est <- resid_sq2 / max(n_resid2, 1)
		if (sigma2_est < 1e-6) sigma2_est <- 1.0

		# merge into the tv_meta we may have already populated (cv, path)
		tv_meta <- c(tv_meta %||% list(), list(
			lambda_used = lambda_numeric,
			lambda_requested = lambda,
			mu = mu,
			tv_converged = tv_fit$converged,
			tv_n_iter = tv_fit$n_iter,
			tv_objective = tv_fit$objective,
			tv_max_iter = tv_max_iter,
			tv_tol_obj = tv_tol_obj,
			stationarity_residual = tv_fit$stationarity_residual
		))
		if (!is.null(pilot)) {
			tv_meta$pilot <- list(
				lambda_raw = pilot$lambda_raw,
				lambda_corrected = pilot$lambda_common,
				lambda_A = pilot$lambda_A,
				lambda_B = pilot$lambda_B,
				sigma2 = pilot$sigma2,
				tau_A2_raw = pilot$tau_A2_raw,
				tau_B2_raw = pilot$tau_B2_raw,
				tau_A2_corrected = pilot$tau_A2_corrected,
				tau_B2_corrected = pilot$tau_B2_corrected,
				tau_A2_meas_trace = pilot$tau_A2_meas_trace,
				tau_B2_meas_trace = pilot$tau_B2_meas_trace,
				correction_factor = pilot$correction_factor,
				boundary_hit = pilot$lambda_boundary_hit
			)
		}
	}

	# reshape for output: each element in list is [n_row, n_row, n_time_keep]
	A_store <- list(array(NA, dim = c(n_row, n_row, length(time_keep))))
	if (use_tv) {
		A_store[[1]] <- A_est  # already n_row x n_row x Tt
	} else {
		A_store[[1]][, , ] <- A_est
	}

	B_store <- list(array(NA, dim = c(n_col, n_col, length(time_keep))))
	if (use_tv) {
		B_store[[1]] <- B_est
	} else {
		B_store[[1]][, , ] <- B_est
	}

	M_store <- array(NA, dim = c(n_row, n_col, p, n_keep))
	if (use_tv) {
		M_store[, , , 1] <- M_est
	} else {
		M_store[, , , 1] <- M_est
	}

	Theta_store <- array(NA, dim = c(n_row, n_col, p, length(time_keep), n_keep))
	Theta_store[, , , , 1] <- Theta_est

	sigma2_store <- c(sigma2_est)
	tau_A2_store <- c(NA_real_)
	tau_B2_store <- c(NA_real_)
	g2_store <- c(NA_real_)

	# build output list (matching dbn_dynamic structure). mirror the
	# symmetric flag to a top-level slot so downstream code can branch on
	# isTRUE(fit$symmetric) on every fit type, not just the MCMC paths.
	out <- list(
		model = "dynamic",
		symmetric = isTRUE(symmetric),
		family = family,
		Y = Y,
		dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
					is_bipartite = (n_row != n_col), is_symmetric = symmetric),
		settings = list(
			nscan = 0,
			burn = 0,
			odens = 1,
			time_thin = 1,
			store_z = FALSE,
			draws = n_keep,
			ridge = ridge,
			max_iter = max_iter,
			tol = tol
		),
		meta = list(
			family = family,
			dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
						is_bipartite = (n_row != n_col), is_symmetric = symmetric),
			draws = n_keep,
			settings = list(nscan = 0, burn = 0, odens = 1, time_thin = 1),
			sampler_requested = "als",
			sampler_used = if (use_tv) "als_tv" else "als",
			uncertainty_available = FALSE,
			approximation_note = if (use_tv)
				sprintf("RW-smoothed time-varying ALS point estimate (lambda=%.4g, mu=%.4g); no credible intervals", lambda_numeric, mu)
				else "Static ALS point estimate; no credible intervals",
			als_solve_summary = als_solve_summary,
			converged = if (use_tv) tv_meta$tv_converged else converged,
			tv = tv_meta
		),
		params = data.frame(sigma2 = sigma2_store),
		draws = list(
			pars = data.frame(sigma2 = sigma2_store),
			misc = list(A = A_store, B = B_store, M = M_store)
		),
		time_kept = time_keep,
		A = A_store,
		B = B_store,
		M = M_store,
		Theta = Theta_store,
		sigma2 = sigma2_store,
		tau_A2 = tau_A2_store,
		tau_B2 = tau_B2_store,
		g2 = g2_store,
		n_iter = n_iter_run,
		burn = 0,
		thin = 1,
		time_thin = 1
	)

	class(out) <- "dbn"

	# optional bootstrap CIs
	if (bootstrap > 0L && identical(out$meta$sampler_used, "als_tv")) {
		boot_type_tv <- switch(
			bootstrap_type,
			"block" = "residual_block",
			"parametric" = "fixed",
			bootstrap_type
		)
		boot <- .dbn_als_tv_bootstrap(out, R = bootstrap, type = boot_type_tv,
			seed = seed, verbose = isTRUE(verbose), parallel = parallel)
		out <- .dbn_attach_bootstrap(out, boot,
			label = sprintf("TV-ALS point estimate plus %d-replicate %s bootstrap CIs (lambda=%.4g).",
				boot$n_valid, boot_type_tv, lambda_numeric),
			source = "bootstrap")
		return(out)
	}

	if (bootstrap > 0L) {
		boot <- .dbn_als_bootstrap(out, R = bootstrap, type = bootstrap_type,
			seed = seed, verbose = isTRUE(verbose), parallel = parallel)
		out <- .dbn_attach_bootstrap(out, boot,
			label = sprintf("Point estimate plus %d-replicate %s bootstrap CIs.",
				boot$n_valid, bootstrap_type),
			source = "bootstrap")
	}

	out
}

####
#' Regularized least-squares solver for an ALS sub-problem
#'
#' @description
#' Solves a single ALS row/column least-squares sub-problem with a numerically
#' stable solver ladder: column equilibration, then a QR solve for
#' well-conditioned full-rank designs, or ridge (Tikhonov) regularization for
#' rank-deficient or ill-conditioned designs. Ridge is the posterior mode under
#' a Gaussian prior, keeping the fallback coherent with the package's Bayesian
#' framing. This replaces a normal-equation solve (which squared the design
#' condition number) and a Moore-Penrose pseudoinverse fallback (which collapsed
#' non-estimable directions to zero with no record).
#'
#' @param X Design matrix.
#' @param y Response vector.
#' @param kappa_qr_max Condition-number ceiling for the QR path (default 1e8).
#' @param kappa_target Target condition number for the ridge path (default 1e6).
#' @return List with `beta` (coefficients on the original scale), `method`
#'   (`"qr"` or `"ridge"`), `lambda` (ridge penalty; 0 for QR), `rank`, `p`,
#'   and `kappa_x` (design condition number).
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}.
#' @keywords internal
#' @noRd
als_solve <- function(X, y, kappa_qr_max = 1e8, kappa_target = 1e6) {
	p <- ncol(X)
	# column equilibration improves conditioning without changing the fit
	col_scale <- sqrt(colSums(X^2))
	col_scale[!is.finite(col_scale) | col_scale == 0] <- 1
	Xs <- sweep(X, 2L, col_scale, "/")

	s <- svd(Xs, nu = 0L, nv = 0L)$d
	smax <- if (length(s)) max(s) else 0
	tol <- max(dim(Xs)) * .Machine$double.eps * max(smax, 1)
	rank <- sum(s > tol)
	smin <- if (rank > 0L) min(s[s > tol]) else 0
	kappa_x <- if (rank < p || smin == 0) Inf else smax / smin

	if (rank == p && is.finite(kappa_x) && kappa_x <= kappa_qr_max) {
		# well-conditioned, full rank: QR avoids squaring the condition number
		beta_s <- qr.solve(Xs, y)
		method <- "qr"
		lambda <- 0
	} else {
		# rank-deficient or ill-conditioned: ridge toward a target condition number
		smin2 <- if (rank == p) smin^2 else 0
		lambda <- max(0, (smax^2 - kappa_target * smin2) / (kappa_target - 1))
		if (!is.finite(lambda) || lambda <= 0) {
			lambda <- max(1e-8 * mean(s^2), 1e-10 * max(smax^2, 1))
		}
		beta_s <- solve(crossprod(Xs) + lambda * diag(p), crossprod(Xs, y))
		method <- "ridge"
	}

	list(beta = as.numeric(beta_s / col_scale),
	     method = method, lambda = lambda,
	     rank = rank, p = p, kappa_x = kappa_x)
}

####
#' Ridge-regularized least-squares solver for an ALS sub-problem
#'
#' @description
#' A companion to [als_solve()] that applies a fixed, absolute ridge penalty
#' \eqn{\lambda} supplied by the caller, rather than ridging only as a
#' rank-deficiency fallback. The asymmetric bilinear ALS sub-problems are
#' routinely ill-conditioned -- the lagged state \eqn{Z_{t-1}} is typically
#' near low rank, so \eqn{A} and \eqn{B} are not separately identified -- and a
#' shared, fixed penalty makes every sub-problem solve a single coherent
#' penalized objective \eqn{\sum_t \|Z_t - A Z_{t-1} B^\top\|^2 + \lambda
#' (\|A\|^2 + \|B\|^2)}. That objective has a unique, bounded minimizer, so the
#' operator estimates cannot inflate or drift along weakly identified
#' directions. The penalty is the posterior mode under a Gaussian prior on the
#' operator row, coherent with the dynamic model's own prior on \eqn{A_t} and
#' \eqn{B_t}.
#'
#' @param X Design matrix.
#' @param y Response vector.
#' @param lambda Absolute ridge penalty (the same value must be used across all
#'   sub-problems of one fit for the penalized objective to be coherent).
#' @return List with the same fields as [als_solve()]; `method` is always
#'   `"ridge"`, and `rank` reports the genuine numerical rank (independent of
#'   the deliberate ridge) so rank-deficiency diagnostics remain meaningful.
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}.
#' @keywords internal
#' @noRd
als_solve_ridged <- function(X, y, lambda) {
	p <- ncol(X)
	# the penalty is on the raw coefficients, so solve the raw normal
	# equations (column equilibration would rescale the penalty per column
	# and break the coherence of the shared penalized objective)
	XtX <- crossprod(X)
	# genuine numerical rank, independent of the deliberate ridge
	ev <- eigen(XtX, symmetric = TRUE, only.values = TRUE)$values
	emax <- max(ev, 0)
	tol_r <- length(ev) * .Machine$double.eps * max(emax, 1)
	rank <- sum(ev > tol_r)
	emin <- if (rank > 0L) min(ev[ev > tol_r]) else 0
	if (!is.finite(lambda) || lambda < 0) lambda <- 0
	beta_s <- solve(XtX + lambda * diag(p), crossprod(X, y))
	list(beta = as.numeric(beta_s),
	     method = "ridge", lambda = lambda, rank = rank, p = p,
	     kappa_x = if (rank < p || emin <= 0) Inf else sqrt(emax / emin))
}

####
#' Summarize ALS linear-solve diagnostics
#'
#' @param solve_log List of `als_solve()` results from one ALS iteration.
#' @return List of aggregate solver diagnostics.
#' @keywords internal
#' @noRd
summarize_als_solves <- function(solve_log) {
	if (length(solve_log) == 0L) {
		return(list(n_solves = 0L, n_qr = 0L, n_ridge = 0L,
		            n_rank_deficient = 0L, worst_kappa = NA_real_,
		            median_lambda = 0))
	}
	methods <- vapply(solve_log, function(s) s$method, character(1))
	lambdas <- vapply(solve_log, function(s) s$lambda, numeric(1))
	kappas  <- vapply(solve_log, function(s) s$kappa_x, numeric(1))
	deficient <- vapply(solve_log, function(s) s$rank < s$p, logical(1))
	list(
		n_solves = length(solve_log),
		n_qr = sum(methods == "qr"),
		n_ridge = sum(methods == "ridge"),
		n_rank_deficient = sum(deficient),
		worst_kappa = max(kappas),
		median_lambda = if (any(methods == "ridge")) stats::median(lambdas[methods == "ridge"]) else 0
	)
}

####
#' Core ALS estimator for the bilinear operator (internal)
#'
#' @description
#' Shared estimation core for [dbn_als()] and [als_refit_warmstart()]. Fits
#' the bilinear transition \eqn{Z_t \approx A Z_{t-1} B^\top} (or, for the
#' symmetric model, \eqn{Z_t \approx A Z_{t-1} A^\top}) by alternating least
#' squares.
#'
#' The bilinear model has the exact invariance \eqn{(A, B) \to (cA, c^{-1}B)}.
#' Plain ALS lets this scale random-walk so the operators explode or collapse,
#' and the symmetric linearisation \eqn{A_{k+1} = \arg\min_A \|Z_t - A Z_{t-1}
#' A_k^\top\|} has a spurious fixed point at \eqn{A = 0}. Two devices remove
#' both pathologies:
#'
#' - Symmetric operator: monotone line-search descent on the true quartic
#'   objective \eqn{f(A) = \sum_t \|Z_t - A Z_{t-1} A^\top\|_F^2}. Each step
#'   line-searches the linearised ALS direction and the negative-gradient
#'   direction and accepts only a strict decrease. Since \eqn{f(0) > f(I)} and
#'   descent is monotone from the identity warm start, the iterate cannot
#'   collapse to the \eqn{A = 0} fixed point.
#' - Asymmetric operator: penalized alternating least squares. With \eqn{A} and
#'   \eqn{B} both free the lagged state is near low rank, so the two factors are
#'   not separately identified and plain ALS drifts. A fixed ridge penalty
#'   \eqn{\lambda(\|A\|^2 + \|B\|^2)} makes every sweep solve one coherent
#'   penalized objective with a unique, bounded minimizer; the penalty also
#'   pins the \eqn{(A, B) \to (cA, c^{-1}B)} gauge to the balanced minimum-norm
#'   representative. Each sweep is followed by a line search that extrapolates
#'   the ALS direction, which accelerates the otherwise slow descent down the
#'   shallow valley of this biquadratic problem.
#'
#' This is the block coordinate descent estimator of Minhas and Hoff (2025),
#' specialised to the dynamic bilinear network and hardened against the scale
#' invariance.
#'
#' Note: the asymmetric penalized ALS converges only linearly, with a rate set
#' by the (typically severe) conditioning of the lagged-state design. The
#' estimate stabilises quickly to a useful point estimate, but very strict
#' tolerances may not be reached within `max_iter`; `converged` reports this
#' honestly.
#'
#' @param Z_list,Z_prev_list Latent-state lists indexed by time.
#' @param Tt Number of time points.
#' @param n_row,n_col Operator dimensions.
#' @param symmetric Logical; if TRUE enforce B = A.
#' @param A_init,B_init Warm-start operators.
#' @param max_iter Maximum iterations.
#' @param tol Relative-objective convergence tolerance.
#' @param ridge_frac Asymmetric ridge penalty, as a fraction of the mean lagged
#'   -state energy; ignored for the symmetric operator.
#' @param verbose Logical; print per-iteration progress.
#' @return List with `A`, `B`, `converged`, `solve_log`, and `objective`.
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}.
#' @keywords internal
#' @noRd
als_core_fit <- function(Z_list, Z_prev_list, Tt, n_row, n_col, symmetric,
                         A_init, B_init, max_iter = 50L, tol = 1e-4,
                         ridge_frac = 3e-3, verbose = FALSE) {

	# joint one-step reconstruction objective sum_t || Z_t - A Z_{t-1} B^T ||_F^2
	obj_fn <- function(A, B) {
		s <- 0
		for (t in 2:Tt) {
			R <- Z_list[[t]] - A %*% Z_prev_list[[t]] %*% t(B)
			s <- s + sum(R * R)
		}
		s
	}

	# --- Symmetric operator: monotone line-search descent on f(A) ----------
	if (symmetric) {
		A_est <- A_init
		f_cur <- obj_fn(A_est, A_est)
		converged <- FALSE
		solve_log <- list()
		alpha_grid <- c(0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1, 1.3, 1.7, 2.5)

		for (iter in seq_len(max_iter)) {

			# direction 1: linearised ALS step. Hold the right factor at
			# A_est, solve the left factor row by row, then symmetrise.
			A_lin <- A_est
			iter_log <- list()
			for (i in seq_len(n_row)) {
				X_lhs <- vector("list", Tt - 1L)
				y_rhs <- numeric(0)
				for (t in 2:Tt) {
					X_lhs[[t - 1L]] <- A_est %*% t(Z_prev_list[[t]])
					y_rhs <- c(y_rhs, Z_list[[t]][i, ])
				}
				X_design <- do.call(rbind, X_lhs)
				if (nrow(X_design) > 0 && ncol(X_design) > 0) {
					sol <- als_solve(X_design, y_rhs)
					A_lin[i, ] <- sol$beta
					iter_log[[length(iter_log) + 1L]] <- sol
				}
			}
			A_lin <- 0.5 * (A_lin + t(A_lin))
			d_lin <- A_lin - A_est

			# direction 2: 0.5 x the symmetric-projected negative gradient of
			# f, a descent direction. grad f(A) = -2 sum_t (R_t A Z_{t-1}^T +
			# R_t^T A Z_{t-1}).
			G <- matrix(0, n_row, n_row)
			for (t in 2:Tt) {
				Zt1 <- Z_prev_list[[t]]
				R <- Z_list[[t]] - A_est %*% Zt1 %*% t(A_est)
				G <- G + R %*% A_est %*% t(Zt1) + t(R) %*% A_est %*% Zt1
			}
			d_grad <- 0.5 * (G + t(G))
			# rescale the gradient direction to the operator's own scale so
			# the alpha grid spans a sensible range of step lengths
			nA <- norm(A_est, "F")
			ng <- norm(d_grad, "F")
			if (ng > 0 && is.finite(ng)) {
				d_grad <- d_grad * (max(nA, 1e-6) / ng)
			}

			# line search: accept the (direction, step) with the lowest
			# objective; keep the current iterate if nothing improves.
			best_f <- f_cur
			best_A <- A_est
			for (d in list(d_lin, d_grad)) {
				if (!all(is.finite(d)) || sum(d * d) <= 0) next
				for (a in alpha_grid) {
					A_try <- A_est + a * d
					f_try <- obj_fn(A_try, A_try)
					if (is.finite(f_try) && f_try < best_f) {
						best_f <- f_try
						best_A <- A_try
					}
				}
			}

			rel <- (f_cur - best_f) / (f_cur + 1e-12)
			A_est <- best_A
			f_cur <- best_f
			solve_log <- iter_log

			if (verbose && (iter %% 5 == 0 || iter == 1)) {
				cli::cli_text("Iter {iter}: objective = {sprintf('%.6g', f_cur)}, rel = {sprintf('%.2e', rel)}")
			}
			if (rel < tol) {
				converged <- TRUE
				if (verbose) cli::cli_text("Converged at iteration {iter}.")
				break
			}
		}

		A_est <- 0.5 * (A_est + t(A_est))
		return(list(A = A_est, B = A_est, converged = converged,
		            solve_log = solve_log, objective = f_cur,
		            iterations = iter))
	}

	# --- Asymmetric operator: extrapolated penalized alternating LS -------
	# every sub-problem solve and the line search use one shared, fixed ridge
	# penalty lambda, so the whole iteration minimises a single coherent
	# penalized objective P(A, B) = g(A, B) + lambda (||A||^2 + ||B||^2). That
	# objective has a unique, bounded minimizer, which removes the operator
	# drift; the penalty also selects the balanced minimum-norm representative
	# of the (A, B) -> (cA, c^-1 B) gauge, so no separate balancing step is
	# needed. Each ALS sweep is followed by a line search that extrapolates
	# the ALS direction to accelerate the otherwise slow valley descent.
	D0 <- mean(vapply(2:Tt, function(t) sum(Z_prev_list[[t]]^2), numeric(1)))
	if (!is.finite(D0) || D0 <= 0) D0 <- 1
	lambda <- ridge_frac * D0
	pen_fn <- function(A, B) obj_fn(A, B) + lambda * (sum(A^2) + sum(B^2))

	A_est <- A_init
	B_est <- B_init
	p_cur <- pen_fn(A_est, B_est)
	best_p <- p_cur
	best_A <- A_est
	best_B <- B_est
	converged <- FALSE
	solve_log <- list()
	alpha_grid <- c(0.5, 1, 1.3, 1.7, 2.2, 3, 4)

	for (iter in seq_len(max_iter)) {
		iter_log <- list()
		A_prev <- A_est
		B_prev <- B_est

		# A-step: fix B, solve each penalized row of A. Row i design is
		# B Z_{t-1}^T; the shared lambda penalizes ||A||^2.
		A_als <- A_est
		for (i in seq_len(n_row)) {
			X_lhs <- vector("list", Tt - 1L)
			y_rhs <- numeric(0)
			for (t in 2:Tt) {
				X_lhs[[t - 1L]] <- B_est %*% t(Z_prev_list[[t]])
				y_rhs <- c(y_rhs, Z_list[[t]][i, ])
			}
			X_design <- do.call(rbind, X_lhs)
			if (nrow(X_design) > 0 && ncol(X_design) > 0) {
				sol <- als_solve_ridged(X_design, y_rhs, lambda)
				A_als[i, ] <- sol$beta
				iter_log[[length(iter_log) + 1L]] <- sol
			}
		}

		# B-step: fix A at the A-step value, solve each penalized row of B.
		# column j of Z_t depends on row j of B (Z = A Z_{t-1} B^T).
		B_als <- B_est
		for (j in seq_len(n_col)) {
			X_lhs <- vector("list", Tt - 1L)
			y_rhs <- numeric(0)
			for (t in 2:Tt) {
				X_lhs[[t - 1L]] <- A_als %*% Z_prev_list[[t]]
				y_rhs <- c(y_rhs, Z_list[[t]][, j])
			}
			X_design <- do.call(rbind, X_lhs)
			if (nrow(X_design) > 0 && ncol(X_design) > 0) {
				sol <- als_solve_ridged(X_design, y_rhs, lambda)
				B_als[j, ] <- sol$beta
				iter_log[[length(iter_log) + 1L]] <- sol
			}
		}

		# line search: extrapolate the ALS direction against the penalized
		# objective. alpha = 1 is the plain penalized ALS step, which already
		# decreases P, so a monotone decrease is guaranteed; alpha > 1 follows
		# the valley faster.
		d_A <- A_als - A_prev
		d_B <- B_als - B_prev
		best_step_p <- Inf
		A_est <- A_prev
		B_est <- B_prev
		for (a in alpha_grid) {
			A_try <- A_prev + a * d_A
			B_try <- B_prev + a * d_B
			p_try <- pen_fn(A_try, B_try)
			if (is.finite(p_try) && p_try < best_step_p) {
				best_step_p <- p_try
				A_est <- A_try
				B_est <- B_try
			}
		}

		p_new <- pen_fn(A_est, B_est)
		rel <- abs(p_cur - p_new) / (p_cur + 1e-12)
		p_cur <- p_new
		solve_log <- iter_log
		if (is.finite(p_new) && p_new < best_p) {
			best_p <- p_new
			best_A <- A_est
			best_B <- B_est
		}

		if (verbose && (iter %% 5 == 0 || iter == 1)) {
			cli::cli_text("Iter {iter}: objective = {sprintf('%.6g', obj_fn(A_est, B_est))}, rel = {sprintf('%.2e', rel)}")
		}
		if (rel < tol) {
			converged <- TRUE
			if (verbose) cli::cli_text("Converged at iteration {iter}.")
			break
		}
	}

	list(A = best_A, B = best_B, converged = converged,
	     solve_log = solve_log, objective = obj_fn(best_A, best_B),
	     iterations = iter)
}

####
#' Internal ALS Refitting with Warm Start
#'
#' @description
#' ALS fitting with optional warm initialization from a prior fit.
#' Used by bootstrap to stabilize convergence on resampled data.
#' Implements the iterative block coordinate descent estimator of Minhas and
#' Hoff (2025) with warm-start initialization.
#'
#' @param Z Latent scores matrix (n_row x n_col)
#' @param Z_prev List of previous Z matrices for AR structure
#' @param A_init Optional initial A matrix (n_row x n_row)
#' @param B_init Optional initial B matrix (n_col x n_col)
#' @param symmetric Logical, enforce B = A
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param ridge_frac Asymmetric ridge penalty fraction; must match the original
#'   fit so the bootstrap replicates the same estimator.
#'
#' @return List with fitted (A, B) matrices
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}.
#' @keywords internal
als_refit_warmstart <- function(Z_list, Z_prev_list, A_init, B_init,
                                 symmetric, max_iter, tol, n_row, n_col, Tt,
                                 ridge_frac = 3e-3) {
	# warm-start initialisations, falling back to the identity operator
	if (is.null(A_init)) A_init <- diag(n_row)
	if (is.null(B_init)) B_init <- if (symmetric) A_init else diag(n_col)
	if (symmetric) B_init <- A_init

	# delegate to the shared, scale-hardened ALS core
	fit <- als_core_fit(Z_list, Z_prev_list, Tt, n_row, n_col, symmetric,
	                    A_init = A_init, B_init = B_init,
	                    max_iter = max_iter, tol = tol,
	                    ridge_frac = ridge_frac, verbose = FALSE)

	list(A = fit$A, B = fit$B, converged = fit$converged)
}

####
#' Internal helper: bootstrap CIs for an ALS fit
#'
#' Two bootstrap strategies:
#'  - Block bootstrap: resample time periods with replacement. Preserves
#'    within-period dependence. Recommended when Tt >= 10.
#'  - Parametric bootstrap: simulate new outcomes from the fitted model.
#'    Better when Tt is small and the model is well-specified.
#'
#' Each replicate refits the full ALS model (warm-started from the point
#' estimate). Replicates that fail to converge are dropped and reported.
#' Standard errors are column standard deviations of the successful
#' replicates; CIs use the percentile method.
#'
#' Follows the resampling-based inference scheme of the Social Influence
#' Regression model of Minhas and Hoff (2025).
#'
#' @param als_fit A dbn object from \code{\link{dbn_als}()} (i.e. with
#'   \code{meta$sampler_used == "als"}).
#' @param R Integer number of bootstrap replicates.
#' @param type \code{"block"} or \code{"parametric"}.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical.
#' @return A list with class \code{dbn_boot} carrying bootstrap coefficients,
#'   standard errors, percentile CIs, and metadata.
#' @keywords internal
#' @noRd
NULL

#' Attach a bootstrap result to a fit, downgrading uncertainty flags when empty
#'
#' If the bootstrap returned `n_valid < 2`, every CI cell is NA and the
#' downstream `predict()` / `dyad_path()` would silently propagate the NAs
#' while reporting `meta$uncertainty_available = TRUE`. Flip the flag off
#' and surface a loud warning so the user knows their CIs are unusable.
#' @keywords internal
#' @noRd
.dbn_attach_bootstrap <- function(out, boot, label, source) {
	n_valid <- boot$n_valid %||% 0L
	if (is.finite(n_valid) && n_valid >= 2L) {
		out$bootstrap <- boot
		out$meta$uncertainty_available <- TRUE
		out$meta$uncertainty_source <- source
		out$meta$approximation_note <- label
	} else {
		cli::cli_warn(c(
			"Bootstrap returned {n_valid} valid replicate{?s}; not enough to compute intervals.",
			"i" = "Marking the fit as point-estimate-only ({.code meta$uncertainty_available = FALSE}); the empty bootstrap slot is dropped.",
			"i" = "Check the per-replicate error output above for the underlying refit failure."
		))
		out$bootstrap <- NULL
		out$meta$uncertainty_available <- FALSE
		out$meta$uncertainty_source <- NA_character_
		out$meta$approximation_note <- sprintf(
			"Bootstrap requested but no valid replicates returned (n_valid = %d).", n_valid)
	}
	out
}

.dbn_als_bootstrap <- function(als_fit, R = 200,
                               type = c("block", "parametric"),
                               seed = NULL, verbose = TRUE,
                               parallel = TRUE) {

	if (!inherits(als_fit, "dbn")) {
		cli::cli_abort("{.arg als_fit} must be a dbn object from {.fun dbn_als}.")
	}
	if (!identical(als_fit$meta$sampler_used, "als")) {
		cli::cli_abort(c(
			"{.arg als_fit} must be a point-estimate fit from {.fun dbn_als}.",
			"x" = "Got a fit from {.fun dbn} (model {.val {als_fit$model}}).",
			"i" = "The bootstrap resamples an ALS fit; for MCMC uncertainty use {.fun dbn} directly."
		))
	}
	type <- match.arg(type)
	if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R < 2 || R != round(R)) {
		cli::cli_abort(c(
			"{.arg R} must be a single whole number of at least 2.",
			"i" = "Bootstrap standard errors and intervals need several replicates; 200 or more is recommended."
		))
	}
	R <- as.integer(R)
	# the top-level `dbn_als(bootstrap = N)` already fired the small-N warning
	# at the user's call site; firing it again here would surface the same
	# message twice (once tagged `bootstrap`, once tagged `R`). Skip when the
	# outer caller has already fired by setting dbn.in_bootstrap.
	if (R < 30 && !isTRUE(getOption("dbn.in_bootstrap", FALSE))) {
		cli::cli_warn(c(
			"{.arg R} = {R} is small for bootstrap inference.",
			"i" = "Standard errors and intervals from fewer than ~30 replicates are unstable; 200 or more is recommended."
		))
	}

	# suppress per-rep diagonal warnings inside the bootstrap loop -- the
	# top-level call already fired once if applicable.
	.dbn_prev_in_boot <- getOption("dbn.in_bootstrap", FALSE)
	options(dbn.in_bootstrap = TRUE)
	on.exit(options(dbn.in_bootstrap = .dbn_prev_in_boot), add = TRUE)

	# `type` already match.arg'd by the caller; do not double-match
	.dbn_restore_seed <- .use_seed_locally(seed)
	on.exit(.dbn_restore_seed(), add = TRUE)

	# note: parametric bootstrap not recommended for asymmetric models due to scale ambiguity
	if (type == "parametric" && !als_fit$dims$is_symmetric && verbose) {
		cli::cli_warn("Parametric bootstrap for asymmetric models: scale ambiguity may cause inflated SE for B. Consider block bootstrap instead.")
	}

	Y <- als_fit$Y
	family <- als_fit$family
	symmetric <- als_fit$dims$is_symmetric
	n_row <- als_fit$dims$n_row
	n_col <- als_fit$dims$n_col
	p <- als_fit$dims$p
	Tt <- als_fit$dims$Tt

	# original point estimates
	A_orig <- als_fit$A[[1]][, , 1]  # n_row x n_row
	B_orig <- als_fit$B[[1]][, , 1]  # n_col x n_col
	M_orig <- als_fit$M[, , , 1]     # n_row x n_col x p

	# ridge penalty of the original fit, so each bootstrap refit replicates
	# the same estimator (older fits without a stored value use the default)
	ridge_boot <- als_fit$settings$ridge
	if (is.null(ridge_boot) || !is.finite(ridge_boot)) ridge_boot <- 3e-3

	n_params_A <- n_row * n_row
	n_params_B <- n_col * n_col
	n_params_M <- n_row * n_col * p

	boot_coefs_A <- matrix(NA, R, n_params_A)
	boot_coefs_B <- matrix(NA, R, n_params_B)
	boot_coefs_M <- matrix(NA, R, n_params_M)
	# capture per-rep failures so they don't disappear into NA rows.
	boot_rep_errors <- character(R)

	if (verbose) {
		cli::cli_h3("Bootstrap for ALS Fit")
		cli::cli_text("Type: {.field {type}} | Replicates: {.val {R}}")
		cli::cli_text("Network: {n_row} x {n_col}, {p} relation{?s}, {Tt} time point{?s}")
	}

	# per-rep closure: returns a list with {A, B, M, err}. extracted from the
	# original for-loop so we can route through lapply OR future.apply.
	do_static_rep <- function(b) {
		if (verbose && b %% 10 == 0) {
			cli::cli_inform("Bootstrap {.val {b}}/{.val {R}}")
		}
		err_msg <- ""

		if (type == "block") {
			# block bootstrap: resample time indices with replacement
			t_idx <- sample(1:Tt, Tt, replace = TRUE)
			Y_b <- Y[, , , t_idx, drop = FALSE]
		} else {
			# parametric bootstrap: simulate from fitted model
			# first, extract Z for prediction
			pre <- shared_preprocess(Y, family = family)
			Z <- pre$Z
			M <- pre$M

			# compute predicted latent state for all times
			Theta_pred <- array(NA, dim = c(n_row, n_col, p, Tt))
			Z_list <- list()
			for (t in 1:Tt) {
				Z_t_all <- Z[, , , t]
				Z_list[[t]] <- apply(Z_t_all, c(1, 2), mean)
			}
			for (t in 1:Tt) {
				Theta_t <- A_orig %*% Z_list[[t]] %*% t(B_orig)
				for (r in 1:p) {
					Theta_pred[, , r, t] <- Theta_t
				}
			}

			# simulate new Y from fitted model
			if (family == "gaussian") {
				sigma <- sqrt(als_fit$sigma2[1])
				Y_b <- array(rnorm(length(Theta_pred), mean = Theta_pred, sd = sigma),
							 dim = dim(Y))
			} else if (family == "ordinal") {
				# ordinal: simulate latent Z, then discretize to the observed
				# category range (do not assume a fixed 1-5 scale)
				ord_rng <- range(Y, na.rm = TRUE)
				Y_b <- array(NA, dim = dim(Y))
				for (t in 1:Tt) {
					for (r in 1:p) {
						Z_sim <- rnorm(n_row * n_col, mean = c(Theta_pred[, , r, t]), sd = sqrt(als_fit$sigma2[1]))
						Y_b[, , r, t] <- matrix(pmin(pmax(round(Z_sim), ord_rng[1]), ord_rng[2]), n_row, n_col)
					}
				}
			} else if (family == "binary") {
				# binary: Bernoulli from the probit link, matching the model's
				# probit-with-data-augmentation likelihood (not a logit link)
				Y_b <- array(NA, dim = dim(Y))
				for (t in 1:Tt) {
					for (r in 1:p) {
						prob <- pnorm(Theta_pred[, , r, t])
						Y_b[, , r, t] <- matrix(rbinom(n_row * n_col, 1, prob), n_row, n_col)
					}
				}
			}
			# preserve diagonal NA structure
			for (tt in 1:Tt) {
				for (rr in 1:p) {
					diag(Y_b[, , rr, tt]) <- NA
				}
			}
		}

		# refit ALS on bootstrap sample using warm-start from original solution
		# this prevents converging to different modes/local optima
		out <- tryCatch({
			# preprocess bootstrap data
			pre_b <- shared_preprocess(Y_b, family = family)
			Z_b <- pre_b$Z
			M_b <- pre_b$M

			# extract latent Z
			Z_list_b <- list()
			Z_prev_list_b <- list()
			for (t in 1:Tt) {
				Z_t_all <- Z_b[, , , t]
				Z_list_b[[t]] <- apply(Z_t_all, c(1, 2), mean)
				if (t > 1) {
					Z_prev_list_b[[t]] <- Z_list_b[[t-1]]
				}
			}

			# warm-start ALS refinement from original estimates
			als_result <- als_refit_warmstart(
				Z_list_b, Z_prev_list_b,
				A_init = A_orig, B_init = B_orig,  # Use original as starting point
				symmetric = symmetric,
				max_iter = 30,  # Fewer iterations needed with warm start
				tol = 1e-5,
				n_row = n_row, n_col = n_col, Tt = Tt,
				ridge_frac = ridge_boot  # same penalty as the original fit
			)

			list(A = als_result$A, B = als_result$B, M = M_b, err = "")
		}, error = function(e) {
			list(A = NULL, B = NULL, M = NULL, err = conditionMessage(e))
		})
		out
	}

	# decide between sequential and future_lapply. parallelise only when a
	# non-sequential plan is active; future.packages = "dbn" ensures workers
	# have the package loaded for the internal helpers.
	use_parallel_static <- isTRUE(parallel) &&
		requireNamespace("future.apply", quietly = TRUE) &&
		requireNamespace("future", quietly = TRUE)
	if (use_parallel_static) {
		plan_cls <- class(future::plan())
		if ("sequential" %in% plan_cls || "uniprocess" %in% plan_cls) {
			use_parallel_static <- FALSE
		}
	}
	# devtools::load_all caveat (see als_tv_bootstrap.R for context). Workers
	# would library() the installed dbn from .libPaths() rather than the dev
	# tree; refuse parallel and fall back to sequential to keep behaviour
	# consistent with the code being audited.
	if (use_parallel_static &&
	    requireNamespace("pkgload", quietly = TRUE) &&
	    pkgload::is_dev_package("dbn")) {
		cli::cli_warn(c(
			"Detected {.fn devtools::load_all} state for {.pkg dbn}.",
			"i" = "Parallel bootstrap workers would load the INSTALLED {.pkg dbn} from {.path .libPaths()}, not the dev tree.",
			"i" = "Falling back to sequential. For real parallel speedup, install the package first."
		))
		use_parallel_static <- FALSE
	}
	reps_out <- if (use_parallel_static) {
		future.apply::future_lapply(seq_len(R), do_static_rep,
			future.seed = if (!is.null(seed)) seed else TRUE,
			future.packages = "dbn")
	} else {
		lapply(seq_len(R), do_static_rep)
	}
	for (b in seq_len(R)) {
		ro <- reps_out[[b]]
		if (nzchar(ro$err)) {
			boot_rep_errors[b] <- ro$err
			next
		}
		boot_coefs_A[b, ] <- c(ro$A)
		boot_coefs_B[b, ] <- c(ro$B)
		boot_coefs_M[b, ] <- c(ro$M)
	}

	# surface per-rep errors loudly so a low n_valid is explained by the
	# top unique error messages rather than left to the caller to diagnose.
	failed_mask <- nzchar(boot_rep_errors)
	if (any(failed_mask)) {
		err_table <- sort(table(boot_rep_errors[failed_mask]), decreasing = TRUE)
		n_shown <- min(3L, length(err_table))
		cli::cli_warn(c(
			"{sum(failed_mask)} of {R} static-bootstrap replicates failed to refit.",
			"i" = "Top unique error{?s} (showing {n_shown}):",
			setNames(paste0("(", as.integer(err_table[1:n_shown]), "x) ",
				substr(names(err_table[1:n_shown]), 1, 80)),
				rep("x", n_shown))
		))
	}

	# compute statistics from successful replicates
	valid_A <- apply(boot_coefs_A, 1, function(r) !any(is.na(r)))
	valid_B <- apply(boot_coefs_B, 1, function(r) !any(is.na(r)))
	valid_M <- apply(boot_coefs_M, 1, function(r) !any(is.na(r)))

	# use intersection of all valid replicates
	valid <- valid_A & valid_B & valid_M
	n_valid <- sum(valid)

	# quality control: flag replicates that are outliers (likely different local optima)
	# by checking Frobenius norm of A relative to others
	if (n_valid > 3) {
		A_norms <- sqrt(rowSums(boot_coefs_A[valid, ]^2))
		A_norm_median <- median(A_norms)
		A_norm_mad <- mad(A_norms)
		# flag replicates with very different norms (>3 MAD from median)
		outlier_threshold <- A_norm_median + 3 * A_norm_mad
		outliers <- which(A_norms > outlier_threshold)

		if (length(outliers) > 0 && verbose) {
			n_outliers <- length(outliers)
			cli::cli_warn("Bootstrap: {n_outliers} replicate{?s} with large coefficient norms (possible local optima). Consider larger R or MCMC.")
		}
	}

	if (n_valid < 10 && verbose) {
		cli::cli_warn("Only {.val {n_valid}} valid bootstrap replicates (of {.val {R}}). Results unreliable.")
	}

	# compute standard errors and CIs
	boot_se_A <- apply(boot_coefs_A[valid, , drop = FALSE], 2, sd)
	boot_se_B <- apply(boot_coefs_B[valid, , drop = FALSE], 2, sd)
	boot_se_M <- apply(boot_coefs_M[valid, , drop = FALSE], 2, sd)

	boot_ci_A <- apply(boot_coefs_A[valid, , drop = FALSE], 2,
					   quantile, probs = c(0.025, 0.975))
	boot_ci_B <- apply(boot_coefs_B[valid, , drop = FALSE], 2,
					   quantile, probs = c(0.025, 0.975))
	boot_ci_M <- apply(boot_coefs_M[valid, , drop = FALSE], 2,
					   quantile, probs = c(0.025, 0.975))

	# basic (reverse-percentile) intervals: [2 theta_hat - q_{1-a/2},
	# 2 theta_hat - q_{a/2}]. better coverage than percentile when the
	# bootstrap distribution is skewed or the estimator is biased
	A_orig_v <- as.numeric(A_orig)
	B_orig_v <- as.numeric(B_orig)
	M_orig_v <- as.numeric(M_orig)
	basic_A_lo <- 2 * A_orig_v - boot_ci_A[2, ]
	basic_A_hi <- 2 * A_orig_v - boot_ci_A[1, ]
	basic_B_lo <- 2 * B_orig_v - boot_ci_B[2, ]
	basic_B_hi <- 2 * B_orig_v - boot_ci_B[1, ]
	basic_M_lo <- 2 * M_orig_v - boot_ci_M[2, ]
	basic_M_hi <- 2 * M_orig_v - boot_ci_M[1, ]

	# reshape flat bootstrap-CI vectors into the matrix shape that matches
	# the operator: static is [n_row, n_row], TV is [n_row, n_row, T-1]
	.reshape_A <- function(v) { d <- v; dim(d) <- c(n_row, n_row); d }
	.reshape_B <- function(v) { d <- v; dim(d) <- c(n_col, n_col); d }
	.reshape_M <- function(v) { d <- v; dim(d) <- c(n_row, n_col, p); d }

	result <- list(
		coefs_A = boot_coefs_A,
		coefs_B = boot_coefs_B,
		coefs_M = boot_coefs_M,
		se_A = .reshape_A(boot_se_A),
		se_B = .reshape_B(boot_se_B),
		se_M = .reshape_M(boot_se_M),
		ci_A_lo = .reshape_A(boot_ci_A[1, ]),
		ci_A_hi = .reshape_A(boot_ci_A[2, ]),
		ci_B_lo = .reshape_B(boot_ci_B[1, ]),
		ci_B_hi = .reshape_B(boot_ci_B[2, ]),
		ci_M_lo = .reshape_M(boot_ci_M[1, ]),
		ci_M_hi = .reshape_M(boot_ci_M[2, ]),
		basic_A_lo = .reshape_A(basic_A_lo),
		basic_A_hi = .reshape_A(basic_A_hi),
		basic_B_lo = .reshape_B(basic_B_lo),
		basic_B_hi = .reshape_B(basic_B_hi),
		basic_M_lo = .reshape_M(basic_M_lo),
		basic_M_hi = .reshape_M(basic_M_hi),
		point_est_A = .reshape_A(c(A_orig)),
		point_est_B = .reshape_B(c(B_orig)),
		point_est_M = .reshape_M(c(M_orig)),
		n_valid = n_valid,
		n_total = R,
		type = type,
		family = family,
		symmetric = symmetric,
		dims = list(n_row = n_row, n_col = n_col, p = p, Tt = Tt)
	)

	class(result) <- "dbn_boot"
	result
}

####
#' Print Bootstrap Results for ALS
#'
#' @param x A `dbn_boot` object from dbn_als(bootstrap = R).
#' @param digits Number of significant digits.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the `dbn_boot` object.
#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cli::cli_text("\n")
	cli::cli_h3("Bootstrap Results for ALS Fit")
	cli::cli_text("Type: {.field {x$type}} | Replicates: {.val {x$n_valid}}/{.val {x$n_total}} valid")
	cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$dims$p} relation{?s}, {x$dims$Tt} time point{?s}")
	cli::cli_text("Family: {.field {x$family}} {ifelse(x$symmetric, '(symmetric)', '')}")

	cli::cli_rule()

	# summarize A matrix
	cli::cli_text("{.strong Operator A ({x$dims$n_row} x {x$dims$n_row})}")
	A_mean_se <- mean(x$se_A)
	A_mean_est <- mean(x$point_est_A)
	cli::cli_text("  Mean estimate: {sprintf('%.4f', A_mean_est)}, Mean SE: {sprintf('%.4f', A_mean_se)}")
	cli::cli_text("  Diagonal mean: {sprintf('%.4f', mean(matrix(x$point_est_A, x$dims$n_row, x$dims$n_row)[cbind(1:x$dims$n_row, 1:x$dims$n_row)]))}")

	# summarize B matrix
	cli::cli_text("{.strong Operator B ({x$dims$n_col} x {x$dims$n_col})}")
	B_mean_se <- mean(x$se_B)
	B_mean_est <- mean(x$point_est_B)
	cli::cli_text("  Mean estimate: {sprintf('%.4f', B_mean_est)}, Mean SE: {sprintf('%.4f', B_mean_se)}")
	cli::cli_text("  Diagonal mean: {sprintf('%.4f', mean(matrix(x$point_est_B, x$dims$n_col, x$dims$n_col)[cbind(1:x$dims$n_col, 1:x$dims$n_col)]))}")

	# summarize M baseline
	cli::cli_text("{.strong Baseline M ({x$dims$n_row} x {x$dims$n_col} x {x$dims$p})}")
	M_mean_se <- mean(x$se_M)
	M_mean_est <- mean(x$point_est_M)
	cli::cli_text("  Mean estimate: {sprintf('%.4f', M_mean_est)}, Mean SE: {sprintf('%.4f', M_mean_se)}")

	cli::cli_text("")
	invisible(x)
}

####
#' Summary of Bootstrap Results for ALS
#'
#' @param object A `dbn_boot` object from dbn_als(bootstrap = R).
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the `dbn_boot` object.
#' @author Tosin Salau and Shahryar Minhas
#' @exportS3Method summary dbn_boot
#' @name summary.dbn_boot
NULL

####
#' Side-by-side percentile and basic bootstrap intervals
#'
#' @description Returns both interval types so users can compare their
#'   widths. Percentile intervals \eqn{[q_{\alpha/2}, q_{1-\alpha/2}]}
#'   are easy to interpret but undercover when the bootstrap
#'   distribution is skewed. Basic (reverse-percentile) intervals
#'   \eqn{[2 \hat\theta - q_{1-\alpha/2},\, 2 \hat\theta - q_{\alpha/2}]}
#'   reduce that skew bias by reflecting the quantiles about the point
#'   estimate.
#' @param object A \code{dbn_boot} bootstrap object (from
#'   \code{dbn_als(..., bootstrap = N)}) OR a \code{dbn} fit whose
#'   \code{$bootstrap} component is a \code{dbn_boot}.
#' @param parm Which operator to return: \code{"A"} (default),
#'   \code{"B"}, or \code{"M"}.
#' @param method One of \code{"both"} (default), \code{"percentile"},
#'   \code{"basic"}.
#' @param level Coverage level. Default \code{0.95} returns the stored
#'   95\% intervals; any other value triggers re-extraction from the
#'   stored \code{coefs_*} draws at the requested level. Requires the
#'   bootstrap to have stored raw coefficient draws (always true for a
#'   default \code{dbn_als(..., bootstrap = N)} fit).
#' @return A list with named \code{percentile} and/or \code{basic}
#'   entries, each a list with \code{lo} and \code{hi} arrays shaped to
#'   match the operator.
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_boot_intervals <- function(object,
                                 parm = c("A", "B", "M"),
                                 method = c("both", "percentile", "basic"),
                                 level = 0.95) {
	parm <- match.arg(parm)
	method <- match.arg(method)
	boot <- if (inherits(object, "dbn_boot")) object else object$bootstrap
	if (is.null(boot) || !inherits(boot, "dbn_boot"))
		cli::cli_abort(c(
			"{.arg object} must be a {.cls dbn_boot} or a {.cls dbn} fit with a stored bootstrap.",
			"i" = "Run {.code dbn_als(data, bootstrap = N)} to populate {.code fit$bootstrap}."
		))
	if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
	    level <= 0 || level >= 1)
		cli::cli_abort("{.arg level} must be a single value in (0, 1).")
	# re-extract at the requested level if it differs from the stored 0.95
	if (abs(level - 0.95) > 1e-6) {
		coefs <- boot[[paste0("coefs_", parm)]]
		if (is.null(coefs))
			cli::cli_abort(c(
				"Cannot re-extract at level {.val {level}}: raw {.code coefs_{parm}} draws not stored.",
				"i" = "Refit with the current package version which retains draws by default."
			))
		valid <- apply(coefs, 1, function(r) !any(is.na(r)))
		alpha <- (1 - level) / 2
		probs <- c(alpha, 1 - alpha)
		qmat <- apply(coefs[valid, , drop = FALSE], 2,
			stats::quantile, probs = probs, na.rm = TRUE)
		point <- as.numeric(boot[[paste0("point_est_", parm)]])
		# reshape flat vectors back to operator shape
		reshape_arr <- function(v) {
			arr <- v
			dim(arr) <- dim(boot[[paste0("point_est_", parm)]])
			arr
		}
		pct_lo <- reshape_arr(qmat[1, ])
		pct_hi <- reshape_arr(qmat[2, ])
		basic_lo <- reshape_arr(2 * point - qmat[2, ])
		basic_hi <- reshape_arr(2 * point - qmat[1, ])
	} else {
		pct_lo <- boot[[paste0("ci_", parm, "_lo")]]
		pct_hi <- boot[[paste0("ci_", parm, "_hi")]]
		basic_lo <- boot[[paste0("basic_", parm, "_lo")]]
		basic_hi <- boot[[paste0("basic_", parm, "_hi")]]
	}
	out <- list()
	if (method %in% c("both", "percentile"))
		out$percentile <- list(lo = pct_lo, hi = pct_hi, level = level)
	if (method %in% c("both", "basic"))
		out$basic <- list(lo = basic_lo, hi = basic_hi, level = level)
	out
}

#' @rdname summary.dbn_boot
#' @author Tosin Salau and Shahryar Minhas
#' @export
summary.dbn_boot <- function(object, ...) {
	cli::cli_h1("Bootstrap Results for ALS Fit")
	cli::cli_text("Type: {.field {object$type}} | Family: {.field {object$family}}")
	cli::cli_text("Replicates: {.val {object$n_valid}} valid of {.val {object$n_total}} total ({.val {round(100 * object$n_valid / object$n_total, 1)}}%)")
	cli::cli_text("Network: {object$dims$n_row} x {object$dims$n_col}, {object$dims$p} relation{?s}, {object$dims$Tt} time point{?s}")

	cli::cli_rule()

	# A matrix summary
	cli::cli_text("{.strong Operator A}")
	A_mat <- matrix(object$point_est_A, object$dims$n_row, object$dims$n_row)
	A_se_mat <- matrix(object$se_A, object$dims$n_row, object$dims$n_row)
	cli::cli_text("Point estimate (diagonal): {paste(sprintf('%.4f', diag(A_mat)), collapse=', ')}")
	cli::cli_text("Bootstrap SE (diagonal): {paste(sprintf('%.4f', diag(A_se_mat)), collapse=', ')}")

	# B matrix summary
	cli::cli_text("{.strong Operator B}")
	B_mat <- matrix(object$point_est_B, object$dims$n_col, object$dims$n_col)
	B_se_mat <- matrix(object$se_B, object$dims$n_col, object$dims$n_col)
	cli::cli_text("Point estimate (diagonal): {paste(sprintf('%.4f', diag(B_mat)), collapse=', ')}")
	cli::cli_text("Bootstrap SE (diagonal): {paste(sprintf('%.4f', diag(B_se_mat)), collapse=', ')}")

	# coverage check: what percentage of params have CIs that include the point est
	A_covers <- (object$ci_A_lo <= object$point_est_A) & (object$point_est_A <= object$ci_A_hi)
	B_covers <- (object$ci_B_lo <= object$point_est_B) & (object$point_est_B <= object$ci_B_hi)
	M_covers <- (object$ci_M_lo <= object$point_est_M) & (object$point_est_M <= object$ci_M_hi)

	cli::cli_rule()
	cli::cli_text("{.strong CI Coverage}")
	cli::cli_text("A: {sum(A_covers)}/{length(A_covers)} params have point est in 95% CI")
	cli::cli_text("B: {sum(B_covers)}/{length(B_covers)} params have point est in 95% CI")
	cli::cli_text("M: {sum(M_covers)}/{length(M_covers)} params have point est in 95% CI")

	invisible(object)
}

#' Attach actor / relation / horizon dimnames to a forecast array
#'
#' Sets dimnames on the array returned by `simulate_dynamic` /
#' `simulate_piecewise` so downstream code (and the user) can read
#' off "actor i sending to actor j on horizon h" without juggling
#' indices. Idempotent: if any axis already has names it is left alone.
#'
#' @param fc Numeric array. Either 4D `[n_row, n_col, p, H]` (when
#'   the simulator collapsed draws via `summary = "mean"`) or 5D
#'   `[n_row, n_col, p, H, draws]`.
#' @param object The `dbn` fit (for actor labels).
#' @param H The forecast horizon (integer).
#' @return The same array with dimnames set on every unnamed axis.
#' @keywords internal
#' @noRd
.dbn_decorate_forecast <- function(fc, object, H) {
	d <- dim(fc)
	if (is.null(d) || length(d) < 4L) return(fc)
	dn <- dimnames(fc)
	if (is.null(dn)) dn <- vector("list", length(d))
	# axis 1: senders (n_row)
	if (is.null(dn[[1]])) {
		nm <- if (!is.null(object$Y)) dimnames(object$Y)[[1]] else NULL
		if (is.null(nm) || length(nm) != d[1]) {
			nm <- paste0("a", seq_len(d[1]))
		}
		dn[[1]] <- nm
	}
	# axis 2: receivers (n_col)
	if (is.null(dn[[2]])) {
		nm <- if (!is.null(object$Y)) dimnames(object$Y)[[2]] else NULL
		if (is.null(nm) || length(nm) != d[2]) {
			nm <- paste0("b", seq_len(d[2]))
		}
		dn[[2]] <- nm
	}
	# axis 3: relations
	if (is.null(dn[[3]])) {
		nm <- if (!is.null(object$Y) && length(dim(object$Y)) == 4L) dimnames(object$Y)[[3]] else NULL
		if (is.null(nm) || length(nm) != d[3]) {
			nm <- paste0("rel", seq_len(d[3]))
		}
		dn[[3]] <- nm
	}
	# axis 4: horizon
	if (is.null(dn[[4]])) {
		dn[[4]] <- paste0("h+", seq_len(d[4]))
	}
	# axis 5: draws (only when present)
	if (length(d) >= 5L && is.null(dn[[5]])) {
		dn[[5]] <- paste0("draw", seq_len(d[5]))
	}
	dimnames(fc) <- dn
	fc
}

#' Print a dbn forecast object with intervals
#' @param x A `dbn_forecast_ci` list of mean / lower / upper arrays.
#' @param ... Unused.
#' @return Invisibly, `x`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_forecast_ci <- function(x, ...) {
	H <- attr(x, "H") %||% NA_integer_
	mod <- attr(x, "model") %||% "?"
	fam <- attr(x, "family") %||% "?"
	prb <- attr(x, "probs") %||% c(0.025, 0.975)
	src <- attr(x, "source") %||% "mcmc"
	d <- dim(x$mean)
	int_label <- if (identical(src, "bootstrap")) "bootstrap intervals" else "credible intervals"
	cli::cli_h2("dbn forecast (with {int_label})")
	cli::cli_bullets(c(
		" " = "model: {.val {mod}}, family: {.val {fam}}",
		" " = "horizon H = {H}",
		" " = "shape: {paste(d, collapse = ' x ')}, {int_label} at {sprintf('%g%%-%g%%', 100 * prb[1], 100 * prb[2])}"
	))
	cli::cli_text("Components: {.var x$mean}, {.var x$lower}, {.var x$upper}.")
	invisible(x)
}

#' Summary for a dbn forecast object with credible intervals
#' @param object A `dbn_forecast_ci` list.
#' @param ... Unused.
#' @return Data frame with one row per horizon: mean / lower / upper aggregated across cells.
#' @author Tosin Salau and Shahryar Minhas
#' @export
summary.dbn_forecast_ci <- function(object, ...) {
	d <- dim(object$mean)
	h_axis <- length(d)  # horizon is last after collapsing draws
	per_h <- vapply(seq_len(d[h_axis]), function(h) {
		idx <- as.list(rep(TRUE, length(d)))
		idx[[h_axis]] <- h
		mn <- do.call(`[`, c(list(object$mean), idx, list(drop = FALSE)))
		lo <- do.call(`[`, c(list(object$lower), idx, list(drop = FALSE)))
		hi <- do.call(`[`, c(list(object$upper), idx, list(drop = FALSE)))
		c(mean  = mean(mn, na.rm = TRUE),
		  lower = mean(lo, na.rm = TRUE),
		  upper = mean(hi, na.rm = TRUE))
	}, numeric(3))
	out <- as.data.frame(t(per_h))
	out$h <- seq_len(d[h_axis])
	out[, c("h", "mean", "lower", "upper")]
}

#' Print a dbn forecast object
#' @param x A `dbn_forecast` array, returned by [predict.dbn()] / [forecast.dbn()] with `H` set.
#' @param ... Unused.
#' @return Invisibly, `x`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_forecast <- function(x, ...) {
	d <- dim(x)
	H <- attr(x, "H") %||% NA_integer_
	mod <- attr(x, "model") %||% "?"
	fam <- attr(x, "family") %||% "?"
	kind <- attr(x, "summary_kind") %||% "draws"
	cli::cli_h2("dbn forecast")
	cli::cli_bullets(c(
		" " = "model: {.val {mod}}, family: {.val {fam}}",
		" " = "horizon H = {H}, summary = {.val {kind}}",
		" " = "shape: {paste(d, collapse = ' x ')}"
	))
	cli::cli_text("Use {.code as.array(fc)} or subset with {.code fc[, , h]} to access values.")
	invisible(x)
}

#' Summary for a dbn forecast object
#' @param object A `dbn_forecast` array.
#' @param ... Unused.
#' @return A data frame summarising each horizon (mean / sd / range across draws and cells).
#' @author Tosin Salau and Shahryar Minhas
#' @export
summary.dbn_forecast <- function(object, ...) {
	d <- dim(object)
	H <- attr(object, "H") %||% (d[length(d)])
	# locate the horizon axis: it is the named "H" dim if present, otherwise
	# the penultimate dim for [n_row, n_col, p, H, draws] outputs from the
	# simulator path, or the last dim if `summary="mean"` collapsed draws.
	h_axis <- if (!is.na(H) && length(d) >= 1L) {
		match(H, d, nomatch = length(d))
	} else length(d)
	per_h <- vapply(seq_len(d[h_axis]), function(h) {
		idx <- as.list(rep(TRUE, length(d)))
		idx[[h_axis]] <- h
		slc <- do.call(`[`, c(list(object), idx, list(drop = FALSE)))
		c(mean = mean(slc, na.rm = TRUE),
		  sd   = stats::sd(as.numeric(slc), na.rm = TRUE),
		  min  = min(slc, na.rm = TRUE),
		  max  = max(slc, na.rm = TRUE))
	}, numeric(4))
	out <- as.data.frame(t(per_h))
	out$h <- seq_len(d[h_axis])
	out <- out[, c("h", "mean", "sd", "min", "max")]
	out
}

#' Reshape a forecast or forecast-CI object into a long data frame
#'
#' Internal helper used by [plot.dbn_forecast()] / [autoplot.dbn_forecast()].
#' Handles both the 4D / 5D `dbn_forecast` array and the
#' `dbn_forecast_ci` list of three 4D arrays.
#' @keywords internal
#' @noRd
.dbn_forecast_long <- function(x) {
	# detect ci-form (list) vs plain forecast array
	if (inherits(x, "dbn_forecast_ci")) {
		arr_m <- x$mean; arr_l <- x$lower; arr_u <- x$upper
		d <- dim(arr_m)
		dn <- dimnames(arr_m)
		grid <- expand.grid(
			sender   = dn[[1]] %||% seq_len(d[1]),
			receiver = dn[[2]] %||% seq_len(d[2]),
			relation = dn[[3]] %||% seq_len(d[3]),
			h        = dn[[4]] %||% paste0("h+", seq_len(d[4])),
			KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
		)
		grid$mean  <- as.vector(arr_m)
		grid$lower <- as.vector(arr_l)
		grid$upper <- as.vector(arr_u)
		return(grid)
	}
	d <- dim(x)
	dn <- dimnames(x)
	# 4D: [sender, receiver, relation, h] (summary = "mean")
	# 5D: [sender, receiver, relation, h, draws] (summary = "none" / "draws")
	if (length(d) == 4L) {
		grid <- expand.grid(
			sender   = dn[[1]] %||% seq_len(d[1]),
			receiver = dn[[2]] %||% seq_len(d[2]),
			relation = dn[[3]] %||% seq_len(d[3]),
			h        = dn[[4]] %||% paste0("h+", seq_len(d[4])),
			KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
		)
		grid$value <- as.vector(x)
		grid
	} else if (length(d) == 5L) {
		grid <- expand.grid(
			sender   = dn[[1]] %||% seq_len(d[1]),
			receiver = dn[[2]] %||% seq_len(d[2]),
			relation = dn[[3]] %||% seq_len(d[3]),
			h        = dn[[4]] %||% paste0("h+", seq_len(d[4])),
			draw     = dn[[5]] %||% seq_len(d[5]),
			KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
		)
		grid$value <- as.vector(x)
		grid
	} else {
		cli::cli_abort("Unexpected forecast array shape: dim = {paste(d, collapse = 'x')}.")
	}
}

#' Plot method for dbn_forecast
#'
#' @description Renders the forecast cells as a faceted line / ribbon plot:
#'   horizon on x, value on y, faceted by (sender, receiver). Uses ggplot2
#'   if available; falls back to a base-R plot of mean-across-cells per h
#'   when not.
#' @param x A `dbn_forecast` array.
#' @param max_cells Cap on the number of (sender, receiver) facets to draw
#'   (default 16). Use a larger value or pre-subset for full networks.
#' @param ... Unused.
#' @return Invisibly, the ggplot object (or NULL for the base fallback).
#' @author Tosin Salau and Shahryar Minhas
#' @export
plot.dbn_forecast <- function(x, max_cells = 16L, ...) {
	if (!inherits(x, "dbn_forecast")) {
		cli::cli_abort(c(
			"{.fn plot.dbn_forecast} requires a {.cls dbn_forecast} object.",
			"x" = "Got {.cls {class(x)[1]}}.",
			"i" = "Produce one with {.code predict(fit, H = h)} or {.fn forecast.dbn}."
		))
	}
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_warn("Install {.pkg ggplot2} for the faceted forecast plot; falling back to base R.")
		d <- dim(x)
		mn <- apply(x, length(d) - if (length(d) == 5L) 1L else 0L, mean, na.rm = TRUE)
		# crude: mean across cells per horizon
		val_per_h <- if (length(d) == 4L) apply(x, 4, mean, na.rm = TRUE)
			else apply(x, 4, mean, na.rm = TRUE)
		plot(seq_along(val_per_h), val_per_h, type = "b",
			xlab = "horizon", ylab = "mean forecast across cells",
			main = "dbn forecast (base fallback)")
		return(invisible(NULL))
	}
	long <- .dbn_forecast_long(x)
	# cap facets to keep the plot legible
	dyads <- unique(long[, c("sender", "receiver")])
	if (nrow(dyads) > max_cells) {
		dyads <- dyads[seq_len(max_cells), , drop = FALSE]
		long <- merge(long, dyads, by = c("sender", "receiver"))
		cli::cli_inform("Displaying {.val {max_cells}} of {.val {nrow(unique(long[, c('sender', 'receiver')]))}} dyads. Pass {.arg max_cells} to raise the cap.")
	}
	p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$h, y = .data$value)) +
		ggplot2::geom_line(ggplot2::aes(group = if ("draw" %in% names(long)) .data$draw else 1),
			alpha = if ("draw" %in% names(long)) 0.3 else 1) +
		ggplot2::facet_grid(rows = ggplot2::vars(.data$sender),
			cols = ggplot2::vars(.data$receiver)) +
		ggplot2::labs(x = "Horizon", y = "Forecast value", title = "dbn forecast") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", colour = "black"),
			strip.text = ggplot2::element_text(colour = "white"),
			legend.position = "top",
			axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
	print(p)
	invisible(p)
}

#' Plot method for dbn_forecast_ci
#'
#' Identical to [plot.dbn_forecast()] but renders the lower / upper
#' bands as a ribbon around the mean.
#' @param x A `dbn_forecast_ci` list.
#' @param max_cells Maximum number of (sender, receiver) facets.
#' @param ... Unused.
#' @return Invisibly, the ggplot object.
#' @author Tosin Salau and Shahryar Minhas
#' @export
plot.dbn_forecast_ci <- function(x, max_cells = 16L, ...) {
	if (!inherits(x, "dbn_forecast_ci")) {
		cli::cli_abort(c(
			"{.fn plot.dbn_forecast_ci} requires a {.cls dbn_forecast_ci} object.",
			"x" = "Got {.cls {class(x)[1]}}.",
			"i" = "Produce one with {.code predict(fit, H = h, summary = \"ci\")}."
		))
	}
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_warn("Install {.pkg ggplot2} for the faceted forecast plot; falling back to base R.")
		return(invisible(NULL))
	}
	long <- .dbn_forecast_long(x)
	dyads <- unique(long[, c("sender", "receiver")])
	if (nrow(dyads) > max_cells) {
		dyads <- dyads[seq_len(max_cells), , drop = FALSE]
		long <- merge(long, dyads, by = c("sender", "receiver"))
		cli::cli_inform("Displaying {.val {max_cells}} of {.val {nrow(unique(long[, c('sender', 'receiver')]))}} dyads.")
	}
	p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$h, y = .data$mean, group = 1)) +
		ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
			fill = "grey80", alpha = 0.5) +
		ggplot2::geom_line() +
		ggplot2::facet_grid(rows = ggplot2::vars(.data$sender),
			cols = ggplot2::vars(.data$receiver)) +
		ggplot2::labs(x = "Horizon", y = "Forecast", title = "dbn forecast with credible intervals") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", colour = "black"),
			strip.text = ggplot2::element_text(colour = "white"),
			legend.position = "top",
			axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
	print(p)
	invisible(p)
}

#' autoplot methods for dbn_forecast objects
#'
#' Aliases for [plot.dbn_forecast()] / [plot.dbn_forecast_ci()] so users
#' coming from the `forecast` / `ggfortify` ecosystem can call
#' `ggplot2::autoplot(fc)` directly. Registered conditionally on
#' `ggplot2::autoplot` at package load so the `ggplot2` dependency stays
#' in `Suggests`, not `Imports`.
#' @param object A `dbn_forecast` or `dbn_forecast_ci`.
#' @param ... Forwarded to the corresponding `plot` method.
#' @return A ggplot object (invisibly, via the underlying plot method).
#' @keywords internal
autoplot.dbn_forecast <- function(object, ...) {
	plot.dbn_forecast(object, ...)
}

#' @rdname autoplot.dbn_forecast
#' @keywords internal
autoplot.dbn_forecast_ci <- function(object, ...) {
	plot.dbn_forecast_ci(object, ...)
}
