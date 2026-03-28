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
#' @seealso [dbn()], [plot_trace()], [check_convergence()], [summary_dbn()]
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' plot(fit)
#' }
#' @export
#' @method plot dbn
plot.dbn <- function(x, ...) {
	if (!inherits(x, "dbn")) cli::cli_abort("{.arg x} must be a {.cls dbn} object.")
	plot_fun <- switch(x$model,
		static = "plot_static",
		dynamic = "plot_dynamic",
		lowrank = "plot_lowrank",
		lowrank_accurate = "plot_lowrank",
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
#' @seealso [dbn()], [param_summary()], [check_convergence()], [plot_dbn()]
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' summary(fit)
#' }
#' @export
#' @method summary dbn
summary.dbn <- function(object, ...) {
	if (!inherits(object, "dbn")) cli::cli_abort("{.arg object} must be a {.cls dbn} object.")
	summary_fun <- switch(object$model,
		static = "summary_static",
		dynamic = "summary_dynamic",
		lowrank = "summary_lowrank",
		lowrank_accurate = "summary_lowrank",
		hmm = "summary_hmm",
		piecewise = "summary_piecewise",
		cli::cli_abort("Unknown model: {.val {object$model}}.")
	)
	do.call(summary_fun, list(object, ...))
}

#' Predict from a Fitted DBN Model
#'
#' @description Generates predictions from a fitted DBN model. With no
#'   arguments (or PPD-only arguments), returns posterior predictive
#'   replications of the observed data for model checking. With forecasting
#'   arguments (`H`, `S`, `summary`), propagates the estimated dynamics
#'   forward to produce multi-step-ahead forecasts.
#'
#' @details
#' **Posterior predictive distribution (default):** Simulates new datasets
#' from the fitted model to compare against the observed data. Use
#' [plot_ppc_ecdf()] or [plot_ppc_density()] for visual checks.
#'
#' **Forecasting:** Pass `H` (horizon), `S` or `draws` (number of
#' simulations), and optionally `summary = "mean"` to average across draws.
#' The forecast propagates \eqn{A_t}, \eqn{B_t}, and process noise forward
#' from the final observed time point.
#'
#' @param object A fitted `dbn` object returned by [dbn()].
#' @param ... Model-specific arguments:
#'   \describe{
#'     \item{H}{Forecast horizon (number of steps ahead). Triggers
#'       forecasting mode.}
#'     \item{S, draws}{Number of posterior draws to use for forecasting.}
#'     \item{summary}{If `"mean"`, return the pointwise mean across draws.
#'       If `NULL` (default), return the full array of simulated
#'       trajectories.}
#'     \item{ndraws}{Number of draws for posterior predictive checks
#'       (default: 100).}
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
#' fc <- predict(fit, H = 5, S = 50, summary = "mean")
#' dim(fc)
#'
#' # posterior predictive check
#' ppd <- predict(fit, ndraws = 10)
#' }
#' @export
#' @method predict dbn
predict.dbn <- function(object, ...) {
	if (!inherits(object, "dbn")) cli::cli_abort("{.arg object} must be a {.cls dbn} object.")
	# HMM and lowrank have dedicated forecast functions with H/draws/summary args
	if (object$model == "hmm") return(predict_hmm(object, ...))
	if (object$model %in% c("lowrank", "lowrank_accurate")) return(predict_lowrank(object, ...))
	# check whether the caller is requesting a forecast (H or S args)
	dots <- list(...)
	is_forecast <- any(c("H", "S", "summary") %in% names(dots))
	if (!is_forecast && length(dots) > 0) {
		known_ppd <- c("ndraws", "seed", "draws")
		unknown <- setdiff(names(dots), known_ppd)
		if (length(unknown) > 0) {
			cli::cli_warn(c(
				"Unrecognized argument{?s}: {.arg {unknown}}.",
				"i" = "For forecasting, use {.arg H}, {.arg S}, or {.arg summary}.",
				"i" = "Falling back to posterior predictive distribution."
			))
		}
	}
	if (is_forecast) {
		# route to model-specific simulation/forecast function
		predict_fun <- switch(object$model,
			static = "simulate_static",
			dynamic = "simulate_dynamic",
			piecewise = "simulate_piecewise",
			cli::cli_abort("Unknown model: {.val {object$model}}.")
		)
		return(do.call(predict_fun, list(object, ...)))
	}
	# default to posterior predictive distribution (requires family info)
	# only pass recognized PPD args (ndraws, seed, draws)
	ppd_args <- dots[intersect(names(dots), c("ndraws", "seed", "draws"))]
	fam <- tryCatch(get_family(object), error = function(e) NULL)
	if (!is.null(fam)) {
		return(do.call(posterior_predict_dbn, c(list(object), ppd_args)))
	}
	# no family info available, try simulation instead
	predict_fun <- switch(object$model,
		static = "simulate_static",
		dynamic = "simulate_dynamic",
		piecewise = "simulate_piecewise",
		cli::cli_abort("Unknown model: {.val {object$model}}.")
	)
	do.call(predict_fun, c(list(object), ppd_args))
}

#' Plot method for DBN objects
#'
#' @description Wrapper function for plot.dbn method
#' @param x A dbn object
#' @param ... Additional arguments passed to model-specific plot functions
#' @return Invisibly returns the plot object
#' @export
plot_dbn <- function(x, ...) plot.dbn(x, ...)

#' Summary method for DBN objects
#'
#' @description Wrapper function for summary.dbn method
#' @param object A dbn object
#' @param ... Additional arguments passed to model-specific summary functions
#' @return A summary object specific to the model type
#' @export
summary_dbn <- function(object, ...) summary.dbn(object, ...)

####
# print method for dbn objects
####

#' @export
#' @method print dbn
print.dbn <- function(x, ...) {
	if (!inherits(x, "dbn")) cli::cli_abort("{.arg x} must be a {.cls dbn} object.")

	# header
	cat("Dynamic Bilinear Network (DBN) Model\n")
	cat(rep("-", 40), "\n", sep = "")

	# model type
	cat("Model Type: ", toupper(x$model), "\n", sep = "")

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

	# MCMC settings
	if (!is.null(x$settings)) {
		cat("\nMCMC Settings:\n")
		cat("  Iterations: ", x$settings$nscan, "\n", sep = "")
		cat("  Burn-in: ", x$settings$burn, "\n", sep = "")
		cat("  Thinning: ", x$settings$odens, "\n", sep = "")
		cat("  Saved draws: ", x$settings$draws, "\n", sep = "")
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

	invisible(x)
}
