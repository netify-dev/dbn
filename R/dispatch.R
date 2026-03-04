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
		cli::cli_abort("Unknown model: {.val {x$model}}.")
	)
	do.call(plot_fun, list(x, ...))
}

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
		cli::cli_abort("Unknown model: {.val {object$model}}.")
	)
	do.call(summary_fun, list(object, ...))
}

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
