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
    if (!inherits(x, "dbn")) stop("x must be a dbn object")

    # dispatch to model-specific plot function
    plot_fun <- switch(x$model,
        static = "plot_static",
        dynamic = "plot_dynamic",
        lowrank = "plot_lowrank",
        hmm = "plot_hmm",
        stop("Unknown model: ", x$model)
    )

    # call the function
    do.call(plot_fun, list(x, ...))
}

#' @export
#' @method summary dbn
summary.dbn <- function(object, ...) {
    if (!inherits(object, "dbn")) stop("object must be a dbn object")

    # dispatch to model-specific summary function
    summary_fun <- switch(object$model,
        static = "summary_static",
        dynamic = "summary_dynamic",
        lowrank = "summary_lowrank",
        hmm = "summary_hmm",
        stop("Unknown model: ", object$model)
    )

    # call the function
    do.call(summary_fun, list(object, ...))
}

#' @export
#' @method predict dbn
predict.dbn <- function(object, ...) {
    if (!inherits(object, "dbn")) stop("object must be a dbn object")

    # dispatch to model-specific predict/simulate function
    predict_fun <- switch(object$model,
        static = "simulate_static",
        dynamic = "simulate_dynamic",
        lowrank = "predict_lowrank",
        hmm = "predict_hmm",
        stop("Unknown model: ", object$model)
    )

    # call the function
    do.call(predict_fun, c(list(object), list(...)))
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

#' @export
#' @method print dbn
print.dbn <- function(x, ...) {
    if (!inherits(x, "dbn")) stop("x must be a dbn object")
    
    # print header
    cat("Dynamic Bilinear Network (DBN) Model\n")
    cat(rep("-", 40), "\n", sep = "")
    
    # model type
    cat("Model Type: ", toupper(x$model), "\n", sep = "")
    
    # family
    if (!is.null(x$family)) {
        cat("Family: ", x$family, "\n", sep = "")
    }
    
    # data dimensions
    if (!is.null(x$dims)) {
        cat("\nData Dimensions:\n")
        cat("  Nodes: ", x$dims$m, "\n", sep = "")
        cat("  Relations: ", x$dims$p, "\n", sep = "")
        if (!is.null(x$dims$n)) {
            cat("  Time points: ", x$dims$n, "\n", sep = "")
        } else if (!is.null(x$dims$Tt)) {
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
    
    # model-specific information
    if (x$model == "hmm" && !is.null(x$R)) {
        cat("\nRegimes: ", x$R, "\n", sep = "")
    } else if (x$model == "lowrank" && !is.null(x$rank)) {
        cat("\nRank: ", x$rank, "\n", sep = "")
    }
    
    # components
    cat("\nModel Components:\n")
    component_names <- names(x)
    # exclude meta fields
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
    
    # diagnostics if available
    if (!is.null(x$diagnostics) && !is.null(x$diagnostics$dic)) {
        cat("\nModel Diagnostics:\n")
        cat("  DIC: ", round(x$diagnostics$dic, 2), "\n", sep = "")
    }
    
    invisible(x)
}
