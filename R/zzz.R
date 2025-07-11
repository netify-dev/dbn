#' Package initialization
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
    # Ensure Rcpp is loaded for proper ABI
    requireNamespace("Rcpp", quietly = TRUE)
    requireNamespace("RcppArmadillo", quietly = TRUE)

    # Load the compiled code - handle development vs installed
    tryCatch(
        {
            library.dynam("dbn", pkgname, libname)
        },
        error = function(e) {
            # During development, the .so might be in src/
            if (file.exists("src/dbn.so")) {
                dyn.load("src/dbn.so")
            }
        }
    )

    # Enable C++ optimizations by default
    op <- options()
    op.dbn <- list(
        dbn.use_cpp_stability = TRUE,
        dbn.use_ffbs_dlm_cpp = TRUE,
        dbn.use_ffbs_cpp = TRUE,
        dbn.use_cpp_ranklik = TRUE,
        dbn.use_batch_ffbs = FALSE, # Disabled due to dimension issues
        dbn.n_threads = 1L # Default to single thread for safety
    )
    toset <- !(names(op.dbn) %in% names(op))
    if (any(toset)) options(op.dbn[toset])

    invisible()
}
