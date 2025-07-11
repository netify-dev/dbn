#' Get the current number of threads used by dbn
#'
#' @return Integer number of threads currently set for parallel computation
#' @export
#' @examples
#' # Get current thread count
#' get_dbn_threads()
get_dbn_threads <- function() {
    getOption("dbn.n_threads", 1L)
}

#' Set the number of threads used by dbn
#'
#' @param n_threads Integer number of threads to use for parallel computation.
#'   Use NULL to reset to default (1 thread).
#' @return The previous value of n_threads (invisibly)
#' @export
#' @examples
#' # Set to use 4 threads
#' set_dbn_threads(4)
#' 
#' # Set to use half of available cores
#' set_dbn_threads(parallel::detectCores() / 2)
#' 
#' # Reset to default (1 thread)
#' set_dbn_threads(NULL)
set_dbn_threads <- function(n_threads = NULL) {
    old <- get_dbn_threads()
    
    if (is.null(n_threads)) {
        n_threads <- 1L
    } else {
        n_threads <- as.integer(n_threads)
        if (n_threads < 1) {
            stop("n_threads must be at least 1")
        }
        max_threads <- parallel::detectCores()
        if (n_threads > max_threads) {
            warning(sprintf("n_threads (%d) exceeds available cores (%d), using %d", 
                          n_threads, max_threads, max_threads))
            n_threads <- max_threads
        }
    }
    
    options(dbn.n_threads = n_threads)
    invisible(old)
}