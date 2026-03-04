####
# OpenMP thread management for parallel MCMC
####

#' Get the current number of threads used by dbn
#'
#' @description Returns the number of OpenMP threads currently configured for
#'   parallel MCMC computation. Defaults to 1 (single-threaded).
#' @return Integer number of threads currently set for parallel computation
#' @seealso \code{\link{set_dbn_threads}}
#' @export
#' @examples
#' get_dbn_threads()
get_dbn_threads <- function() {
	getOption("dbn.n_threads", 1L)
}

#' Set the number of threads used by dbn
#'
#' @description Controls the number of OpenMP threads used for parallel MCMC
#'   updates in the dynamic model. Parallelization applies to the row-wise A
#'   and B FFBS updates and variance computations. The setting persists for
#'   the current R session until changed.
#'
#'   For best performance, set to the number of physical (not logical) cores.
#'   Speedup scales with network size: expect 2--4x improvement with 4--8
#'   cores for networks with 15+ actors. For small networks the overhead of
#'   thread management can exceed the benefit.
#'
#'   Requires OpenMP support at compile time (standard on Linux and Windows;
#'   on macOS install via \code{brew install libomp}). Without OpenMP the
#'   package runs single-threaded regardless of this setting.
#' @param n_threads Integer number of threads to use for parallel computation.
#'   Use NULL to reset to default (1 thread).
#' @return The previous value of n_threads (invisibly)
#' @seealso \code{\link{get_dbn_threads}}
#' @export
#' @examples
#' set_dbn_threads(4)
#' set_dbn_threads(parallel::detectCores(logical = FALSE))
#' set_dbn_threads(NULL)
set_dbn_threads <- function(n_threads = NULL) {
	old <- get_dbn_threads()

	if (is.null(n_threads)) {
		n_threads <- 1L
	} else {
		n_threads <- as.integer(n_threads)
		if (n_threads < 1) {
			cli::cli_abort("{.arg n_threads} must be at least 1.")
		}
		max_threads <- parallel::detectCores()
		if (n_threads > max_threads) {
			warning(sprintf(
				"n_threads (%d) exceeds available cores (%d), using %d",
				n_threads, max_threads, max_threads
			))
			n_threads <- max_threads
		}
	}

	options(dbn.n_threads = n_threads)
	invisible(old)
}
