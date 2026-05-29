####
# openmp thread management for parallel mcmc
####

#' Get the current number of threads used by dbn
#'
#' @description Returns the number of OpenMP threads currently configured for
#'   parallel MCMC computation. Defaults to 1 (single-threaded).
#' @return Integer number of threads currently set for parallel computation
#' @seealso \code{\link{set_dbn_threads}}
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' get_dbn_threads()
get_dbn_threads <- function() {
	getOption("dbn.n_threads", 1L)
}

#' Set the number of threads used by dbn
#'
#' @description Controls the number of OpenMP threads used for parallel
#'   computation. The setting persists for the current R session until changed.
#'
#'   \strong{Scope of the speedup.} Within a single dynamic-model fit, the MCMC
#'   iteration is dominated by the forward-filtering/backward-sampling (FFBS)
#'   pass over the latent state, a Kalman recursion over time in which each
#'   step depends on the previous one and so cannot be parallelized across
#'   time. Thread-level parallelism covers only the row-wise A and B
#'   influence-matrix updates, a minority of per-iteration cost. On the
#'   dynamic Gaussian model, raising \code{n_threads} therefore yields little
#'   or no within-fit speedup and can be marginally slower, because of
#'   per-iteration thread-spawn overhead and BLAS oversubscription (each
#'   OpenMP thread spawning its own multi-threaded BLAS calls).
#'
#'   \strong{Recommended use of parallelism.} Research workflows usually fit
#'   several independent models (different cases, specifications, or seeds).
#'   That is the level at which parallelism pays off: run each fit as a
#'   separate single-threaded process (or via \code{parallel::mclapply}),
#'   keeping \code{set_dbn_threads(1)} and the BLAS at one thread per process
#'   (\code{RhpcBLASctl::blas_set_num_threads(1)} or environment variable
#'   \code{OMP_NUM_THREADS=1}). Six independent fits across six processes give
#'   a near-linear wall-clock speedup; sixteen threads inside one fit do not.
#'
#'   The \code{lowrank} and \code{hmm} variants have larger parallelizable
#'   kernels and can benefit more from \code{n_threads > 1}.
#'
#'   Requires OpenMP support at compile time (standard on Linux and Windows;
#'   on macOS install via \code{brew install libomp}). Without OpenMP the
#'   package runs single-threaded regardless of this setting.
#' @param n_threads Integer number of threads to use for parallel computation.
#'   Use NULL to reset to default (1 thread).
#' @param quiet Logical. If FALSE (default), print a one-time note about where
#'   threading does and does not help when \code{n_threads > 1} is requested.
#' @return The previous value of n_threads (invisibly)
#' @seealso \code{\link{get_dbn_threads}}
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' set_dbn_threads(4)
#' set_dbn_threads(NULL)
set_dbn_threads <- function(n_threads = NULL, quiet = FALSE) {
	old <- get_dbn_threads()

	if (is.null(n_threads)) {
		n_threads <- 1L
	} else {
		# validate upfront so non-numeric / NA / decimal inputs give a
		# directive error instead of leaking base R coercion warnings.
		if (length(n_threads) != 1L || !is.numeric(n_threads) || !is.finite(n_threads)) {
			cli::cli_abort(c(
				"{.arg n_threads} must be a single positive whole number.",
				"x" = "Got {.cls {class(n_threads)[1]}} of length {length(n_threads)} = {.val {n_threads}}."
			))
		}
		if (n_threads != round(n_threads)) {
			cli::cli_abort(c(
				"{.arg n_threads} must be a whole number; got {.val {n_threads}}."
			))
		}
		n_threads <- as.integer(n_threads)
		if (n_threads < 1) {
			cli::cli_abort("{.arg n_threads} must be at least 1.")
		}
		max_threads <- parallel::detectCores()
		if (n_threads > max_threads) {
			cli::cli_warn(c(
				"{.arg n_threads} ({n_threads}) exceeds available cores ({max_threads}); clamping to {max_threads}."
			))
			n_threads <- max_threads
		}
	}

	if (!quiet && n_threads > 1L && !isTRUE(getOption("dbn.threads_note_shown"))) {
		cli::cli_inform(c(
			"i" = "{.code set_dbn_threads({n_threads})}: within-fit threading gives limited speedup on the dynamic model.",
			"i" = "The per-iteration FFBS pass is a sequential Kalman recursion over time and cannot be parallelized across time steps.",
			"i" = "For multiple fits, run them as separate single-threaded processes instead; see {.code ?set_dbn_threads}."
		))
		options(dbn.threads_note_shown = TRUE)
	}

	options(dbn.n_threads = n_threads)
	invisible(old)
}
