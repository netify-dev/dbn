####
#' Fit a DBN Model to Several Networks in Parallel
#'
#' @description Fits [dbn()] to each network in a list, running the fits as
#'   independent parallel processes. This is the supported way to use multiple
#'   cores: the dynamic sampler's per-iteration cost is dominated by a
#'   forward-filtering/backward-sampling pass that is sequential over time and
#'   cannot be parallelized within a single fit, so [set_dbn_threads()] gives
#'   little within-fit speedup. Independent fits, by contrast, are
#'   embarrassingly parallel; six fits across six cores finish in about the
#'   wall-clock time of one.
#'
#'   Each worker pins its BLAS and dbn thread counts to one, because the
#'   parallelism is across fits and nested threading only oversubscribes
#'   cores. A fit that errors does not abort the batch: its slot in the
#'   returned list holds a `dbn_fit_error` object instead.
#'
#' @param Y_list A named list of network data arrays, one per fit. Unnamed
#'   lists are given names `fit1`, `fit2`, and so on.
#' @param ... Arguments forwarded to [dbn()] for every fit, for example
#'   `model`, `family`, `nscan`, `burn`, `odens`, `symmetric`, and `keep`.
#'   Do not pass `seed` here; use `seeds`.
#' @param mc.cores Number of parallel workers. Defaults to the smaller of the
#'   number of fits and one less than the detected core count. Forking is
#'   unavailable on Windows, where the fits run sequentially.
#' @param seeds Optional integer vector of per-fit random seeds, one per fit.
#'   Defaults to a distinct, reproducible seed per fit.
#' @param save_dir Optional directory. When set, each completed fit is saved
#'   there as `<name>.rda` (xz-compressed) as soon as it finishes, giving
#'   crash-safety and a visible progress signal for long batches.
#' @param verbose Whether to print batch-level progress (per-fit output is
#'   always suppressed to avoid interleaved logs).
#' @return A named list parallel to `Y_list`. Each element is a fitted `dbn`
#'   object, or a `dbn_fit_error` object if that fit failed.
#' @seealso [dbn()], [set_dbn_threads()], [estimate_memory()]
#' @examples
#' \donttest{
#' sims <- list(
#'   net_a = simulate_dynamic_dbn(n = 6, time = 8, seed = 1)$Y,
#'   net_b = simulate_dynamic_dbn(n = 6, time = 8, seed = 2)$Y
#' )
#' fits <- dbn_fit_many(sims, model = "dynamic", family = "gaussian",
#'                      nscan = 200, burn = 100, mc.cores = 2)
#' }
#' @export
dbn_fit_many <- function(Y_list, ..., mc.cores = NULL, seeds = NULL,
						 save_dir = NULL, verbose = TRUE) {
	if (!is.list(Y_list) || length(Y_list) == 0L) {
		cli::cli_abort("{.arg Y_list} must be a non-empty list of network data arrays.")
	}
	n_fit <- length(Y_list)

	dots_names <- names(list(...))
	if ("seed" %in% dots_names) {
		cli::cli_abort(c(
			"Do not pass {.arg seed} through {.fn dbn_fit_many}.",
			"i" = "Use the {.arg seeds} argument to set one seed per fit."
		))
	}

	fit_names <- names(Y_list)
	if (is.null(fit_names) || any(fit_names == "" | is.na(fit_names))) {
		fit_names <- paste0("fit", seq_len(n_fit))
	}

	if (is.null(seeds)) {
		seeds <- 6886L + seq_len(n_fit) - 1L
	} else if (length(seeds) != n_fit) {
		cli::cli_abort("{.arg seeds} must have one entry per fit ({n_fit}).")
	}

	avail <- parallel::detectCores()
	if (is.na(avail) || avail < 1L) avail <- 1L
	if (is.null(mc.cores)) {
		mc.cores <- min(n_fit, max(1L, avail - 1L))
	}
	mc.cores <- max(1L, min(as.integer(mc.cores), n_fit))

	# mclapply relies on forking, which is unavailable on Windows
	if (.Platform$OS.type != "unix" && mc.cores > 1L) {
		if (verbose) {
			cli::cli_alert_info("Forking is unavailable on Windows; running the {n_fit} fit{?s} sequentially.")
		}
		mc.cores <- 1L
	}

	if (!is.null(save_dir) && !dir.exists(save_dir)) {
		dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
	}

	if (verbose) {
		cli::cli_h2("dbn_fit_many: {n_fit} fit{?s} on {mc.cores} worker{?s}")
		cli::cli_bullets(c(
			"*" = "One single-threaded process per fit (BLAS and dbn threads pinned to 1).",
			"i" = "Per-fit progress is suppressed; this returns when all fits finish."
		))
	}

	fit_one <- function(i) {
		# the speedup is across fits, not within them: pin nested threading
		# off so workers do not oversubscribe the cores.
		if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
			try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
		}
		try(set_dbn_threads(1L, quiet = TRUE), silent = TRUE)
		res <- tryCatch(
			dbn(Y_list[[i]], ..., seed = seeds[i], verbose = FALSE),
			error = function(e) structure(
				list(error = conditionMessage(e), fit_name = fit_names[i]),
				class = "dbn_fit_error"
			)
		)
		if (!is.null(save_dir) && !inherits(res, "dbn_fit_error")) {
			fit <- res
			try(save(fit, compress = "xz",
					 file = file.path(save_dir, paste0(fit_names[i], ".rda"))),
				silent = TRUE)
		}
		res
	}

	t0 <- proc.time()[[3]]
	results <- if (mc.cores > 1L) {
		parallel::mclapply(seq_len(n_fit), fit_one, mc.cores = mc.cores)
	} else {
		lapply(seq_len(n_fit), fit_one)
	}
	names(results) <- fit_names

	# a worker that died is reported by mclapply as a try-error element
	for (i in seq_len(n_fit)) {
		if (inherits(results[[i]], "try-error") || is.null(results[[i]])) {
			detail <- paste(as.character(results[[i]]), collapse = " ")
			if (!nzchar(detail)) detail <- "process returned no result"
			results[[i]] <- structure(
				list(error = paste("worker failed:", detail),
					 fit_name = fit_names[i]),
				class = "dbn_fit_error"
			)
		}
	}

	failed <- vapply(results, inherits, logical(1), "dbn_fit_error")
	if (verbose) {
		wall <- .dbn_fmt_duration(proc.time()[[3]] - t0)
		if (any(failed)) {
			cli::cli_alert_warning("{sum(failed)} of {n_fit} fit{?s} failed: {.val {fit_names[failed]}}.")
		}
		cli::cli_alert_success("dbn_fit_many complete in {wall}: {sum(!failed)}/{n_fit} succeeded.")
	}
	results
}
####

####
#' Print a Failed Fit from dbn_fit_many
#'
#' @description Print method for a `dbn_fit_error`, the placeholder object
#'   [dbn_fit_many()] returns in place of a fit that errored.
#' @param x A `dbn_fit_error` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.dbn_fit_error <- function(x, ...) {
	cli::cli_alert_danger("dbn fit {.val {x$fit_name}} failed: {x$error}")
	invisible(x)
}
####
