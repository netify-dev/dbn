####
# goodness-of-fit by network statistic, per time slice.
# for each named statistic and each time index it reports the observed
# value, the posterior-predictive distribution, and a per-time coverage
# indicator -- the ergm::gof() analog for dynamic fits.
####

#' Goodness-of-fit by network statistic over time
#'
#' For each requested network statistic and each time slice, compares
#' the observed value to the posterior-predictive distribution and
#' returns a tidy data frame with observed/posterior_mean/ci_lo/ci_hi/
#' in_ci columns. The default statistics are network-summary scalars
#' (density, reciprocity, transitivity, mean in-/out-degree); custom
#' statistics may be passed via the `stat_funs` argument as a named
#' list of functions.
#'
#' This is the analog of `ergm::gof()` for dynamic DBN fits: rather
#' than asking "did the model recover the global structural statistics
#' of the observed network?", it asks "for each time period, does the
#' observed statistic fall inside the model's posterior-predictive
#' interval?".
#'
#' @param fit A `dbn` fit object (dynamic / static; HMM and low-rank
#'   are not yet supported).
#' @param stats Character vector of built-in statistics to include.
#'   Subset of: `"density"`, `"reciprocity"`, `"transitivity"`,
#'   `"mean_in_degree"`, `"mean_out_degree"`. Default: all five.
#' @param stat_funs Optional named list of functions taking an
#'   `n_row x n_col` matrix and returning a scalar. If supplied,
#'   overrides `stats`.
#' @param draws Number of posterior-predictive draws to use (default 50).
#'   For posterior-package convention, `nsim` and `n_draws` are accepted
#'   aliases via `...`.
#' @param level Coverage level for the per-time PI (default 0.95).
#' @param threshold Edge threshold for stats that need to binarise
#'   continuous data (default `NULL`; passed to the underlying
#'   `stat_*` helpers).
#' @param seed Optional RNG seed for reproducibility of the PPD draws.
#' @param ... Additional arguments. `nsim` and `n_draws` are accepted as
#'   aliases for `draws` (matching the `posterior` package's convention).
#' @return A data frame with columns `stat`, `t`, `observed`,
#'   `posterior_mean`, `ci_lo`, `ci_hi`, `in_ci`. Class
#'   `c("dbn_gof", "data.frame")`.
#' @seealso [posterior_predict_dbn()], [stat_density()],
#'   [stat_reciprocity()], [stat_transitivity()]
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 8, seed = 1)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'            nscan = 200, burn = 100, verbose = FALSE)
#' g <- gof(fit, stats = c("density", "reciprocity"), draws = 30)
#' head(g)
#' }
gof.dbn <- function(fit,
                     stats = c("density", "reciprocity", "transitivity",
                               "mean_in_degree", "mean_out_degree"),
                     stat_funs = NULL,
                     draws = 50L, level = 0.95,
                     threshold = NULL, seed = NULL, ...) {
	if (!inherits(fit, "dbn")) cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	# alias the conventional posterior-package `nsim` to gof.dbn's `draws`,
	# rather than letting it fall through to `...` and producing a cryptic
	# "unused argument" error. mirror the same alias for `n_draws`.
	dots <- list(...)
	for (alias in c("nsim", "n_draws")) {
		if (alias %in% names(dots) && missing(draws)) {
			draws <- dots[[alias]]
			dots[[alias]] <- NULL
		} else if (alias %in% names(dots)) {
			cli::cli_warn(c(
				"Both {.arg draws} and {.arg {alias}} were specified; using {.arg draws}.",
				"i" = "Pass only one."
			))
			dots[[alias]] <- NULL
		}
	}
	if (length(dots) > 0L) {
		cli::cli_warn(c(
			"{.fun gof.dbn} ignoring unused argument{?s} {.arg {names(dots)}}.",
			"i" = "Pass {.arg stats}, {.arg stat_funs}, {.arg draws}, {.arg threshold}, or {.arg seed}."
		))
	}
	if (level <= 0 || level >= 1) cli::cli_abort("{.arg level} must be in (0, 1).")
	alpha <- (1 - level) / 2

	# build the function list from the named-stat default or a user override.
	# stat_reciprocity takes only X; stat_transitivity also accepts threshold.
	if (is.null(stat_funs)) {
		stats <- match.arg(stats, several.ok = TRUE)
		thr <- if (is.null(threshold)) 0 else threshold
		stat_funs <- list(
			density         = function(X) stat_density(X),
			reciprocity     = function(X) stat_reciprocity(X),
			transitivity    = function(X) stat_transitivity(X, threshold = thr),
			mean_in_degree  = function(X) mean(stat_in_degree(X), na.rm = TRUE),
			mean_out_degree = function(X) mean(stat_out_degree(X), na.rm = TRUE)
		)[stats]
	} else {
		if (!is.list(stat_funs) || is.null(names(stat_funs)) || any(names(stat_funs) == "")) {
			cli::cli_abort("{.arg stat_funs} must be a named list of functions.")
		}
		if (!all(vapply(stat_funs, is.function, logical(1L)))) {
			cli::cli_abort("All elements of {.arg stat_funs} must be functions.")
		}
	}

	Y <- fit$Y
	if (is.null(Y)) {
		cli::cli_abort("{.code fit$Y} not present; cannot compute GOF.")
	}
	dims <- fit$dims
	n_row <- dims$n_row; n_col <- dims$n_col; p <- dims$p %||% 1L; Tt <- dims$Tt

	# posterior-predictive replicates
	ppd <- posterior_predict_dbn(fit, draws = draws, seed = seed)
	if (is.list(ppd) && !is.null(ppd$replicates)) ppd <- ppd$replicates
	if (!is.list(ppd) || length(ppd) == 0L) {
		cli::cli_abort("Could not extract posterior-predictive replicates from {.fun posterior_predict_dbn}.")
	}

	# normalise Y_t and Yrep_t onto the family's observation scale so the
	# stat functions get apples-to-apples inputs; binary gets thresholded at
	# zero (the probit cutpoint).
	fam_name <- if (is.list(fit$family)) fit$family$name else fit$family
	to_obs_scale <- function(M) {
		if (identical(fam_name, "binary")) {
			out <- (M > 0) * 1
			out[is.na(M)] <- NA
			return(out)
		}
		M
	}
	# surface the first cell-level error per statistic so users see it
	# instead of silently getting NaN coverage
	err_seen <- character(0L)
	report_err <- function(sn, e) {
		if (!sn %in% err_seen) {
			err_seen[[length(err_seen) + 1L]] <<- sn
			cli::cli_warn(c(
				"Statistic {.field {sn}} errored on a cell -- masking with NA.",
				"x" = "{conditionMessage(e)}",
				"i" = "Subsequent errors on the same statistic are silenced."
			))
		}
		NA_real_
	}

	# for each (statistic, time slice), compute observed and per-replicate values
	rows <- list()
	for (sn in names(stat_funs)) {
		fn <- stat_funs[[sn]]
		for (t in seq_len(Tt)) {
			Y_t <- to_obs_scale(Y[, , 1L, t])
			obs <- tryCatch(as.numeric(fn(Y_t)), error = function(e) report_err(sn, e))
			rep_vals <- vapply(ppd, function(Yrep) {
				Yrep_t <- if (length(dim(Yrep)) == 4L) Yrep[, , 1L, t] else
					if (length(dim(Yrep)) == 3L) Yrep[, , t] else Yrep
				Yrep_t <- to_obs_scale(Yrep_t)
				tryCatch(as.numeric(fn(Yrep_t)), error = function(e) report_err(sn, e))
			}, numeric(1L))
			finite <- is.finite(rep_vals)
			pm <- if (any(finite)) mean(rep_vals[finite]) else NA_real_
			lo <- if (sum(finite) >= 2L) as.numeric(quantile(rep_vals[finite], alpha, na.rm = TRUE)) else NA_real_
			hi <- if (sum(finite) >= 2L) as.numeric(quantile(rep_vals[finite], 1 - alpha, na.rm = TRUE)) else NA_real_
			in_ci <- if (is.finite(obs) && is.finite(lo) && is.finite(hi)) (obs >= lo && obs <= hi) else NA
			rows[[length(rows) + 1L]] <- data.frame(
				stat = sn, t = t, observed = obs,
				posterior_mean = pm, ci_lo = lo, ci_hi = hi, in_ci = in_ci,
				stringsAsFactors = FALSE
			)
		}
	}
	out <- do.call(rbind, rows)
	rownames(out) <- NULL
	attr(out, "level") <- level
	attr(out, "n_draws") <- length(ppd)
	class(out) <- c("dbn_gof", "data.frame")
	out
}

#' Print method for dbn_gof
#'
#' @param x A `dbn_gof` object.
#' @param ... Unused.
#' @export
#' @keywords internal
print.dbn_gof <- function(x, ...) {
	level <- attr(x, "level") %||% 0.95
	n_draws <- attr(x, "n_draws") %||% NA
	cli::cli_h2("DBN goodness-of-fit summary")
	cli::cli_text("{nrow(x)} (stat, time) rows; {length(unique(x$stat))} statistic{?s}; {length(unique(x$t))} time slice{?s}; {n_draws} posterior-predictive draws; nominal coverage {sprintf('%.0f%%', 100 * level)}")
	# per-stat empirical coverage with masked-na count
	by_stat <- split(x, x$stat)
	for (sn in names(by_stat)) {
		sub <- by_stat[[sn]]
		n_total <- nrow(sub)
		n_used <- sum(!is.na(sub$in_ci))
		n_masked <- n_total - n_used
		cov <- mean(sub$in_ci, na.rm = TRUE)
		msg <- sprintf(
			"%s: empirical coverage = %.0f%% (%d/%d time slices in PI",
			sn, 100 * cov, sum(sub$in_ci, na.rm = TRUE), n_used
		)
		if (n_masked > 0L)
			msg <- paste0(msg, sprintf("; %d slice%s masked as NA", n_masked, if (n_masked == 1L) "" else "s"))
		msg <- paste0(msg, ")")
		cli::cli_li(msg)
	}
	cli::cli_text("")
	print.data.frame(head(x, 10L))
	if (nrow(x) > 10L) cli::cli_text("... ({nrow(x) - 10L} more rows)")
	invisible(x)
}

#' Goodness-of-fit by network statistic (generic)
#'
#' Compute a network statistic on the observed data and on posterior /
#' bootstrap replicates, then compare. Dispatches to [gof.dbn()] for
#' single-chain fits and [gof.dbn_multichain()] for multichain fits.
#'
#' @param fit A `dbn` or `dbn_multichain` fit.
#' @param ... Forwarded to the method (e.g. `stats`, `draws`, `seed`).
#'   See [gof.dbn()] for the full argument list.
#' @return A long data frame with one row per (statistic, time, replicate
#'   comparison). See [gof.dbn()] for the precise schema.
#' @seealso [gof.dbn()], [gof.dbn_multichain()],
#'   [stat_density()], [stat_reciprocity()], [stat_transitivity()].
#' @author Tosin Salau and Shahryar Minhas
#' @export
gof <- function(fit, ...) UseMethod("gof")
