####
# Run-experience helpers: honest status reporting for long MCMC fits.
#
# A dynamic DBN fit can run for hours. These helpers make the package
# communicate that run honestly. A pre-flight panel states the plan and cost
# before the loop; a log-safe heartbeat reports measured progress and ETA
# during the loop (surviving non-interactive or redirected output, where an
# animated cli progress bar does not render); a closing summary reports the
# actual wall time and object size.
####

#' Format a duration in seconds as a compact human-readable string
#' @keywords internal
#' @noRd
.dbn_fmt_duration <- function(seconds) {
	if (length(seconds) != 1L || !is.finite(seconds) || seconds < 0) {
		return("--")
	}
	seconds <- round(seconds)
	h <- seconds %/% 3600
	m <- (seconds %% 3600) %/% 60
	s <- seconds %% 60
	if (h > 0) {
		sprintf("%dh%02dm", h, m)
	} else if (m > 0) {
		sprintf("%dm%02ds", m, s)
	} else {
		sprintf("%ds", s)
	}
}

#' Pre-flight panel printed once before the MCMC loop
#'
#' States dimensions, the iteration plan, and the estimated output object size
#' and peak RAM. Runtime is deliberately not predicted here: a number guessed
#' before any iteration runs would be unreliable, and the heartbeat reports a
#' measured ETA within the first few iterations instead.
#' @keywords internal
#' @noRd
.dbn_preflight <- function(model, n_row, n_col, p, Tt,
						   burn, nscan, odens, n_keep,
						   object_gb = NA_real_, peak_gb = NA_real_,
						   notes = character()) {
	n_iter <- burn + nscan
	cli::cli_h2("dbn {model} fit")
	bullets <- c(
		"*" = "Network {n_row} x {n_col}, {p} relation{?s}, {Tt} time point{?s}.",
		"*" = "MCMC {n_iter} iterations ({burn} burn-in + {nscan} sampling); {n_keep} draw{?s} kept (odens {odens})."
	)
	if (is.finite(object_gb)) {
		size_line <- paste0("Output object ~", sprintf("%.2f", object_gb),
							" GB; peak RAM ~", sprintf("%.2f", peak_gb), " GB.")
		bullets <- c(bullets, "*" = size_line)
	}
	bullets <- c(bullets,
				 "*" = "Runtime calibrating; first measured ETA at iteration 10.")
	cli::cli_bullets(bullets)
	for (nt in notes) cli::cli_alert_warning(nt)
	invisible()
}

#' Log-safe progress heartbeat emitted from within the MCMC loop
#'
#' Prints a plain status line (iteration, percent, measured elapsed time, and
#' measured ETA) on an early-milestone schedule and then every `every`
#' iterations. Unlike an animated progress bar, this survives non-interactive
#' and redirected output, which is where long batch fits actually run.
#' @keywords internal
#' @noRd
.dbn_heartbeat <- function(iter, n_iter, t_start, every) {
	milestones <- c(10L, 25L, 50L, 100L)
	periodic <- is.numeric(every) && every > 0 && (iter %% every == 0)
	if (!(periodic || iter %in% milestones || iter == n_iter)) {
		return(invisible())
	}
	elapsed <- proc.time()[[3]] - t_start
	frac <- iter / n_iter
	eta  <- if (frac > 0) elapsed * (1 - frac) / frac else NA_real_
	pct  <- round(100 * frac)
	el   <- .dbn_fmt_duration(elapsed)
	et   <- .dbn_fmt_duration(eta)
	cli::cli_inform("dbn | iter {iter}/{n_iter} ({pct}%) | elapsed {el} | ETA {et}")
	invisible()
}

#' Closing summary printed once after the MCMC loop
#' @keywords internal
#' @noRd
.dbn_closing <- function(model, t_start, object_gb = NA_real_, conv = NULL) {
	wall <- .dbn_fmt_duration(proc.time()[[3]] - t_start)
	msg <- c("v" = "dbn {model} fit complete in {wall}.")
	if (is.finite(object_gb)) {
		msg <- c(msg, "*" = "Output object ~{sprintf('%.2f', object_gb)} GB; save with {.code compress = \"xz\"} to shrink the file.")
	}
	if (!is.null(conv)) {
		msg <- c(msg, "i" = conv)
	}
	cli::cli_bullets(msg)
	invisible()
}

#' Heuristic: does a square network's data look undirected (symmetric)?
#'
#' Used to nudge a user who set `symmetric = FALSE` on data that is in fact
#' symmetric. Returns FALSE (no nudge) for anything it cannot cleanly check.
#' @keywords internal
#' @noRd
.dbn_data_looks_undirected <- function(Y, n_row, n_col) {
	if (n_row != n_col || is.null(Y) || !is.numeric(Y)) return(FALSE)
	d <- dim(Y)
	if (is.null(d) || length(d) < 2L || d[1] != n_row || d[2] != n_col) {
		return(FALSE)
	}
	n_slice <- if (length(d) <= 2L) 1L else as.integer(prod(d[-(1:2)]))
	if (is.na(n_slice) || n_slice < 1L) return(FALSE)
	dim(Y) <- c(n_row, n_col, n_slice)
	take <- unique(round(seq(1, n_slice, length.out = min(6L, n_slice))))
	for (k in take) {
		M <- Y[, , k]
		off <- suppressWarnings(max(abs(M - t(M)), na.rm = TRUE))
		if (!is.finite(off) || off > 1e-8) return(FALSE)
	}
	TRUE
}

#' One-line mixing note from a kept variance trace (cheap, dependency-free)
#'
#' Reports the lag-1 autocorrelation of a scalar parameter trace and flags it
#' when high. Not a substitute for full convergence diagnostics.
#' @keywords internal
#' @noRd
.dbn_mixing_note <- function(trace, label = "sigma2") {
	trace <- trace[is.finite(trace)]
	if (length(trace) < 10L || stats::sd(trace) == 0) return(NULL)
	ac1 <- suppressWarnings(stats::cor(trace[-1], trace[-length(trace)]))
	if (!is.finite(ac1)) return(NULL)
	tag <- if (ac1 > 0.9) " (high; consider more iterations or thinning)" else ""
	paste0(label, " lag-1 autocorrelation ", sprintf("%.2f", ac1), tag,
		   "; verify convergence before interpreting.")
}
