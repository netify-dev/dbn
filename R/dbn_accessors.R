####
# standard model-accessor S3 methods for `dbn` fits.
#
# users coming from lm/glm reach for coef(fit), vcov(fit), confint(fit),
# AIC(fit), tidy(fit), simulate(fit) first; each method below either
# returns the right substitute or refuses with a directive cli error.
#
# `tidy` is registered on generics::tidy so plain tidy(fit) and
# broom::tidy(fit) both dispatch. static MCMC fits store B as a
# 3-element Tucker list and A = NULL; the accessors branch on
# fit$model == "static" and map the sender Tucker factor B[[1]] into
# the conventional A slot. ALS+bootstrap fits store CIs at t = 2..T (no
# transition at t = 1); tidy() aligns coef [n,n,T] and CI [n,n,T-1] by
# NA-padding t = 1. vcov.dbn errors with a directive message since the
# parameters are array-valued, not a coefficient vector; dbn_se() is
# the exported entry-wise SD accessor. confint.dbn honours parm,
# computes M intervals on MCMC fits, and propagates a non-default
# level to M.
####

#' @title Re-exported broom and forecast generics
#'
#' @description The following objects are re-exported from other packages
#'   so that `dbn::tidy(fit)` and friends work directly on a `dbn` fit
#'   without forcing the user to `library(broom)` / `library(forecast)`.
#'   Method implementations live in this package; the canonical generic
#'   definitions live in the package named in the link.
#'
#' \describe{
#'   \item{[generics::tidy()]}{One row per parameter long-form data frame.}
#'   \item{[generics::glance()]}{One-row model-level summary.}
#'   \item{[generics::augment()]}{Observation-level data frame with `.fitted` / `.resid`.}
#'   \item{[generics::forecast()]}{H-step-ahead forecast wrapper.}
#' }
#'
#' @name reexports
#' @keywords internal
#' @aliases tidy glance augment forecast
NULL

#' @importFrom generics tidy
#' @rdname reexports
#' @author Tosin Salau and Shahryar Minhas
#' @export
generics::tidy

#' @importFrom generics glance
#' @rdname reexports
#' @author Tosin Salau and Shahryar Minhas
#' @export
generics::glance

#' @importFrom generics augment
#' @rdname reexports
#' @author Tosin Salau and Shahryar Minhas
#' @export
generics::augment

#' @importFrom generics forecast
#' @rdname reexports
#' @author Tosin Salau and Shahryar Minhas
#' @export
generics::forecast

# -- internal helpers --------------------------------------------------------

#' @keywords internal
#' @noRd
.dbn_is_static <- function(object) {
	identical(object$model, "static")
}

#' @keywords internal
#' @noRd
.dbn_is_als <- function(object) {
	sused <- object$meta$sampler_used %||% ""
	isTRUE(sused %in% c("als", "als_tv", "als_piecewise"))
}

#' @keywords internal
#' @noRd
.collapse_draws <- function(lst, fun) {
	# stack a list of arrays along a new trailing dim, then reduce with `fun`
	# over that dim. Returns NULL for empty lists.
	if (is.null(lst) || length(lst) == 0L) return(NULL)
	d <- dim(lst[[1]])
	if (is.null(d)) {
		v <- vapply(lst, identity, numeric(length(lst[[1]])))
		return(apply(v, 1, fun))
	}
	arr <- do.call(abind::abind, c(lst, list(along = length(d) + 1L)))
	apply(arr, seq_len(length(d)), fun)
}

#' @keywords internal
#' @noRd
.collapse_draws_pair <- function(lst, probs) {
	if (is.null(lst) || length(lst) == 0L) return(NULL)
	d <- dim(lst[[1]])
	arr <- do.call(abind::abind, c(lst, list(along = length(d) + 1L)))
	lo <- apply(arr, seq_len(length(d)), quantile, probs = probs[1], na.rm = TRUE)
	hi <- apply(arr, seq_len(length(d)), quantile, probs = probs[2], na.rm = TRUE)
	mn <- apply(arr, seq_len(length(d)), mean)
	list(lo = lo, hi = hi, mean = mn)
}

#' @keywords internal
#' @noRd
.dbn_static_is_trivial <- function(arr) {
	# returns TRUE if the Tucker factor is structurally the identity across all
	# kept draws (i.e. the static MCMC sampler does not update this slot).
	# suppresses the misleading "B = I" output from coef.dbn(static_fit).
	if (is.null(arr)) return(TRUE)
	d <- dim(arr)
	if (is.null(d) || length(d) < 3L) return(FALSE)
	core <- seq_len(length(d) - 1L)
	sds <- apply(arr, core, sd)
	if (!all(sds < 1e-12, na.rm = TRUE)) return(FALSE)
	# variance is zero across draws; check if the constant value is the identity.
	first <- if (length(d) == 3L) arr[, , 1] else arr[, , , 1]
	if (!is.matrix(first) || nrow(first) != ncol(first)) return(FALSE)
	isTRUE(all.equal(first, diag(nrow(first)), tolerance = 1e-10))
}

#' @keywords internal
#' @noRd
.dbn_static_means <- function(object) {
	# map Tucker factors into A/B slots for cross-method conceptual consistency:
	#   A = B[[1]]  (sender factor, n_row x n_row)
	#   B = B[[2]]  (receiver factor, n_col x n_col)
	# when B[[2]] is structurally identity (sampler doesn't update it), report
	# NULL so users don't read it as "no receiver effect estimated".
	# M is averaged across stored draws when available.
	arr_mean <- function(arr) {
		if (is.null(arr)) return(NULL)
		d <- dim(arr)
		if (is.null(d) || length(d) < 3L) return(arr)
		apply(arr, seq_len(length(d) - 1L), mean)
	}
	A_mean <- if (is.list(object$B) && length(object$B) >= 1L) arr_mean(object$B[[1]]) else NULL
	B_raw <- if (is.list(object$B) && length(object$B) >= 2L) object$B[[2]] else NULL
	B_mean <- if (!is.null(B_raw) && !.dbn_static_is_trivial(B_raw)) arr_mean(B_raw) else NULL
	Msave  <- object$draws$misc$M
	if (is.list(Msave) && length(Msave) >= 1L && requireNamespace("abind", quietly = TRUE)) {
		M_arr <- do.call(abind::abind, c(Msave, list(along = length(dim(Msave[[1]])) + 1L)))
		M_mean <- apply(M_arr, seq_len(length(dim(M_arr)) - 1L), mean)
	} else {
		M_mean <- object$M
	}
	list(A = A_mean, B = B_mean, M = M_mean)
}

#' @keywords internal
#' @noRd
.dbn_static_sds <- function(object) {
	arr_sd <- function(arr) {
		if (is.null(arr)) return(NULL)
		d <- dim(arr)
		if (is.null(d) || length(d) < 3L) return(NULL)
		apply(arr, seq_len(length(d) - 1L), sd)
	}
	A_sd <- if (is.list(object$B) && length(object$B) >= 1L) arr_sd(object$B[[1]]) else NULL
	B_raw <- if (is.list(object$B) && length(object$B) >= 2L) object$B[[2]] else NULL
	B_sd <- if (!is.null(B_raw) && !.dbn_static_is_trivial(B_raw)) arr_sd(B_raw) else NULL
	Msave <- object$draws$misc$M
	if (is.list(Msave) && length(Msave) >= 2L && requireNamespace("abind", quietly = TRUE)) {
		M_arr <- do.call(abind::abind, c(Msave, list(along = length(dim(Msave[[1]])) + 1L)))
		M_sd <- apply(M_arr, seq_len(length(dim(M_arr)) - 1L), sd)
	} else {
		M_sd <- NULL
	}
	list(A = A_sd, B = B_sd, M = M_sd)
}

#' @keywords internal
#' @noRd
.dbn_static_qs <- function(object, probs) {
	arr_q <- function(arr) {
		if (is.null(arr)) return(NULL)
		d <- dim(arr)
		if (is.null(d) || length(d) < 3L) return(NULL)
		core <- seq_len(length(d) - 1L)
		list(
			lo = apply(arr, core, quantile, probs = probs[1], na.rm = TRUE),
			hi = apply(arr, core, quantile, probs = probs[2], na.rm = TRUE),
			mean = apply(arr, core, mean)
		)
	}
	A_q <- if (is.list(object$B) && length(object$B) >= 1L) arr_q(object$B[[1]]) else NULL
	B_raw <- if (is.list(object$B) && length(object$B) >= 2L) object$B[[2]] else NULL
	B_q <- if (!is.null(B_raw) && !.dbn_static_is_trivial(B_raw)) arr_q(B_raw) else NULL
	Msave <- object$draws$misc$M
	M_q <- NULL
	if (is.list(Msave) && length(Msave) >= 2L && requireNamespace("abind", quietly = TRUE)) {
		M_arr <- do.call(abind::abind, c(Msave, list(along = length(dim(Msave[[1]])) + 1L)))
		core <- seq_len(length(dim(M_arr)) - 1L)
		M_q <- list(
			lo = apply(M_arr, core, quantile, probs = probs[1], na.rm = TRUE),
			hi = apply(M_arr, core, quantile, probs = probs[2], na.rm = TRUE),
			mean = apply(M_arr, core, mean)
		)
	}
	list(A = A_q, B = B_q, M = M_q)
}

#' @keywords internal
#' @noRd
.dbn_resolve_parm <- function(parm, what) {
	# accept `parm = "A"/"B"/"M"` as an alias for `what`; warn and fall
	# back to `what` if `parm` is anything else.
	if (is.null(parm)) return(what)
	if (length(parm) == 1L && is.character(parm) && parm %in% c("all", "A", "B", "M")) {
		return(parm)
	}
	cli::cli_warn(c(
		"{.fun confint.dbn}: argument {.arg parm} = {.val {parm}} is not one of {.val all}, {.val A}, {.val B}, {.val M}.",
		"i" = "Use {.arg what} to select the quantity. Ignoring {.arg parm}."
	))
	what
}

# -- coef.dbn ---------------------------------------------------------------

#' Extract estimated operators and baseline mean from a dbn fit
#'
#' For MCMC fits, returns posterior MEANS of the operator trajectories.
#' For ALS fits (static or TV), returns the point estimate.
#'
#' Static MCMC fits have no transition operator; their `A` slot returns the
#' Tucker sender factor (the n_row x n_row matrix `fit$B[[1]]` averaged across
#' draws) and `B` returns the Tucker receiver factor.
#'
#' @param object A `dbn` fit object.
#' @param what Character. Which quantity to extract. One of \code{"all"}
#'   (default; named list with A, B, M), \code{"A"}, \code{"B"}, \code{"M"}.
#' @param ... Unused.
#' @return If `what = "all"`, a named list with elements `A`, `B`, `M`.
#'   Otherwise a single array. Operator arrays are
#'   `[n_row, n_row, T]` (A) and `[n_col, n_col, T]` (B); the baseline `M`
#'   is `[n_row, n_col, p]`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
coef.dbn <- function(object, what = c("all", "A", "B", "M"), ...) {
	what <- match.arg(what)
	if (.dbn_is_static(object)) {
		out <- .dbn_static_means(object)
		return(switch(what, all = out, A = out$A, B = out$B, M = out$M))
	}
	# dynamic path: A/B are either lists of cubes (one per draw) or single arrays.
	if (is.list(object$A)) {
		A_mean <- Reduce(`+`, object$A) / length(object$A)
	} else if (is.array(object$A)) {
		A_mean <- if (length(dim(object$A)) >= 4L)
			apply(object$A, seq_len(length(dim(object$A)) - 1L), mean) else object$A
	} else {
		A_mean <- NULL
	}
	if (is.list(object$B)) {
		B_mean <- Reduce(`+`, object$B) / length(object$B)
	} else if (is.array(object$B)) {
		B_mean <- if (length(dim(object$B)) >= 4L)
			apply(object$B, seq_len(length(dim(object$B)) - 1L), mean) else object$B
	} else {
		B_mean <- NULL
	}
	M_mean <- if (is.array(object$M) && length(dim(object$M)) == 4L)
		apply(object$M, 1:3, mean) else object$M
	switch(what,
		all = list(A = A_mean, B = B_mean, M = M_mean),
		A   = A_mean,
		B   = B_mean,
		M   = M_mean)
}

# -- dbn_se -----------------------------------------------------------------

#' Entry-wise standard deviations for a dbn fit
#'
#' Returns posterior (MCMC) or bootstrap (ALS) standard deviations as arrays
#' shaped like the operators / baseline themselves. Use this when you want
#' "how uncertain is each entry of A_t" -- this is NOT a covariance matrix.
#' Base R's `vcov()` generic expects a covariance matrix; calling
#' `vcov(fit)` instead errors directively and points you here.
#'
#' @param object A `dbn` fit.
#' @param what Character: which quantity ("all" / "A" / "B" / "M"). Default "all".
#' @param ... Unused.
#' @return Named list (or array) of standard deviations with the same shape
#'   as the corresponding operator / baseline.
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_se <- function(object, what = c("all", "A", "B", "M"), ...) {
	what <- match.arg(what)
	if (.dbn_is_static(object)) {
		out <- .dbn_static_sds(object)
		return(switch(what, all = out, A = out$A, B = out$B, M = out$M))
	}
	# bootstrap-equipped fit: use stored SEs.
	if (!is.null(object$bootstrap) && inherits(object$bootstrap, "dbn_boot")) {
		out <- list(
			A = if (!is.null(object$bootstrap$se_A)) object$bootstrap$se_A else NULL,
			B = if (!is.null(object$bootstrap$se_B)) object$bootstrap$se_B else NULL,
			M = if (!is.null(object$bootstrap$se_M)) object$bootstrap$se_M else NULL
		)
		return(switch(what, all = out, A = out$A, B = out$B, M = out$M))
	}
	# MCMC: posterior SDs.
	if (is.list(object$A) && length(object$A) >= 2L) {
		if (!requireNamespace("abind", quietly = TRUE)) {
			cli::cli_warn(c(
				"{.fun dbn_se} on an MCMC fit needs the {.pkg abind} package.",
				"i" = "Install with {.code install.packages(\"abind\")}."
			))
			return(NULL)
		}
		A_sd <- .collapse_draws(object$A, sd)
		B_sd <- .collapse_draws(object$B, sd)
		M_sd <- if (is.array(object$M) && length(dim(object$M)) == 4L)
			apply(object$M, 1:3, sd) else NULL
		out <- list(A = A_sd, B = B_sd, M = M_sd)
		return(switch(what, all = out, A = out$A, B = out$B, M = out$M))
	}
	cli::cli_inform(c(
		"i" = "No uncertainty available for this fit ({.code dbn_se} returns {.code NULL}).",
		"i" = "Refit with {.code dbn(..., model = \"dynamic\")} for posterior draws, or {.code dbn_als(..., bootstrap = N)} for bootstrap SEs."
	))
	NULL
}

# -- vcov.dbn ---------------------------------------------------------------

#' Variance-covariance method for a dbn fit (directive error)
#'
#' The base R `vcov()` generic returns a *covariance matrix* over a flat
#' coefficient vector. A dbn fit's parameters are high-dimensional arrays
#' (A_t, B_t, M); a coefficient-vector covariance is not the right summary.
#' This method errors with a directive message pointing at the right
#' function:
#'
#'   - `dbn_se(fit)`  -- entry-wise posterior / bootstrap standard deviations
#'   - `confint(fit)` -- per-entry credible / bootstrap intervals
#'   - `vcov(fit$pars)` if you want the covariance over scalar parameters
#'     (sigma^2, tau_A^2, ...)
#'
#' @param object A `dbn` fit.
#' @param ... Unused.
#' @return Does not return; calls `cli::cli_abort`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
vcov.dbn <- function(object, ...) {
	cli::cli_abort(c(
		"{.fun vcov} on a {.cls dbn} fit is not defined: operators are array-valued, not a flat coefficient vector.",
		"i" = "For entry-wise standard deviations, use {.code dbn_se(fit)}.",
		"i" = "For credible / bootstrap intervals, use {.code confint(fit)}."
	))
}

# -- confint.dbn ------------------------------------------------------------

#' Credible / confidence intervals for a dbn fit
#'
#' For MCMC fits: posterior quantiles across draws (including M). For
#' ALS-with-bootstrap fits: bootstrap percentile intervals. For
#' ALS-without-bootstrap: refuses with a hint about running the bootstrap.
#'
#' @param object A `dbn` fit object.
#' @param parm Optional. If supplied with one of \code{"A"}, \code{"B"},
#'   \code{"M"}, or \code{"all"}, overrides `what`. Otherwise warns and is
#'   ignored.
#' @param level Coverage level (default 0.95).
#' @param what Character: "all" / "A" / "B" / "M".
#' @param ... Unused.
#' @return Named list with `lo`/`hi` arrays per quantity (and `mean`).
#' @author Tosin Salau and Shahryar Minhas
#' @export
confint.dbn <- function(object, parm = NULL, level = 0.95,
                         what = c("all", "A", "B", "M"), ...) {
	what <- match.arg(what)
	what <- .dbn_resolve_parm(parm, what)
	if (level <= 0 || level >= 1) {
		cli::cli_abort("{.arg level} must be in (0, 1); got {.val {level}}.")
	}
	alpha <- (1 - level) / 2
	probs <- c(alpha, 1 - alpha)

	if (.dbn_is_static(object)) {
		out <- .dbn_static_qs(object, probs)
		# gaussian static fits use a closed-form M update (no MCMC noise on M),
		# so the posterior on M is a point mass and intervals are degenerate.
		# warn so the user does not read this as a real "tight" interval.
		if (!is.null(out$M) && all(abs(out$M$hi - out$M$lo) < 1e-12, na.rm = TRUE)) {
			fam <- if (is.list(object$family)) object$family$name else object$family
			cli::cli_warn(c(
				"M posterior is degenerate (lo == hi) for this static fit.",
				"i" = if (identical(fam, "gaussian"))
					"For a static Gaussian fit, M is updated deterministically (sample-mean of Z); there is no MCMC variability on M. The interval is a point estimate."
					else
					"All stored draws of M are identical -- check that {.code odens} > 0 and that {.code nscan} produced multiple kept iterations."
			))
		}
		return(switch(what, all = out, A = out$A, B = out$B, M = out$M))
	}

	# bootstrap fit
	if (!is.null(object$bootstrap) && inherits(object$bootstrap, "dbn_boot")) {
		.quant <- function(arr, q) {
			if (is.null(arr)) return(NULL)
			d <- dim(arr)
			if (is.null(d)) return(quantile(arr, probs = q, na.rm = TRUE))
			if (length(d) == 2L) {
				apply(arr, 2, quantile, probs = q, na.rm = TRUE)
			} else if (length(d) >= 3L) {
				apply(arr, 2:length(d), quantile, probs = q, na.rm = TRUE)
			} else NULL
		}
		if (abs(level - 0.95) < 1e-6) {
			out <- list(
				A = list(lo = object$bootstrap$ci_A_lo, hi = object$bootstrap$ci_A_hi,
				         mean = object$bootstrap$point_est_A),
				B = list(lo = object$bootstrap$ci_B_lo, hi = object$bootstrap$ci_B_hi,
				         mean = object$bootstrap$point_est_B),
				M = if (!is.null(object$bootstrap$ci_M_lo))
					list(lo = object$bootstrap$ci_M_lo, hi = object$bootstrap$ci_M_hi,
					     mean = object$bootstrap$point_est_M) else NULL
			)
		} else {
			A_lo <- .quant(object$bootstrap$coefs_A, probs[1])
			A_hi <- .quant(object$bootstrap$coefs_A, probs[2])
			B_lo <- .quant(object$bootstrap$coefs_B, probs[1])
			B_hi <- .quant(object$bootstrap$coefs_B, probs[2])
			M_lo <- .quant(object$bootstrap$coefs_M, probs[1])
			M_hi <- .quant(object$bootstrap$coefs_M, probs[2])
			out <- list(
				A = list(lo = A_lo, hi = A_hi, mean = object$bootstrap$point_est_A),
				B = list(lo = B_lo, hi = B_hi, mean = object$bootstrap$point_est_B),
				M = if (!is.null(M_lo))
					list(lo = M_lo, hi = M_hi, mean = object$bootstrap$point_est_M) else NULL
			)
		}
	} else if (is.list(object$A) && length(object$A) >= 2L) {
		# MCMC posterior quantiles
		if (!requireNamespace("abind", quietly = TRUE)) {
			cli::cli_warn("{.fun confint.dbn} on an MCMC fit needs {.pkg abind}; install it for posterior quantiles.")
			return(NULL)
		}
		A_q <- .collapse_draws_pair(object$A, probs)
		B_q <- .collapse_draws_pair(object$B, probs)
		# M is stored as a 4D array [n_row, n_col, p, n_keep]
		M_q <- NULL
		if (is.array(object$M) && length(dim(object$M)) == 4L) {
			M_lo <- apply(object$M, 1:3, quantile, probs = probs[1], na.rm = TRUE)
			M_hi <- apply(object$M, 1:3, quantile, probs = probs[2], na.rm = TRUE)
			M_mn <- apply(object$M, 1:3, mean)
			M_q <- list(lo = M_lo, hi = M_hi, mean = M_mn)
		}
		out <- list(A = A_q, B = B_q, M = M_q)
	} else if (.dbn_is_als(object)) {
		# ALS fit without a stored bootstrap: run one on demand so the user
		# gets actual intervals instead of an abort. Default to 200 reps
		# (the recommended count for stable percentile intervals); user can
		# pre-empt with `dbn_als(..., bootstrap = N)` to control N.
		boot_n <- getOption("dbn.confint_auto_bootstrap", 200L)
		cli::cli_inform(c(
			"i" = "ALS fit has no stored uncertainty; running {boot_n}-replicate bootstrap to compute intervals.",
			"i" = "Pre-bootstrap with {.code dbn_als(..., bootstrap = {boot_n})} to skip this step, or set {.code options(dbn.confint_auto_bootstrap = N)} to change the default."
		))
		su <- object$meta$sampler_used %||% ""
		boot <- if (su %in% c("als_tv", "als_piecewise")) {
			.dbn_als_tv_bootstrap(object, R = boot_n, type = "fixed",
				verbose = FALSE)
		} else {
			.dbn_als_bootstrap(object, R = boot_n, type = "block",
				verbose = FALSE)
		}
		# attach to the fit-local object and recurse via a fresh confint call
		# (which will hit the bootstrap-fit branch above). route through the
		# shared attacher so confint also degrades cleanly if the auto
		# bootstrap returned no valid replicates.
		object <- .dbn_attach_bootstrap(object, boot,
			label = sprintf("Point estimate plus %d-replicate auto-bootstrap CIs.", boot$n_valid %||% 0L),
			source = "bootstrap")
		if (!isTRUE(object$meta$uncertainty_available)) {
			cli::cli_abort(c(
				"Auto-bootstrap returned no valid replicates; cannot compute confidence intervals.",
				"i" = "Pre-bootstrap with {.code dbn_als(..., bootstrap = N)} and inspect per-rep errors."
			))
		}
		return(confint(object, parm = parm, level = level, what = what, ...))
	} else {
		cli::cli_abort(c(
			"This fit has no uncertainty quantification.",
			"i" = "Refit with {.code dbn(..., model = \"dynamic\")} for MCMC posterior intervals, or {.code dbn_als(..., bootstrap = N)} for bootstrap intervals."
		))
	}
	switch(what, all = out, A = out$A, B = out$B, M = out$M)
}

# -- logLik.dbn -------------------------------------------------------------

#' Log-likelihood (approximate) for a dbn fit
#'
#' Returns the Gaussian log-likelihood of the latent state under the fitted
#' operators, evaluated at the point estimate (MCMC: posterior mean; ALS:
#' point estimate). This is a *coarse* model-comparison object: it ignores
#' MCMC uncertainty on operator parameters and uses the model's process
#' variance, not an observation-noise variance.
#'
#' @param object A `dbn` fit object.
#' @param ... Unused.
#' @return Object of class `logLik` with the usual `df` and `nobs` attrs.
#' @author Tosin Salau and Shahryar Minhas
#' @export
logLik.dbn <- function(object, ...) {
	# short-circuit static fits before the A/B null check: coef(static)$B
	# is NULL when the second Tucker factor is structurally identity, and
	# static fits have no transition operator anyway. point users at the
	# DIC stored on the fit.
	if (.dbn_is_static(object)) {
		cli::cli_warn(c(
			"{.fun logLik.dbn} is defined only for dynamic / TV / piecewise fits.",
			"i" = "Static fits have no transition operator; use the DIC stored in {.code fit$diagnostics$dic}."
		))
		return(structure(NA_real_, df = NA, nobs = NA, class = "logLik"))
	}
	# HMM fits: A is per-regime (n x n x R), not per-time, and the
	# log-likelihood needs the regime sequence S_t to know which A to
	# apply at each t. Return NA with a directive warning rather than
	# silently mis-computing (which is what happened when the dynamic
	# code path applied A[min(t, R)] regardless of regime).
	if (identical(object$model, "hmm")) {
		cli::cli_warn(c(
			"{.fun logLik.dbn} is not yet implemented for HMM fits.",
			"i" = "HMM models have per-regime (not per-time) operators; computing logLik requires marginalizing over the regime path.",
			"i" = "Use {.code fit$diagnostics$dic} if available, or {.fun compute_waic_dbn} for model comparison."
		))
		return(structure(NA_real_, df = NA, nobs = NA, class = "logLik"))
	}
	pe <- coef(object, "all")
	A <- pe$A; B <- pe$B; M <- pe$M
	if (is.null(A) || is.null(B)) {
		cli::cli_warn(c(
			"{.fun logLik.dbn} could not extract A/B point estimates from this fit.",
			"i" = "Returning NA. Inspect {.code coef(fit)} to see what's missing."
		))
		return(structure(NA_real_, df = NA, nobs = NA, class = "logLik"))
	}
	dims <- object$dims
	n_row <- dims$n_row; n_col <- dims$n_col; p <- dims$p; Tt <- dims$Tt
	Y <- object$Y
	if (is.null(Y)) cli::cli_abort("{.fun logLik.dbn} needs {.code fit$Y} to be present.")
	Y_f <- Y; Y_f[is.na(Y_f)] <- 0
	sigma2 <- if (is.numeric(object$sigma2)) mean(object$sigma2, na.rm = TRUE)
		else 1
	Omega <- array(1, dim = c(n_row, n_col, Tt))
	if (n_row == n_col) for (t in 1:Tt) diag(Omega[, , t]) <- 0
	for (r in 1:p) for (t in 1:Tt) Omega[, , t][is.na(Y[, , r, t])] <- 0

	n_obs_total <- 0; sse <- 0
	for (r in 1:p) {
		Phi_r <- array(0, c(n_row, n_col, Tt))
		for (t in 1:Tt) {
			Phi_r[, , t] <- Y_f[, , r, t] - M[, , r]
			Phi_r[, , t][Omega[, , t] == 0] <- 0
		}
		get_A <- function(t) if (length(dim(A)) == 3L) A[, , min(t, dim(A)[3])] else A
		get_B <- function(t) if (length(dim(B)) == 3L) B[, , min(t, dim(B)[3])] else B
		for (t in 2:Tt) {
			pred <- get_A(t) %*% Phi_r[, , t - 1L] %*% t(get_B(t))
			resid <- Omega[, , t] * (Phi_r[, , t] - pred)
			sse <- sse + sum(resid^2)
			n_obs_total <- n_obs_total + sum(Omega[, , t])
		}
	}
	ll <- -0.5 * n_obs_total * log(2 * pi * sigma2) - sse / (2 * sigma2)
	# degrees of freedom: count actual parameter slices, not the static-shape
	# default. For dynamic / TV fits we count T_eff = (Tt - 1) slices of A and
	# B (t=1 is unestimated identity), plus M (one [n_row, n_col, p] slice)
	# and the 3 variance scalars (sigma2, tauA2, tauB2). For piecewise fits we
	# count K slices of A and B (one per block).
	get_n_slices <- function(arr) {
		if (is.null(arr)) return(1L)
		d <- dim(arr)
		if (is.null(d) || length(d) < 3L) return(1L)
		d[length(d)]
	}
	su <- object$meta$sampler_used %||% ""
	T_eff <- if (identical(su, "als_piecewise")) {
		# K block count: try to recover from coef shape
		get_n_slices(coef(object, "A"))
	} else if (object$model %in% c("dynamic")) {
		max(1L, get_n_slices(coef(object, "A")) - 1L)
	} else {
		1L
	}
	# under the symmetric specification A_t = B_t and A_t is upper-triangle
	# only (n*(n+1)/2 free entries), so the asymmetric formula
	# T_eff*(n_row^2 + n_col^2) would over-count by ~4x and bias AIC/BIC
	# against symmetric fits.
	is_sym <- isTRUE(object$dims$is_symmetric)
	per_slice_AB <- if (is_sym) {
		n_row * (n_row + 1L) / 2L
	} else {
		n_row * n_row + n_col * n_col
	}
	df_approx <- T_eff * per_slice_AB +
	             n_row * n_col * p + 3L
	structure(ll, df = df_approx, nobs = n_obs_total, class = "logLik")
}

####
#' Pointwise log-likelihood matrix for use with `loo::loo()`
#'
#' @description Computes the per-draw, per-observation Gaussian
#'   log-likelihood matrix expected by [loo::loo()] (PSIS-LOO) and
#'   [loo::waic()]. Returned matrix has shape `[n_draws, n_obs]` where
#'   `n_obs` is the count of non-NA dyad-time observations the model
#'   was fit on (off-diagonal cells only on unipartite slices).
#'
#'   Only implemented for Gaussian dynamic/piecewise fits where draws of
#'   the operator are stored. Returns NA with a directive on other
#'   families / models, matching the [logLik.dbn()] policy.
#'
#' @param object A `dbn` fit.
#' @param ... Unused.
#' @return A numeric matrix with `n_draws` rows and `n_obs` columns,
#'   suitable as input to `loo::loo(log_lik(fit))`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
log_lik.dbn <- function(object, ...) {
	if (!inherits(object, "dbn"))
		cli::cli_abort("{.arg object} must be a {.cls dbn} object.")
	fam <- if (is.list(object$family)) object$family$name else object$family
	if (!identical(fam, "gaussian")) {
		cli::cli_warn(c(
			"{.fun log_lik.dbn} is implemented for the Gaussian family only.",
			"i" = "Got family {.val {fam}}; returning NA matrix.",
			"i" = "For non-Gaussian fits, use {.fun compute_waic_dbn} which already implements the appropriate working-likelihood approximation."
		))
		return(matrix(NA_real_, nrow = 1L, ncol = 1L))
	}
	if (!object$model %in% c("dynamic", "piecewise"))
		cli::cli_abort(c(
			"{.fun log_lik.dbn} requires {.code model %in% c(\"dynamic\", \"piecewise\")}.",
			"x" = "Got {.val {object$model}}."
		))
	A_list <- object$A; B_list <- object$B
	if (!is.list(A_list) || !is.list(B_list) || length(A_list) == 0L)
		cli::cli_abort("Could not find posterior draws of {.code A} and {.code B}.")
	# observation noise: prefer sigma2_obs (Gaussian latent likelihood),
	# fall back to sigma2 (process noise) only if sigma2_obs absent.
	sigma2_obs_vec <- if (is.numeric(object$sigma2_obs)) as.numeric(object$sigma2_obs)
	                  else if (is.numeric(object$sigma2)) as.numeric(object$sigma2)
	                  else cli::cli_abort("Neither {.code sigma2_obs} nor {.code sigma2} is numeric on this fit.")
	M_arr <- if (is.array(object$M) && length(dim(object$M)) == 4L) {
		apply(object$M, 1:3, mean)
	} else {
		object$M
	}
	if (is.null(M_arr))
		cli::cli_abort("Could not extract baseline mean {.code M} from this fit.")
	dims <- object$dims
	n_row <- dims$n_row; n_col <- dims$n_col; p <- dims$p; Tt <- dims$Tt
	Y <- object$Y
	if (is.null(Y))
		cli::cli_abort("{.fun log_lik.dbn} requires {.code fit$Y} to be present.")
	# build observation mask: off-diagonal on unipartite, all cells on
	# bipartite. Replicate the convention used by logLik.dbn.
	mask <- array(TRUE, dim = c(n_row, n_col, p, Tt))
	if (n_row == n_col)
		for (r in 1:p) for (t in 1:Tt) diag(mask[, , r, t]) <- FALSE
	for (r in 1:p) for (t in 1:Tt) mask[, , r, t][is.na(Y[, , r, t])] <- FALSE
	# only t >= 2 contributes (t = 1 has no predictor for Gaussian latent)
	if (Tt >= 2L) for (r in 1:p) mask[, , r, 1L] <- FALSE
	obs_idx <- which(mask)
	n_obs <- length(obs_idx)
	if (n_obs == 0L)
		cli::cli_abort("No usable observations found; cannot compute log-likelihood.")
	n_draws <- length(A_list)
	# recycle sigma2_obs to match draws
	sig2 <- rep(sigma2_obs_vec, length.out = n_draws)
	ll <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
	# helper: extract A_m at time t (supports dynamic 3D cube and
	# piecewise 2D-per-block-list).
	get_At <- function(Am, t) {
		if (is.array(Am) && length(dim(Am)) == 3L) Am[, , min(t, dim(Am)[3L])] else Am
	}
	for (m in seq_len(n_draws)) {
		Am <- A_list[[m]]; Bm <- B_list[[m]]
		# predicted Y at each (i, j, r, t): M + A_t (Y_{t-1} - M) B_t^T
		pred <- array(NA_real_, dim = c(n_row, n_col, p, Tt))
		Y_filled <- Y; Y_filled[is.na(Y_filled)] <- 0
		for (r in 1:p) {
			for (t in 2:Tt) {
				At <- get_At(Am, t); Bt <- get_At(Bm, t)
				dev_prev <- Y_filled[, , r, t - 1L] - M_arr[, , r]
				pred[, , r, t] <- M_arr[, , r] + At %*% dev_prev %*% t(Bt)
			}
		}
		resid <- as.vector(Y - pred)
		# Gaussian log-density at masked entries only
		ll[m, ] <- stats::dnorm(as.vector(Y)[obs_idx],
		                        mean = as.vector(pred)[obs_idx],
		                        sd = sqrt(sig2[m]), log = TRUE)
	}
	ll
}

# -- tidy.dbn ---------------------------------------------------------------

#' Tidy data-frame view of a dbn fit
#'
#' Returns a long data frame with columns `parameter` (A / B / M), `t`, `i`,
#' `j`, `estimate`, `std.error`, `conf.low`, `conf.high`. Column names follow
#' the broom convention so downstream tooling (`ggdist`, `dotwhisker`,
#' `modelsummary`) dispatches correctly.
#'
#' Each row is one entry of A_{t,ij}, B_{t,ij}, or M_{ij}. For M, `t = NA`
#' (the baseline is time-invariant).
#'
#' `std.error` and CI columns are NA-padded at t=1 for ALS-with-bootstrap fits
#' (the bootstrap does not estimate the transition at t=1; only t=2..T).
#'
#' @param x A `dbn` fit object.
#' @param level Coverage level for intervals (default 0.95).
#' @param ... Unused.
#' @return A data frame with broom-conventional columns.
#' @author Tosin Salau and Shahryar Minhas
#' @export
tidy.dbn <- function(x, level = 0.95, ...) {
	pe <- coef(x, "all")
	# decide whether to ask for uncertainty up front, so unexpected errors
	# below are real bugs (not silently swallowed). For ALS fits without
	# bootstrap we explicitly know there is none -- skip the call and fill NAs.
	has_uncertainty <- isTRUE(x$meta$uncertainty_available) ||
	                   (is.list(x$A) && length(x$A) >= 2L) ||
	                   !is.null(x$bootstrap)
	if (has_uncertainty) {
		se <- dbn_se(x, "all")
		ci <- confint(x, level = level)
	} else {
		se <- NULL
		ci <- NULL
	}

	# align an uncertainty array (se / lo / hi) to the estimate shape:
	# - ALS-bootstrap case: coef is [n,n,T] but uncertainty is [n,n,T-1]
	#   (no transition at t = 1), so NA-pad t = 1.
	# - static-SE / tiled-coef case: uncertainty is [n_row, n_col] but the
	#   estimate is [n_row, n_col, T], so tile the SE across time rather
	#   than dropping it.
	align_to <- function(unc_arr, est_arr, n_row_df) {
		if (is.null(unc_arr)) return(rep(NA_real_, n_row_df))
		d_est <- dim(est_arr)
		d_unc <- dim(unc_arr)
		# length already matches: simple as.vector.
		if (length(as.vector(unc_arr)) == n_row_df) {
			return(as.vector(unc_arr))
		}
		# 3D shape mismatch with T-1 vs T trailing dim: pad with NA at t=1.
		if (!is.null(d_est) && !is.null(d_unc) &&
		    length(d_est) == 3L && length(d_unc) == 3L &&
		    d_est[1] == d_unc[1] && d_est[2] == d_unc[2] &&
		    d_unc[3] == d_est[3] - 1L) {
			pad <- array(NA_real_, dim = c(d_unc[1], d_unc[2], 1L))
			padded <- array(c(as.vector(pad), as.vector(unc_arr)),
			                dim = c(d_unc[1], d_unc[2], d_unc[3] + 1L))
			return(as.vector(padded))
		}
		# 2D unc / 3D estimate: broadcast (tile) the SE across time. This
		# fires for static-ALS bootstraps where coef is tiled-across-time
		# but the bootstrap SE is the (time-invariant) static SE.
		if (!is.null(d_est) && !is.null(d_unc) &&
		    length(d_est) == 3L && length(d_unc) == 2L &&
		    d_est[1] == d_unc[1] && d_est[2] == d_unc[2]) {
			tiled <- array(rep(as.vector(unc_arr), d_est[3]),
			               dim = c(d_unc[1], d_unc[2], d_est[3]))
			return(as.vector(tiled))
		}
		# anything else: bail with NAs.
		rep(NA_real_, n_row_df)
	}

	make_long <- function(estimate, se_arr, ci_pair, label) {
		if (is.null(estimate)) return(NULL)
		d <- dim(estimate)
		# handle 3D M as [n_row, n_col, p] by stacking relations.
		if (length(d) == 3L && label == "M") {
			df <- expand.grid(i = seq_len(d[1]), j = seq_len(d[2]),
			                  rel = seq_len(d[3]), KEEP.OUT.ATTRS = FALSE)
			df$t <- NA_integer_
			df$estimate <- as.vector(estimate)
		} else if (length(d) == 2L) {
			df <- expand.grid(i = seq_len(d[1]), j = seq_len(d[2]),
			                  KEEP.OUT.ATTRS = FALSE)
			df$t <- NA_integer_
			df$rel <- NA_integer_
			df$estimate <- as.vector(estimate)
		} else if (length(d) == 3L) {
			df <- expand.grid(i = seq_len(d[1]), j = seq_len(d[2]),
			                  t = seq_len(d[3]), KEEP.OUT.ATTRS = FALSE)
			df$rel <- NA_integer_
			df$estimate <- as.vector(estimate)
		} else return(NULL)
		df$std.error <- align_to(se_arr, estimate, nrow(df))
		df$conf.low  <- if (!is.null(ci_pair)) align_to(ci_pair$lo, estimate, nrow(df)) else NA_real_
		df$conf.high <- if (!is.null(ci_pair)) align_to(ci_pair$hi, estimate, nrow(df)) else NA_real_
		df$parameter <- label
		df[, c("parameter", "t", "rel", "i", "j", "estimate",
		       "std.error", "conf.low", "conf.high")]
	}

	df_A <- make_long(pe$A, se$A, ci$A, "A")
	df_B <- make_long(pe$B, se$B, ci$B, "B")
	df_M <- make_long(pe$M, se$M, ci$M, "M")
	# drop the t = 1 anchor rows from A and B for dynamic / piecewise fits:
	# A_0 = B_0 = I is the structural identity anchor, not an estimated
	# parameter, so the rows show estimate = 1, std.error = 0, ci = (1, 1)
	# which reads as a numeric bug.
	dyn_like <- identical(x$model, "dynamic") || identical(x$model, "piecewise") ||
		isTRUE(x$meta$piecewise_expanded) || isTRUE(x$meta$bootstrap_expanded)
	if (dyn_like) {
		if (!is.null(df_A)) df_A <- df_A[is.na(df_A$t) | df_A$t > 1L, , drop = FALSE]
		if (!is.null(df_B)) df_B <- df_B[is.na(df_B$t) | df_B$t > 1L, , drop = FALSE]
	}
	out <- do.call(rbind, list(df_A, df_B, df_M))
	rownames(out) <- NULL
	# propagate time-axis labels (e.g. years) when fit$Y carries them.
	# tidy()'s `t` column is a post-thin integer 1..length(time_kept); the
	# matching label is dimnames(fit$Y)[[4]][time_kept[t]]. M rows have t = NA.
	if (!is.null(x$Y) && length(dim(x$Y)) == 4L) {
		t_names <- dimnames(x$Y)[[4]]
		if (!is.null(t_names)) {
			tk <- x$time_kept
			t_idx <- out$t
			if (!is.null(tk) && length(tk) >= max(t_idx, na.rm = TRUE)) {
				out$time_label <- ifelse(is.na(t_idx), NA_character_, t_names[tk[t_idx]])
			} else if (length(t_names) >= max(t_idx, na.rm = TRUE)) {
				out$time_label <- ifelse(is.na(t_idx), NA_character_, t_names[t_idx])
			}
		}
	}
	out
}

# -- actor_embedding -------------------------------------------------------

#' Per-actor coupling-profile embedding for a dbn fit
#'
#' Extracts an actor-by-time matrix of scalar coupling scores summarising
#' how each row-actor's sender (or receiver) role evolves over the fit's
#' time horizon. The default scalar is the per-actor row norm of the
#' posterior-mean transition operator: `||A_t[i, ]||` for `what =
#' "send"`, `||B_t[i, ]||` for `what = "receive"`. These are the
#' coupling-strength scores that paper with `dbn_coupling_rank_probs()`
#' and `compute_irf()`.
#'
#' The intended use is "give me a tidy matrix I can hand to xgboost,
#' write to CSV for Python, or cluster with kmeans" for graph-mining
#' and downstream-ML pipelines.
#'
#' @param fit A `dbn` fit object (dynamic / piecewise / HMM). Static
#'   fits return a single-column matrix with the static operator's row
#'   norm.
#' @param what Character: which role to summarise. `"send"` (default)
#'   pulls `A_t`; `"receive"` pulls `B_t`. Ignored on symmetric fits
#'   (B = A by construction).
#' @param summary Character: `"mean"` (default) returns posterior-mean
#'   row norms; `"median"` returns posterior-median row norms (more
#'   robust under heavy-tailed posteriors).
#' @param fun Optional user-supplied scalar function of an
#'   `n_row x n_col` matrix slice, e.g. `function(M) max(abs(M))`. When
#'   provided, overrides the default row-norm reduction. Should return
#'   a numeric vector of length `n_row`.
#' @return Numeric matrix. For static fits a `n_row x 1` matrix (the
#'   row norm of the single operator). For dynamic / piecewise fits a
#'   `n_row x (Tt - 1)` matrix: the `t = 1` column is dropped because
#'   `A_0` is pinned at the identity anchor (every row norm is exactly
#'   1, which isn't an estimated parameter). Under auto-time-thinning,
#'   columns are named from `dimnames(fit$Y)[[4]][fit$time_kept[-1]]`
#'   rather than `[-1]` of the full time axis, so the labels match the
#'   stored slices.
#' @seealso [coef.dbn()], [dbn_coupling_rank_probs()], [compute_irf()]
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 8, seed = 1)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'            nscan = 200, burn = 100, verbose = FALSE)
#' emb <- actor_embedding(fit, what = "send")
#' dim(emb)  # n_row x (Tt - 1) -- t=1 anchor column is dropped
#' }
actor_embedding <- function(fit, what = c("send", "receive"),
                             summary = c("mean", "median"), fun = NULL) {
	if (!inherits(fit, "dbn")) {
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	}
	what <- match.arg(what)
	summary <- match.arg(summary)
	pe <- coef(fit, "all")
	# pick the operator array to summarise.
	op <- if (identical(what, "send")) pe$A else pe$B
	if (is.null(op)) {
		# symmetric or static fits may store only A; reuse it.
		op <- pe$A
		if (is.null(op)) {
			cli::cli_abort(c(
				"This fit has neither {.code A} nor {.code B} operators to extract embeddings from.",
				"i" = "Got model {.val {fit$model}} with sampler {.val {fit$meta$sampler_used %||% NA}}."
			))
		}
	}
	d <- dim(op)
	if (is.null(d) || !(length(d) %in% c(2L, 3L))) {
		cli::cli_abort(c(
			"Unexpected operator shape from {.fun coef}.",
			"x" = "Got dim {.val {paste(d, collapse = 'x')}}; expected [n_row, n_col] or [n_row, n_col, Tt]."
		))
	}
	# optional summary-of-posterior-summaries: if the user asked for the
	# median we recompute from raw draws when available (otherwise fall back
	# to mean = posterior-mean operator).
	if (identical(summary, "median") && is.list(fit$A) && length(fit$A) >= 2L) {
		# stack draws and compute per-cell median to get a median operator.
		dr <- fit$A
		if (identical(what, "receive") && is.list(fit$B) && length(fit$B) >= 2L) dr <- fit$B
		op_arr <- array(unlist(dr), dim = c(dim(dr[[1]]), length(dr)))
		op <- apply(op_arr, seq_len(length(dim(op_arr)) - 1L), median)
		d <- dim(op)
	}
	# reducer: default is the L2 row norm; user may override.
	if (is.null(fun)) {
		fun <- function(M) sqrt(rowSums(M^2, na.rm = TRUE))
	}
	# build the embedding [n_row, Tt].
	if (length(d) == 2L) {
		emb <- matrix(fun(op), ncol = 1L)
		colnames(emb) <- "t=1"
	} else {
		emb <- vapply(seq_len(d[3]), function(t) as.numeric(fun(op[, , t])),
		              numeric(d[1]))
		colnames(emb) <- paste0("t=", seq_len(d[3]))
		# drop the t = 1 column on dynamic / piecewise fits: A_0 is pinned at
		# the identity anchor so the row norm is exactly 1 for every actor,
		# which looks degenerate and isn't an estimated parameter. static
		# fits keep the single column they have.
		if (d[3] >= 2L && identical(fit$model, "dynamic") ||
		    identical(fit$model, "piecewise") ||
		    !is.null(fit$meta$piecewise_expanded)) {
			emb <- emb[, -1L, drop = FALSE]
		}
	}
	# preserve actor row names. coef() strips operator dimnames by
	# default, so fall back to fit$Y row names for labeled-subset
	# workflows like emb[c("a", "b"), ].
	op_names <- if (!is.null(dimnames(op))) dimnames(op)[[1]] else NULL
	if (is.null(op_names) && !is.null(fit$Y)) {
		dn <- dimnames(fit$Y)
		if (!is.null(dn) && length(dn) >= 1L) op_names <- dn[[1]]
	}
	if (!is.null(op_names) && length(op_names) == nrow(emb)) {
		rownames(emb) <- op_names
	}
	# propagate time-axis column names from the input if present.
	# under auto-time-thinning, fit$time_kept indexes which of the original
	# Tt time points are stored; pull labels from that subset (skipping t=1
	# for dynamic / piecewise where the embedding drops the identity anchor).
	if (!is.null(fit$Y) && length(dim(fit$Y)) == 4L) {
		t_names <- dimnames(fit$Y)[[4]]
		drop_anchor <- d[3] >= 2L && identical(fit$model, "dynamic") ||
		    identical(fit$model, "piecewise") ||
		    !is.null(fit$meta$piecewise_expanded)
		t_labels <- if (!is.null(fit$time_kept)) {
			if (isTRUE(drop_anchor) && length(fit$time_kept) >= 2L)
				t_names[fit$time_kept[-1L]] else t_names[fit$time_kept]
		} else {
			if (isTRUE(drop_anchor)) t_names[-1L] else t_names
		}
		if (!is.null(t_labels) && length(t_labels) == ncol(emb)) {
			colnames(emb) <- t_labels
		} else if (!is.null(t_names) && ncol(emb) == length(t_names) - 1L) {
			colnames(emb) <- t_names[-1L]
		} else if (!is.null(t_names) && ncol(emb) == length(t_names)) {
			colnames(emb) <- t_names
		}
	}
	emb
}

####
#' Latent node positions from a low-rank dbn fit
#'
#' @description Returns the posterior-mean `n_row x r` node-position
#'   matrix from a `model = "lowrank"` fit, with column signs aligned
#'   across draws via the Procrustes-style trick of fixing the first
#'   draw and flipping subsequent draws' column signs to maximise their
#'   inner product with the reference. Unlike [actor_embedding()],
#'   which returns a scalar (the row-norm) per actor per time, this
#'   accessor returns the full latent position so the geometry is
#'   inspectable.
#'
#' @param fit A `dbn` fit with `model = "lowrank"`. Other model variants
#'   abort with a directive message.
#' @param align Logical; if `TRUE` (default), align column signs across
#'   draws before averaging. Set to `FALSE` for the raw posterior mean
#'   (which will likely cancel toward zero column-by-column).
#' @return An `n_row x r` matrix with `rownames` set from
#'   `dimnames(fit$Y)[[1]]` when available.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 10, r = 2, seed = 1)
#' fit <- dbn(sim$Z, model = "lowrank", family = "gaussian", r = 2,
#'            nscan = 200, burn = 100, verbose = FALSE)
#' U <- node_embedding(fit)
#' dim(U)  # n_row x r
#' }
node_embedding <- function(fit, align = TRUE) {
	if (!inherits(fit, "dbn"))
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	if (!identical(fit$model, "lowrank"))
		cli::cli_abort(c(
			"{.fun node_embedding} is defined only for {.code model = \"lowrank\"} fits.",
			"i" = "For other model variants, use {.fun actor_embedding} (row-norm summary)."
		))
	U_draws <- fit$U %||% fit$draws$U
	if (!is.list(U_draws) || length(U_draws) == 0L)
		cli::cli_abort("Could not find posterior draws of {.code U} in the fit object.")
	# all draws share the same n x r shape
	d <- dim(U_draws[[1L]])
	if (length(d) != 2L)
		cli::cli_abort(c(
			"Expected each {.code U} draw to be a 2D {.code n x r} matrix.",
			"x" = "Got dim {.val {paste(d, collapse = 'x')}}."
		))
	n_row <- d[1L]; r <- d[2L]
	if (isTRUE(align)) {
		ref <- U_draws[[1L]]
		acc <- ref
		# for each subsequent draw, flip column signs to maximise inner
		# product with the reference, then add to the running sum.
		for (m in seq_along(U_draws)[-1L]) {
			Um <- U_draws[[m]]
			for (k in seq_len(r)) {
				if (sum(Um[, k] * ref[, k]) < 0) Um[, k] <- -Um[, k]
			}
			acc <- acc + Um
		}
		out <- acc / length(U_draws)
	} else {
		out <- Reduce("+", U_draws) / length(U_draws)
	}
	# preserve actor labels.
	if (!is.null(dimnames(fit$Y)[[1L]]) && length(dimnames(fit$Y)[[1L]]) == n_row)
		rownames(out) <- dimnames(fit$Y)[[1L]]
	colnames(out) <- paste0("dim", seq_len(r))
	out
}

# -- simulate.dbn -----------------------------------------------------------

#' One-row summary of a dbn fit (broom-style `glance`)
#'
#' Returns a single-row data frame with model-level statistics: log-likelihood,
#' AIC, BIC, number of observations, sigma^2 (process variance), and metadata
#' (model type, sampler, family, dimensions). Designed for use with
#' `modelsummary` and `broom`-aware comparison tools.
#'
#' @param x A `dbn` fit object.
#' @param ... Unused.
#' @return A one-row data frame.
#' @author Tosin Salau and Shahryar Minhas
#' @export
glance.dbn <- function(x, ...) {
	# capture logLik / AIC / BIC failures with a single-fire informative
	# warning per session so silent NA columns don't hide a real dispatch bug
	warn_once <- function(quantity, msg) {
		opt <- paste0("dbn.glance_", quantity, "_warned")
		if (!isTRUE(getOption(opt, FALSE))) {
			cli::cli_warn(c(
				"{.code glance(fit)$ {quantity}} returned NA: {msg}",
				"i" = "Suppress via {.code options({opt} = TRUE)}."
			))
			options(structure(list(TRUE), names = opt))
		}
	}
	ll <- tryCatch(logLik(x), error = function(e) {
		warn_once("logLik", conditionMessage(e))
		structure(NA_real_, df = NA, nobs = NA)
	})
	dims <- x$dims %||% x$meta$dims
	fam <- if (is.list(x$family)) x$family$name else x$family
	su <- x$meta$sampler_used %||% NA_character_
	sigma2 <- if (is.numeric(x$sigma2)) mean(x$sigma2, na.rm = TRUE) else NA_real_
	n_draws <- if (is.list(x$A)) length(x$A) else (x$meta$draws %||% NA_integer_)
	aic <- tryCatch(AIC(x), error = function(e) {
		warn_once("AIC", conditionMessage(e))
		NA_real_
	})
	bic <- tryCatch(BIC(x), error = function(e) {
		warn_once("BIC", conditionMessage(e))
		NA_real_
	})
	data.frame(
		model        = x$model %||% NA_character_,
		sampler_used = su,
		family       = fam %||% NA_character_,
		n_row        = dims$n_row %||% NA_integer_,
		n_col        = dims$n_col %||% NA_integer_,
		p            = dims$p %||% NA_integer_,
		Tt           = dims$Tt %||% NA_integer_,
		n_draws      = n_draws,
		sigma2       = sigma2,
		logLik       = as.numeric(ll),
		df           = attr(ll, "df") %||% NA_integer_,
		AIC          = aic,
		BIC          = bic,
		nobs         = attr(ll, "nobs") %||% NA_integer_,
		stringsAsFactors = FALSE
	)
}

#' Augment data with fitted values and residuals (broom-style `augment`)
#'
#' Returns a long data frame with one row per (i, j, rel, t) cell of the
#' observed array, containing `.observed`, `.fitted`, `.resid`. Includes
#' cells with `NA` observed values (the diagonal and missing entries) so the
#' output aligns with the original 4D shape after `tidyr::pivot_*`.
#'
#' @param x A `dbn` fit object.
#' @param ... Unused.
#' @return A long data frame.
#' @author Tosin Salau and Shahryar Minhas
#' @export
augment.dbn <- function(x, ...) {
	if (is.null(x$Y)) cli::cli_abort("{.fun augment.dbn} needs {.code fit$Y} to be present.")
	fitted_arr <- fitted(x)
	resid_arr  <- residuals(x)
	dims <- x$dims %||% x$meta$dims
	grid <- expand.grid(
		i   = seq_len(dims$n_row),
		j   = seq_len(dims$n_col),
		rel = seq_len(dims$p),
		t   = seq_len(dims$Tt),
		KEEP.OUT.ATTRS = FALSE
	)
	grid$.observed <- as.vector(x$Y)
	grid$.fitted   <- as.vector(fitted_arr)
	grid$.resid    <- as.vector(resid_arr)
	grid
}

#' Coda compatibility shim for dbn fits
#'
#' Provides an `as.mcmc` method so users coming from a coda / BUGS
#' background don't silently get the raw fit list back when they call
#' `coda::as.mcmc(fit)`. The dbn posterior-package interop is
#' [as_draws.dbn()] / [as_draws.dbn_multichain()]; this method delegates
#' to it and converts via [posterior::as_draws_matrix()] when coda is
#' loaded, or errors directively with a pointer otherwise.
#'
#' @param x A `dbn` fit.
#' @param ... Forwarded to [as_draws.dbn()].
#' @return A `coda::mcmc` object covering the scalar variance traces.
#' @keywords internal
as.mcmc.dbn <- function(x, ...) {
	if (!requireNamespace("coda", quietly = TRUE)) {
		cli::cli_abort(c(
			"{.fn as.mcmc.dbn} requires the {.pkg coda} package.",
			"i" = "Install {.pkg coda}, or use the canonical {.fn as_draws.dbn} for the {.pkg posterior} / {.pkg bayesplot} workflow."
		))
	}
	if (!requireNamespace("posterior", quietly = TRUE)) {
		cli::cli_abort(c(
			"{.fn as.mcmc.dbn} relies on {.pkg posterior} to extract scalar traces.",
			"i" = "Install {.pkg posterior}, or call {.fn as_draws.dbn} directly once it is available."
		))
	}
	da <- as_draws.dbn(x, ...)
	mat <- posterior::as_draws_matrix(da)
	coda::as.mcmc(mat)
}

#' Simulate from a fitted dbn model
#'
#' Posterior-predictive simulation. For an MCMC fit, calls
#' [posterior_predict_dbn]. For ALS fits without bootstrap, refuses (no
#' uncertainty available).
#'
#' @param object A `dbn` fit.
#' @param nsim Number of replicated datasets (default 1).
#' @param seed Optional integer seed.
#' @param ... Passed to `posterior_predict_dbn`.
#' @return A list of replicated arrays (one per `nsim`), each shaped like
#'   `fit$Y`.
#' @importFrom stats simulate coef confint fitted logLik residuals vcov setNames
#' @author Tosin Salau and Shahryar Minhas
#' @export
simulate.dbn <- function(object, nsim = 1L, seed = NULL, ...) {
	if (length(nsim) != 1L || !is.numeric(nsim) || !is.finite(nsim) ||
	    nsim < 1 || nsim != round(nsim)) {
		cli::cli_abort(c(
			"{.arg nsim} must be a single positive integer.",
			"x" = "Got {.val {nsim}}."
		))
	}
	if (!is.null(seed)) {
		if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed)) {
			cli::cli_abort(c(
				"{.arg seed} must be a single finite number (or {.code NULL}).",
				"x" = "Got {.val {seed}}."
			))
		}
		set.seed(seed)
	}
	posterior_predict_dbn(object, draws = as.integer(nsim), seed = seed, ...)
}

# -- fitted.dbn / residuals.dbn ---------------------------------------------

#' Predicted values from a dbn fit
#'
#' Returns the posterior-mean (MCMC) or point-estimate (ALS) prediction
#' `M + A_t Phi_{t-1} B_t^T`. Use [predict.dbn] for forecasts (`H = N`)
#' or [posterior_predict_dbn] for replicated draws.
#'
#' At t=1 the model has no transition, so `fitted` returns `Y[..., 1]`
#' unchanged and `residuals` is identically zero.
#'
#' @param object A `dbn` fit.
#' @param ... Unused.
#' @return Array shaped like `fit$Y` with predicted values.
#' @author Tosin Salau and Shahryar Minhas
#' @export
fitted.dbn <- function(object, ...) {
	# HMM fits: A is per-regime, need regime sequence S_t to construct the
	# per-time A_{s_t}. Refuse with a directive warning rather than silently
	# applying A[min(t, R)] (which was the dynamic-code-path bug).
	if (identical(object$model, "hmm")) {
		cli::cli_warn(c(
			"{.fun fitted.dbn} is not yet implemented for HMM fits.",
			"i" = "HMM models have per-regime operators applied conditional on the regime sequence; constructing fitted values requires the per-draw regime path.",
			"i" = "Use {.fun posterior_predict_dbn} for replicated data or inspect {.code fit$S} for per-draw regime paths."
		))
		return(NULL)
	}
	pe <- coef(object, "all")
	A <- pe$A; B <- pe$B; M <- pe$M
	dims <- object$dims
	n_row <- dims$n_row; n_col <- dims$n_col; p <- dims$p; Tt <- dims$Tt
	Y <- object$Y
	if (is.null(Y)) cli::cli_abort("{.fun fitted.dbn} needs {.code fit$Y} to be present.")
	Y_f <- Y; Y_f[is.na(Y_f)] <- 0
	out <- array(NA, dim = c(n_row, n_col, p, Tt))
	if (.dbn_is_static(object)) {
		# static prediction: M (constant in time)
		# handle the empty-M edge case (HMM may have fit$M of length 0; for
		# static, M is always [n_row, n_col, p]).
		if (is.null(M) || length(M) == 0L) {
			cli::cli_warn("{.fun fitted.dbn}: {.code coef(fit)$M} is empty; returning Y unchanged.")
			return(Y)
		}
		for (r in 1:p) for (t in 1:Tt) out[, , r, t] <- M[, , r]
		return(out)
	}
	get_A <- function(t) if (length(dim(A)) == 3L) A[, , min(t, dim(A)[3])] else A
	get_B <- function(t) if (length(dim(B)) == 3L) B[, , min(t, dim(B)[3])] else B
	for (r in 1:p) {
		out[, , r, 1] <- Y_f[, , r, 1]
		for (t in 2:Tt) {
			Phi_prev <- Y_f[, , r, t - 1L] - M[, , r]
			out[, , r, t] <- M[, , r] + get_A(t) %*% Phi_prev %*% t(get_B(t))
		}
	}
	out
}

#' Residuals from a dbn fit
#'
#' `Y - fitted(Y)`. NAs in `Y` propagate (cells without observations).
#' At t=1 residuals are identically zero for dynamic fits (no transition).
#'
#' @param object A `dbn` fit.
#' @param ... Unused.
#' @return Array shaped like `fit$Y` with residuals.
#' @author Tosin Salau and Shahryar Minhas
#' @export
residuals.dbn <- function(object, ...) {
	object$Y - fitted(object)
}

####
#' Posterior operator draws on block or time scale
#'
#' @description Returns posterior operator draws. For piecewise fits the
#'   natural scale is per-block (one operator per regime block); time-scale
#'   expansion repeats each block operator across the time points that
#'   belong to it. For dynamic and lowrank fits draws are already per-time;
#'   block scale aborts because there are no blocks. Static fits abort
#'   either way.
#'
#' @param fit A `dbn` fit object.
#' @param scale Either `"block"` (default, piecewise only) or `"time"`.
#' @return A list with `A` and `B`, each a 4D array shaped
#'   `[n_draws, n_row, n_col, K_or_T]`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
operator_draws <- function(fit, scale = c("block", "time")) {
	scale <- match.arg(scale)
	if (!inherits(fit, "dbn"))
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	mdl <- fit$model %||% ""
	if (identical(mdl, "piecewise")) {
		A4 <- fit$draws$A_blocks
		B4 <- fit$draws$B_blocks
		if (is.null(A4))
			cli::cli_abort("Piecewise fit is missing {.code fit$draws$A_blocks}; refit.")
		if (identical(scale, "block")) return(list(A = A4, B = B4))
		block_info <- fit$blocks
		if (is.null(block_info$block_indices))
			cli::cli_abort("Cannot find block-to-time mapping in {.code fit$blocks$block_indices}.")
		Tt <- max(unlist(block_info$block_indices))
		block_id <- integer(Tt)
		for (k in seq_along(block_info$block_indices))
			block_id[block_info$block_indices[[k]]] <- k
		d <- dim(A4)
		n_draws <- d[1]; n_row <- d[3]
		n_col <- dim(B4)[3]
		A_t <- array(0, dim = c(n_draws, n_row, n_row, Tt))
		B_t <- array(0, dim = c(n_draws, n_col, n_col, Tt))
		for (t in seq_len(Tt)) {
			A_t[, , , t] <- A4[, block_id[t], , ]
			B_t[, , , t] <- B4[, block_id[t], , ]
		}
		return(list(A = A_t, B = B_t))
	}
	if (mdl %in% c("dynamic", "lowrank", "hmm")) {
		if (identical(scale, "block"))
			cli::cli_abort(c(
				"{.code scale = \"block\"} is only defined for piecewise fits.",
				"i" = "Use {.code scale = \"time\"} for {.val {mdl}} fits."
			))
		if (!is.list(fit$A) || length(fit$A) == 0L)
			cli::cli_abort("Could not find per-draw operators on the fit.")
		d1 <- dim(fit$A[[1L]])
		n_draws <- length(fit$A)
		A_t <- array(unlist(fit$A), dim = c(d1, n_draws))
		B_t <- array(unlist(fit$B), dim = c(dim(fit$B[[1L]]), n_draws))
		A_t <- aperm(A_t, c(4, 1, 2, 3))
		B_t <- aperm(B_t, c(4, 1, 2, 3))
		return(list(A = A_t, B = B_t))
	}
	cli::cli_abort("{.fun operator_draws} not defined for {.code model = \"{mdl}\"}.")
}
