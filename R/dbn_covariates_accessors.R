####
# s3 methods for dbn_covariates_fit
#
# extends coef / tidy / summary / print for the covariate fit class to
# expose the beta posterior alongside the standard A / B / M operators.
####

####
#' Coefficient accessor for dbn_covariates_fit
#'
#' @description Returns the posterior mean of the requested quantity. The
#'   covariate fit augments the base \code{coef.dbn} method with a
#'   \code{"beta"} option that returns the posterior mean of the
#'   covariate coefficient vector.
#'
#' @param object A \code{dbn_covariates_fit} object.
#' @param what One of \code{"all"} (default), \code{"beta"}, \code{"A"},
#'   \code{"B"}, \code{"M"}. With \code{"all"}, returns a named list
#'   including the beta vector.
#' @param ... Unused.
#' @return Posterior means; either a numeric vector (\code{"beta"}), an
#'   array (\code{"A"} / \code{"B"} / \code{"M"}), or a named list
#'   (\code{"all"}).
#' @author Tosin Salau and Shahryar Minhas
#' @export
coef.dbn_covariates_fit <- function(object,
                                     what = c("all", "beta", "A", "B", "M"),
                                     ...) {
	what <- match.arg(what)
	if (what == "beta") {
		return(.dbn_coef_beta_mean(object))
	}
	base <- NextMethod()
	if (what == "all") {
		base$beta <- .dbn_coef_beta_mean(object)
	}
	base
}

#' @keywords internal
#' @noRd
.dbn_coef_beta_mean <- function(object) {
	if (is.null(object$beta))
		cli::cli_abort("Fit does not contain a {.code beta} matrix; was it fit with {.arg covariates}?")
	if (length(dim(object$beta)) == 3L) {
		# time-varying: posterior mean is K x Tt
		return(apply(object$beta, c(2, 3), mean))
	}
	colMeans(object$beta)
}

####
#' Tidy data-frame view of a dbn_covariates_fit
#'
#' @description Appends per-coefficient rows for \code{beta} to the
#'   long-format frame returned by \code{tidy.dbn}. Each beta row has
#'   \code{parameter = "beta"}, \code{i = NA}, \code{j = NA},
#'   \code{t = NA}, and the term name in column \code{term}.
#'
#' @param x A \code{dbn_covariates_fit} object.
#' @param level Coverage level for credible intervals (default 0.95).
#' @param ... Unused.
#' @return Data frame with broom-conventional columns plus a \code{term}
#'   column carrying the covariate name (NA for the A / B / M rows).
#' @author Tosin Salau and Shahryar Minhas
#' @export
tidy.dbn_covariates_fit <- function(x, level = 0.95, ...) {
	base <- NextMethod()
	if (!"term" %in% names(base)) base$term <- NA_character_
	beta_df <- .dbn_tidy_beta(x, level = level)
	pieces <- list(base, beta_df)
	if (!is.null(x$u_row))
		pieces[[length(pieces) + 1L]] <- .dbn_tidy_actor(x$u_row, level, "u_row")
	if (!is.null(x$u_col))
		pieces[[length(pieces) + 1L]] <- .dbn_tidy_actor(x$u_col, level, "u_col")
	# align columns across pieces (NA-fill missing columns)
	all_cols <- unique(unlist(lapply(pieces, names)))
	pieces <- lapply(pieces, function(df) {
		for (col in setdiff(all_cols, names(df))) df[[col]] <- NA
		df[, all_cols]
	})
	do.call(rbind, pieces)
}

#' @keywords internal
#' @noRd
.dbn_tidy_actor <- function(u_mat, level, label) {
	n <- ncol(u_mat)
	mu <- colMeans(u_mat)
	se <- apply(u_mat, 2, stats::sd)
	alpha <- (1 - level) / 2
	lo <- apply(u_mat, 2, stats::quantile, probs = alpha)
	hi <- apply(u_mat, 2, stats::quantile, probs = 1 - alpha)
	data.frame(
		parameter = label,
		term      = paste0("actor_", seq_len(n)),
		t         = NA_integer_,
		rel       = NA_integer_,
		i         = seq_len(n),
		j         = NA_integer_,
		estimate  = as.numeric(mu),
		std.error = as.numeric(se),
		conf.low  = as.numeric(lo),
		conf.high = as.numeric(hi),
		stringsAsFactors = FALSE
	)
}

#' @keywords internal
#' @noRd
.dbn_tidy_beta <- function(x, level = 0.95) {
	beta_arr <- x$beta
	if (is.null(beta_arr))
		cli::cli_abort("Fit does not contain a {.code beta} matrix.")
	alpha <- (1 - level) / 2
	if (length(dim(beta_arr)) == 3L) {
		# time-varying: emit one row per (term, t)
		d <- dim(beta_arr)
		K <- d[2]; Tt <- d[3]
		terms <- dimnames(beta_arr)[[2]]
		if (is.null(terms)) terms <- paste0("beta", seq_len(K))
		mu <- apply(beta_arr, c(2, 3), mean)
		se <- apply(beta_arr, c(2, 3), stats::sd)
		lo <- apply(beta_arr, c(2, 3), stats::quantile, probs = alpha)
		hi <- apply(beta_arr, c(2, 3), stats::quantile, probs = 1 - alpha)
		grid <- expand.grid(k = seq_len(K), t = seq_len(Tt),
			KEEP.OUT.ATTRS = FALSE)
		return(data.frame(
			parameter = "beta",
			term      = terms[grid$k],
			t         = as.integer(grid$t),
			rel       = NA_integer_,
			i         = NA_integer_,
			j         = NA_integer_,
			estimate  = as.numeric(mu),
			std.error = as.numeric(se),
			conf.low  = as.numeric(lo),
			conf.high = as.numeric(hi),
			stringsAsFactors = FALSE
		))
	}
	# static beta: matrix view
	K <- ncol(beta_arr)
	terms <- colnames(beta_arr)
	if (is.null(terms)) terms <- paste0("beta", seq_len(K))
	mu <- colMeans(beta_arr)
	se <- apply(beta_arr, 2, stats::sd)
	lo <- apply(beta_arr, 2, stats::quantile, probs = alpha)
	hi <- apply(beta_arr, 2, stats::quantile, probs = 1 - alpha)
	data.frame(
		parameter = "beta",
		term      = terms,
		t         = NA_integer_,
		rel       = NA_integer_,
		i         = NA_integer_,
		j         = NA_integer_,
		estimate  = as.numeric(mu),
		std.error = as.numeric(se),
		conf.low  = as.numeric(lo),
		conf.high = as.numeric(hi),
		stringsAsFactors = FALSE
	)
}

####
#' Summary method for a dbn_covariates_fit
#'
#' @description Returns the same summary structure as \code{summary.dbn}
#'   plus a \code{beta} data frame with posterior mean / SD / 95% CI for
#'   each covariate term.
#' @param object A \code{dbn_covariates_fit} object.
#' @param level Coverage for credible intervals (default 0.95).
#' @param ... Passed to \code{summary.dbn}.
#' @author Tosin Salau and Shahryar Minhas
#' @export
summary.dbn_covariates_fit <- function(object, level = 0.95, ...) {
	out <- NextMethod()
	if (is.null(out)) out <- list()
	out$beta <- .dbn_tidy_beta(object, level = level)
	class(out) <- c("summary.dbn_covariates_fit", class(out))
	out
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
print.summary.dbn_covariates_fit <- function(x, ...) {
	# strip both summary.dbn_covariates_fit AND dbn_covariates_fit from the
	# class chain so dispatch lands on the base summary printer (not
	# print.dbn_covariates_fit, which assumes $beta is the draws matrix
	# rather than the tidy data frame summary() stores)
	cls <- class(x)
	x_base <- x
	class(x_base) <- setdiff(cls,
		c("summary.dbn_covariates_fit", "dbn_covariates_fit"))
	# avoid double-printing the beta table from the base printer
	x_base$beta <- NULL
	print(x_base, ...)
	beta_df <- x$beta
	if (!is.null(beta_df) && is.data.frame(beta_df) && nrow(beta_df) > 0L) {
		cli::cli_h3("Covariate coefficients (posterior summary)")
		print(beta_df[, c("term", "estimate", "std.error", "conf.low", "conf.high")],
			row.names = FALSE)
	}
	invisible(x)
}

####
#' Print method for dbn_covariates_fit
#'
#' @param x A \code{dbn_covariates_fit} object.
#' @param ... Unused.
#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_covariates_fit <- function(x, ...) {
	NextMethod()
	if (!is.null(x$beta)) {
		if (length(dim(x$beta)) == 3L) {
			# time-varying: print average-over-t mean per term
			mu_t <- apply(x$beta, c(2, 3), mean)  # K x Tt
			mu_avg <- rowMeans(mu_t)
			terms <- dimnames(x$beta)[[2]]
			if (is.null(terms)) terms <- paste0("beta", seq_along(mu_avg))
			cli::cli_text("Covariate posterior means (averaged over t): {.field {paste(sprintf('%s = %.3f', terms, mu_avg), collapse = ', ')}}")
		} else {
			mu <- colMeans(x$beta)
			se <- apply(x$beta, 2, stats::sd)
			terms <- colnames(x$beta)
			if (is.null(terms)) terms <- paste0("beta", seq_along(mu))
			cli::cli_text("Covariate posterior means: {.field {paste(sprintf('%s = %.3f (%.3f)', terms, mu, se), collapse = ', ')}}")
		}
	}
	invisible(x)
}

####
#' Variance-covariance matrix of the covariate coefficient vector
#'
#' @description Returns the empirical posterior covariance of the beta
#'   draws stored in \code{x$beta}. The result has the standard
#'   \code{vcov} layout: a K x K symmetric matrix.
#' @param object A \code{dbn_covariates_fit} object.
#' @param ... Unused.
#' @return A K x K covariance matrix.
#' @author Tosin Salau and Shahryar Minhas
#' @export
vcov.dbn_covariates_fit <- function(object, ...) {
	if (is.null(object$beta))
		cli::cli_abort("Fit does not contain a {.code beta} matrix; was it fit with {.arg covariates}?")
	if (length(dim(object$beta)) == 3L) {
		# time-varying: average per-t covariance across time
		d <- dim(object$beta)
		K <- d[2]; Tt <- d[3]
		V_sum <- matrix(0, K, K)
		for (t in seq_len(Tt)) {
			slice <- object$beta[, , t]
			# preserve matrix shape when K == 1
			if (is.null(dim(slice))) slice <- matrix(slice, ncol = 1L)
			V_sum <- V_sum + stats::cov(slice)
		}
		V <- V_sum / Tt
		terms <- dimnames(object$beta)[[2]]
		if (is.null(terms)) terms <- paste0("beta", seq_len(K))
		dimnames(V) <- list(terms, terms)
		return(V)
	}
	V <- stats::cov(object$beta)
	K <- ncol(object$beta)
	terms <- colnames(object$beta)
	if (is.null(terms)) terms <- paste0("beta", seq_len(K))
	dimnames(V) <- list(terms, terms)
	V
}

####
#' Credible intervals for a dbn_covariates_fit
#'
#' @description Extends \code{confint.dbn} to add a \code{beta} component.
#'   When \code{what = "beta"} (or \code{parm = "beta"}), returns a K x 2
#'   quantile matrix in the CRAN \code{confint} style. Otherwise returns
#'   the parent list structure with an extra \code{$beta} matrix appended.
#' @param object A \code{dbn_covariates_fit} object.
#' @param parm Optional; pass \code{"beta"} to short-circuit to the K x 2
#'   matrix view.
#' @param level Coverage level (default 0.95).
#' @param what One of \code{"all"} (default), \code{"beta"}, \code{"A"},
#'   \code{"B"}, \code{"M"}.
#' @param ... Unused.
#' @return Either a list (parent contract, with \code{$beta} appended) or
#'   a K x 2 matrix when only beta is requested.
#' @author Tosin Salau and Shahryar Minhas
#' @export
confint.dbn_covariates_fit <- function(object, parm = NULL, level = 0.95,
                                        what = c("all", "beta", "A", "B", "M"),
                                        ...) {
	what <- match.arg(what)
	if (identical(parm, "beta")) what <- "beta"
	if (what == "beta") {
		return(.dbn_beta_qmat(object, level = level))
	}
	# delegate to confint.dbn for A / B / M; append beta to the returned list
	# when what = "all"
	cls <- class(object)
	on.exit(class(object) <- cls, add = TRUE)
	class(object) <- setdiff(class(object), "dbn_covariates_fit")
	out <- stats::confint(object, parm = parm, level = level, what = what)
	if (what == "all" && is.list(out)) {
		out$beta <- .dbn_beta_qmat(object, level = level)
	}
	out
}

#' @keywords internal
#' @noRd
.dbn_beta_qmat <- function(object, level = 0.95) {
	if (is.null(object$beta))
		cli::cli_abort("Fit does not contain a {.code beta} matrix.")
	alpha <- (1 - level) / 2
	if (length(dim(object$beta)) == 3L) {
		# time-varying: return K x Tt x 2 array of {lo, hi}
		lo <- apply(object$beta, c(2, 3), stats::quantile, probs = alpha)
		hi <- apply(object$beta, c(2, 3), stats::quantile, probs = 1 - alpha)
		out <- array(NA_real_, dim = c(dim(lo), 2L),
			dimnames = list(dimnames(object$beta)[[2]], NULL,
				c(sprintf("%g%%", 100 * alpha), sprintf("%g%%", 100 * (1 - alpha)))))
		out[, , 1] <- lo
		out[, , 2] <- hi
		return(out)
	}
	terms <- colnames(object$beta)
	if (is.null(terms)) terms <- paste0("beta", seq_len(ncol(object$beta)))
	lo <- apply(object$beta, 2, stats::quantile, probs = alpha)
	hi <- apply(object$beta, 2, stats::quantile, probs = 1 - alpha)
	CI <- cbind(lo, hi)
	colnames(CI) <- c(sprintf("%g%%", 100 * alpha),
		sprintf("%g%%", 100 * (1 - alpha)))
	rownames(CI) <- terms
	CI
}
