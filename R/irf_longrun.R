####
# long-run / asymptotic impulse response
#
# the bilinear lag operator at a fixed (A, B) maps a shock S to its
# cumulative long-run response by stacking geometric series:
#   Delta_inf = vec^{-1} ( (I - B kron A)^{-1} vec(S) )
# the inverse exists iff the spectral radius of B kron A = rho(A) rho(B)
# is strictly below 1. when the gain >= 1 we return NA with a directive
# message.
####

#' Long-run / asymptotic impulse response
#'
#' @description Returns the cumulative response of the latent state to a
#'   one-time shock, summed over all future horizons. Defined only when
#'   the lag operator is contractive, i.e. the spectral radius of the
#'   Kronecker operator \eqn{B^\top \otimes A} is below 1 (equivalently,
#'   \eqn{\rho(A)\rho(B) < 1}). The closed form is
#'   \deqn{\mathrm{vec}(\Delta_\infty) = (I - B \otimes A)^{-1} \, \mathrm{vec}(S).}
#'
#' @param fit A fitted \code{dbn} object (static, dynamic, piecewise, or hmm).
#' @param shock A shock matrix (\eqn{n_{row} \times n_{col}}), or one
#'   built via \code{\link{build_shock}}.
#' @param t For dynamic / piecewise / hmm fits, the time index (or block /
#'   regime label) whose operator pair (A_t, B_t) to use. Default
#'   \code{NULL} uses the last observed time.
#' @param n_draws Number of posterior draws to average over (default 50);
#'   set \code{Inf} or \code{NA} to use all stored draws.
#' @param prob Credible-interval coverage (default 0.95).
#' @param seed Optional RNG seed for draw subsampling.
#' @return A list with \code{Delta_inf} (posterior mean response, n_row x
#'   n_col matrix), \code{lower} and \code{upper} bounds, \code{prob},
#'   and \code{gain_mean} (posterior mean of \eqn{\rho(A)\rho(B)} —
#'   diagnostic; \code{Delta_inf} entries return \code{NA} for draws
#'   where this exceeds 1).
#' @author Tosin Salau and Shahryar Minhas
#' @export
irf_longrun <- function(fit, shock, t = NULL, n_draws = 50L,
                          prob = 0.95, seed = NULL) {
	if (!inherits(fit, "dbn"))
		cli::cli_abort("{.arg fit} must be a {.cls dbn} fit.")
	if (!is.matrix(shock))
		cli::cli_abort("{.arg shock} must be a numeric matrix.")
	if (!is.null(seed)) set.seed(seed)
	dims <- fit$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	if (!identical(dim(shock), c(as.integer(n_row), as.integer(n_col))))
		cli::cli_abort(c(
			"{.arg shock} dimensions don't match the fit.",
			"x" = "Got {.val {paste(dim(shock), collapse = ' x ')}}; expected {.val {n_row} x {n_col}}."
		))

	# select per-draw operators for the requested t
	pick_AB <- .irf_longrun_pick_operators(fit, t)
	A_list <- pick_AB$A_list
	B_list <- pick_AB$B_list

	# subsample draws
	n_avail <- length(A_list)
	if (is.na(n_draws) || !is.finite(n_draws) || n_draws >= n_avail) {
		idx <- seq_len(n_avail)
	} else {
		idx <- sample(n_avail, as.integer(n_draws), replace = FALSE)
	}

	# per-draw long-run computation
	delta_draws <- array(NA_real_, dim = c(n_row, n_col, length(idx)))
	gains <- numeric(length(idx))
	vec_s <- as.numeric(shock)
	I_nc <- diag(n_row * n_col)
	for (k in seq_along(idx)) {
		A_k <- A_list[[idx[k]]]
		B_k <- B_list[[idx[k]]]
		ra <- tryCatch(max(Mod(eigen(A_k, only.values = TRUE)$values)),
			error = function(e) NA_real_)
		rb <- tryCatch(max(Mod(eigen(B_k, only.values = TRUE)$values)),
			error = function(e) NA_real_)
		gains[k] <- ra * rb
		if (!is.finite(gains[k]) || gains[k] >= 1) next
		# (I - B kron A)^{-1} vec(S), where vec(A X B') = (B kron A) vec(X)
		K_op <- kronecker(B_k, A_k)
		v <- tryCatch(solve(I_nc - K_op, vec_s),
			error = function(e) rep(NA_real_, n_row * n_col))
		delta_draws[, , k] <- matrix(v, n_row, n_col)
	}

	# warn if every draw was non-contractive
	n_contractive <- sum(is.finite(delta_draws[1, 1, ]))
	if (n_contractive == 0L)
		cli::cli_warn(c(
			"Every selected posterior draw has gain >= 1 at this t; the long-run IRF is undefined.",
			"i" = "Posterior mean gain: {sprintf('%.3f', mean(gains, na.rm = TRUE))}.",
			"i" = "Use {.fn compute_irf} for finite-horizon (decaying or growing) responses instead."
		))
	else if (n_contractive < length(idx))
		cli::cli_inform("Used {n_contractive} of {length(idx)} draws (others non-contractive).")

	alpha <- (1 - prob) / 2
	Delta_mean <- apply(delta_draws, c(1, 2), mean, na.rm = TRUE)
	Delta_lo   <- apply(delta_draws, c(1, 2), stats::quantile, probs = alpha, na.rm = TRUE)
	Delta_hi   <- apply(delta_draws, c(1, 2), stats::quantile, probs = 1 - alpha, na.rm = TRUE)
	out <- list(
		Delta_inf = Delta_mean,
		lower     = Delta_lo,
		upper     = Delta_hi,
		prob      = prob,
		gain_mean = mean(gains, na.rm = TRUE),
		n_contractive = n_contractive,
		n_total       = length(idx)
	)
	class(out) <- c("dbn_irf_longrun", "list")
	out
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_irf_longrun <- function(x, digits = 3L, ...) {
	cli::cli_h2("Long-run impulse response")
	cli::cli_text("Contractive draws: {x$n_contractive} / {x$n_total}; mean gain = {sprintf('%.3f', x$gain_mean)}")
	cli::cli_text("{format(100 * x$prob, nsmall = 0)}% credible-interval bounds available in $lower / $upper.")
	cli::cli_text("Posterior mean Delta_inf:")
	print(round(x$Delta_inf, digits))
	invisible(x)
}

####
# pick per-draw (A, B) operators for the requested time / block / regime
.irf_longrun_pick_operators <- function(fit, t = NULL) {
	model <- fit$model %||% "static"
	if (identical(model, "static")) {
		# fit$B is a list of three Tucker factors; per-draw sender = B[[1]],
		# receiver = B[[2]] (length(fit$B) == 3, with each element a cube)
		if (!is.list(fit$B) || length(fit$B) < 2L)
			cli::cli_abort("Static fit must store {.code B[[1]]} (sender) and {.code B[[2]]} (receiver) Tucker factors.")
		n_draws <- dim(fit$B[[1]])[3]
		A_list <- lapply(seq_len(n_draws), function(d) fit$B[[1]][, , d])
		B_list <- lapply(seq_len(n_draws), function(d) fit$B[[2]][, , d])
	} else if (identical(model, "dynamic")) {
		if (is.null(t)) t <- dim(fit$A[[1]])[3]
		t <- as.integer(t)
		A_list <- lapply(fit$A, function(a) a[, , t])
		B_list <- lapply(fit$B, function(b) b[, , t])
	} else if (identical(model, "piecewise")) {
		# pick by block index
		if (is.null(t)) {
			# default: last block
			n_blocks <- dim(fit$A[[1]])[3]
			t <- n_blocks
		}
		t <- as.integer(t)
		A_list <- lapply(fit$A, function(a) a[, , t])
		B_list <- lapply(fit$B, function(b) b[, , t])
	} else if (identical(model, "hmm")) {
		# pick by regime index (default = regime 1)
		if (is.null(t)) t <- 1L
		t <- as.integer(t)
		# HMM A: list of cubes [n_row, n_row, R]
		A_list <- lapply(fit$A, function(a) a[, , t])
		B_list <- lapply(fit$B, function(b) b[, , t])
	} else {
		cli::cli_abort(c(
			"Long-run IRF not implemented for model {.val {model}}.",
			"i" = "Supported: static, dynamic, piecewise, hmm."
		))
	}
	list(A_list = A_list, B_list = B_list)
}
