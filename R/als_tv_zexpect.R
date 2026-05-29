####
# deterministic conditional-mean z imputation for binary / ordinal
# families inside the tv-als outer loop. closed-form truncated-normal
# mean (the e-step of the ecm), distinct from the mcmc sampler
# `update_Z_optimized` which returns draws.
####

#' Conditional mean of truncated normal for binary probit
#'
#' E[Z | Y=1, theta] = theta + phi(theta) / Phi(theta)         (left-truncated at 0, Z > 0)
#' E[Z | Y=0, theta] = theta - phi(theta) / Phi(-theta)        (right-truncated at 0, Z < 0)
#' Vectorized.
#'
#' @param theta predicted mean array (same dims as Y)
#' @param Y binary array, same dims as theta
#' @return Z (same dims) with conditional means; preserves NA
#' @keywords internal
#' @noRd
.expected_Z_binary <- function(theta, Y) {
	Z <- theta
	nonna <- !is.na(Y)
	idx1 <- nonna & Y == 1
	idx0 <- nonna & Y == 0
	# inverse Mills ratios
	# guard against numerical underflow/overflow: clamp theta to [-8, 8]
	t1 <- pmin(pmax(theta[idx1], -8), 8)
	t0 <- pmin(pmax(theta[idx0], -8), 8)
	Z[idx1] <- t1 + dnorm(t1) / pmax(pnorm(t1), 1e-12)
	Z[idx0] <- t0 - dnorm(t0) / pmax(pnorm(-t0), 1e-12)
	Z
}

#' Conditional mean of truncated normal for ordinal/probit
#'
#' For observed category c in {1, ..., K} with thresholds
#'   gamma_0 = -Inf < gamma_1 < gamma_2 < ... < gamma_{K-1} < gamma_K = +Inf
#' Z | Y = c, theta ~ N(theta, 1) truncated to (gamma_{c-1}, gamma_c).
#' E[Z | Y = c, theta] = theta + (phi(gamma_{c-1} - theta) - phi(gamma_c - theta)) /
#'                                (Phi(gamma_c - theta) - Phi(gamma_{c-1} - theta))
#'
#' @param theta predicted mean array
#' @param Y ordinal array with integer categories 1..K
#' @param thresholds numeric vector of length K-1 (interior cut points). If
#'   NULL, computed from the empirical quantiles of theta.
#' @keywords internal
#' @noRd
.expected_Z_ordinal <- function(theta, Y, thresholds = NULL) {
	Z <- theta
	nonna <- !is.na(Y)
	Y_int <- as.integer(round(Y))
	cats <- sort(unique(Y_int[nonna]))
	K <- length(cats)
	if (K < 2) {
		# degenerate: only one category observed; nothing to update
		return(Z)
	}
	# construct thresholds if not provided. Use the same convention as the
	# package's MCMC rank-likelihood preprocessing (cuts at quantiles of theta).
	if (is.null(thresholds)) {
		probs <- seq(1 / K, 1 - 1 / K, length.out = K - 1)
		thresholds <- quantile(theta[nonna], probs, na.rm = TRUE)
	}
	gammas <- c(-Inf, thresholds, Inf)
	# index gammas by position (1..K), not by category value -- categories
	# can start at 0, and gammas[0] would be numeric(0).
	for (k in seq_along(cats)) {
		c <- cats[k]
		idx <- which(nonna & Y_int == c)
		if (!length(idx)) next
		theta_c <- pmin(pmax(theta[idx], -8), 8)
		lo <- gammas[k] - theta_c
		hi <- gammas[k + 1L] - theta_c
		# phi(lo) - phi(hi) all finite; pnorm(hi) - pnorm(lo) > 0
		denom <- pnorm(hi) - pnorm(lo)
		denom <- pmax(denom, 1e-12)
		Z[idx] <- theta_c + (dnorm(lo) - dnorm(hi)) / denom
	}
	Z
}

#' Conditional mean Z update (deterministic ECM)
#'
#' @param Y observed array (same dims as theta)
#' @param theta current predicted mean array
#' @param family one of "binary", "ordinal"
#' @return Z array, same dims, with NAs preserved
#' @keywords internal
#' @noRd
.expected_Z <- function(Y, theta, family) {
	if (family == "binary") return(.expected_Z_binary(theta, Y))
	if (family == "ordinal") return(.expected_Z_ordinal(theta, Y))
	stop("expected_Z does not apply to family: ", family)
}
