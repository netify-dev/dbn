####
#' Utility Functions for DBN Package
#'
#' @description Core utilities shared across all models
#' @name utils
#' @keywords internal
NULL
####

####
#' Null Coalescing Operator
#'
#' @name grapes-or-or-grapes
#' @param a First value
#' @param b Default value if a is NULL
#' @return a if not NULL, otherwise b
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
####

####
#' Compute Theta for Static Model
#'
#' @description On-demand theta computation for static DBN model
#' @param B Sender/receiver effect matrix (m x m)
#' @param Z Latent positions array (m x m x p x Tt) or specific slice
#' @param M Baseline mean array (m x m x p) or specific slice
#' @return Theta array with same dimensions as Z
#' @keywords internal
compute_theta_static <- function(B, Z, M) {
	if (length(dim(Z)) == 4) {
		dims <- dim(Z)
		m <- dims[1]
		p <- dims[3]
		Tt <- dims[4]

		Theta <- array(NA, dim = dims)
		for (r in 1:p) {
			for (t in 1:Tt) {
				deviation <- Z[, , r, t] - M[, , r]
				# missing entries (e.g. the NA diagonal) contribute zero
				# deviation; without this the matrix products below
				# propagate NA across the whole reconstruction
				deviation[!is.finite(deviation)] <- 0
				Theta[, , r, t] <- M[, , r] + B %*% deviation %*% t(B)
			}
		}
		return(Theta)
	} else if (length(dim(Z)) == 2) {
		deviation <- Z - M
		deviation[!is.finite(deviation)] <- 0
		return(M + B %*% deviation %*% t(B))
	} else {
		cli::cli_abort("{.arg Z} must be either a 4D array or 2D matrix.")
	}
}
####

####
#' Estimate Memory Usage for Dynamic DBN
#'
#' @description Estimates RAM (in GB) needed for a dynamic model fit: both the
#'   size of the returned object and the peak working memory during sampling.
#' @param n_row Number of sender nodes, or a data array: a 3D
#'   `[actors, actors, time]` or 4D `[actors, actors, relations, time]` array
#'   passed here has its dimensions unpacked automatically.
#' @param n_col Number of receiver nodes (default: n_row)
#' @param p Number of relations
#' @param Tt Number of time points
#' @param nscan Number of MCMC iterations after burn-in
#' @param burn Number of burn-in iterations
#' @param odens Thinning interval
#' @param time_thin Time-thinning interval
#' @param family Data family ("ordinal", "gaussian", "binary")
#' @param store_z Logical, or `NULL` (default) to match `dbn()`'s behavior:
#'   `TRUE` for ordinal/binary fits, `FALSE` for gaussian.
#' @param keep Character vector of draw components to retain, e.g.
#'   \code{c("A", "B", "M")} to drop the Theta draws. \code{NULL} keeps all.
#' @param quiet Suppress printed output
#' @return A named numeric vector with \code{object_gb} (size of the fitted
#'   object) and \code{peak_gb} (peak memory during sampling), invisibly.
#' @seealso \code{\link{dbn}}, \code{\link{dbn_dynamic}}
#' @examples
#' # Estimate memory for a moderate network
#' estimate_memory(n_row = 20, Tt = 30, nscan = 5000)
#'
#' # Bipartite network with time thinning
#' estimate_memory(n_row = 15, n_col = 25, Tt = 50,
#'                 nscan = 10000, time_thin = 5)
#'
#' # Suppress printed output, just get the value
#' gb <- estimate_memory(n_row = 10, Tt = 20, nscan = 2000, quiet = TRUE)
#' gb
#' @author Tosin Salau and Shahryar Minhas
#' @export
estimate_memory <- function(n_row, n_col = n_row, p = 1, Tt = 50,
							nscan = 5000, burn = 1000, odens = 1,
							time_thin = NULL, family = "ordinal",
							store_z = NULL, keep = NULL,
							quiet = FALSE) {
	# convenience: allow estimate_memory(Y) on a data array, not just scalar
	# dimensions -- a 3D [actors, actors, time] or 4D [actors, actors,
	# relations, time] array passed as the first argument is unpacked here
	if (is.array(n_row) && length(dim(n_row)) %in% c(3L, 4L)) {
		d <- dim(n_row)
		if (length(d) == 4L) {
			n_col <- d[2]; p <- d[3]; Tt <- d[4]; n_row <- d[1]
		} else {
			n_col <- d[2]; p <- 1L;  Tt <- d[3]; n_row <- d[1]
		}
	}
	# match dbn()'s actual default: time_thin = NULL maps to
	# max(1, Tt %/% 20), so this estimate matches the fit that will run.
	if (is.null(time_thin)) time_thin <- max(1L, Tt %/% 20L)
	# validate dimensions and MCMC settings -- without this a negative or
	# non-numeric input silently yields a nonsensical (e.g. negative) estimate
	pos <- list(n_row = n_row, n_col = n_col, p = p, Tt = Tt,
				nscan = nscan, odens = odens, time_thin = time_thin)
	for (nm in names(pos)) {
		v <- pos[[nm]]
		if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v < 1) {
			cli::cli_abort("{.arg {nm}} must be a single positive number, got {.val {v}}.")
		}
	}
	if (length(burn) != 1L || !is.numeric(burn) || !is.finite(burn) || burn < 0) {
		cli::cli_abort("{.arg burn} must be a single non-negative number, got {.val {burn}}.")
	}

	# match dbn()'s default: ordinal/binary fits store the latent Z, gaussian
	# fits do not. A fixed store_z = FALSE under-predicts the ordinal/binary
	# workflow (the common case) by ~25%.
	if (is.null(store_z)) store_z <- family %in% c("ordinal", "binary")
	nc <- n_row * n_col
	n_keep <- floor(nscan / odens)
	n_time_keep <- length(seq(1, Tt, by = time_thin))
	bpd <- 8

	# per-draw-array byte sizes (one representation)
	theta_bytes <- nc * p * n_time_keep * n_keep * bpd
	a_bytes <- n_row * n_row * n_time_keep * n_keep * bpd
	b_bytes <- n_col * n_col * n_time_keep * n_keep * bpd
	m_bytes <- nc * p * n_keep * bpd
	z_bytes <- if (isTRUE(store_z) && family %in% c("ordinal", "binary")) theta_bytes else 0

	# the fitted object holds Theta, A, B, M (and Z) in two representations
	# each -- a list form under $draws and an array accessor at the top level --
	# so the serialized object counts every draw array twice
	keep_theta <- is.null(keep) || "Theta" %in% keep
	keep_a <- is.null(keep) || "A" %in% keep
	keep_b <- is.null(keep) || "B" %in% keep
	keep_m <- is.null(keep) || "M" %in% keep

	draws_bytes <- 2 * (
		(if (keep_theta) theta_bytes else 0) +
		(if (keep_a) a_bytes else 0) +
		(if (keep_b) b_bytes else 0) +
		(if (keep_m) m_bytes else 0) +
		z_bytes)

	# input data (stored once) and scalar parameter draws. The draw arrays
	# dominate; list/metadata overhead is negligible (well under 1%), so the
	# byte sum here matches the serialized object size to within a percent.
	y_bytes <- nc * p * Tt * bpd
	scalar_bytes <- 6 * n_keep * bpd
	object_bytes <- draws_bytes + y_bytes + scalar_bytes

	# peak working memory: object plus full-resolution current-state arrays
	# (Theta, Z, A, B, M) held during the MCMC sweep, with FFBS scratch
	state_bytes <- (2 * nc * p * Tt + n_row * n_row * Tt +
	                n_col * n_col * Tt + nc * p) * bpd
	peak_bytes <- object_bytes + 4 * state_bytes

	object_gb <- object_bytes / (1024^3)
	peak_gb <- peak_bytes / (1024^3)

	if (!quiet) {
		cat(sprintf("Dynamic DBN memory estimate:\n"))
		cat(sprintf("  Network:   %d x %d, %d relation(s), %d time points\n", n_row, n_col, p, Tt))
		cat(sprintf("  MCMC:      %d draws (nscan=%d, odens=%d)\n", n_keep, nscan, odens))
		cat(sprintf("  Time thin: %d (keeping %d of %d time points)\n", time_thin, n_time_keep, Tt))
		if (!is.null(keep)) {
			cat(sprintf("  Keeping:   %s (other draws dropped)\n", paste(keep, collapse = ", ")))
		}
		cat(sprintf("  --------------------------------\n"))
		# auto-scale to KB / MB / GB
		fmt_mem <- function(gb) {
			if (is.na(gb) || !is.finite(gb)) return(sprintf("%6s", "NA"))
			if (gb < 1e-3)   sprintf("%6.2f KB", gb * 1e6)
			else if (gb < 1) sprintf("%6.2f MB", gb * 1024)
			else             sprintf("%6.2f GB", gb)
		}
		cat(sprintf("  Object:    %s\n", fmt_mem(object_gb)))
		cat(sprintf("  Peak:      %s\n", fmt_mem(peak_gb)))

		if (peak_gb > 4) {
			cat(sprintf("\n  WARNING: Estimated peak memory exceeds 4 GB.\n"))
			cat(sprintf("  Increase time_thin or odens, pass keep = c(\"A\", \"B\", \"M\"),\n"))
			cat(sprintf("  or use model = 'lowrank'.\n"))
		} else if (peak_gb > 2) {
			cat(sprintf("\n  Note: Estimated peak memory exceeds 2 GB; may be tight on 8GB-RAM laptops.\n"))
			cat(sprintf("  Consider time_thin or keep = c(\"A\", \"B\", \"M\") if needed.\n"))
		}
	}

	invisible(c(object_gb = object_gb, peak_gb = peak_gb))
}

####
#' Gibbs draw of the operator mean Abar (symmetric, lambda_diag > 0)
#'
#' @description
#' Full-conditional draw of `Abar`, the stationary mean of the AR(1) operator
#' process `A_t = (1 - phi) Abar + phi A_{t-1} + eta_t`, with per-entry
#' innovations `eta_t ~ N(0, tauA2)` and a per-entry Gaussian prior
#' `Abar ~ N(0, kappa_Abar2)`. The full conditional is conjugate Gaussian;
#' entries are drawn on the upper triangle and mirrored so `Abar` stays
#' symmetric. Only the AR(1) transitions (t = 2..Tt) contribute, matching the
#' innovation accounting used for the `tauA2` update.
#'
#' @param Aarray Array `[n, n, Tt]` of symmetric operator draws.
#' @param phi_diag AR(1) persistence scalar.
#' @param tauA2 Per-entry innovation variance scalar.
#' @param kappa_Abar2 Per-entry prior variance scalar.
#' @return Symmetric `n x n` matrix.
#' @keywords internal
#' @noRd
sample_Abar_internal <- function(Aarray, phi_diag, tauA2, kappa_Abar2) {
	d <- dim(Aarray)
	n <- d[1]
	Tt <- d[3]
	Abar <- matrix(0, n, n)
	prec_prior <- 1 / kappa_Abar2

	if (Tt < 2) {
		# no AR(1) information: draw from the prior
		post_sd <- sqrt(1 / prec_prior)
		for (i in seq_len(n)) {
			for (j in i:n) {
				Abar[i, j] <- Abar[j, i] <- stats::rnorm(1, 0, post_sd)
			}
		}
		return(Abar)
	}

	one_minus_phi <- 1 - phi_diag
	# per-entry sufficient statistic: sum_{t=2}^{Tt} (A_t - phi A_{t-1})
	d_sum <- matrix(0, n, n)
	for (t in 2:Tt) {
		d_sum <- d_sum + (Aarray[, , t] - phi_diag * Aarray[, , t - 1])
	}

	prec_post <- prec_prior + (Tt - 1) * one_minus_phi^2 / tauA2
	post_sd <- sqrt(1 / prec_post)
	mean_post <- (one_minus_phi / tauA2) * d_sum / prec_post

	for (i in seq_len(n)) {
		for (j in i:n) {
			val <- stats::rnorm(1, mean_post[i, j], post_sd)
			Abar[i, j] <- val
			Abar[j, i] <- val
		}
	}
	Abar
}
####
