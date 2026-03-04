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
				Theta[, , r, t] <- M[, , r] + B %*% deviation %*% t(B)
			}
		}
		return(Theta)
	} else if (length(dim(Z)) == 2) {
		deviation <- Z - M
		return(M + B %*% deviation %*% t(B))
	} else {
		cli::cli_abort("{.arg Z} must be either a 4D array or 2D matrix.")
	}
}
####

####
#' Estimate Memory Usage for Dynamic DBN
#'
#' @description Estimates RAM (in GB) needed for dynamic model output storage
#' @param n_row Number of sender nodes
#' @param n_col Number of receiver nodes (default: n_row)
#' @param p Number of relations
#' @param Tt Number of time points
#' @param nscan Number of MCMC iterations after burn-in
#' @param burn Number of burn-in iterations
#' @param odens Thinning interval
#' @param time_thin Time-thinning interval
#' @param family Data family ("ordinal", "gaussian", "binary")
#' @param quiet Suppress printed output
#' @return Estimated memory in GB (invisibly)
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
#' @export
estimate_memory <- function(n_row, n_col = n_row, p = 1, Tt = 50,
							nscan = 5000, burn = 1000, odens = 1,
							time_thin = 1, family = "ordinal",
							quiet = FALSE) {
	nc <- n_row * n_col
	n_keep <- floor(nscan / odens)
	n_time_keep <- length(seq(1, Tt, by = time_thin))
	bytes_per_double <- 8

	theta_bytes <- nc * p * n_time_keep * n_keep * bytes_per_double
	z_bytes <- if (family %in% c("ordinal", "binary")) theta_bytes else 0
	a_bytes <- n_row * n_row * n_time_keep * n_keep * bytes_per_double
	b_bytes <- n_col * n_col * n_time_keep * n_keep * bytes_per_double
	m_bytes <- nc * p * n_keep * bytes_per_double
	scalar_bytes <- 4 * n_keep * bytes_per_double

	total_bytes <- theta_bytes + z_bytes + a_bytes + b_bytes + m_bytes + scalar_bytes
	total_gb <- total_bytes / (1024^3)

	if (!quiet) {
		cat(sprintf("Dynamic DBN memory estimate:\n"))
		cat(sprintf("  Network:   %d x %d, %d relation(s), %d time points\n", n_row, n_col, p, Tt))
		cat(sprintf("  MCMC:      %d draws (nscan=%d, odens=%d)\n", n_keep, nscan, odens))
		cat(sprintf("  Time thin: %d (keeping %d of %d time points)\n", time_thin, n_time_keep, Tt))
		cat(sprintf("  --------------------------------\n"))
		cat(sprintf("  Theta:  %6.2f GB\n", theta_bytes / 1024^3))
		if (z_bytes > 0) cat(sprintf("  Z:      %6.2f GB\n", z_bytes / 1024^3))
		cat(sprintf("  A:      %6.2f GB\n", a_bytes / 1024^3))
		cat(sprintf("  B:      %6.2f GB\n", b_bytes / 1024^3))
		cat(sprintf("  M:      %6.2f GB\n", m_bytes / 1024^3))
		cat(sprintf("  --------------------------------\n"))
		cat(sprintf("  TOTAL:  %6.2f GB\n", total_gb))

		if (total_gb > 4) {
			cat(sprintf("\n  WARNING: Estimated memory exceeds 4 GB.\n"))
			cat(sprintf("  Consider increasing time_thin or odens, or using model='lowrank'.\n"))
		}
	}

	invisible(total_gb)
}
####
