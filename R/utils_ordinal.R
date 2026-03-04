####
#' Ordinal Utilities
#'
#' @description Helper functions for ordinal data handling
#' @keywords internal
#' @noRd
NULL
####

####
#' Check if Ordinal Data Needs Gaussian Approximation
#'
#' @description Ordinal data with many unique values (>10) is slow with exact
#'   truncated normal sampling. Returns TRUE when Gaussian approximation is faster.
#' @param R Rank array (3D or 4D)
#' @param threshold Maximum unique values before switching to approximation
#' @return Logical
#' @keywords internal
#' @noRd
should_use_gaussian_approximation <- function(R, threshold = 10) {
	n_unique <- numeric()

	if (length(dim(R)) == 4) {
		p <- dim(R)[3]
		for (j in 1:p) {
			R_j <- R[, , j, , drop = FALSE]
			unique_vals <- unique(as.vector(R_j))
			unique_vals <- unique_vals[!is.na(unique_vals)]
			n_unique[j] <- length(unique_vals)
		}
	} else if (length(dim(R)) == 3) {
		unique_vals <- unique(as.vector(R))
		unique_vals <- unique_vals[!is.na(unique_vals)]
		n_unique <- length(unique_vals)
	}

	return(any(n_unique > threshold))
}
####

####
#' Gaussian Approximation for Many-Level Ordinal Data
#'
#' @description Fast approximation for ordinal data with many levels.
#'   Avoids expensive truncated normal sampling.
#' @param R Rank array
#' @param Z Current latent Z values
#' @param EZ Expected Z values
#' @param sigma Noise scale
#' @return Updated Z array
#' @keywords internal
#' @noRd
rz_gaussian_approx <- function(R, Z, EZ, sigma = 1.0) {
	Z_new <- Z
	valid_idx <- which(!is.na(R))

	if (length(valid_idx) > 0) {
		R_vec <- as.vector(R)
		EZ_vec <- as.vector(EZ)

		R_valid <- R_vec[valid_idx]
		R_min <- min(R_valid, na.rm = TRUE)
		R_max <- max(R_valid, na.rm = TRUE)

		if (R_max > R_min) {
			R_norm <- (R_valid - R_min) / (R_max - R_min)
			for (i in seq_along(valid_idx)) {
				idx <- valid_idx[i]
				mean_val <- EZ_vec[idx] + 0.1 * (R_norm[i] - 0.5)
				Z_new[idx] <- rnorm(1, mean_val, sigma)
			}
		} else {
			Z_new[valid_idx] <- rnorm(length(valid_idx), EZ_vec[valid_idx], sigma)
		}
	}

	return(Z_new)
}
####
