####
#' Z Update Utilities
#'
#' @description Latent Z updates for ordinal and binary data
#' @name utils_z_update
#' @keywords internal
NULL
####

####
#' Update Z for Ordinal or Binary Data
#'
#' @description Dispatches Z updates using either Gaussian approximation or
#'   exact truncated normal sampling based on data characteristics
#' @param R Rank data array
#' @param Z Current latent values
#' @param Theta Current theta values
#' @param M Mean array
#' @param IR Rank indices (only needed for exact sampling)
#' @param family Data family
#' @return Updated Z array
#' @keywords internal
update_Z_optimized <- function(R, Z, Theta, M, IR = NULL, family = "ordinal") {
	if (!family %in% c("ordinal", "binary")) {
		return(Z)
	}

	dims <- dim(R)

	# binary: truncated normal
	if (family == "binary") {
		if (!requireNamespace("truncnorm", quietly = TRUE)) {
			cli::cli_abort(c(
				"Package {.pkg truncnorm} is required for binary outcomes.",
				"i" = "Install with {.code install.packages(\"truncnorm\")}"
			))
		}

		# broadcast M across time
		if (length(dim(M)) == 3 && length(dims) == 4) {
			EZ <- Theta
			for (j in 1:dims[3]) {
				for (t in 1:dims[4]) {
					EZ[, , j, t] <- Theta[, , j, t] + M[, , j]
				}
			}
		} else {
			EZ <- Theta + M
		}

		for (j in 1:dims[3]) {
			for (t in 1:dims[4]) {
				eta <- EZ[, , j, t]
				Y_jt <- R[, , j, t]
				pos <- which(Y_jt == 1)
				neg <- which(Y_jt == 0)
				if (length(pos) > 0) {
					Z[, , j, t][pos] <- truncnorm::rtruncnorm(length(pos), a = 0, b = Inf,
						mean = eta[pos], sd = 1)
				}
				if (length(neg) > 0) {
					Z[, , j, t][neg] <- truncnorm::rtruncnorm(length(neg), a = -Inf, b = 0,
						mean = eta[neg], sd = 1)
				}
			}
		}
		return(Z)
	}
	####

	# ordinal: exact or approximate
	use_approx <- should_use_gaussian_approximation(R) ||
		(prod(dims) > 5000)

	# broadcast M across time
	if (length(dim(M)) == 3 && length(dims) == 4) {
		EZ <- Theta
		for (j in 1:dims[3]) {
			for (t in 1:dims[4]) {
				EZ[, , j, t] <- Theta[, , j, t] + M[, , j]
			}
		}
	} else {
		EZ <- Theta + M
	}

	n_row <- dims[1]
	n_col <- dims[2]
	p <- dims[3]
	Tt <- dims[4]

	R_flat <- array(R, c(n_row, n_col, p * Tt))
	Z_flat <- array(Z, c(n_row, n_col, p * Tt))
	EZ_flat <- array(EZ, c(n_row, n_col, p * Tt))

	if (use_approx) {
		if (exists("rz_gaussian_approx_cpp", mode = "function")) {
			Z_flat <- rz_gaussian_approx_cpp(R_flat, Z_flat, EZ_flat)
		} else {
			Z_flat <- rz_gaussian_approx(R_flat, Z_flat, EZ_flat)
		}
	} else {
		if (!is.null(IR)) {
			Z_flat <- rz_fc_batch(R_flat, Z_flat, EZ_flat, IR, n_row, n_col, p, Tt)
		} else {
			Z_flat <- rz_gaussian_approx(R_flat, Z_flat, EZ_flat)
		}
	}

	array(Z_flat, dims)
}
####

####
#' Compute Gaussian Observation Residuals
#'
#' @description Sum of squared residuals for sigma2_obs update
#' @param Z Observation array (for Gaussian, Z = Y)
#' @param Theta Latent mean array
#' @param M Baseline mean array
#' @return Sum of squared residuals
#' @keywords internal
compute_gaussian_obs_residuals <- function(Z, Theta, M) {
	dims <- dim(Z)
	n_row <- dims[1]
	n_col <- dims[2]
	p <- dims[3]
	Tt <- dims[4]

	if (exists("compute_gaussian_obs_residuals_4d_cpp", mode = "function")) {
		Z_flat <- array(Z, c(n_row, n_col, p * Tt))
		Theta_flat <- array(Theta, c(n_row, n_col, p * Tt))
		return(compute_gaussian_obs_residuals_4d_cpp(Z_flat, Theta_flat, M, n_row, p, Tt))
	}

	resid_obs <- 0
	if (length(dim(M)) == 3 && length(dims) == 4) {
		for (j in 1:p) {
			for (t in 1:Tt) {
				diff <- Z[, , j, t] - Theta[, , j, t] - M[, , j]
				resid_obs <- resid_obs + sum(diff^2)
			}
		}
	} else {
		diff <- Z - Theta - M
		resid_obs <- sum(diff^2)
	}
	resid_obs
}
####

####
#' Vectorized Gaussian Residual Computation
#'
#' @description Vectorized version for flattened arrays
#' @param Z_flat Flattened observation array
#' @param Theta_flat Flattened theta array
#' @param M_flat Flattened M array
#' @param m Number of nodes
#' @param p Number of relations
#' @param Tt Number of time points
#' @return Sum of squared residuals
#' @keywords internal
compute_gaussian_obs_residuals_flat <- function(Z_flat, Theta_flat, M_flat, m, p, Tt) {
	resid_obs <- 0
	for (j in 1:p) {
		offset <- (j - 1) * Tt
		for (t in 1:Tt) {
			diff <- Z_flat[, , offset + t] - Theta_flat[, , offset + t] - M_flat[, , j]
			resid_obs <- resid_obs + sum(diff^2)
		}
	}
	resid_obs
}
####
