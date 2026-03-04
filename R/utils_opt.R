####
#' Performance Utilities
#'
#' @description Optimized update routines using batch C++ implementations
#' @name utils_opt
#' @keywords internal
NULL
####

####
#' Update A and B with Pre-allocated Memory
#'
#' @description A/B updates using batch C++ implementation
#' @param Theta_all Full Theta array (n_row x n_col x p x Tt)
#' @param Aarray Current A array (n_row x n_row x Tt)
#' @param Barray Current B array (n_col x n_col x Tt)
#' @param sigma2 Observation variance
#' @param tauA2 A innovation variance
#' @param tauB2 B innovation variance
#' @param ar1 Use AR(1) dynamics
#' @param rhoA AR(1) coefficient for A
#' @param rhoB AR(1) coefficient for B
#' @return List with updated A and B arrays
#' @keywords internal
update_AB_optimized <- function(Theta_all, Aarray, Barray,
								sigma2, tauA2, tauB2,
								ar1 = FALSE, rhoA = 0, rhoB = 0) {
	dims <- dim(Theta_all)
	n_row <- dims[1]
	n_col <- dims[2]
	p <- dims[3]
	Tt <- dims[4]

	# reshape 4D array into list of 3D cubes
	Theta_cubes <- list(
		Theta_all[, , 1, , drop = FALSE],
		if (p >= 2) Theta_all[, , 2, , drop = FALSE] else array(0, c(n_row, n_col, Tt)),
		if (p >= 3) Theta_all[, , 3, , drop = FALSE] else array(0, c(n_row, n_col, Tt)),
		if (p >= 4) Theta_all[, , 4, , drop = FALSE] else array(0, c(n_row, n_col, Tt))
	)
	for (i in 1:4) {
		dim(Theta_cubes[[i]]) <- c(n_row, n_col, Tt)
	}

	A_result <- update_A_batch(
		Theta_cubes[[1]], Theta_cubes[[2]],
		Theta_cubes[[3]], Theta_cubes[[4]],
		Aarray, Barray, sigma2, tauA2, ar1, rhoA, p
	)

	B_result <- update_B_batch(
		Theta_cubes[[1]], Theta_cubes[[2]],
		Theta_cubes[[3]], Theta_cubes[[4]],
		A_result$A, Barray, sigma2, tauB2, ar1, rhoB, p
	)

	list(A = A_result$A, B = B_result$B)
}
####

####
#' Compute Residual Sum of Squares
#'
#' @description Efficient RSS computation for variance updates
#' @param Z Latent values array (n_row x n_col x p x Tt)
#' @param M Mean array (n_row x n_col x p)
#' @return Residual sum of squares
#' @keywords internal
compute_rss_fast <- function(Z, M) {
	dims <- dim(Z)
	n_row <- dims[1]
	n_col <- dims[2]
	p <- dims[3]
	Tt <- dims[4]

	Z_field <- list()
	for (j in 1:p) {
		Z_field[[j]] <- Z[, , j, , drop = FALSE]
		dim(Z_field[[j]]) <- c(n_row, n_col, Tt)
	}

	compute_residual_sum_squares(Z_field, M, p)
}
####

####
#' Update M with Optimal Memory Usage
#'
#' @description Computes posterior mean of M efficiently
#' @param Z Latent values array (n_row x n_col x p x Tt)
#' @param g2 Prior variance for M
#' @return Updated M array
#' @keywords internal
update_M_fast <- function(Z, g2) {
	dims <- dim(Z)
	n_row <- dims[1]
	n_col <- dims[2]
	p <- dims[3]
	Tt <- dims[4]

	Z_field <- list()
	for (j in 1:p) {
		Z_field[[j]] <- Z[, , j, , drop = FALSE]
		dim(Z_field[[j]]) <- c(n_row, n_col, Tt)
	}

	compute_M_update(Z_field, g2, n_row, p)
}
####
