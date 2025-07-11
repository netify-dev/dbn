#' Rank Likelihood Utilities
#'
#' @description Shared functions for Gaussian rank likelihood
#' @name utils_ranklik
#' @keywords internal
NULL

#' Z-scores from Ranks
#'
#' @description Computes normal z-scores from ordinal ranks
#' @param y Ordinal data matrix
#' @param ties.method Method for handling ties
#' @return Matrix of z-scores
#' @keywords internal
zscores <- function(y, ties.method = "average") {
    # General version that preserves array dimensions
    z <- qnorm(rank(c(y),
        na.last = "keep",
        ties.method = ties.method
    ) / (sum(!is.na(y)) + 1))
    array(z, dim = dim(y))
}

#' Sample Z given Rank Constraints
#'
#' @description Samples latent Z values given ordinal constraints via truncated normal
#' @param R Ordinal ranks
#' @param Z Current Z values
#' @param EZ Expected Z values
#' @param IR List of indices for each rank
#' @return Updated Z matrix
#' @keywords internal
rz_fc <- function(R, Z, EZ, IR) {
    # Convert to vectors for C++
    R_vec <- c(R)
    Z_vec <- c(Z)
    EZ_vec <- c(EZ)

    # Call C++ implementation
    Z_new <- rz_fc_cpp(R_vec, Z_vec, EZ_vec, IR)

    # Reshape back to original dimensions
    dim(Z_new) <- dim(Z)
    return(Z_new)
}

#' Precompute Rank Indices
#'
#' @description Precomputes indices for each rank level to speed up MCMC
#' @param Y Ordinal data array
#' @return List of index lists by relation
#' @keywords internal
precompute_ranks <- function(Y) {
    p <- dim(Y)[3]
    IR <- vector("list", p)

    for (j in 1:p) {
        # Extract slice as matrix (actors*actors x time)
        Y_slice <- Y[, , j, ]  # This gives us m x m x Tt array
        Y_mat <- matrix(Y_slice, nrow = dim(Y)[1] * dim(Y)[2], ncol = dim(Y)[4])
        IR[[j]] <- build_rank_indices(Y_mat)
    }
    IR
}
