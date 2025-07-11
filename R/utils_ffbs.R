#' FFBS Utilities
#'
#' @description Forward-Filter Backward-Sample algorithms shared across models
#' @name utils_ffbs
#' @keywords internal
NULL

# Move the existing ffbs_dlm and ffbs_theta from ffbs.R to here
# This file will contain all shared FFBS implementations


#' Forward-Filter Backward-Sample for DLM
#'
#' @description FFBS algorithm for dynamic linear model with optional AR(1) dynamics
#' @param y List of observation vectors
#' @param Flist List of design matrices
#' @param V Observation variance matrix
#' @param W State innovation variance matrix (for AR(1), this is innovation variance)
#' @param m0 Prior mean vector
#' @param C0 Prior covariance matrix
#' @param ar1 Logical: use AR(1) dynamics instead of random walk
#' @param rho AR(1) coefficient (used if ar1=TRUE)
#' @return Matrix of sampled state vectors (columns are time points)
#' @keywords internal
ffbs_dlm <- function(y, Flist, V, W, m0, C0, ar1 = FALSE, rho = 0) {
    return(ffbs_dlm_cpp(y, Flist, V, W, m0, C0, ar1, rho))
}

#' FFBS for Theta Matrix
#'
#' @description FFBS for the latent Theta process in bilinear model
#' @param Z Observed data (m x m x Tt)
#' @param mu Baseline mean (m x m)
#' @param Aarray Time-varying A matrices (m x m x Tt)
#' @param Barray Time-varying B matrices (m x m x Tt)
#' @param sigma2 Innovation variance
#' @return Array of sampled Theta matrices (m x m x Tt)
#' @keywords internal
ffbs_theta <- function(Z, mu, Aarray, Barray, sigma2) {
    # ensure we have finite values
    if (any(!is.finite(Z))) {
        Z[!is.finite(Z)] <- 0
    }
    if (any(!is.finite(mu))) {
        mu[!is.finite(mu)] <- 0
    }
    if (!is.finite(sigma2) || sigma2 <= 0) {
        sigma2 <- 1e-6
    }

    # f-ing avoid kron products
    return(ffbs_bilinear(Z, mu, Aarray, Barray, sigma2, 1.0))
}
