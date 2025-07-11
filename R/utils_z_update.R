#' Z update for ordinal and binary data
#' 
#' @description handles z updates for ordinal data using either gaussian
#' approximation or exact truncated normal sampling based on data characteristics
#' 
#' @param R rank data array
#' @param Z current latent values
#' @param Theta current theta values
#' @param M mean array
#' @param IR rank indices (only needed for exact sampling)
#' @param family data family
#' @return updated Z array
#' @keywords internal
update_Z_optimized <- function(R, Z, Theta, M, IR = NULL, family = "ordinal") {
    if (!family %in% c("ordinal", "binary")) {
        return(Z)  # only for ordinal and binary data
    }
    
    dims <- dim(R)
    
    # for binary, always use truncated normal (no approximation)
    if (family == "binary") {
        # check if truncnorm is available
        if (!requireNamespace("truncnorm", quietly = TRUE)) {
            stop("Package 'truncnorm' is required for binary outcomes. Please install it.")
        }
        
        # compute EZ = Theta + M
        if (length(dim(M)) == 3 && length(dims) == 4) {
            # broadcast M across time
            EZ <- Theta
            for (j in 1:dims[3]) {
                for (t in 1:dims[4]) {
                    EZ[,,j,t] <- Theta[,,j,t] + M[,,j]
                }
            }
        } else {
            EZ <- Theta + M
        }
        
        # update Z using truncated normal
        for (j in 1:dims[3]) {
            for (t in 1:dims[4]) {
                eta <- EZ[,,j,t]
                Y_jt <- R[,,j,t]
                
                # positions of 1s and 0s
                pos <- which(Y_jt == 1)
                neg <- which(Y_jt == 0)
                
                if (length(pos) > 0) {
                    Z[,,j,t][pos] <- truncnorm::rtruncnorm(length(pos), a = 0, b = Inf, 
                                                          mean = eta[pos], sd = 1)
                }
                if (length(neg) > 0) {
                    Z[,,j,t][neg] <- truncnorm::rtruncnorm(length(neg), a = -Inf, b = 0, 
                                                          mean = eta[neg], sd = 1)
                }
            }
        }
        return(Z)
    }
    
    # for ordinal, check if we should use gaussian approximation
    use_approx <- should_use_gaussian_approximation(R) || 
                 (prod(dims) > 5000)  # also use for large problems
    
    # compute EZ = Theta + M
    if (length(dim(M)) == 3 && length(dims) == 4) {
        # broadcast M across time
        EZ <- Theta
        for (j in 1:dims[3]) {
            for (t in 1:dims[4]) {
                EZ[,,j,t] <- Theta[,,j,t] + M[,,j]
            }
        }
    } else {
        EZ <- Theta + M
    }
    
    # flatten to 3D for processing
    m <- dims[1]
    p <- dims[3]
    Tt <- dims[4]
    
    R_flat <- array(R, c(m, m, p * Tt))
    Z_flat <- array(Z, c(m, m, p * Tt))
    EZ_flat <- array(EZ, c(m, m, p * Tt))
    
    if (use_approx) {
        # use fast gaussian approximation
        if (exists("rz_gaussian_approx_cpp", mode = "function")) {
            Z_flat <- rz_gaussian_approx_cpp(R_flat, Z_flat, EZ_flat)
        } else {
            Z_flat <- rz_gaussian_approx(R_flat, Z_flat, EZ_flat)
        }
    } else {
        # use exact truncated normal sampling
        if (!is.null(IR)) {
            Z_flat <- rz_fc_batch(R_flat, Z_flat, EZ_flat, IR, m, p, Tt)
        } else {
            # fallback to simple gaussian if no rank indices
            Z_flat <- rz_gaussian_approx(R_flat, Z_flat, EZ_flat)
        }
    }
    
    # reshape back to 4D
    array(Z_flat, dims)
}

#' Compute observation residuals for gaussian family
#' 
#' @description efficiently computes sum of squared residuals for sigma2_obs update
#' 
#' @param Z observation array (for gaussian, Z = Y)
#' @param Theta latent mean array
#' @param M baseline mean array
#' @return sum of squared residuals
#' @keywords internal
compute_gaussian_obs_residuals <- function(Z, Theta, M) {
    dims <- dim(Z)
    m <- dims[1]
    p <- dims[3]
    Tt <- dims[4]
    
    # check for c++ function first
    if (exists("compute_gaussian_obs_residuals_4d_cpp", mode = "function")) {
        # flatten arrays for c++ processing
        Z_flat <- array(Z, c(m, m, p * Tt))
        Theta_flat <- array(Theta, c(m, m, p * Tt))
        # M might be 3D (m x m x p) - keep as is
        return(compute_gaussian_obs_residuals_4d_cpp(Z_flat, Theta_flat, M, m, p, Tt))
    }
    
    # r implementation
    resid_obs <- 0
    
    # handle dimension mismatch between M and Z/Theta
    if (length(dim(M)) == 3 && length(dims) == 4) {
        # M is m x m x p, need to broadcast across time
        for (j in 1:p) {
            for (t in 1:Tt) {
                diff <- Z[, , j, t] - Theta[, , j, t] - M[, , j]
                resid_obs <- resid_obs + sum(diff^2)
            }
        }
    } else {
        # same dimensions
        diff <- Z - Theta - M
        resid_obs <- sum(diff^2)
    }
    
    resid_obs
}

#' Vectorized gaussian residual computation
#' 
#' @description vectorized version for flattened arrays
#' 
#' @param Z_flat flattened observation array
#' @param Theta_flat flattened theta array  
#' @param M_flat flattened M array
#' @param m number of nodes
#' @param p number of relations
#' @param Tt number of time points
#' @return sum of squared residuals
#' @keywords internal
compute_gaussian_obs_residuals_flat <- function(Z_flat, Theta_flat, M_flat, m, p, Tt) {
    resid_obs <- 0
    
    for (j in 1:p) {
        offset <- (j-1) * Tt
        for (t in 1:Tt) {
            diff <- Z_flat[, , offset + t] - Theta_flat[, , offset + t] - M_flat[, , j]
            resid_obs <- resid_obs + sum(diff^2)
        }
    }
    
    resid_obs
}