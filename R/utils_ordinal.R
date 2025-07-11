#' Check if ordinal data has too many unique values for exact sampling
#' 
#' @description for ordinal data with many unique values (>10), exact truncated
#' normal sampling becomes very slow due to the need to maintain ordering 
#' constraints. in these cases, a gaussian approximation is much faster and
#' still provides reasonable results.
#' 
#' @keywords internal
#' @noRd
should_use_gaussian_approximation <- function(R, threshold = 10) {
    # Count unique non-NA values per relation
    n_unique <- numeric()
    
    if (length(dim(R)) == 4) {
        # 4D array
        p <- dim(R)[3]
        for (j in 1:p) {
            R_j <- R[, , j, , drop = FALSE]
            unique_vals <- unique(as.vector(R_j))
            unique_vals <- unique_vals[!is.na(unique_vals)]
            n_unique[j] <- length(unique_vals)
        }
    } else if (length(dim(R)) == 3) {
        # 3D array (single relation)
        unique_vals <- unique(as.vector(R))
        unique_vals <- unique_vals[!is.na(unique_vals)]
        n_unique <- length(unique_vals)
    }
    
    # Use Gaussian approximation if any relation has many unique values
    return(any(n_unique > threshold))
}

#' Fast Gaussian approximation for ordinal data with many levels
#' @keywords internal
#' @noRd
rz_gaussian_approx <- function(R, Z, EZ, sigma = 1.0) {
    # For ordinal data with many levels, use Gaussian approximation
    # This avoids the expensive truncated normal sampling
    
    Z_new <- Z
    
    # Only update non-NA entries
    valid_idx <- which(!is.na(R))
    
    if (length(valid_idx) > 0) {
        # Add small noise based on rank differences
        # Higher ranks get slightly positive perturbation
        R_vec <- as.vector(R)
        EZ_vec <- as.vector(EZ)
        
        # Normalize ranks to [0,1]
        R_valid <- R_vec[valid_idx]
        R_min <- min(R_valid, na.rm = TRUE)
        R_max <- max(R_valid, na.rm = TRUE)
        
        if (R_max > R_min) {
            R_norm <- (R_valid - R_min) / (R_max - R_min)
            
            # Sample with slight bias based on normalized rank
            # This maintains rank ordering in expectation
            for (i in seq_along(valid_idx)) {
                idx <- valid_idx[i]
                mean_val <- EZ_vec[idx] + 0.1 * (R_norm[i] - 0.5)
                Z_new[idx] <- rnorm(1, mean_val, sigma)
            }
        } else {
            # All same rank - just sample from prior
            Z_new[valid_idx] <- rnorm(length(valid_idx), EZ_vec[valid_idx], sigma)
        }
    }
    
    return(Z_new)
}