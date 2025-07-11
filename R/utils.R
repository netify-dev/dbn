#' Utility Functions for DBN Package
#'
#' @description Legacy utility functions - most have been moved to utils_*.R files
#' @name utils
#' @keywords internal
NULL


#' Null coalescing operator
#'
#' @name grapes-or-or-grapes
#' @param a First value
#' @param b Default value if a is NULL
#' @return a if not NULL, otherwise b
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Compute theta for static model
#'
#' @description Compute theta values on-demand for static DBN model
#' @param B Sender/receiver effect matrix (m x m)
#' @param Z Latent positions array (m x m x p x Tt) or specific slice
#' @param M Baseline mean array (m x m x p) or specific slice
#' @return Theta array with same dimensions as Z
#' @keywords internal
compute_theta_static <- function(B, Z, M) {
    # Handle both full arrays and specific slices
    if (length(dim(Z)) == 4) {
        # Full array case
        dims <- dim(Z)
        m <- dims[1]
        p <- dims[3]
        Tt <- dims[4]
        
        Theta <- array(NA, dim = dims)
        
        for (r in 1:p) {
            for (t in 1:Tt) {
                # For static model: Theta = M + B * (Z - M) * B'
                deviation <- Z[, , r, t] - M[, , r]
                Theta[, , r, t] <- M[, , r] + B %*% deviation %*% t(B)
            }
        }
        
        return(Theta)
    } else if (length(dim(Z)) == 2) {
        # Single slice case (m x m matrix)
        # Assumes M is also m x m matrix
        deviation <- Z - M
        return(M + B %*% deviation %*% t(B))
    } else {
        stop("Z must be either a 4D array or 2D matrix")
    }
}
