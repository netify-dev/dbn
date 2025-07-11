#' Tensor Operations Utilities
#'
#' @description Shared tensor and matrix operations for all DBN models
#' @name utils_tensor
#' @keywords internal
NULL

#' Matrix Unfold Operation
#'
#' @description Unfolds an array along a specified mode
#' @param X Array to unfold
#' @param mode Mode along which to unfold (1, 2, or 3)
#' @return Matrix with mode as rows
#' @keywords internal
mat <- function(X, mode) {
    d <- dim(X)

    # Check for empty array
    if (any(d == 0)) {
        stop("Cannot unfold empty array (dimension 0)")
    }

    X <- aperm(X, c(mode, (1:length(d))[-mode]))
    dim(X) <- c(d[mode], prod(d[-mode]))
    X
}

#' Array-Matrix Product
#'
#' @description Multiply an array by a matrix along a given mode
#' @param A Array
#' @param M Matrix
#' @param k Mode
#' @return Transformed array
#' @keywords internal
amprod <- function(A, M, k) {
    K <- length(dim(A))
    AM <- M %*% mat(A, k)
    AMA <- array(AM, dim = c(dim(M)[1], dim(A)[-k]))
    aperm(AMA, match(1:K, c(k, (1:K)[-k])))
}

#' Tensor-Matrix Product
#'
#' @description Computes tensor product with matrices along specified modes
#' @param A Input tensor
#' @param B List of matrices
#' @param modes Modes for multiplication (default: along each dimension)
#' @return Transformed tensor
#' @keywords internal
tprod <- function(A, B, modes = 1:length(B)) {
    X <- A
    for (k in modes) {
        X <- amprod(X, B[[k]], k)
    }
    X
}

#' Random Standard Array Normal
#'
#' @description Generates array of standard normal random variables
#' @param d Dimensions of array
#' @return Array of standard normal values
#' @keywords internal
rsan <- function(d) {
    array(rnorm(prod(d)), d)
}

# Note: qr.Q is already provided by base R, no need to redefine
# Just use base::qr.Q(base::qr(X)) directly in code

#' Log-sum-exp trick
#'
#' @description Numerically stable log-sum-exp computation
#' @param x Numeric vector
#' @return log(sum(exp(x)))
#' @keywords internal
log_sum_exp <- function(x) {
    if (length(x) == 0) {
        return(-Inf)
    }
    m <- max(x[is.finite(x)])
    if (is.infinite(m)) {
        return(-Inf)
    }
    m + log(sum(exp(x - m)))
}

#' Kronecker Product with Sparse Option
#'
#' @description Computes Kronecker product with optional sparse matrices
#' @param A First matrix
#' @param B Second matrix
#' @param use_sparse Use sparse matrices for large dimensions
#' @return Kronecker product A ⊗ B
#' @keywords internal
kron_mult <- function(A, B, use_sparse = FALSE) {
    if (use_sparse && requireNamespace("Matrix", quietly = TRUE)) {
        Matrix::kronecker(A, B)
    } else {
        kronecker(A, B)
    }
}
