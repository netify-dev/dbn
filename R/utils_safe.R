#' Safe Statistical Utilities
#'
#' @description Numerically stable implementations of common statistical operations
#' @name utils_safe
#' @keywords internal
NULL

#' Draw from an inverse-Gamma in log-scale (avoids under-/overflow)
#'
#' @description Safe inverse-gamma sampling for variance parameters
#' @param shape Shape parameter
#' @param rate Rate parameter
#' @param floor Minimum value (default 1e-8)
#' @param ceiling Maximum value (default 1e8)
#' @return Sample from inverse-gamma distribution
#' @keywords internal
safe_rinv_gamma <- function(shape, rate, floor = 1e-8, ceiling = 1e8) {
    # clamp arguments
    shape <- max(shape, 1e-10)
    rate <- max(rate, 1e-10)

    # draw on the log-scale => always finite
    log_x <- (log(rate) - log(stats::rgamma(1, shape = shape, rate = 1))) # inv-gamma
    x <- exp(log_x)

    # hard floor / ceiling
    x <- min(max(x, floor), ceiling)
    return(x)
}

#' Sparse-aware matrix multiply
#'
#' @description Efficiently multiply matrices that may be sparse
#' @param A First matrix (regular or sparse)
#' @param B Second matrix (regular or sparse)
#' @return Product A %*% B
#' @keywords internal
sparse_mult <- function(A, B) {
    # Check if either matrix is sparse
    is_sparse_A <- inherits(A, "Matrix") || inherits(A, "sparseMatrix")
    is_sparse_B <- inherits(B, "Matrix") || inherits(B, "sparseMatrix")

    if (is_sparse_A || is_sparse_B) {
        # Use Matrix package multiplication
        requireNamespace("Matrix", quietly = TRUE)
        return(A %*% B)
    } else {
        # Regular matrix multiplication
        return(A %*% B)
    }
}

#' Sparse-aware quadratic form
#'
#' @description Compute t(A) %*% B %*% A efficiently
#' @param A Matrix
#' @param B Matrix (often diagonal/sparse)
#' @return Quadratic form result
#' @keywords internal
sparse_quad <- function(A, B) {
    # Check if B is diagonal
    if (inherits(B, "ddiMatrix") || inherits(B, "diagonalMatrix") ||
        (is.matrix(B) && all(B[row(B) != col(B)] == 0))) {
        # Efficient computation for diagonal B
        diag_B <- if (inherits(B, "ddiMatrix") || inherits(B, "diagonalMatrix")) {
            Matrix::diag(B)
        } else {
            diag(B)
        }
        return(t(A * sqrt(diag_B)) %*% (A * sqrt(diag_B)))
    } else {
        # General case
        return(sparse_mult(t(A), sparse_mult(B, A)))
    }
}

#' Create sparse diagonal matrix
#'
#' @description Create diagonal matrix, using sparse format if large
#' @param x Diagonal values or dimension
#' @param sparse_threshold Size threshold for using sparse format
#' @return Diagonal matrix
#' @keywords internal
sparse_diag <- function(x, sparse_threshold = 100) {
    n <- if (length(x) == 1) x else length(x)

    if (n >= sparse_threshold && requireNamespace("Matrix", quietly = TRUE)) {
        if (length(x) == 1) {
            return(Matrix::Diagonal(n))
        } else {
            return(Matrix::Diagonal(x = x))
        }
    } else {
        return(diag(x))
    }
}

#' Simulate low-rank factors
#'
#' @description Generate orthogonal factors for low-rank decomposition
#' @param m Dimension
#' @param r Rank
#' @param sparse Use sparse matrices if beneficial
#' @return List with U and V matrices
#' @keywords internal
simulate_lowrank_factors <- function(m, r, sparse = FALSE) {
    if (r > m) stop("Rank cannot exceed dimension")

    # Generate random matrices
    U_raw <- matrix(rnorm(m * r), m, r)
    V_raw <- matrix(rnorm(m * r), m, r)

    # Orthogonalize via QR
    U <- qr.Q(qr(U_raw))
    V <- qr.Q(qr(V_raw))

    # Convert to sparse if requested and beneficial
    if (sparse && m >= 100 && r < m / 2 && requireNamespace("Matrix", quietly = TRUE)) {
        U <- Matrix::Matrix(U, sparse = TRUE)
        V <- Matrix::Matrix(V, sparse = TRUE)
    }

    list(U = U, V = V)
}

#' Initialize HMM states
#'
#' @description Initialize hidden states for HMM model
#' @param Tt Number of time points
#' @param K Number of states
#' @param init_probs Initial state probabilities
#' @param trans_probs Transition probability matrix
#' @return Vector of state assignments
#' @keywords internal
initialize_hmm_states <- function(Tt, K, init_probs = NULL, trans_probs = NULL) {
    # Default uniform initial probabilities
    if (is.null(init_probs)) {
        init_probs <- rep(1 / K, K)
    }

    # Default sticky transition matrix
    if (is.null(trans_probs)) {
        trans_probs <- matrix(0.1 / (K - 1), K, K)
        diag(trans_probs) <- 0.9
    }

    # Initialize states
    states <- integer(Tt)
    states[1] <- sample(K, 1, prob = init_probs)

    for (t in 2:Tt) {
        states[t] <- sample(K, 1, prob = trans_probs[states[t - 1], ])
    }

    states
}

#' Sparse matrix multiplication operator
#'
#' @description %*% that keeps sparse objects sparse
#' @param A First matrix
#' @param B Second matrix
#' @return Matrix product
#' @keywords internal
`%sp%` <- function(A, B) {
    if (inherits(A, "sparseMatrix") || inherits(B, "sparseMatrix")) {
        return(Matrix::crossprod(Matrix::t(A), B)) # avoids dense promote
    }
    A %*% B
}

#' Bilinear multiplication step
#'
#' @description Left and right bilinear multiply used everywhere
#' @param A Left matrix
#' @param Theta_prev Middle matrix
#' @param B Right matrix
#' @return A %*% Theta_prev %*% t(B)
#' @keywords internal
bilinear_step <- function(A, Theta_prev, B) {
    A %sp% Theta_prev %sp% Matrix::t(B)
}

#' Remove diagonal from sparse matrix
#'
#' @description Efficiently remove diagonal elements from a sparse matrix
#' @param mat Sparse or regular matrix
#' @return Matrix with diagonal elements removed (set to 0)
#' @keywords internal
drop_diag <- function(mat) {
    if (inherits(mat, "sparseMatrix")) {
        # For sparse matrices, set diagonal to 0 and drop zeros
        diag(mat) <- 0
        return(Matrix::drop0(mat))
    } else {
        # For regular matrices, just set diagonal to 0
        diag(mat) <- 0
        return(mat)
    }
}

#' Regularize matrix for safe Cholesky decomposition
#'
#' @description Add small ridge to diagonal for numerical stability
#' @param mat Matrix to regularize
#' @param eps Regularization strength (default: 5e-8 * mean diagonal)
#' @return Regularized matrix
#' @keywords internal
regularize_matrix <- function(mat, eps = NULL) {
    if (is.null(eps)) {
        eps <- 5e-8 * mean(diag(mat), na.rm = TRUE)
    }
    mat + eps * diag(nrow(mat))
}
