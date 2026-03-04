####
#' Safe Statistical Utilities
#'
#' @description Numerically stable implementations of common statistical operations
#' @name utils_safe
#' @keywords internal
NULL
####

####
#' Safe Inverse-Gamma Sampling
#'
#' @description Draws from inverse-Gamma on log-scale to avoid under/overflow
#' @param shape Shape parameter
#' @param rate Rate parameter
#' @param floor Minimum return value
#' @param ceiling Maximum return value
#' @return Sample from inverse-gamma distribution
#' @keywords internal
safe_rinv_gamma <- function(shape, rate, floor = 1e-8, ceiling = 1e8) {
	shape <- max(shape, 1e-10)
	rate <- max(rate, 1e-10)

	log_x <- (log(rate) - log(stats::rgamma(1, shape = shape, rate = 1)))
	x <- exp(log_x)
	x <- min(max(x, floor), ceiling)
	return(x)
}
####

####
#' Sparse-Aware Matrix Multiply
#'
#' @description Multiplies matrices that may be sparse
#' @param A First matrix
#' @param B Second matrix
#' @return Product A %*% B
#' @keywords internal
sparse_mult <- function(A, B) {
	is_sparse_A <- inherits(A, "Matrix") || inherits(A, "sparseMatrix")
	is_sparse_B <- inherits(B, "Matrix") || inherits(B, "sparseMatrix")

	if (is_sparse_A || is_sparse_B) {
		requireNamespace("Matrix", quietly = TRUE)
		return(A %*% B)
	} else {
		return(A %*% B)
	}
}
####

####
#' Sparse-Aware Quadratic Form
#'
#' @description Computes t(A) %*% B %*% A efficiently for diagonal B
#' @param A Matrix
#' @param B Matrix (often diagonal/sparse)
#' @return Quadratic form result
#' @keywords internal
sparse_quad <- function(A, B) {
	if (inherits(B, "ddiMatrix") || inherits(B, "diagonalMatrix") ||
		(is.matrix(B) && all(B[row(B) != col(B)] == 0))) {
		diag_B <- if (inherits(B, "ddiMatrix") || inherits(B, "diagonalMatrix")) {
			Matrix::diag(B)
		} else {
			diag(B)
		}
		return(t(A * sqrt(diag_B)) %*% (A * sqrt(diag_B)))
	} else {
		return(sparse_mult(t(A), sparse_mult(B, A)))
	}
}
####

####
#' Sparse Diagonal Matrix
#'
#' @description Creates diagonal matrix using sparse format for large dimensions
#' @param x Diagonal values or dimension
#' @param sparse_threshold Size threshold for sparse format
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
####

####
#' Simulate Low-Rank Factors
#'
#' @description Generates orthogonal U and V factors for low-rank decomposition
#' @param m Dimension
#' @param r Rank
#' @param sparse Use sparse matrices if beneficial
#' @return List with U and V matrices
#' @keywords internal
simulate_lowrank_factors <- function(m, r, sparse = FALSE) {
	if (r > m) cli::cli_abort("Rank cannot exceed dimension.")

	U_raw <- matrix(rnorm(m * r), m, r)
	V_raw <- matrix(rnorm(m * r), m, r)
	U <- qr.Q(qr(U_raw))
	V <- qr.Q(qr(V_raw))

	if (sparse && m >= 100 && r < m / 2 && requireNamespace("Matrix", quietly = TRUE)) {
		U <- Matrix::Matrix(U, sparse = TRUE)
		V <- Matrix::Matrix(V, sparse = TRUE)
	}

	list(U = U, V = V)
}
####

####
#' Initialize HMM States
#'
#' @description Generates initial hidden state sequence for HMM model
#' @param Tt Number of time points
#' @param K Number of states
#' @param init_probs Initial state probabilities
#' @param trans_probs Transition probability matrix
#' @return Integer vector of state assignments
#' @keywords internal
initialize_hmm_states <- function(Tt, K, init_probs = NULL, trans_probs = NULL) {
	if (is.null(init_probs)) {
		init_probs <- rep(1 / K, K)
	}
	if (is.null(trans_probs)) {
		trans_probs <- matrix(0.1 / (K - 1), K, K)
		diag(trans_probs) <- 0.9
	}

	states <- integer(Tt)
	states[1] <- sample(K, 1, prob = init_probs)
	for (t in 2:Tt) {
		states[t] <- sample(K, 1, prob = trans_probs[states[t - 1], ])
	}
	states
}
####

####
#' Sparse Matrix Multiply Operator
#'
#' @description Infix operator that keeps sparse objects sparse
#' @param A First matrix
#' @param B Second matrix
#' @return Matrix product
#' @keywords internal
`%sp%` <- function(A, B) {
	if (inherits(A, "sparseMatrix") || inherits(B, "sparseMatrix")) {
		return(Matrix::crossprod(Matrix::t(A), B))
	}
	A %*% B
}
####

####
#' Bilinear Multiplication Step
#'
#' @description Computes A %*% Theta %*% t(B)
#' @param A Left matrix
#' @param Theta_prev Middle matrix
#' @param B Right matrix
#' @return Bilinear product
#' @keywords internal
bilinear_step <- function(A, Theta_prev, B) {
	A %sp% Theta_prev %sp% Matrix::t(B)
}
####

####
#' Remove Diagonal from Matrix
#'
#' @description Sets diagonal elements to zero, handling sparse matrices
#' @param mat Sparse or regular matrix
#' @return Matrix with zero diagonal
#' @keywords internal
drop_diag <- function(mat) {
	if (inherits(mat, "sparseMatrix")) {
		diag(mat) <- 0
		return(Matrix::drop0(mat))
	} else {
		diag(mat) <- 0
		return(mat)
	}
}
####

####
#' Regularize Matrix for Cholesky
#'
#' @description Adds small ridge to diagonal for numerical stability
#' @param mat Matrix to regularize
#' @param eps Regularization strength
#' @return Regularized matrix
#' @keywords internal
regularize_matrix <- function(mat, eps = NULL) {
	if (is.null(eps)) {
		eps <- 5e-8 * mean(diag(mat), na.rm = TRUE)
	}
	mat + eps * diag(nrow(mat))
}
####
