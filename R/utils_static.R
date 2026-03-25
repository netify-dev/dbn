####
#' Static Model Utilities
#'
#' @description Tucker product based update functions for static DBN models
#' @name utils_static
#' @keywords internal
NULL
####

####
#' Update B Matrices Using Tucker Products
#'
#' @description updates all K influence matrices using Tucker product formulation
#' @param Y response array (n_row x n_col x p x Tt)
#' @param X predictor array (n_row x n_col x p x Tt)
#' @param B list of K influence matrices
#' @param s2 observation variance
#' @param t2 prior precision for B
#' @param K number of modes (default 3)
#' @return list with updated B matrices and sse
#' @keywords internal
update_B_tucker <- function(Y, X, B, s2, t2, K = length(B)) {
	d <- dim(Y)[1:K]
	sse <- 0

	for (k in 1:K) {
		# prior mean for B[[k]] (diagonal with shrinkage)
		Ek <- diag(d[k]) * rnorm(1,
			(mean(diag(B[[k]])) * (d[k] / t2) + 1) / (d[k] / t2 + 1),
			sqrt(1 / (d[k] / t2 + 1))
		)

		# matricize Y along mode k
		Yk <- mat(Y, k)

		# apply B to all modes except k, then matricize
		other_modes <- setdiff(1:K, k)
		Xk <- mat(tprod(X, B, modes = other_modes), k)

		# sufficient statistics
		YXk <- tcrossprod(Yk, Xk)
		XXk <- tcrossprod(Xk)

		# posterior precision and mean
		Vk <- tryCatch(
			solve(XXk / s2 + diag(d[k]) / t2),
			error = function(e) {
				# regularize if singular
				solve(XXk / s2 + diag(d[k]) / t2 + diag(d[k]) * 1e-6)
			}
		)

		# sample from posterior
		B[[k]] <- (YXk / s2 + Ek / t2) %*% Vk + rsan(c(d[k], d[k])) %*% chol(Vk)

		# accumulate sse for t2 update
		sse <- sse + sum((B[[k]] - Ek)^2)
	}

	list(B = B, sse = sse)
}
####

####
#' Update Theta Using Tucker Products
#'
#' @description updates theta array using B matrices via Tucker products
#' @param Z centered data array (n_row x n_col x p x Tt)
#' @param M baseline mean (n_row x n_col x p)
#' @param B list of K influence matrices
#' @param s2 observation variance
#' @return list with updated Y (theta deviations) and X (lagged theta)
#' @keywords internal
update_theta_tucker <- function(Z, M, B, s2) {
	dims <- dim(Z)
	n_row <- dims[1]
	n_col <- dims[2]
	p <- dims[3]
	Tt <- dims[4]
	K <- 3

	# center data
	Y <- sweep(Z, c(1, 2, 3), M, "-")

	# lagged predictor (X_t = Y_{t-1})
	X <- abind::abind(Y[, , , 1, drop = FALSE] + rsan(c(n_row, n_col, p, 1)),
					  Y[, , , -Tt, drop = FALSE], along = 4)

	# precompute SVD of B matrices for efficient sampling
	V <- W <- list()
	for (k in 1:K) {
		sk <- svd(B[[k]])
		V[[k]] <- sk$v
		W[[k]] <- sk$d^2
	}

	# eigenvalue tensor (outer product of singular values squared)
	W_tensor <- outer(outer(W[[1]], W[[2]]), W[[3]])

	# transpose of B matrices
	Bt <- lapply(B, "t")

	# sample Y for each time point
	for (t in sample(1:Tt)) {
		# likelihood contribution from current time as response
		L <- (Z[, , , t] - M) + tprod(X[, , , t], B) / s2

		# likelihood contribution from next time as predictor
		if (t < Tt) {
			L <- L + tprod(Y[, , , t + 1], Bt) / s2
		}

		# posterior variance (depends on whether t is predictor for t+1)
		S2 <- 1 / (1 + 1 / s2 + W_tensor * (t < Tt) / s2)

		# sample Y[,,,t] from posterior
		Y[, , , t] <- tprod(S2 * tprod(L, lapply(V, "t")) + sqrt(S2) * rsan(c(n_row, n_col, p)), V)

		# update X for next time point
		if (t < Tt) {
			X[, , , t + 1] <- Y[, , , t]
		}
	}

	# sample X[,,,1] (initial predictor with uninformative prior)
	L <- tprod(Y[, , , 1], B) / s2
	S2 <- 1 / (1 / n_row + W_tensor / s2)
	X[, , , 1] <- tprod(S2 * tprod(L, lapply(V, "t")) + sqrt(S2) * rsan(c(n_row, n_col, p)), V)

	list(Y = Y, X = X)
}
####

####
#' Static Model MCMC Using Tucker Products
#'
#' @description alternative static model implementation using full Tucker products
#' @param Y_obs observed data array (n_row x n_col x p x Tt)
#' @param Z latent continuous array (n_row x n_col x p x Tt)
#' @param M baseline mean (n_row x n_col x p)
#' @param B list of K influence matrices
#' @param s2 observation variance
#' @param t2 prior precision for B
#' @param g2 prior variance for M
#' @return list with updated B, M, s2, t2, g2, Y_theta, X_lag
#' @keywords internal
static_gibbs_step_tucker <- function(Y_obs, Z, M, B, s2, t2, g2) {
	dims <- dim(Z)
	n_row <- dims[1]
	n_col <- dims[2]
	p <- dims[3]
	Tt <- dims[4]
	K <- 3
	d <- c(n_row, n_col, p)

	# update theta deviations
	theta_result <- update_theta_tucker(Z, M, B, s2)
	Y_theta <- theta_result$Y
	X_lag <- theta_result$X

	# update B matrices
	B_result <- update_B_tucker(Y_theta, X_lag, B, s2, t2, K)
	B <- B_result$B
	sse_B <- B_result$sse

	# update t2 (B precision)
	t2 <- 1 / rgamma(1, (1 + sum(d^2)) / 2, (1 + sse_B) / 2)

	# update M (baseline mean)
	M_num <- apply(Z - Y_theta, c(1, 2, 3), sum)
	M <- M_num / (Tt + 1 / g2) + rsan(c(n_row, n_col, p)) / sqrt(Tt + 1 / g2)

	# update g2 (M prior variance)
	g2 <- 1 / rgamma(1, (1 + n_row * n_col * p) / 2, (1 + sum(M^2)) / 2)

	# update s2 (observation variance)
	XB <- tprod(X_lag, B)
	rss <- sum((Y_theta - XB)^2)
	s2 <- 1 / rgamma(1, (1 + length(Y_theta)) / 2, (1 + rss) / 2)

	list(
		B = B,
		M = M,
		s2 = s2,
		t2 = t2,
		g2 = g2,
		Y_theta = Y_theta,
		X_lag = X_lag
	)
}
####
