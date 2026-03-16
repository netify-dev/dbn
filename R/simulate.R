####
#' Simulate Data from DBN Models
#'
#' @description Functions to simulate data from various DBN model types
#' @name simulate
NULL
####

####
#' Simulate from Static DBN Model
#'
#' @description Generates data from a static DBN with fixed A and B matrices
#' @param n Number of row actors / senders
#' @param n_col Number of column actors / receivers (default: n)
#' @param p Number of relation types
#' @param time Number of time points
#' @param sigma2 Innovation variance
#' @param tau2 Variance for A/B deviations from identity
#' @param K Number of ordinal categories
#' @param return_truth Return true parameters in a truth sub-list
#' @param seed Random seed
#' @param symmetric Logical. If TRUE, set B = A for symmetric/undirected networks. Default: FALSE.
#' @return List containing simulated data and true parameters
#' @seealso \code{\link{dbn}} for model fitting, \code{\link{simulate_test_data}} for quick test data
#' @examples
#' sim <- simulate_static_dbn(n = 8, time = 5, seed = 42)
#' str(sim$Y)
#' @export
simulate_static_dbn <- function(n = 30, n_col = n, p = 2, time = 50,
								sigma2 = 0.5, tau2 = 0.1, K = 5,
								return_truth = TRUE, seed = NULL,
								symmetric = FALSE) {
	if (!is.null(seed)) set.seed(seed)
	n_row <- n
	is_bipartite <- (n_row != n_col)
	if (symmetric && is_bipartite) cli::cli_abort("Symmetric networks require {.code n_row == n_col}.")
	nc <- n_row * n_col

	# true parameters
	A <- diag(n_row) + matrix(rnorm(n_row^2, 0, sqrt(tau2)), n_row, n_row)
	B <- diag(n_col) + matrix(rnorm(n_col^2, 0, sqrt(tau2)), n_col, n_col)

	max_eig_A <- max(abs(eigen(A)$values))
	max_eig_B <- max(abs(eigen(B)$values))
	if (max_eig_A > 1e-10) A <- A / (max_eig_A * 1.01) else A <- diag(n_row)
	if (max_eig_B > 1e-10) B <- B / (max_eig_B * 1.01) else B <- diag(n_col)
	if (symmetric) B <- A

	M <- array(rnorm(nc * p, 0, 1), dim = c(n_row, n_col, p))
	####

	# latent Z
	Z <- array(NA, dim = c(n_row, n_col, p, time))
	for (r in 1:p) {
		Z[, , r, 1] <- M[, , r] + matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
	}
	if (time > 1) {
		for (t in 2:time) {
			for (r in 1:p) {
				deviation <- Z[, , r, t - 1] - M[, , r]
				mean_t <- M[, , r] + A %*% deviation %*% t(B)
				Z[, , r, t] <- mean_t + matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
			}
		}
	}
	####

	# theta
	Theta <- array(NA, dim = c(n_row, n_col, p, time))
	Theta[, , , 1] <- M
	if (time > 1) {
		for (t in 2:time) {
			for (r in 1:p) {
				deviation <- Z[, , r, t - 1] - M[, , r]
				Theta[, , r, t] <- M[, , r] + A %*% deviation %*% t(B)
			}
		}
	}
	####

	# ordinal Y
	Y <- array(NA, dim = c(n_row, n_col, p, time))
	cuts <- vector("list", p)
	for (r in 1:p) {
		z_flat <- c(Z[, , r, ])
		probs <- seq(0, 1, length.out = K + 1)
		cuts[[r]] <- quantile(z_flat, probs = probs, na.rm = TRUE)
		cuts[[r]][1] <- -Inf
		cuts[[r]][K + 1] <- Inf
		for (t in 1:time) {
			Y[, , r, t] <- cut(Z[, , r, t], breaks = cuts[[r]], labels = 1:K)
		}
	}
	####

	# mask self-loops in unipartite networks
	if (!is_bipartite) {
		for (t in 1:time) {
			for (r in 1:p) {
				diag(Y[, , r, t]) <- NA
			}
		}
	}
	####

	out <- list(
		Y = Y, Z = Z, Theta = Theta, A = A, B = B, M = M,
		cuts = cuts, sigma2 = sigma2, tau2 = tau2, K = K,
		n_row = n_row, n_col = n_col, is_bipartite = is_bipartite,
		is_symmetric = symmetric
	)
	if (return_truth) {
		out$truth <- list(A = A, B = B, Theta = Theta, cuts = cuts)
	}
	return(out)
}
####

####
#' Simulate from Dynamic DBN Model
#'
#' @description Generates data from a dynamic DBN with time-varying A and B
#' @param n Number of row actors / senders
#' @param n_col Number of column actors / receivers (default: n)
#' @param p Number of relation types
#' @param time Number of time points
#' @param sigma2 Innovation variance
#' @param tauA2 Variance for A innovations
#' @param tauB2 Variance for B innovations
#' @param ar1 Use AR(1) dynamics instead of random walk
#' @param rhoA AR coefficient for A
#' @param rhoB AR coefficient for B
#' @param K Number of ordinal categories
#' @param return_truth Return true parameters in a truth sub-list
#' @param seed Random seed
#' @param symmetric Logical. If TRUE, set B = A at each time point. Default: FALSE.
#' @return List containing simulated data and true parameters
#' @seealso \code{\link{dbn}} for model fitting, \code{\link{simulate_test_data}} for quick test data
#' @examples
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 42)
#' str(sim$Y)
#' @export
simulate_dynamic_dbn <- function(n = 30, n_col = n, p = 2, time = 50,
								 sigma2 = 0.5, tauA2 = 0.05, tauB2 = 0.05,
								 ar1 = FALSE, rhoA = 0.9, rhoB = 0.9,
								 K = 5, return_truth = TRUE, seed = NULL,
								 symmetric = FALSE) {
	if (!is.null(seed)) set.seed(seed)
	n_row <- n
	is_bipartite <- (n_row != n_col)
	if (symmetric && is_bipartite) cli::cli_abort("Symmetric networks require {.code n_row == n_col}.")
	nc <- n_row * n_col

	# time-varying A and B matrices
	Aarray <- array(0, dim = c(n_row, n_row, time))
	Barray <- array(0, dim = c(n_col, n_col, time))
	Aarray[, , 1] <- diag(n_row)
	Barray[, , 1] <- diag(n_col)

	if (time > 1) {
		for (t in 2:time) {
			if (ar1) {
				innov_a <- matrix(rnorm(n_row^2, 0, sqrt(tauA2)), n_row, n_row)
				innov_b <- matrix(rnorm(n_col^2, 0, sqrt(tauB2)), n_col, n_col)
				Aarray[, , t] <- rhoA * Aarray[, , t - 1] + (1 - rhoA) * diag(n_row) + innov_a
				Barray[, , t] <- rhoB * Barray[, , t - 1] + (1 - rhoB) * diag(n_col) + innov_b
			} else {
				Aarray[, , t] <- Aarray[, , t - 1] + matrix(rnorm(n_row^2, 0, sqrt(tauA2)), n_row, n_row)
				Barray[, , t] <- Barray[, , t - 1] + matrix(rnorm(n_col^2, 0, sqrt(tauB2)), n_col, n_col)
			}
			max_eig_A <- max(abs(eigen(Aarray[, , t])$values))
			max_eig_B <- max(abs(eigen(Barray[, , t])$values))
			if (max_eig_A > 1e-10) Aarray[, , t] <- Aarray[, , t] / (max_eig_A * 1.01)
			if (max_eig_B > 1e-10) Barray[, , t] <- Barray[, , t] / (max_eig_B * 1.01)
			if (symmetric) Barray[, , t] <- Aarray[, , t]
		}
	}
	####

	# baseline mean M and latent Z
	M <- array(rnorm(nc * p, 0, 1), dim = c(n_row, n_col, p))
	Z <- array(NA, dim = c(n_row, n_col, p, time))
	for (r in 1:p) {
		Z[, , r, 1] <- M[, , r] + matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
	}
	if (time > 1) {
		for (t in 2:time) {
			for (r in 1:p) {
				deviation <- Z[, , r, t - 1] - M[, , r]
				mean_t <- M[, , r] + Aarray[, , t] %*% deviation %*% t(Barray[, , t])
				Z[, , r, t] <- mean_t + matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
			}
		}
	}
	####

	# theta
	Theta <- array(NA, dim = c(n_row, n_col, p, time))
	Theta[, , , 1] <- M
	if (time > 1) {
		for (t in 2:time) {
			for (r in 1:p) {
				deviation <- Z[, , r, t - 1] - M[, , r]
				Theta[, , r, t] <- M[, , r] + Aarray[, , t] %*% deviation %*% t(Barray[, , t])
			}
		}
	}
	####

	# ordinal Y
	Y <- array(NA, dim = c(n_row, n_col, p, time))
	cuts <- vector("list", p)
	for (r in 1:p) {
		z_flat <- c(Z[, , r, ])
		probs <- seq(0, 1, length.out = K + 1)
		cuts[[r]] <- quantile(z_flat, probs = probs, na.rm = TRUE)
		cuts[[r]][1] <- -Inf
		cuts[[r]][K + 1] <- Inf
		for (t in 1:time) {
			Y[, , r, t] <- cut(Z[, , r, t], breaks = cuts[[r]], labels = 1:K)
		}
	}
	####

	# mask self-loops in unipartite networks
	if (!is_bipartite) {
		for (t in 1:time) {
			for (r in 1:p) {
				diag(Y[, , r, t]) <- NA
			}
		}
	}
	####

	out <- list(
		Y = Y, Z = Z, Theta = Theta, A = Aarray, B = Barray, M = M,
		cuts = cuts, sigma2 = sigma2, tauA2 = tauA2, tauB2 = tauB2,
		ar1 = ar1, rhoA = if (ar1) rhoA else NULL, rhoB = if (ar1) rhoB else NULL,
		K = K, n_row = n_row, n_col = n_col, is_bipartite = is_bipartite,
		is_symmetric = symmetric
	)
	if (return_truth) {
		out$truth <- list(A = Aarray, B = Barray, Theta = Theta, cuts = cuts)
	}
	return(out)
}
####

####
#' Simulate from Low-Rank DBN Model
#'
#' @description Generates data from a low-rank DBN with factored A = U diag(alpha) U'
#' @param n Number of sender actors
#' @param n_col Number of receiver actors (default: n)
#' @param p Number of relation types
#' @param time Number of time points
#' @param r Rank
#' @param sigma2 Innovation variance
#' @param tau_alpha2 Variance for alpha innovations
#' @param tauB2 Variance for B innovations
#' @param ar1_alpha Use AR(1) for alpha dynamics
#' @param rho_alpha AR coefficient for alpha
#' @param seed Random seed
#' @param return_truth Return true latent factors and parameters
#' @return List containing simulated data and true parameters
#' @seealso \code{\link{dbn}} for model fitting, \code{\link{simulate_test_data}} for quick test data
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_lowrank_dbn(n = 25, p = 1, time = 12, r = 1,
#'     sigma2 = 0.2, tau_alpha2 = 0.02, seed = 42)
#' fit <- dbn_lowrank(sim$Y, r = 1, n_iter = 600, burn = 200, thin = 2)
#' }
simulate_lowrank_dbn <- function(n = 30, n_col = n, p = 2, time = 50,
								 r = 3, sigma2 = 0.5, tau_alpha2 = 0.1,
								 tauB2 = 0.05, ar1_alpha = TRUE,
								 rho_alpha = 0.9, seed = NULL,
								 return_truth = TRUE) {
	if (!is.null(seed)) set.seed(seed)
	n_row <- n
	is_bipartite <- (n_row != n_col)
	nc <- n_row * n_col

	# latent factor paths
	U_init <- matrix(rnorm(n_row * r), n_row, r)
	U <- qr.Q(qr(U_init))

	alpha <- matrix(0, r, time)
	alpha[, 1] <- rnorm(r, 0, sqrt(tau_alpha2))
	if (time > 1) {
		for (t in 2:time) {
			innov <- rnorm(r, 0, sqrt(tau_alpha2))
			alpha[, t] <- if (ar1_alpha) rho_alpha * alpha[, t - 1] + innov else alpha[, t - 1] + innov
		}
	}
	####

	# A_t and B_t
	Aarray <- array(0, c(n_row, n_row, time))
	for (t in 1:time) {
		if (r == 1) {
			Aarray[, , t] <- as.numeric(alpha[, t]) * (U %*% t(U))
		} else {
			Aarray[, , t] <- U %*% diag(alpha[, t]) %*% t(U)
		}
	}

	Barray <- array(0, c(n_col, n_col, time))
	Barray[, , 1] <- diag(n_col)
	if (time > 1) {
		for (t in 2:time) {
			Barray[, , t] <- Barray[, , t - 1] +
				matrix(rnorm(n_col^2, 0, sqrt(tauB2)), n_col, n_col)
			Barray[, , t] <- Barray[, , t] /
				(max(abs(eigen(Barray[, , t])$values)) * 1.01)
		}
	}
	####

	# baseline mean M and state-space dynamics
	M <- array(rnorm(nc * p, 0, 1), dim = c(n_row, n_col, p))
	Theta <- array(0, c(n_row, n_col, p, time))
	Z <- array(NA, dim = c(n_row, n_col, p, time))

	for (rel in 1:p) {
		Theta[, , rel, 1] <- M[, , rel]
		Z[, , rel, 1] <- Theta[, , rel, 1] + matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
	}
	if (time > 1) {
		for (t in 2:time) {
			for (rel in 1:p) {
				deviation <- Z[, , rel, t - 1] - M[, , rel]
				Theta[, , rel, t] <- M[, , rel] + Aarray[, , t] %*% deviation %*% t(Barray[, , t])
				Z[, , rel, t] <- Theta[, , rel, t] + matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
			}
		}
	}
	####

	# ordinal Y
	cuts <- vector("list", p)
	Y <- array(NA_integer_, dim = c(n_row, n_col, p, time))
	for (rel in 1:p) {
		cuts[[rel]] <- quantile(c(Z[, , rel, ]), probs = seq(0, 1, length = 6))
		cuts[[rel]][1] <- -Inf
		cuts[[rel]][length(cuts[[rel]])] <- Inf
		for (t in 1:time) {
			Y[, , rel, t] <- findInterval(Z[, , rel, t], cuts[[rel]][-1]) + 1
		}
	}
	####

	# mask self-loops in unipartite networks
	if (!is_bipartite) {
		for (rel in 1:p) {
			for (t in 1:time) {
				diag(Y[, , rel, t]) <- NA
			}
		}
	}
	####

	out <- list(
		Y = Y, Z = Z, Theta = Theta, U = U, alpha = alpha,
		A = Aarray, B = Barray, M = M, cuts = cuts, r = r,
		sigma2 = sigma2, tau_alpha2 = tau_alpha2, tauB2 = tauB2,
		ar1_alpha = ar1_alpha, rho_alpha = if (ar1_alpha) rho_alpha else NULL,
		n_row = n_row, n_col = n_col, is_bipartite = is_bipartite
	)
	if (return_truth) {
		out$truth <- list(U = U, alpha = alpha, A = Aarray,
			B = Barray, M = M, Theta = Theta, cuts = cuts)
	}
	out
}
####

####
#' Simulate from HMM DBN Model
#'
#' @description Generates data from a regime-switching HMM DBN
#' @param n Number of sender actors
#' @param n_col Number of receiver actors (default: n)
#' @param p Number of relation types
#' @param time Number of time points
#' @param R Number of regimes
#' @param sigma2 Innovation variance
#' @param tau_A2 Prior variance for regime A matrices
#' @param tau_B2 Prior variance for regime B matrices
#' @param transition_prob Diagonal transition probability
#' @param stickiness Deprecated alias for transition_prob
#' @param seed Random seed
#' @param return_truth Return true latent states and parameters
#' @param symmetric Logical. If TRUE, set B = A for each regime. Default: FALSE.
#' @return List containing simulated data and true parameters
#' @seealso \code{\link{dbn}} for model fitting, \code{\link{simulate_test_data}} for quick test data
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_hmm_dbn(n = 20, p = 1, time = 30, R = 2,
#'     transition_prob = 0.9, sigma2 = 0.1, seed = 123)
#' fit <- dbn_hmm(sim$Y, R = 2, n_iter = 800, burn = 200, thin = 2)
#' }
simulate_hmm_dbn <- function(n = 30, n_col = n, p = 2, time = 50,
							 R = 3, sigma2 = 0.5, tau_A2 = 0.2,
							 tau_B2 = 0.2, transition_prob = 0.8,
							 stickiness = NULL, seed = NULL,
							 return_truth = TRUE, symmetric = FALSE) {
	if (!is.null(stickiness)) transition_prob <- stickiness
	if (!is.null(seed)) set.seed(seed)
	n_row <- n
	is_bipartite <- (n_row != n_col)
	if (symmetric && is_bipartite) cli::cli_abort("Symmetric networks require {.code n_row == n_col}.")
	nc <- n_row * n_col

	# regime-specific A and B
	A_list <- B_list <- vector("list", R)
	for (r in 1:R) {
		A_list[[r]] <- diag(n_row) + matrix(rnorm(n_row^2, 0, sqrt(tau_A2)), n_row)
		B_list[[r]] <- diag(n_col) + matrix(rnorm(n_col^2, 0, sqrt(tau_B2)), n_col)
		if (r == 1) {
			A_list[[r]] <- 1.2 * A_list[[r]]
		} else if (r == 2) {
			B_list[[r]] <- 1.2 * B_list[[r]]
		}
		A_list[[r]] <- A_list[[r]] / (max(abs(eigen(A_list[[r]])$values)) * 1.01)
		B_list[[r]] <- B_list[[r]] / (max(abs(eigen(B_list[[r]])$values)) * 1.01)
	}
	if (symmetric) for (r in 1:R) B_list[[r]] <- A_list[[r]]
	####

	# Markov chain for regime sequence
	Pi <- matrix((1 - transition_prob) / (R - 1), R, R)
	diag(Pi) <- transition_prob
	S <- integer(time)
	S[1] <- sample.int(R, 1)
	if (time > 1) {
		for (t in 2:time) S[t] <- sample.int(R, 1, prob = Pi[S[t - 1], ])
	}
	####

	# baseline mean M and state-space dynamics
	M <- array(rnorm(nc * p, 0, 1), dim = c(n_row, n_col, p))
	Theta <- array(0, c(n_row, n_col, p, time))
	Z <- array(NA, dim = c(n_row, n_col, p, time))

	for (rel in 1:p) {
		Theta[, , rel, 1] <- M[, , rel]
		Z[, , rel, 1] <- Theta[, , rel, 1] + matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
	}
	if (time > 1) {
		for (t in 2:time) {
			for (rel in 1:p) {
				deviation <- Z[, , rel, t - 1] - M[, , rel]
				Theta[, , rel, t] <- M[, , rel] + A_list[[S[t]]] %*% deviation %*% t(B_list[[S[t]]])
				Z[, , rel, t] <- Theta[, , rel, t] + matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
			}
		}
	}
	####

	# ordinal Y
	cuts <- vector("list", p)
	Y <- array(NA_integer_, c(n_row, n_col, p, time))
	for (rel in 1:p) {
		cuts[[rel]] <- quantile(c(Z[, , rel, ]), probs = seq(0, 1, length = 6))
		cuts[[rel]][1] <- -Inf
		cuts[[rel]][length(cuts[[rel]])] <- Inf
		for (t in 1:time) {
			Y[, , rel, t] <- findInterval(Z[, , rel, t], cuts[[rel]][-1]) + 1
		}
	}
	####

	# mask self-loops in unipartite networks
	if (!is_bipartite) {
		for (rel in 1:p) {
			for (t in 1:time) {
				diag(Y[, , rel, t]) <- NA
			}
		}
	}
	####

	out <- list(
		Y = Y, Z = Z, Theta = Theta, S = S, M = M,
		A_list = A_list, B_list = B_list, Pi = Pi, cuts = cuts,
		R = R, sigma2 = sigma2, tau_A2 = tau_A2, tau_B2 = tau_B2,
		transition_prob = transition_prob,
		n_row = n_row, n_col = n_col, is_bipartite = is_bipartite,
		is_symmetric = symmetric
	)
	if (return_truth) {
		out$truth <- list(S = S, A = A_list, B = B_list,
			M = M, Theta = Theta, cuts = cuts, Pi = Pi)
	}
	out
}
####

####
#' Simulate Simple Test Data
#'
#' @description Generates simple ordinal network data for quick testing
#' @param n Number of row actors / senders
#' @param n_col Number of column actors / receivers (default: n)
#' @param p Number of relations
#' @param time Number of time points
#' @param seed Random seed
#' @return Array of ordinal network data
#' @seealso \code{\link{dbn}} for model fitting, \code{\link{simulate_static_dbn}} for full simulation with true parameters
#' @examples
#' Y <- simulate_test_data(n = 10, time = 5, seed = 42)
#' str(Y)
#' @export
simulate_test_data <- function(n = 10, n_col = n, p = 2, time = 20, seed = NULL) {
	if (!is.null(seed)) set.seed(seed)
	n_row <- n
	is_bipartite <- (n_row != n_col)
	nc <- n_row * n_col

	Y <- array(NA, dim = c(n_row, n_col, p, time))
	for (t in 1:time) {
		for (r in 1:p) {
			if (r == 1) {
				Y[, , r, t] <- matrix(sample(3:5, nc, replace = TRUE, prob = c(0.3, 0.5, 0.2)), n_row, n_col)
			} else {
				Y[, , r, t] <- matrix(sample(1:3, nc, replace = TRUE, prob = c(0.5, 0.3, 0.2)), n_row, n_col)
			}
		}
	}

	if (time > 1) {
		for (t in 2:time) {
			for (r in 1:p) {
				persist_idx <- sample(1:nc, round(0.7 * nc))
				Y[, , r, t][persist_idx] <- Y[, , r, t - 1][persist_idx]
			}
		}
	}

	if (!is_bipartite) {
		for (t in 1:time) {
			for (r in 1:p) {
				diag(Y[, , r, t]) <- NA
			}
		}
	}

	return(Y)
}
####
