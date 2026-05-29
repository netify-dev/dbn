####
#' Simulate Data from DBN Models
#'
#' @description Functions to simulate data from various DBN model types
#' @name simulate
NULL
####

####
#' Validate common simulate_*_dbn arguments
#'
#' @description Checks the network-dimension and variance arguments shared by
#'   the `simulate_*_dbn()` functions, so a bad value produces a clear early
#'   error instead of a cryptic failure deep in array construction.
#' @param n,n_col,p,time Network dimensions; each must be a positive whole
#'   number, with at least 2 actors on each side.
#' @param variances Optional named list of variance arguments; each must be a
#'   single non-negative finite number.
#' @return Invisibly `NULL`; called for its side effect (aborts on bad input).
#' @keywords internal
#' @noRd
validate_sim_args <- function(n, n_col, p, time, variances = NULL,
							  R = NULL, r = NULL) {
	dims <- list(n = n, n_col = n_col, p = p, time = time)
	if (!is.null(R)) dims$R <- R
	if (!is.null(r)) dims$r <- r
	for (nm in names(dims)) {
		v <- dims[[nm]]
		if (length(v) != 1L || !is.numeric(v) || !is.finite(v) ||
		    v != round(v) || v < 1) {
			cli::cli_abort("{.arg {nm}} must be a single positive whole number, got {.val {v}}.")
		}
	}
	if (n < 2 || n_col < 2) {
		cli::cli_abort("Simulated networks need at least 2 actors on each side.")
	}
	if (!is.null(r) && r > min(n, n_col)) {
		cli::cli_abort(c(
			"{.arg r} ({r}) cannot exceed {.code min(n, n_col)} = {min(n, n_col)}.",
			"i" = "The low-rank factor dimension must fit within the network."
		))
	}
	for (nm in names(variances)) {
		v <- variances[[nm]]
		if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v < 0) {
			cli::cli_abort("{.arg {nm}} must be a single non-negative finite number, got {.val {v}}.")
		}
	}
	invisible(NULL)
}
####

####
#' Simulate from Static DBN Model
#'
#' @description Generates synthetic network data from a static DBN with fixed
#'   A and B influence matrices. Useful for testing model recovery and
#'   understanding the data-generating process.
#' @param n Number of actors (senders). For bipartite networks, this is the
#'   number of senders.
#' @param n_col Number of receivers (default: same as `n` for unipartite)
#' @param p Number of relation types (default: 2)
#' @param time Number of time periods to simulate
#' @param sigma2 Innovation variance of the simulated continuous series `Z`.
#'   This is a property of the generated data, not of the fitted model; judge
#'   recovery on the operators `A` and `B` rather than on `sigma2`.
#' @param tau2 Prior variance for A and B deviations from the identity matrix.
#'   Larger values produce stronger cross-actor influence.
#' @param K Number of ordinal categories for the observed data (default: 5).
#'   The continuous latent values are discretized into 1 through K.
#' @param return_truth If TRUE (default), include the true parameters
#'   in a `$truth` sub-list for validation.
#' @param seed Random seed for reproducibility
#' @param symmetric If TRUE, generate a symmetric (undirected) network: the
#'   operator (`B = A`, symmetric), baseline, and noise are all symmetrized so
#'   `Z`, `Theta`, and `Y` are symmetric.
#' @return A list containing:
#'   \item{Y}{Observed ordinal data array `[n_row, n_col, p, time]`}
#'   \item{Z}{Continuous latent values (use with `family = "gaussian"`)}
#'   \item{Theta}{True latent network state at each time point}
#'   \item{A}{True sender influence matrix}
#'   \item{B}{True receiver influence matrix}
#'   \item{M}{True baseline mean array `[n_row, n_col, p]`}
#'   \item{sigma2, tau2, K}{True parameter values used in simulation}
#' @seealso [dbn()] for model fitting, [simulate_dynamic_dbn()] for
#'   time-varying version, [simulate_test_data()] for quick test data
#' @examples
#' sim <- simulate_static_dbn(n = 8, time = 5, seed = 42)
#' dim(sim$Y)    # observed ordinal data
#' dim(sim$Z)    # continuous latent (for gaussian family)
#' dim(sim$A)    # true sender influence matrix
#' @author Tosin Salau and Shahryar Minhas
#' @export
simulate_static_dbn <- function(n = 30, n_col = n, p = 2, time = 50,
								sigma2 = 0.5, tau2 = 0.1, K = 5,
								return_truth = TRUE, seed = NULL,
								symmetric = FALSE) {
	validate_sim_args(n, n_col, p, time, list(sigma2 = sigma2, tau2 = tau2))
	if (!is.null(seed)) set.seed(seed)
	n_row <- n
	is_bipartite <- (n_row != n_col)
	if (symmetric && is_bipartite) cli::cli_abort("Symmetric networks require {.code n_row == n_col}.")
	if (symmetric && p > 1) cli::cli_abort(c(
		"Symmetric simulation requires a single relation ({.code p = 1}).",
		"i" = "The symmetric {.fun dbn} specification accepts only unipartite, single-relation data."
	))
	nc <- n_row * n_col

	# true parameters
	A <- diag(n_row) + matrix(rnorm(n_row^2, 0, sqrt(tau2)), n_row, n_row)
	B <- diag(n_col) + matrix(rnorm(n_col^2, 0, sqrt(tau2)), n_col, n_col)

	# symmetrize the operator under the symmetric specification so the
	# generated network is genuinely undirected
	if (symmetric) A <- 0.5 * (A + t(A))
	max_eig_A <- max(abs(eigen(A)$values))
	max_eig_B <- max(abs(eigen(B)$values))
	if (max_eig_A > 1e-10) A <- A / (max_eig_A * 1.01) else A <- diag(n_row)
	if (max_eig_B > 1e-10) B <- B / (max_eig_B * 1.01) else B <- diag(n_col)
	if (symmetric) B <- A

	# baseline and noise are symmetrized too, so Z, Theta, and Y come out
	# symmetric at every time point under the symmetric specification
	M <- array(rnorm(nc * p, 0, 1), dim = c(n_row, n_col, p))
	if (symmetric) for (r in seq_len(p)) M[, , r] <- 0.5 * (M[, , r] + t(M[, , r]))
	####

	# latent Z
	Z <- array(NA, dim = c(n_row, n_col, p, time))
	for (r in 1:p) {
		noise1 <- matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
		if (symmetric) noise1 <- 0.5 * (noise1 + t(noise1))
		Z[, , r, 1] <- M[, , r] + noise1
	}
	if (time > 1) {
		for (t in 2:time) {
			for (r in 1:p) {
				deviation <- Z[, , r, t - 1] - M[, , r]
				mean_t <- M[, , r] + A %*% deviation %*% t(B)
				noise_t <- matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
				if (symmetric) noise_t <- 0.5 * (noise_t + t(noise_t))
				Z[, , r, t] <- mean_t + noise_t
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

	# mask self-loops in unipartite networks: the model does not model
	# self-ties, so the diagonal is NA across the observed, latent, and
	# state arrays alike
	if (!is_bipartite) {
		for (t in 1:time) {
			for (r in 1:p) {
				diag(Y[, , r, t]) <- NA
				diag(Z[, , r, t]) <- NA
				diag(Theta[, , r, t]) <- NA
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
		out$truth <- list(A = A, B = B, M = M, Theta = Theta, cuts = cuts,
			sigma2 = sigma2, tau_A2 = tau2, tau_B2 = tau2)
	}
	return(out)
}
####

####
#' Simulate from Dynamic DBN Model
#'
#' @description Generates synthetic network data from a dynamic DBN with
#'   time-varying A and B influence matrices. Each period's influence
#'   structure evolves from the previous period's via a random walk
#'   (default) or AR(1) process.
#' @param n Number of actors (senders)
#' @param n_col Number of receivers (default: same as `n`)
#' @param p Number of relation types (default: 2)
#' @param time Number of time periods to simulate
#' @param sigma2 Innovation variance of the simulated continuous series `Z`:
#'   each period adds Gaussian noise of variance `sigma2` to the one-step mean.
#'   Note this is a property of the generated series, not of the fitted model:
#'   a dynamic [dbn()] fit separates a latent process variance from an
#'   observation variance, so its posterior `sigma2` is a different
#'   parametrization and will not numerically match this value. Judge recovery
#'   on the operators `A`, `B`, and `W_t = A_t B_t^T`, not on `sigma2`.
#' @param tauA2 Innovation variance for A (how fast sender influence changes).
#'   Larger values produce more volatile influence dynamics.
#' @param tauB2 Innovation variance for B (how fast receiver influence changes)
#' @param ar1 If TRUE, use AR(1) dynamics (smooth evolution). If FALSE (default),
#'   use random walk (less smooth).
#' @param rhoA AR(1) persistence parameter for A (only used if `ar1 = TRUE`).
#'   Values near 1 produce slowly changing influence.
#' @param rhoB AR(1) persistence parameter for B
#' @param K Number of ordinal categories for observed data (default: 5)
#' @param return_truth If TRUE (default), include true parameters in output
#' @param seed Random seed for reproducibility
#' @param symmetric If TRUE, generate a symmetric (undirected) network: the
#'   operator (`B_t = A_t`, both symmetric), the baseline, and the noise are
#'   all symmetrized, so `Z`, `Theta`, and `Y` are symmetric at every time
#'   point.
#' @param tau_A2,tau_B2 Underscore-form aliases for `tauA2` / `tauB2`. If
#'   both forms are supplied, the underscore form wins.
#' @param stabilize_operator If `TRUE` (default), rescale per-period
#'   operator draws to keep the bilinear lag operator inside the unit
#'   circle, which keeps simulated trajectories bounded.
#' @return A list containing:
#'   \item{Y}{Observed ordinal data array `[n_row, n_col, p, time]`}
#'   \item{Z}{Continuous observed series (use with `family = "gaussian"`)}
#'   \item{Theta}{Noise-free one-step conditional mean of `Z` at each time
#'     point (`M + A_t (Z_{t-1} - M) B_t^T`). This is the predicted mean, not
#'     a separately-noised latent state; a fitted [dbn()] object's `Theta` is
#'     the model's latent state and is a different quantity.}
#'   \item{A}{True time-varying sender influence `[n_row, n_row, time]`}
#'   \item{B}{True time-varying receiver influence `[n_col, n_col, time]`}
#'   \item{M}{True baseline mean array}
#'   \item{sigma2, tauA2, tauB2, rhoA, rhoB}{True parameter values}
#' @seealso [dbn()] for model fitting, [simulate_static_dbn()] for
#'   fixed-influence version
#' @examples
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 42)
#' dim(sim$Y)    # observed ordinal data
#' dim(sim$A)    # true A matrices [n, n, time]
#' @author Tosin Salau and Shahryar Minhas
#' @export
simulate_dynamic_dbn <- function(n = 30, n_col = n, p = 2, time = 50,
								 sigma2 = 0.5, tauA2 = 0.05, tauB2 = 0.05,
								 ar1 = FALSE, rhoA = 0.9, rhoB = 0.9,
								 K = 5, return_truth = TRUE, seed = NULL,
								 symmetric = FALSE,
								 tau_A2 = NULL, tau_B2 = NULL,
								 stabilize_operator = TRUE) {
	# accept tau_A2 / tau_B2 as aliases for tauA2 / tauB2 (outputs use the
	# underscore form). if both are supplied, the underscore form wins.
	if (!is.null(tau_A2)) tauA2 <- tau_A2
	if (!is.null(tau_B2)) tauB2 <- tau_B2
	# stabilize_operator = TRUE (default) rescales each A_t and B_t to
	# have spectral radius < 1 after each RW increment, so the simulated
	# theta path never explodes. consequence: the realized increments are
	# NOT N(0, tauA2 I) -- they have smaller variance than the user-supplied
	# tauA2 because each step also includes the rescale. for tau-recovery
	# tests use stabilize_operator = FALSE and pick a small tauA2 by hand.
	validate_sim_args(n, n_col, p, time,
		list(sigma2 = sigma2, tauA2 = tauA2, tauB2 = tauB2))
	if (!is.null(seed)) set.seed(seed)
	n_row <- n
	is_bipartite <- (n_row != n_col)
	if (symmetric && is_bipartite) cli::cli_abort("Symmetric networks require {.code n_row == n_col}.")
	if (symmetric && p > 1) cli::cli_abort(c(
		"Symmetric simulation requires a single relation ({.code p = 1}).",
		"i" = "The symmetric {.fun dbn} specification accepts only unipartite, single-relation data."
	))
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
			# symmetrize the operator under the symmetric specification so the
			# generated network is genuinely undirected
			if (symmetric) Aarray[, , t] <- 0.5 * (Aarray[, , t] + t(Aarray[, , t]))
			if (stabilize_operator) {
				max_eig_A <- max(abs(eigen(Aarray[, , t])$values))
				max_eig_B <- max(abs(eigen(Barray[, , t])$values))
				if (max_eig_A > 1e-10) Aarray[, , t] <- Aarray[, , t] / (max_eig_A * 1.01)
				if (max_eig_B > 1e-10) Barray[, , t] <- Barray[, , t] / (max_eig_B * 1.01)
			}
			if (symmetric) Barray[, , t] <- Aarray[, , t]
		}
	}
	####

	# baseline mean M and latent Z. Under the symmetric specification every
	# random component (baseline, operator, noise) is symmetrized, so the
	# generated Z, Theta, and Y are symmetric at every time point.
	M <- array(rnorm(nc * p, 0, 1), dim = c(n_row, n_col, p))
	if (symmetric) for (r in seq_len(p)) M[, , r] <- 0.5 * (M[, , r] + t(M[, , r]))
	Z <- array(NA, dim = c(n_row, n_col, p, time))
	for (r in 1:p) {
		noise1 <- matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
		if (symmetric) noise1 <- 0.5 * (noise1 + t(noise1))
		Z[, , r, 1] <- M[, , r] + noise1
	}
	if (time > 1) {
		for (t in 2:time) {
			for (r in 1:p) {
				deviation <- Z[, , r, t - 1] - M[, , r]
				mean_t <- M[, , r] + Aarray[, , t] %*% deviation %*% t(Barray[, , t])
				noise_t <- matrix(rnorm(nc, 0, sqrt(sigma2)), n_row, n_col)
				if (symmetric) noise_t <- 0.5 * (noise_t + t(noise_t))
				Z[, , r, t] <- mean_t + noise_t
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

	# mask self-loops in unipartite networks: the model does not model
	# self-ties, so the diagonal is NA across the observed, latent, and
	# state arrays alike
	if (!is_bipartite) {
		for (t in 1:time) {
			for (r in 1:p) {
				diag(Y[, , r, t]) <- NA
				diag(Z[, , r, t]) <- NA
				diag(Theta[, , r, t]) <- NA
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
		# standardised $truth payload across all simulate_*_dbn functions:
		# A, B, M, Theta, cuts, sigma2, tau_A2, tau_B2. enables a uniform SBC
		# story (every $truth carries M for baseline-mean coverage).
		out$truth <- list(A = Aarray, B = Barray, M = M, Theta = Theta,
			cuts = cuts, sigma2 = sigma2,
			tau_A2 = tauA2, tau_B2 = tauB2)
	}
	return(out)
}
####

####
#' Simulate a dynamic binary DBN dataset (convenience wrapper)
#'
#' Generates a binary `[n_row, n_col, 1, time]` array by simulating
#' continuous latent dynamics via [simulate_dynamic_dbn()] and
#' thresholding the latent state at zero via the probit link.
#'
#' @param n Number of actors.
#' @param time Number of time periods.
#' @param keep_relation_dim Logical. If `TRUE` (default) the returned
#'   `Y` is a 4D `[n, n, 1, time]` array. If `FALSE`, the singleton
#'   relation axis is dropped and `Y` is a 3D `[n, n, time]` array,
#'   which matches what `dim(Y)[3]` users expect to be `time` rather
#'   than `1`.
#' @param ... Passed to [simulate_dynamic_dbn()] (e.g. `sigma2`,
#'   `tauA2`, `symmetric`, `seed`).
#' @return A list with elements `Y` (binary array; shape depends on
#'   `keep_relation_dim`), `Theta` (the underlying latent state),
#'   and `truth` (the operator trajectory and cuts).
#' @seealso [simulate_dynamic_dbn()], [dbn()] with `family = "binary"`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_binary_dbn(n = 6, time = 10, seed = 1)
#' table(sim$Y, useNA = "always")
#' fit <- dbn(sim$Y, model = "dynamic", family = "binary",
#'            nscan = 200, burn = 100, verbose = FALSE)
#' }
simulate_binary_dbn <- function(n = 12, time = 10, keep_relation_dim = TRUE, ...) {
	dots <- list(...)
	# force p = 1 (binary single-relation) regardless of what the user passed.
	if (!is.null(dots$p) && dots$p != 1L) {
		cli::cli_warn(c(
			"{.fun simulate_binary_dbn} forces {.code p = 1} (single relation).",
			"i" = "Use {.fun simulate_dynamic_dbn} directly for multi-relation simulation."
		))
	}
	dots$p <- 1L
	sim <- do.call(simulate_dynamic_dbn, c(list(n = n, time = time), dots))
	# threshold the underlying continuous Theta to a {0, 1} array.
	Theta <- sim$Theta
	if (is.null(Theta)) {
		Theta <- sim$truth$Theta %||% sim$Z
	}
	Y_bin <- (Theta > 0) * 1L
	# restore the structural NA diagonal for unipartite square networks.
	if (length(dim(Y_bin)) == 4L && dim(Y_bin)[1] == dim(Y_bin)[2]) {
		for (t in seq_len(dim(Y_bin)[4])) {
			diag(Y_bin[, , 1L, t]) <- NA_integer_
		}
	}
	# optionally drop the singleton relation axis so dim(Y) reads
	# [n, n, time] rather than [n, n, 1, time].
	if (!isTRUE(keep_relation_dim) && length(dim(Y_bin)) == 4L && dim(Y_bin)[3] == 1L) {
		dn <- dimnames(Y_bin)
		Y_bin <- array(Y_bin, dim = c(dim(Y_bin)[1], dim(Y_bin)[2], dim(Y_bin)[4]))
		if (!is.null(dn)) dimnames(Y_bin) <- list(dn[[1]], dn[[2]], dn[[4]])
	}
	# attach the true M (baseline mean) to the truth list so SBC-style
	# coverage on M can use it.
	truth <- sim$truth
	if (!is.null(sim$M)) truth$M <- sim$M
	list(
		Y = Y_bin,
		Theta = Theta,
		truth = truth
	)
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
#' @param r Rank of the factorization
#' @param sigma2 Process noise variance
#' @param tau_alpha2 Variance for alpha factor innovations
#' @param tauB2 Variance for B innovations
#' @param ar1_alpha Use AR(1) for alpha dynamics (default TRUE)
#' @param rho_alpha AR(1) persistence for alpha (default 0.9)
#' @param seed Random seed for reproducibility
#' @param return_truth If TRUE (default), include true parameters in output
#' @return A list containing:
#'   \item{Y}{Observed ordinal data array `[n, n, p, time]`}
#'   \item{Z}{Continuous latent values (use with `family = "gaussian"`)}
#'   \item{Theta}{True latent network state at each time point}
#'   \item{U}{True orthogonal factor matrix `[n, r]`}
#'   \item{alpha}{True factor trajectories `[r, time]`}
#'   \item{A}{True time-varying sender influence `[n, n, time]` (reconstructed from U and alpha)}
#'   \item{B}{True time-varying receiver influence `[n_col, n_col, time]`}
#'   \item{M}{True baseline mean array `[n, n_col, p]`}
#'   \item{sigma2, tau_alpha2, tauB2, r}{True parameter values used in simulation}
#' @seealso \code{\link{dbn}} for model fitting, \code{\link{simulate_test_data}} for quick test data
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 6886)
#' dim(sim$Y)
#' dim(sim$alpha)
simulate_lowrank_dbn <- function(n = 30, n_col = n, p = 2, time = 50,
								 r = 3, sigma2 = 0.5, tau_alpha2 = 0.1,
								 tauB2 = 0.05, ar1_alpha = TRUE,
								 rho_alpha = 0.9, seed = NULL,
								 return_truth = TRUE) {
	validate_sim_args(n, n_col, p, time,
		list(sigma2 = sigma2, tau_alpha2 = tau_alpha2, tauB2 = tauB2), r = r)
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

	# mask self-loops in unipartite networks (observed, latent, and state)
	if (!is_bipartite) {
		for (rel in 1:p) {
			for (t in 1:time) {
				diag(Y[, , rel, t]) <- NA
				diag(Z[, , rel, t]) <- NA
				diag(Theta[, , rel, t]) <- NA
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
			B = Barray, M = M, Theta = Theta, cuts = cuts,
			sigma2 = sigma2, tau_alpha2 = tau_alpha2, tau_B2 = tauB2)
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
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' sim <- simulate_hmm_dbn(n = 8, p = 1, time = 10, R = 2, seed = 6886)
#' dim(sim$Y)
#' table(sim$S)
simulate_hmm_dbn <- function(n = 30, n_col = n, p = 2, time = 50,
							 R = 3, sigma2 = 0.5, tau_A2 = 0.2,
							 tau_B2 = 0.2, transition_prob = 0.8,
							 stickiness = NULL, seed = NULL,
							 return_truth = TRUE, symmetric = FALSE) {
	validate_sim_args(n, n_col, p, time,
		list(sigma2 = sigma2, tau_A2 = tau_A2, tau_B2 = tau_B2), R = R)
	if (!is.null(stickiness)) transition_prob <- stickiness
	if (!is.null(seed)) set.seed(seed)
	n_row <- n
	is_bipartite <- (n_row != n_col)
	if (symmetric && is_bipartite) cli::cli_abort("Symmetric networks require {.code n_row == n_col}.")
	if (symmetric && p > 1) cli::cli_abort(c(
		"Symmetric simulation requires a single relation ({.code p = 1}).",
		"i" = "The symmetric {.fun dbn} specification accepts only unipartite, single-relation data."
	))
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

	# markov chain for regime sequence
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

	# mask self-loops in unipartite networks (observed, latent, and state)
	if (!is_bipartite) {
		for (rel in 1:p) {
			for (t in 1:time) {
				diag(Y[, , rel, t]) <- NA
				diag(Z[, , rel, t]) <- NA
				diag(Theta[, , rel, t]) <- NA
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
			M = M, Theta = Theta, cuts = cuts, Pi = Pi,
			sigma2 = sigma2, tau_A2 = tau_A2, tau_B2 = tau_B2)
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
#' @author Tosin Salau and Shahryar Minhas
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
