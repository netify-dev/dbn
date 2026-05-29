####
# Simulated-truth regimes for validating RW-smoothed time-varying ALS.
#
# Each regime is a generator that returns:
#   list(name, A_true, B_true, M_true, sigma2, Y, Z, dims, family, regime_type)
#
# Where:
#   A_true, B_true are arrays [n_row, n_row, T] and [n_col, n_col, T] with
#                  A_true[,,t], B_true[,,t] the truth at time t.
#   M_true         is the baseline mean [n_row, n_col, p].
#   Y              is the observed array [n_row, n_col, p, T].
#   Z              is the latent (= Y for Gaussian).
#   dims, family   match the dbn() convention.
#   regime_type    one of "constant", "smooth_rw", "piecewise", "near_unstable".
#
# These regimes are used by:
#   - tests/validation/test_recovery.R     (truth recovery)
#   - tests/validation/test_cv.R           (lambda selection behavior)
#   - tests/validation/test_bootstrap.R    (CI coverage)
#   - tests/validation/test_mcmc_match.R   (ALS vs MCMC)
####

# Helper: generate a smooth random orthogonal-ish A matrix with spectral
# radius near `rho`. Uses a small random matrix + scaling.
.regime_rand_operator <- function(n, rho = 0.5, seed = NULL) {
	if (!is.null(seed)) set.seed(seed)
	A <- matrix(rnorm(n * n, 0, 1 / sqrt(n)), n, n)
	# scale to target spectral radius
	curr_rho <- max(Mod(eigen(A, only.values = TRUE)$values))
	if (curr_rho > 0) A <- A * (rho / curr_rho)
	A
}

# Simulate Z trajectory under the model
#   Phi_t = A_t Phi_{t-1} B_t^T + eps_t,    eps_t ~ N(0, sigma2)
#   Z_t = M + Phi_t
.regime_simulate_Z <- function(A_arr, B_arr, M, sigma2, n_row, n_col, p, Tt) {
	Z <- array(0, dim = c(n_row, n_col, p, Tt))
	Phi <- array(0, dim = c(n_row, n_col, Tt))
	# Initial state
	Phi[, , 1] <- matrix(rnorm(n_row * n_col, 0, sqrt(sigma2)), n_row, n_col)
	for (t in 2:Tt) {
		mean_t <- A_arr[, , t] %*% Phi[, , t - 1] %*% t(B_arr[, , t])
		Phi[, , t] <- mean_t + matrix(rnorm(n_row * n_col, 0, sqrt(sigma2)), n_row, n_col)
	}
	for (r in 1:p) for (t in 1:Tt) Z[, , r, t] <- M[, , r] + Phi[, , t]
	# diagonal NA for square unipartite networks
	if (n_row == n_col) for (r in 1:p) for (t in 1:Tt) diag(Z[, , r, t]) <- NA
	list(Z = Z, Phi = Phi)
}

# Family-specific Y observation given latent Z
.regime_observe <- function(Z, family) {
	if (family == "gaussian") {
		return(Z)
	} else if (family == "binary") {
		Y <- Z
		Y[!is.na(Z)] <- as.numeric(Z[!is.na(Z)] > 0)
		return(Y)
	} else if (family == "ordinal") {
		Y <- Z
		nonna <- !is.na(Z)
		# 5 categories with cuts at the quintiles of the latent
		cuts <- quantile(Z[nonna], probs = seq(0.2, 0.8, 0.2))
		Y[nonna] <- as.numeric(cut(Z[nonna], breaks = c(-Inf, cuts, Inf), labels = 1:5))
		return(Y)
	}
	stop("unknown family: ", family)
}

#' Regime 1: constant operator
#'
#' A_t = A_0, B_t = B_0 for all t. The static ALS should recover this
#' exactly at any reasonable T; CV should select lambda -> Inf.
regime_constant <- function(n = 8, Tt = 20, p = 1, family = "gaussian",
                            rho_A = 0.5, rho_B = 0.5, sigma2 = 0.5, seed = 1) {
	set.seed(seed)
	n_row <- n_col <- n
	A0 <- .regime_rand_operator(n_row, rho_A, seed = seed + 1)
	B0 <- .regime_rand_operator(n_col, rho_B, seed = seed + 2)
	A_arr <- array(rep(c(A0), Tt), dim = c(n_row, n_row, Tt))
	B_arr <- array(rep(c(B0), Tt), dim = c(n_col, n_col, Tt))
	M <- array(rnorm(n_row * n_col * p, 0, 0.3), dim = c(n_row, n_col, p))
	sim <- .regime_simulate_Z(A_arr, B_arr, M, sigma2, n_row, n_col, p, Tt)
	Y <- .regime_observe(sim$Z, family)
	list(name = "constant", regime_type = "constant",
		A_true = A_arr, B_true = B_arr, M_true = M, sigma2 = sigma2,
		Y = Y, Z = sim$Z, Phi = sim$Phi,
		dims = list(n_row = n_row, n_col = n_col, p = p, Tt = Tt),
		family = family, seed = seed)
}

#' Regime 2: smooth random walk on the operator
#'
#' A_t = A_{t-1} + Gaussian innovation of variance tau2. CV should select
#' a finite, intermediate lambda.
regime_smooth_rw <- function(n = 8, Tt = 20, p = 1, family = "gaussian",
                              rho_A_init = 0.5, tau2 = 0.005, sigma2 = 0.5,
                              seed = 2) {
	set.seed(seed)
	n_row <- n_col <- n
	A0 <- .regime_rand_operator(n_row, rho_A_init, seed = seed + 1)
	B0 <- .regime_rand_operator(n_col, rho_A_init, seed = seed + 2)
	A_arr <- array(0, dim = c(n_row, n_row, Tt))
	B_arr <- array(0, dim = c(n_col, n_col, Tt))
	A_arr[, , 1] <- A0
	B_arr[, , 1] <- B0
	for (t in 2:Tt) {
		A_arr[, , t] <- A_arr[, , t - 1] + matrix(rnorm(n_row^2, 0, sqrt(tau2)), n_row, n_row)
		B_arr[, , t] <- B_arr[, , t - 1] + matrix(rnorm(n_col^2, 0, sqrt(tau2)), n_col, n_col)
	}
	M <- array(rnorm(n_row * n_col * p, 0, 0.3), dim = c(n_row, n_col, p))
	sim <- .regime_simulate_Z(A_arr, B_arr, M, sigma2, n_row, n_col, p, Tt)
	Y <- .regime_observe(sim$Z, family)
	list(name = sprintf("smooth_rw_tau2_%g", tau2), regime_type = "smooth_rw",
		A_true = A_arr, B_true = B_arr, M_true = M, sigma2 = sigma2,
		Y = Y, Z = sim$Z, Phi = sim$Phi,
		dims = list(n_row = n_row, n_col = n_col, p = p, Tt = Tt),
		family = family, seed = seed, tau2 = tau2)
}

#' Regime 3: piecewise-constant operator with one or two break-points
#'
#' The RW penalty will smear these breaks; the residual-burst diagnostic
#' should fire. CV should select small lambda but the trajectory still
#' misfits at the break.
regime_piecewise <- function(n = 8, Tt = 20, p = 1, family = "gaussian",
                              break_points = c(0.4, 0.7),
                              rho_levels = c(0.3, 0.7, 0.5),
                              sigma2 = 0.5, seed = 3) {
	set.seed(seed)
	n_row <- n_col <- n
	# break_points are fractions of T
	break_idx <- round(break_points * Tt)
	break_idx <- pmax(2, pmin(Tt - 1, break_idx))
	regime_start <- c(1, break_idx + 1)
	regime_end <- c(break_idx, Tt)
	A_arr <- array(0, dim = c(n_row, n_row, Tt))
	B_arr <- array(0, dim = c(n_col, n_col, Tt))
	for (r in seq_along(regime_start)) {
		A_r <- .regime_rand_operator(n_row, rho_levels[r], seed = seed + 10 * r)
		B_r <- .regime_rand_operator(n_col, rho_levels[r], seed = seed + 10 * r + 1)
		for (t in regime_start[r]:regime_end[r]) {
			A_arr[, , t] <- A_r
			B_arr[, , t] <- B_r
		}
	}
	M <- array(rnorm(n_row * n_col * p, 0, 0.3), dim = c(n_row, n_col, p))
	sim <- .regime_simulate_Z(A_arr, B_arr, M, sigma2, n_row, n_col, p, Tt)
	Y <- .regime_observe(sim$Z, family)
	list(name = sprintf("piecewise_breaks_%s", paste(break_idx, collapse = "_")),
		regime_type = "piecewise",
		A_true = A_arr, B_true = B_arr, M_true = M, sigma2 = sigma2,
		Y = Y, Z = sim$Z, Phi = sim$Phi,
		dims = list(n_row = n_row, n_col = n_col, p = p, Tt = Tt),
		family = family, seed = seed, break_idx = break_idx)
}

#' Regime 4: near stability boundary
#'
#' rho(A_t) * rho(B_t) ~ 0.95, close to the explosive regime. Recursive
#' bootstrap should be unstable here; the auto-switch should activate.
regime_near_unstable <- function(n = 8, Tt = 20, p = 1, family = "gaussian",
                                  target_rho = 0.95, sigma2 = 0.2, seed = 4) {
	set.seed(seed)
	n_row <- n_col <- n
	# split target_rho across A and B
	rho_each <- sqrt(target_rho)
	A0 <- .regime_rand_operator(n_row, rho_each, seed = seed + 1)
	B0 <- .regime_rand_operator(n_col, rho_each, seed = seed + 2)
	A_arr <- array(rep(c(A0), Tt), dim = c(n_row, n_row, Tt))
	B_arr <- array(rep(c(B0), Tt), dim = c(n_col, n_col, Tt))
	M <- array(rnorm(n_row * n_col * p, 0, 0.3), dim = c(n_row, n_col, p))
	sim <- .regime_simulate_Z(A_arr, B_arr, M, sigma2, n_row, n_col, p, Tt)
	Y <- .regime_observe(sim$Z, family)
	list(name = sprintf("near_unstable_rho_%g", target_rho),
		regime_type = "near_unstable",
		A_true = A_arr, B_true = B_arr, M_true = M, sigma2 = sigma2,
		Y = Y, Z = sim$Z, Phi = sim$Phi,
		dims = list(n_row = n_row, n_col = n_col, p = p, Tt = Tt),
		family = family, seed = seed, target_rho = target_rho)
}

#' Run all regimes once for a given family (smoke test)
all_regimes_smoke <- function(family = "gaussian", n = 8, Tt = 15, seed = 1) {
	regimes <- list(
		regime_constant(n = n, Tt = Tt, family = family, seed = seed),
		regime_smooth_rw(n = n, Tt = Tt, family = family, seed = seed + 1),
		regime_piecewise(n = n, Tt = Tt, family = family, seed = seed + 2),
		regime_near_unstable(n = n, Tt = Tt, family = family, seed = seed + 3)
	)
	names(regimes) <- sapply(regimes, function(r) r$name)
	regimes
}

# Quick sanity check when sourced directly
if (sys.nframe() == 0L) {
	cat("Regime generators loaded. Quick test:\n")
	for (fam in c("gaussian", "binary", "ordinal")) {
		cat("  family:", fam, "\n")
		reg <- regime_constant(n = 6, Tt = 10, family = fam)
		cat("    constant: dim Y =", paste(dim(reg$Y), collapse = "x"),
			"  unique Y values =", length(unique(c(reg$Y))), "\n")
	}
}
