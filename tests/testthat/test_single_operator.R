# symmetric = TRUE path: matrix-free PCG on Theta, sufficient-stat
# row-FFBS on A, triangle-only RSS, conjugate IG on tauA2. enforces
# B = A. requires unipartite square (n_row == n_col, p == 1).

test_that("symmetric = TRUE requires unipartite square", {
	skip_on_cran()
	# bipartite: should error
	Y_bip <- array(rnorm(5 * 6 * 1 * 8), c(5, 6, 1, 8))
	expect_error(
		dbn_dynamic(Y_bip, family = "gaussian", nscan = 10, burn = 5,
			symmetric = TRUE, verbose = FALSE),
		"unipartite square"
	)
	# multilayer (p > 1): should error
	Y_p2 <- array(rnorm(6 * 6 * 2 * 8), c(6, 6, 2, 8))
	expect_error(
		dbn_dynamic(Y_p2, family = "gaussian", nscan = 10, burn = 5,
			symmetric = TRUE, verbose = FALSE),
		"unipartite square"
	)
})

test_that("symmetric = TRUE runs end-to-end on small Gaussian data", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 6, time = 12, sigma2 = 0.3,
		tauA2 = 0.05, symmetric = TRUE, seed = 13,
		return_truth = TRUE)
	Y <- sim$Y[, , 1, , drop = FALSE]

	fit <- dbn_dynamic(Y, family = "gaussian",
		nscan = 80, burn = 40, odens = 5, time_thin = 1,
		symmetric = TRUE, seed = 1, verbose = FALSE)

	expect_s3_class(fit, "dbn")
	# A and B are stored identically under B = A
	expect_true(identical(fit$A[[1]], fit$B[[1]]))
})

test_that("tau_A_fixed skips the conjugate IG update", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 6, time = 10, sigma2 = 0.3,
		tauA2 = 0.05, symmetric = TRUE, seed = 7,
		return_truth = TRUE)
	Y <- sim$Y[, , 1, , drop = FALSE]

	fit <- dbn_dynamic(Y, family = "gaussian",
		nscan = 40, burn = 20, odens = 5, time_thin = 1,
		symmetric = TRUE, tau_A_fixed = 0.01, seed = 2, verbose = FALSE)

	# every stored tau_A2 draw equals the fixed value
	expect_true(all(abs(fit$tau_A2 - 0.01) < 1e-12))
})

test_that("tau_A_fixed validates input", {
	Y <- array(rnorm(6 * 6 * 1 * 8), c(6, 6, 1, 8))
	for (t in 1:8) {
		s <- Y[, , 1, t]; Y[, , 1, t] <- 0.5 * (s + t(s)); diag(Y[, , 1, t]) <- 0
	}
	expect_error(
		dbn_dynamic(Y, family = "gaussian", nscan = 10, burn = 5,
			symmetric = TRUE, tau_A_fixed = -0.1, verbose = FALSE),
		"positive finite"
	)
	expect_error(
		dbn_dynamic(Y, family = "gaussian", nscan = 10, burn = 5,
			symmetric = TRUE, tau_A_fixed = c(0.1, 0.2), verbose = FALSE),
		"single positive finite"
	)
})

# diagonal-penalty (lambda_diag > 0) path: A_t = Abar + Delta_t with
# stationary AR(1) on Delta_t. checks Abar storage shape, the IG update
# for tauA2 doesn't blow up, and warm-start preserves the Abar trajectory
# and the adapted M-H proposal sd (regression for the warm-start reset
# bug fixed in 1.1.0).

test_that("lambda_diag > 0 produces Abar draws and adapted proposal sd", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 6, time = 10, sigma2 = 0.3,
		tauA2 = 0.05, symmetric = TRUE, seed = 27,
		return_truth = TRUE)
	Y <- sim$Y[, , 1, , drop = FALSE]

	fit <- dbn_dynamic(Y, family = "gaussian",
		nscan = 60, burn = 30, odens = 5, time_thin = 1,
		symmetric = TRUE, lambda_diag = 0.5, seed = 11,
		verbose = FALSE)

	# Abar trajectory is stored, one entry per kept draw
	expect_false(is.null(fit$Abar))
	expect_equal(length(fit$Abar), length(fit$sigma2))
	# Abar is symmetric n x n
	expect_equal(dim(fit$Abar[[1]]), c(6L, 6L))
	expect_true(all(abs(fit$Abar[[1]] - t(fit$Abar[[1]])) < 1e-12))
	# adapted MH proposal sd exposed for warm start
	expect_true(is.numeric(fit$mh_proposal_sd))
	expect_true(fit$mh_proposal_sd > 0)
	# tauA2 stays positive and finite
	expect_true(all(is.finite(fit$tau_A2) & fit$tau_A2 > 0))
})

test_that("lambda_diag > 0 warm-start restores Abar_state and mh_proposal_sd", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 6, time = 10, sigma2 = 0.3,
		tauA2 = 0.05, symmetric = TRUE, seed = 27,
		return_truth = TRUE)
	Y <- sim$Y[, , 1, , drop = FALSE]

	fit1 <- dbn_dynamic(Y, family = "gaussian",
		nscan = 80, burn = 40, odens = 5, time_thin = 1,
		symmetric = TRUE, lambda_diag = 0.5, seed = 19,
		verbose = FALSE)

	# continuation should use the last Abar / proposal sd, not reset them
	fit2 <- dbn_dynamic(Y, family = "gaussian",
		nscan = 40, burn = 0, odens = 5, time_thin = 1,
		symmetric = TRUE, lambda_diag = 0.5, seed = 19,
		previous = fit1, verbose = FALSE)

	# proposal sd is restored from the last fit
	expect_equal(fit2$mh_proposal_sd, fit1$mh_proposal_sd, tolerance = 1e-12)
	# Abar starting point is the last Abar of fit1 (warm-start uses it as
	# the initial state, and the first stored draw is one MH+Gibbs sweep
	# away, so we check the magnitude is comparable not exact)
	abar_end_fit1 <- fit1$Abar[[length(fit1$Abar)]]
	abar_start_fit2 <- fit2$Abar[[1]]
	expect_lt(
		max(abs(abar_start_fit2 - abar_end_fit1)),
		2 * sd(as.numeric(abar_end_fit1)) + 1
	)
})
