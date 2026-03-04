####
# symmetric network (B = A constraint)
####

####
# simulation tests
####

test_that("simulate_static_dbn with symmetric=TRUE produces A == B", {
	sim <- simulate_static_dbn(n = 8, p = 1, time = 5, symmetric = TRUE, seed = 100)
	expect_equal(sim$A, sim$B)
	expect_true(sim$is_symmetric)
})

test_that("simulate_dynamic_dbn with symmetric=TRUE produces matching A and B", {
	sim <- simulate_dynamic_dbn(n = 8, p = 1, time = 10, symmetric = TRUE, seed = 101)
	for (t in 1:10) {
		expect_equal(sim$A[, , t], sim$B[, , t])
	}
	expect_true(sim$is_symmetric)
})

test_that("simulate_hmm_dbn with symmetric=TRUE produces matching A and B lists", {
	sim <- simulate_hmm_dbn(n = 8, p = 1, time = 10, R = 2,
		symmetric = TRUE, seed = 102)
	for (r in 1:2) {
		expect_equal(sim$A_list[[r]], sim$B_list[[r]])
	}
	expect_true(sim$is_symmetric)
})

test_that("simulate functions error when symmetric + bipartite", {
	expect_error(
		simulate_static_dbn(n = 10, n_col = 8, symmetric = TRUE),
		"n_row == n_col"
	)
	expect_error(
		simulate_dynamic_dbn(n = 10, n_col = 8, symmetric = TRUE),
		"n_row == n_col"
	)
	expect_error(
		simulate_hmm_dbn(n = 10, n_col = 8, symmetric = TRUE),
		"n_row == n_col"
	)
})

####
# model fitting tests
####

test_that("dbn static + symmetric + ordinal runs and stores flag", {
	sim <- simulate_static_dbn(n = 8, p = 1, time = 5,
		symmetric = TRUE, seed = 110)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 100, burn = 50, verbose = FALSE, symmetric = TRUE)
	expect_s3_class(fit, "dbn")
	expect_true(fit$dims$is_symmetric)
})

test_that("dbn dynamic + symmetric + ordinal enforces A == B in output", {
	sim <- simulate_dynamic_dbn(n = 8, p = 1, time = 10,
		symmetric = TRUE, seed = 111)
	fit <- dbn(sim$Y, model = "dynamic", family = "ordinal",
		nscan = 100, burn = 50, verbose = FALSE, symmetric = TRUE)
	expect_s3_class(fit, "dbn")
	expect_true(fit$dims$is_symmetric)
	# A == B in stored draws
	for (s in seq_along(fit$A)) {
		expect_equal(fit$A[[s]], fit$B[[s]], tolerance = 1e-10)
	}
})

test_that("dbn dynamic + symmetric + gaussian works", {
	set.seed(112)
	Y <- array(rnorm(8 * 8 * 1 * 10), dim = c(8, 8, 1, 10))
	for (t in 1:10) diag(Y[, , 1, t]) <- NA
	fit <- dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 100, burn = 50, verbose = FALSE, symmetric = TRUE)
	expect_s3_class(fit, "dbn")
	for (s in seq_along(fit$A)) {
		expect_equal(fit$A[[s]], fit$B[[s]], tolerance = 1e-10)
	}
})

test_that("dbn dynamic + symmetric + binary works", {
	skip_if_not_installed("truncnorm")
	set.seed(113)
	Y <- array(rbinom(8 * 8 * 1 * 10, 1, 0.3), dim = c(8, 8, 1, 10))
	for (t in 1:10) diag(Y[, , 1, t]) <- NA
	fit <- dbn(Y, model = "dynamic", family = "binary",
		nscan = 100, burn = 50, verbose = FALSE, symmetric = TRUE)
	expect_s3_class(fit, "dbn")
	for (s in seq_along(fit$A)) {
		expect_equal(fit$A[[s]], fit$B[[s]], tolerance = 1e-10)
	}
})

test_that("dbn hmm + symmetric enforces A == B per regime", {
	sim <- simulate_hmm_dbn(n = 8, p = 1, time = 10, R = 2,
		symmetric = TRUE, seed = 114)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
		nscan = 100, burn = 50, verbose = FALSE, R = 2, symmetric = TRUE)
	expect_s3_class(fit, "dbn")
	expect_true(fit$dims$is_symmetric)
	# A == B per regime per draw
	for (s in seq_along(fit$A)) {
		R_dim <- dim(fit$A[[s]])[3]
		for (r in 1:R_dim) {
			expect_equal(fit$A[[s]][, , r], fit$B[[s]][, , r], tolerance = 1e-10)
		}
	}
})

####
# error handling
####

test_that("symmetric + bipartite errors in dbn()", {
	sim <- simulate_dynamic_dbn(n = 10, n_col = 8, p = 1, time = 10, seed = 120)
	expect_error(
		dbn(sim$Y, model = "dynamic", symmetric = TRUE, verbose = FALSE),
		"equal sender and receiver"
	)
})

test_that("symmetric + lowrank errors in dbn()", {
	sim <- simulate_static_dbn(n = 8, p = 1, time = 10, seed = 121)
	expect_error(
		dbn(sim$Y, model = "lowrank", r = 2, symmetric = TRUE, verbose = FALSE),
		"not yet supported"
	)
})

test_that("symmetric tau_A2 == tau_B2 in dynamic model", {
	sim <- simulate_dynamic_dbn(n = 8, p = 1, time = 10,
		symmetric = TRUE, seed = 130)
	fit <- dbn(sim$Y, model = "dynamic", family = "ordinal",
		nscan = 100, burn = 50, verbose = FALSE, symmetric = TRUE)
	expect_equal(fit$tau_A2, fit$tau_B2)
})
