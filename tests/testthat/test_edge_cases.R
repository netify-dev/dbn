####
# edge case tests
####

test_that("burn >= nscan errors", {
	sim <- simulate_static_dbn(n = 5, time = 5, seed = 101)
	expect_error(
		dbn(sim$Y, model = "static", nscan = 50, burn = 60, verbose = FALSE),
		"burn.*must be less than.*nscan"
	)
})

test_that("burn == nscan errors", {
	sim <- simulate_static_dbn(n = 5, time = 5, seed = 102)
	expect_error(
		dbn(sim$Y, model = "static", nscan = 50, burn = 50, verbose = FALSE),
		"burn.*must be less than.*nscan"
	)
})

test_that("odens too large errors", {
	sim <- simulate_static_dbn(n = 5, time = 5, seed = 103)
	expect_error(
		dbn(sim$Y, model = "static", nscan = 50, burn = 20, odens = 100, verbose = FALSE),
		"odens.*too large"
	)
})

test_that("r > n errors for lowrank", {
	sim <- simulate_lowrank_dbn(n = 5, time = 5, r = 2, seed = 104)
	expect_error(
		dbn(sim$Y, model = "lowrank", r = 10, nscan = 30, burn = 10, verbose = FALSE),
		"cannot exceed"
	)
})

test_that("R = 1 HMM works (single regime)", {
	sim <- simulate_dynamic_dbn(n = 5, time = 5, seed = 105)
	fit <- dbn(sim$Y, model = "hmm", R = 1, nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "hmm")
})

test_that("all-constant Y matrix runs for gaussian family", {
	Y <- array(3, dim = c(5, 5, 1, 5))
	for (t in 1:5) diag(Y[,,1,t]) <- NA
	# constant data is valid for gaussian (uninformative)
	fit <- dbn(Y, model = "static", family = "gaussian", nscan = 30, burn = 10, verbose = FALSE)
	expect_s3_class(fit, "dbn")
})

test_that("binary non-0/1 data is rejected", {
	Y <- array(c(0, 1, 2, 3), dim = c(5, 5, 1, 3))
	for (t in 1:3) diag(Y[,,1,t]) <- NA
	expect_error(
		dbn(Y, model = "static", family = "binary", nscan = 30, burn = 10, verbose = FALSE),
		"Binary family requires 0/1"
	)
})

test_that("ordinal non-integer data warns", {
	Y <- array(runif(5 * 5 * 1 * 3), dim = c(5, 5, 1, 3))
	for (t in 1:3) diag(Y[,,1,t]) <- NA
	expect_warning(
		dbn(Y, model = "static", family = "ordinal", nscan = 30, burn = 10, verbose = FALSE),
		"integer-valued"
	)
})

test_that("predict.dbn warns on unrecognized arguments", {
	sim <- simulate_static_dbn(n = 6, time = 5, seed = 106)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 50, burn = 20, verbose = FALSE)
	expect_warning(
		predict(fit, sumary = "mean"),
		"Unrecognized"
	)
})

test_that("negative nscan errors", {
	sim <- simulate_static_dbn(n = 5, time = 5, seed = 107)
	expect_error(
		dbn(sim$Y, model = "static", nscan = -10, verbose = FALSE),
		"nscan.*must be positive"
	)
})

test_that("verbose = TRUE works like verbose = 100", {
	sim <- simulate_static_dbn(n = 5, time = 3, seed = 108)
	# verbose = TRUE should not error
	expect_no_error(
		dbn(sim$Y, model = "static", nscan = 50, burn = 20, verbose = TRUE)
	)
})

####
# bipartite guards for HMM and lowrank
####

test_that("bipartite data errors for HMM model", {
	Y <- array(rnorm(4 * 6 * 1 * 5), dim = c(4, 6, 1, 5))
	expect_error(
		dbn(Y, model = "hmm", R = 2, nscan = 30, burn = 10, verbose = FALSE),
		"bipartite"
	)
})

test_that("bipartite data errors for lowrank model", {
	Y <- array(rnorm(4 * 6 * 1 * 5), dim = c(4, 6, 1, 5))
	expect_error(
		dbn(Y, model = "lowrank", r = 2, nscan = 30, burn = 10, verbose = FALSE),
		"bipartite"
	)
})

test_that("dynamic binary with small n produces warning", {
	set.seed(201)
	n <- 8
	Y <- array(sample(0:1, n * n * 1 * 5, replace = TRUE), dim = c(n, n, 1, 5))
	for (t in 1:5) diag(Y[,,1,t]) <- NA
	expect_warning(
		dbn(Y, model = "dynamic", family = "binary",
			nscan = 30, burn = 10, verbose = FALSE),
		"singularit"
	)
})

####
# additional dimension validation
####

test_that("2D matrix input errors", {
	Y <- matrix(rnorm(25), 5, 5)
	expect_error(
		dbn(Y, model = "static", nscan = 30, burn = 10, verbose = FALSE)
	)
})

test_that("minimal valid run: nscan=2, burn=1, odens=1", {
	sim <- simulate_static_dbn(n = 5, time = 3, seed = 110)
	fit <- dbn(sim$Y, model = "static", nscan = 2, burn = 1, odens = 1, verbose = FALSE)
	expect_s3_class(fit, "dbn")
})

test_that("lowrank with r = 0 errors", {
	sim <- simulate_lowrank_dbn(n = 5, time = 5, r = 2, seed = 111)
	expect_error(
		dbn(sim$Y, model = "lowrank", r = 0, nscan = 30, burn = 10, verbose = FALSE),
		"at least 1"
	)
})
