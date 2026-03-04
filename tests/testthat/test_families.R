#### model x family combinations

# helper to create gaussian data
make_gaussian_data <- function(n = 8, p = 1, time = 10, seed = 100) {
	set.seed(seed)
	Y <- array(rnorm(n * n * p * time), dim = c(n, n, p, time))
	for (t in 1:time) for (r in 1:p) diag(Y[,,r,t]) <- NA
	Y
}

# helper to create binary data
make_binary_data <- function(n = 8, p = 1, time = 10, seed = 100) {
	set.seed(seed)
	Y <- array(rbinom(n * n * p * time, 1, 0.3), dim = c(n, n, p, time))
	for (t in 1:time) for (r in 1:p) diag(Y[,,r,t]) <- NA
	Y
}

#### static

test_that("static x ordinal works", {
	sim <- simulate_static_dbn(n = 8, p = 1, time = 8, seed = 10)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "static")
	expect_equal(fit$family, "ordinal")
})

test_that("static x gaussian works", {
	Y <- make_gaussian_data(n = 8, time = 8, seed = 20)
	fit <- dbn(Y, model = "static", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$family, "gaussian")
})

test_that("static x binary works", {
	Y <- make_binary_data(n = 8, time = 8, seed = 30)
	fit <- dbn(Y, model = "static", family = "binary",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$family, "binary")
})

#### dynamic

test_that("dynamic x ordinal works", {
	sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 40)
	fit <- dbn(sim$Y, model = "dynamic", family = "ordinal",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "dynamic")
	expect_true(length(fit$A) > 0)
})

test_that("dynamic x gaussian works", {
	Y <- make_gaussian_data(n = 6, time = 5, seed = 50)
	fit <- dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_true(!is.null(fit$sigma2_obs))
})

test_that("dynamic x binary works", {
	Y <- make_binary_data(n = 6, time = 5, seed = 60)
	fit <- dbn(Y, model = "dynamic", family = "binary",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_true(!is.null(fit$Z))
})

#### lowrank

test_that("lowrank x ordinal works", {
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 8, r = 2, seed = 70)
	fit <- dbn(sim$Y, model = "lowrank", family = "ordinal",
		nscan = 50, burn = 20, verbose = FALSE, r = 2)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "lowrank")
})

test_that("lowrank x gaussian works", {
	Y <- make_gaussian_data(n = 8, time = 8, seed = 80)
	fit <- dbn(Y, model = "lowrank", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE, r = 2)
	expect_s3_class(fit, "dbn")
})

test_that("lowrank x binary works", {
	Y <- make_binary_data(n = 8, time = 8, seed = 90)
	fit <- dbn(Y, model = "lowrank", family = "binary",
		nscan = 50, burn = 20, verbose = FALSE, r = 2)
	expect_s3_class(fit, "dbn")
})

#### HMM

test_that("hmm x ordinal works", {
	sim <- simulate_hmm_dbn(n = 8, p = 1, time = 10, R = 2, seed = 100)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
		nscan = 50, burn = 20, verbose = FALSE, R = 2)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "hmm")
})

test_that("hmm x gaussian works", {
	Y <- make_gaussian_data(n = 8, time = 10, seed = 110)
	fit <- dbn(Y, model = "hmm", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE, R = 2)
	expect_s3_class(fit, "dbn")
})

test_that("hmm x binary works", {
	Y <- make_binary_data(n = 8, time = 10, seed = 120)
	fit <- dbn(Y, model = "hmm", family = "binary",
		nscan = 50, burn = 20, verbose = FALSE, R = 2)
	expect_s3_class(fit, "dbn")
})
