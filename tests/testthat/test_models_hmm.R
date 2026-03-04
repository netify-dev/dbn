####
# HMM model
####

test_that("simulate_hmm_dbn produces valid output structure", {
	sim <- simulate_hmm_dbn(n = 8, p = 1, time = 15, R = 2, seed = 42)
	expect_true(is.list(sim))
	expect_true("Y" %in% names(sim))
	expect_equal(dim(sim$Y)[1], 8)
	expect_equal(dim(sim$Y)[2], 8)
	expect_equal(dim(sim$Y)[3], 1)
	expect_equal(dim(sim$Y)[4], 15)
	# ordinal: positive integers
	y_vals <- sim$Y[!is.na(sim$Y)]
	expect_true(all(y_vals == round(y_vals)))
	expect_true(all(y_vals > 0))
})

test_that("simulate_hmm_dbn works with multiple relations", {
	sim <- simulate_hmm_dbn(n = 6, p = 2, time = 10, R = 2, seed = 43)
	expect_equal(dim(sim$Y)[3], 2)
})

test_that("dbn_hmm fits ordinal data and returns expected object", {
	sim <- simulate_hmm_dbn(n = 8, p = 1, time = 10, R = 2, seed = 44)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 nscan = 100, burn = 50, verbose = FALSE, R = 2)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "hmm")
	# required components
	expect_true(!is.null(fit$A))
	expect_true(!is.null(fit$B))
	expect_true(!is.null(fit$sigma2))
	expect_true(!is.null(fit$Pi))
	expect_true(!is.null(fit$S))
})

test_that("dbn_hmm fits gaussian data", {
	set.seed(45)
	Y <- array(rnorm(8 * 8 * 1 * 10), dim = c(8, 8, 1, 10))
	for (t in 1:10) diag(Y[,,1,t]) <- NA
	fit <- dbn(Y, model = "hmm", family = "gaussian",
						 nscan = 100, burn = 50, verbose = FALSE, R = 2)
	expect_s3_class(fit, "dbn")
	expect_true(!is.null(fit$draws$pars$sigma2_obs))
})

test_that("dbn_hmm fits binary data", {
	set.seed(46)
	Y <- array(rbinom(8 * 8 * 1 * 10, 1, 0.3), dim = c(8, 8, 1, 10))
	for (t in 1:10) diag(Y[,,1,t]) <- NA
	fit <- dbn(Y, model = "hmm", family = "binary",
						 nscan = 100, burn = 50, verbose = FALSE, R = 2)
	expect_s3_class(fit, "dbn")
})

test_that("posterior_predict works for HMM model", {
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 47)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 nscan = 80, burn = 40, verbose = FALSE, R = 2)
	ppd <- posterior_predict_dbn(fit, ndraws = 5, seed = 1)
	expect_s3_class(ppd, "dbn_ppd")
	expect_equal(length(ppd), 5)
})

test_that("summary and print work for HMM model", {
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 48)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 nscan = 80, burn = 40, verbose = FALSE, R = 2)
	# should not error
	expect_output(print(fit))
	expect_output(summary(fit))
})

test_that("HMM regime probabilities are valid", {
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 10, R = 2, seed = 49)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 nscan = 100, burn = 50, verbose = FALSE, R = 2)
	# Pi: list of transition matrices
	expect_true(is.list(fit$Pi))
	Pi_1 <- fit$Pi[[1]]
	# rows sum to 1
	expect_true(all(abs(rowSums(Pi_1) - 1) < 1e-6))
	# non-negative
	expect_true(all(Pi_1 >= 0))
})

test_that("predict_hmm returns forecasts", {
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 55)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 nscan = 80, burn = 40, verbose = FALSE, R = 2)
	pred <- predict(fit, H = 3, draws = 5)
	expect_true(!is.null(pred))
})

test_that("HMM rejects T < 2", {
	set.seed(50)
	Y <- array(rnorm(8 * 8 * 1 * 1), dim = c(8, 8, 1, 1))
	expect_error(dbn(Y, model = "hmm", family = "ordinal", R = 2),
							 "at least 2 time points")
})
