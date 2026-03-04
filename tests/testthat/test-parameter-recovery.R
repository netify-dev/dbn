####
# parameter recovery tests
####

test_that("static ordinal model recovers baseline mean M", {
	set.seed(1234)
	sim <- simulate_static_dbn(n = 10, p = 1, time = 15, sigma2 = 0.3,
														 tau2 = 0.05, seed = 1234)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 400, burn = 200, odens = 2, verbose = FALSE)

	# extract posterior mean of M
	M_draws <- fit$draws$misc$M
	M_post <- Reduce("+", M_draws) / length(M_draws)

	# check correlation between true and recovered M
	true_M <- sim$M[,,1]
	cor_M <- cor(as.vector(true_M), as.vector(M_post[,,1]))
	expect_gt(cor_M, 0.3, label = "M recovery correlation")
})

test_that("static gaussian model recovers M and estimates sigma2", {
	set.seed(2345)
	# simulate gaussian data directly
	n <- 10; p <- 1; time <- 15
	M_true <- array(rnorm(n * n * p, 0, 1), dim = c(n, n, p))
	Z <- array(NA, dim = c(n, n, p, time))
	for (t in 1:time) Z[,,1,t] <- M_true[,,1] + rnorm(n^2, 0, 0.5)
	for (t in 1:time) diag(Z[,,1,t]) <- NA

	fit <- dbn(Z, model = "static", family = "gaussian",
						 nscan = 400, burn = 200, odens = 2, verbose = FALSE)

	M_draws <- fit$draws$misc$M
	M_post <- Reduce("+", M_draws) / length(M_draws)
	cor_M <- cor(as.vector(M_true[,,1]), as.vector(M_post[,,1]), use = "complete.obs")
	expect_gt(cor_M, 0.5, label = "Gaussian M recovery")
})

test_that("static binary model runs and produces valid output", {
	set.seed(3456)
	n <- 10; p <- 1; time <- 10
	# simulate binary data
	prob_mat <- array(0.3, dim = c(n, n, p, time))
	Y <- array(rbinom(n * n * p * time, 1, prob_mat), dim = c(n, n, p, time))
	for (t in 1:time) diag(Y[,,1,t]) <- NA

	fit <- dbn(Y, model = "static", family = "binary",
						 nscan = 100, burn = 50, odens = 2, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$family, "binary")
})

test_that("dynamic ordinal model recovers M structure", {
	set.seed(4567)
	sim <- simulate_dynamic_dbn(n = 8, p = 1, time = 10, sigma2 = 0.3,
															tauA2 = 0.02, tauB2 = 0.02, seed = 4567)
	fit <- dbn(sim$Y, model = "dynamic", family = "ordinal",
						 nscan = 300, burn = 150, odens = 2, verbose = FALSE)

	# check output structure
	expect_equal(fit$model, "dynamic")
	expect_equal(fit$dims$n_row, 8)
	expect_equal(fit$dims$n_col, 8)
	expect_equal(length(fit$sigma2), 150)

	# check A is stored and has correct dims
	expect_true(length(fit$A) > 0)
	expect_equal(dim(fit$A[[1]])[1], 8)
	expect_equal(dim(fit$A[[1]])[2], 8)
})

test_that("dynamic gaussian model estimates observation variance", {
	set.seed(5678)
	n <- 8; p <- 1; time <- 10
	M_true <- array(rnorm(n * n * p, 0, 0.5), dim = c(n, n, p))
	Z <- array(NA, dim = c(n, n, p, time))
	for (t in 1:time) Z[,,1,t] <- M_true[,,1] + rnorm(n^2, 0, 0.5)
	for (t in 1:time) diag(Z[,,1,t]) <- NA

	fit <- dbn(Z, model = "dynamic", family = "gaussian",
						 nscan = 200, burn = 100, odens = 2, verbose = FALSE)

	expect_true(!is.null(fit$sigma2_obs))
	expect_true(all(fit$sigma2_obs > 0))
	expect_true(all(is.finite(fit$sigma2_obs)))
})
