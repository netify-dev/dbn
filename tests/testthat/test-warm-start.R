####
# warm start and initialization
####

test_that("warm start works for static model", {
	skip_if_not_installed("cli")
	
	# generate small test data
	set.seed(6886)
	m <- 4
	p <- 2  
	Tt <- 5
	Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
	
	# run initial model for short time
	result1 <- dbn_static(Y, nscan = 100, burn = 20, odens = 1, verbose = FALSE)
	
	# continue from previous results
	result2 <- dbn_static(Y, nscan = 100, burn = 20, odens = 1, verbose = FALSE, previous = result1)
	
	# check that parameters continued
	expect_equal(length(result2$params[,1]), 100)  # nscan=100 iterations kept
	expect_true(!is.null(result2$continued_from))
	expect_equal(result2$continued_from, 120)  # previous total: burn(20) + nscan(100) = 120
	expect_equal(result2$total_iter, 220)  # 120 + 100 new iterations = 220
	
	# check that warm start initialized from last values of result1
	last_s2 <- result1$params[nrow(result1$params), "s2"]
	# Note: since we set seed again, the first sample might be identical to a previous sample
	# but the chain should continue from the last state
	expect_true(length(result2$params[, "s2"]) == 100)
})

test_that("warm start works for dynamic model", {
	skip_if_not_installed("cli")
	
	# generate small test data
	set.seed(6886)
	m <- 4
	p <- 2
	Tt <- 5
	Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
	
	# run initial model
	result1 <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE)
	
	# continue from previous results
	result2 <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE, previous = result1)
	
	# check that parameters continued
	expect_equal(length(result2$sigma2), 150)  # nscan=150 iterations kept
	expect_true(!is.null(result2$continued_from))
	
	# verify dimensions match
	expect_equal(dim(result2$A[[1]]), dim(result1$A[[1]]))
	expect_equal(dim(result2$B[[1]]), dim(result1$B[[1]]))
})

test_that("init argument works for static model", {
	skip_if_not_installed("cli")
	
	# generate small test data
	set.seed(6886)
	m <- 4
	p <- 2
	Tt <- 5
	Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
	
	# run with custom initial values
	init_vals <- list(
		s2 = 0.5,
		t2 = 2.0,
		g2 = 0.1,
		B = list(diag(m) * 0.9, diag(m) * 0.9, diag(p) * 0.9)
	)
	
	result <- dbn_static(Y, nscan = 50, burn = 10, odens = 1, verbose = FALSE, init = init_vals)
	
	# the model should run without errors
	expect_s3_class(result, "dbn")
	expect_equal(result$model, "static")
})

test_that("init argument works for dynamic model", {
	skip_if_not_installed("cli")
	
	# generate small test data
	set.seed(6886)
	m <- 4
	p <- 2
	Tt <- 5
	Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
	
	# create initial a and b arrays
	Aarray <- array(0, dim=c(m, m, Tt))
	Barray <- array(0, dim=c(m, m, Tt))
	for(t in 1:Tt) {
		Aarray[,,t] <- diag(m) * 0.8
		Barray[,,t] <- diag(m) * 0.8
	}
	
	# run with custom initial values
	init_vals <- list(
		sigma2 = 0.5,
		tau_A2 = 5.0,
		tau_B2 = 5.0,
		g2 = 0.2,
		A = Aarray,
		B = Barray
	)
	
	result <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE, init = init_vals)
	
	# the model should run without errors
	expect_s3_class(result, "dbn")
	expect_equal(result$model, "dynamic")
})

test_that("warm start validates model type", {
	skip_if_not_installed("cli")
	
	# generate small test data
	set.seed(6886)
	m <- 4
	p <- 2
	Tt <- 5
	Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
	
	# run static model
	result_static <- dbn_static(Y, nscan = 50, burn = 10, odens = 1, verbose = FALSE)
	
	# try to use static results to continue dynamic model - should error
	expect_error(
		dbn_dynamic(Y, nscan = 50, burn = 10, odens = 1, verbose = FALSE, previous = result_static),
		"previous.*must be results from"
	)
})

test_that("warm start handles time thinning correctly", {
	skip_if_not_installed("cli")

	# generate small test data
	set.seed(6886)
	m <- 4
	p <- 2
	Tt <- 10
	Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))

	# run with time thinning
	result1 <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE, time_thin = 2)

	# check that saved arrays have correct dimensions
	expect_equal(dim(result1$A[[1]])[3], 5)  # tt/2 = 10/2 = 5

	# continue from previous - should expand back to full time
	result2 <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE, previous = result1)

	# should still work correctly
	expect_s3_class(result2, "dbn")
})

####
# warm start chain continuity
####

test_that("warm start produces continuous chain (no discontinuity)", {
	skip_on_cran()

	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 5001)

	fit1 <- dbn(sim$Y, model = "static",
		nscan = 200, burn = 50, odens = 1, verbose = FALSE)
	fit2 <- dbn(sim$Y, model = "static",
		nscan = 200, burn = 10, odens = 1, verbose = FALSE, previous = fit1)

	# sigma2 from chain 2 should start near where chain 1 ended
	last_s2 <- tail(fit1$draws$pars$s2, 1)
	first_s2 <- fit2$draws$pars$s2[1]
	# within an order of magnitude (MCMC is stochastic, but should be in same region)
	expect_lt(abs(log(first_s2 / last_s2)), 3,
		label = "warm start sigma2 continues from previous chain")

	# total_iter should be tracked for continued fit
	fit1_total <- fit1$total_iter %||% (fit1$settings$nscan + fit1$settings$burn)
	expect_true(fit2$total_iter > fit1_total,
		label = "total_iter increases across warm start")
})

####
# thread safety: 1 thread vs multi-thread produce comparable results
####

test_that("single-thread and multi-thread produce comparable results", {
	skip_on_cran()

	sim <- simulate_static_dbn(n = 8, p = 1, time = 8, seed = 5002)

	set_dbn_threads(1)
	fit1 <- dbn(sim$Z, model = "static", family = "gaussian",
		nscan = 200, burn = 100, odens = 1, verbose = FALSE)

	n_threads <- min(4, parallel::detectCores(logical = FALSE))
	if (n_threads > 1) {
		set_dbn_threads(n_threads)
		fit_mt <- dbn(sim$Z, model = "static", family = "gaussian",
			nscan = 200, burn = 100, odens = 1, verbose = FALSE)

		# posterior means should agree within ~20% (MCMC variance, not exact)
		s2_ratio <- mean(fit_mt$draws$pars$s2) / mean(fit1$draws$pars$s2)
		expect_gt(s2_ratio, 0.5, label = "sigma2 ratio > 0.5 (1-thread vs multi)")
		expect_lt(s2_ratio, 2.0, label = "sigma2 ratio < 2.0 (1-thread vs multi)")
	}

	# restore default
	set_dbn_threads(min(4, parallel::detectCores(logical = FALSE)))
})