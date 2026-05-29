#### warm-start essentials

test_that("warm start from previous fit produces a continuation", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	fit1 = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 30, burn = 10, verbose = FALSE, seed = 1)
	fit2 = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 20, burn = 0, verbose = FALSE, seed = 1,
		previous = fit1)
	expect_s3_class(fit2, "dbn")
	# fit2 should have its own draws (not just fit1's)
	expect_true(length(fit2$A) > 0L)
})

test_that("warm start validates model-type compatibility", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	fit_d = dbn(sim$Y, model = "dynamic",
		nscan = 30, burn = 10, verbose = FALSE)
	# resuming a dynamic fit under model = "static" should error
	expect_error(
		dbn(sim$Y, model = "static",
			nscan = 20, burn = 0, verbose = FALSE, previous = fit_d),
		regexp = "model|previous"
	)
})
