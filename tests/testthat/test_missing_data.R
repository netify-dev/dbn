#### missing-data handling

test_that("dbn handles NA entries beyond the diagonal", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 6, time = 8, seed = 1)
	Y = sim$Y
	# inject NAs at off-diagonal cells
	Y[1, 3, 1, 2:4] = NA
	Y[5, 2, 1, 6] = NA
	fit = dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 30, burn = 10, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	# the original NAs in fit$Y should still be NA (the sampler imputes
	# Theta but Y stays as observed).
	expect_true(is.na(fit$Y[1, 3, 1, 2]))
})
