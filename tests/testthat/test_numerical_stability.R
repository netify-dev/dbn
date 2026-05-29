# tiny variance does not crash
test_that("tiny variance does not crash", {
	skip_on_cran()
	Y = array(rnorm(5 * 5 * 1 * 5, mean = 0, sd = 1e-6), dim = c(5, 5, 1, 5))
	for (t in 1:5) diag(Y[, , 1, t]) = NA
	fit = suppressWarnings(suppressMessages(
		dbn(
			Y,
			model = "static",
			family = "gaussian",
			nscan = 30,
			burn = 10,
			verbose = FALSE
		)
	))
	expect_s3_class(fit, "dbn")
})

# imbalanced binary data runs without crash
test_that("imbalanced binary data runs", {
	skip_on_cran()
	Y = array(rbinom(8 * 8 * 1 * 8, 1, 0.05), dim = c(8, 8, 1, 8))
	for (t in 1:8) diag(Y[, , 1, t]) = NA
	fit = suppressWarnings(suppressMessages(
		dbn(
			Y,
			model = "dynamic",
			family = "binary",
			nscan = 50,
			burn = 20,
			verbose = FALSE
		)
	))
	expect_s3_class(fit, "dbn")
})
