####
# cross-sectional data tests
####

test_that("static model works with T=1 data", {
	set.seed(100)
	Y = array(sample(1:5, 8 * 8, replace = TRUE), dim = c(8, 8, 1, 1))
	diag(Y[,,1,1]) = NA

	fit = dbn(Y, model = "static", family = "ordinal",
						 nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "static")
})

test_that("dbn() auto-detects T=1 and uses static model", {
	set.seed(200)
	Y = array(sample(1:5, 8 * 8, replace = TRUE), dim = c(8, 8, 1, 1))
	diag(Y[,,1,1]) = NA

	expect_message(
		fit <- dbn(Y, family = "ordinal", nscan = 30, burn = 10, verbose = FALSE),
		"static"
	)
	expect_equal(fit$model, "static")
})

test_that("dynamic model rejects T=1 data", {
	Y = array(sample(1:5, 8 * 8, replace = TRUE), dim = c(8, 8, 1, 1))
	expect_error(
		dbn(Y, model = "dynamic", family = "ordinal", nscan = 10, burn = 5),
		"at least 2"
	)
})

test_that("hmm model rejects T=1 data", {
	Y = array(sample(1:5, 8 * 8, replace = TRUE), dim = c(8, 8, 1, 1))
	expect_error(
		dbn(Y, model = "hmm", family = "ordinal", nscan = 10, burn = 5),
		"at least 2"
	)
})

test_that("lowrank model rejects T=1 data", {
	Y = array(sample(1:5, 8 * 8, replace = TRUE), dim = c(8, 8, 1, 1))
	expect_error(
		dbn(Y, model = "lowrank", family = "ordinal", nscan = 10, burn = 5),
		"at least 2"
	)
})

test_that("static model also works with T>1 data", {
	sim = simulate_static_dbn(n = 8, p = 1, time = 5, seed = 300)
	fit = dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
})
