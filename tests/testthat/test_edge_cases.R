#### input validation and edge cases

# helper: a small, valid Y array
.mk_Y = function(n = 6, time = 8) {
	array(rnorm(n * n * 1 * time), dim = c(n, n, 1, time))
}

test_that("mcmc control arguments are validated", {
	Y = .mk_Y()
	expect_error(dbn(Y, model = "static", nscan = 10, burn = 10, verbose = FALSE),
		regexp = "burn")
	expect_error(dbn(Y, model = "static", nscan = 10, burn = 15, verbose = FALSE),
		regexp = "burn")
	expect_error(dbn(Y, model = "static", nscan = 10, burn = 1, odens = 999, verbose = FALSE),
		regexp = "odens")
	expect_error(dbn(Y, model = "static", nscan = -5, burn = 1, verbose = FALSE),
		regexp = "positive|nscan")
})

test_that("malformed data arrays are rejected", {
	expect_error(dbn(matrix(0, 5, 5), model = "static",
		nscan = 10, burn = 5, verbose = FALSE), regexp = "3D|4D|array")
})

test_that("binary family rejects non-{0,1} data", {
	Y = array(rnorm(6 * 6 * 1 * 5, mean = 2), dim = c(6, 6, 1, 5))
	expect_error(
		dbn(Y, model = "static", family = "binary",
			nscan = 10, burn = 5, verbose = FALSE),
		regexp = "binary|0|1"
	)
})

test_that("ordinal family warns on non-integer values", {
	sim = simulate_static_dbn(n = 5, p = 1, time = 5, seed = 1)
	Y_cont = sim$Z  # continuous; user mistakenly tags as ordinal
	expect_warning(
		dbn(Y_cont, model = "static", family = "ordinal",
			nscan = 10, burn = 5, verbose = FALSE),
		regexp = "ordinal|integer"
	)
})

test_that("predict.dbn warns on unrecognised arguments", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 30, burn = 10, verbose = FALSE)
	expect_warning(predict(fit, banana = TRUE), regexp = "Unrecognized|banana")
})

test_that("verbose = TRUE behaves like verbose = 100", {
	sim = simulate_dynamic_dbn(n = 5, time = 5, seed = 1)
	fit_t = suppressMessages(suppressWarnings(
		dbn(sim$Y, model = "dynamic", family = "gaussian",
			nscan = 10, burn = 5, verbose = TRUE)))
	fit_n = suppressMessages(suppressWarnings(
		dbn(sim$Y, model = "dynamic", family = "gaussian",
			nscan = 10, burn = 5, verbose = 100)))
	expect_s3_class(fit_t, "dbn")
	expect_s3_class(fit_n, "dbn")
})

test_that("minimal valid run (nscan=2, burn=1) succeeds", {
	Y = .mk_Y(n = 4, time = 4)
	fit = suppressWarnings(suppressMessages(
		dbn(Y, model = "static", family = "gaussian",
			nscan = 2, burn = 1, verbose = FALSE)))
	expect_s3_class(fit, "dbn")
})

test_that("R = 1 HMM works (single regime)", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 5, seed = 1)
	fit = dbn(sim$Y, model = "hmm", R = 1, nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "hmm")
})

test_that("dynamic binary with small n produces a warning", {
	skip_on_cran()
	Y = array(rbinom(6 * 6 * 1 * 5, 1, 0.3), dim = c(6, 6, 1, 5))
	expect_warning(
		fit <- dbn(Y, model = "dynamic", family = "binary",
			nscan = 50, burn = 20, verbose = FALSE),
		regexp = "small networks|singularity"
	)
	expect_s3_class(fit, "dbn")
})

test_that("all-constant Y matrix runs without crash", {
	skip_on_cran()
	Y = array(0, dim = c(5, 5, 1, 5))
	for (t in 1:5) diag(Y[, , 1, t]) = NA
	fit = suppressWarnings(suppressMessages(
		dbn(Y, model = "static", family = "gaussian",
			nscan = 30, burn = 10, verbose = FALSE)))
	expect_s3_class(fit, "dbn")
})
