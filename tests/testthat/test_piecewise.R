#### piecewise model: essentials only

test_that("parse_blocks accepts integer K and breakpoint vector", {
	bi_int = parse_blocks(4, Tt = 40, n = 10)
	expect_equal(bi_int$K, 4)
	bi_vec = parse_blocks(c(10, 20, 30), Tt = 40, n = 10)
	expect_equal(bi_vec$K, 4)
})

test_that("parse_blocks rejects invalid input", {
	expect_error(parse_blocks("auto", Tt = 40, n = 10), regexp = "auto")
	expect_error(parse_blocks(0, Tt = 40, n = 10), regexp = "at least 1")
})

test_that("dbn piecewise fits with both integer-K and breakpoint-vector", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 6, time = 20, seed = 1)
	fit_int = dbn(sim$Y, model = "piecewise", blocks = 2,
		nscan = 30, burn = 10, verbose = FALSE)
	fit_vec = dbn(sim$Y, model = "piecewise", blocks = c(10),
		nscan = 30, burn = 10, verbose = FALSE)
	expect_s3_class(fit_int, "dbn")
	expect_s3_class(fit_vec, "dbn")
})

test_that("piecewise requires blocks parameter", {
	sim = simulate_dynamic_dbn(n = 5, time = 10, seed = 1)
	expect_error(
		dbn(sim$Y, model = "piecewise", nscan = 30, burn = 10, verbose = FALSE),
		regexp = "blocks"
	)
})

test_that("piecewise rejects short series", {
	sim = simulate_dynamic_dbn(n = 5, time = 3, seed = 1)
	expect_error(
		dbn(sim$Y, model = "piecewise", blocks = 2,
			nscan = 30, burn = 10, verbose = FALSE),
		regexp = "4 time points|short"
	)
})
