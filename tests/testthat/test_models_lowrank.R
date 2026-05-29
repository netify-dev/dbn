#### low-rank model: essentials only

test_that("simulate_lowrank_dbn produces valid output", {
	sim = simulate_lowrank_dbn(n = 6, p = 1, time = 10, r = 2, seed = 1)
	expect_equal(dim(sim$Y), c(6, 6, 1, 10))
	expect_equal(dim(sim$U)[2], 2)
	expect_equal(dim(sim$alpha), c(2, 10))
})

test_that("dbn_lowrank fits and has correct shape", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 6, p = 1, time = 10, r = 2, seed = 1)
	fit = dbn(sim$Y, model = "lowrank", r = 2, nscan = 40, burn = 15,
		verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "lowrank")
})

test_that("lowrank rejects r > n and r = 0", {
	sim = simulate_dynamic_dbn(n = 5, time = 8, seed = 1)
	expect_error(
		dbn(sim$Y, model = "lowrank", r = 10,
			nscan = 20, burn = 5, verbose = FALSE),
		regexp = "r"
	)
	expect_error(
		dbn(sim$Y, model = "lowrank", r = 0,
			nscan = 20, burn = 5, verbose = FALSE),
		regexp = "r"
	)
})
