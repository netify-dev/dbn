#### bipartite network: essentials only

test_that("simulate_dynamic_dbn produces bipartite data", {
	sim = simulate_dynamic_dbn(n = 8, n_col = 6, p = 1, time = 5, seed = 1)
	expect_equal(dim(sim$Y), c(8, 6, 1, 5))
	expect_true(sim$is_bipartite)
})

test_that("dbn fits bipartite data with dynamic and static models", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 8, n_col = 6, p = 1, time = 8, seed = 1)
	fit_s = dbn(sim$Y, model = "static", family = "gaussian",
		nscan = 40, burn = 15, verbose = FALSE)
	fit_d = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 40, burn = 15, verbose = FALSE)
	expect_s3_class(fit_s, "dbn")
	expect_s3_class(fit_d, "dbn")
	expect_true(fit_d$dims$is_bipartite)
})

test_that("hmm and lowrank reject bipartite data", {
	Y_bip = array(rnorm(8 * 6 * 1 * 5), dim = c(8, 6, 1, 5))
	expect_error(
		dbn(Y_bip, model = "hmm", R = 2,
			nscan = 20, burn = 5, verbose = FALSE),
		regexp = "bipartite|square|unipartite"
	)
	expect_error(
		dbn(Y_bip, model = "lowrank", r = 2,
			nscan = 20, burn = 5, verbose = FALSE),
		regexp = "bipartite|square|unipartite"
	)
})
