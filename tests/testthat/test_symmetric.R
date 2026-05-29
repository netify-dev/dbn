#### symmetric (B = A) constraint: essentials only

test_that("simulate_*_dbn with symmetric=TRUE produces matching A and B", {
	sim_s = simulate_static_dbn(n = 6, p = 1, time = 5, symmetric = TRUE, seed = 1)
	expect_true(all(abs(sim_s$A - sim_s$B) < 1e-10))
	sim_d = simulate_dynamic_dbn(n = 6, p = 1, time = 5, symmetric = TRUE, seed = 1)
	expect_true(all(abs(sim_d$A - sim_d$B) < 1e-10))
})

test_that("simulate functions error when symmetric + bipartite", {
	expect_error(
		simulate_dynamic_dbn(n = 6, n_col = 8, time = 5, symmetric = TRUE),
		regexp = "n_row|n_col|symmetric|square|bipartite"
	)
})

test_that("dbn dynamic + symmetric enforces A == B in output", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 6, p = 1, time = 10, symmetric = TRUE, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian", symmetric = TRUE,
		nscan = 40, burn = 15, verbose = FALSE)
	for (s in seq_along(fit$A)) {
		expect_equal(fit$A[[s]], fit$B[[s]], tolerance = 1e-10)
	}
})

test_that("dbn rejects incompatible symmetric combinations", {
	# bipartite + symmetric: error
	Y_bip = array(rnorm(6 * 8 * 1 * 5), dim = c(6, 8, 1, 5))
	expect_error(
		dbn(Y_bip, model = "dynamic", symmetric = TRUE,
			nscan = 20, burn = 5, verbose = FALSE),
		regexp = "symmetric|bipartite|square"
	)
	# lowrank + symmetric: error
	sim = simulate_dynamic_dbn(n = 6, p = 1, time = 8, seed = 1)
	expect_error(
		dbn(sim$Y, model = "lowrank", r = 2, symmetric = TRUE,
			nscan = 20, burn = 5, verbose = FALSE),
		regexp = "low-?rank|symmetric"
	)
})
