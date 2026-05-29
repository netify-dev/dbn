#### simulator essentials

test_that("simulate_*_dbn generate valid arrays", {
	sim_s = simulate_static_dbn(n = 6, p = 1, time = 5, seed = 1)
	sim_d = simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 1)
	sim_h = simulate_hmm_dbn(n = 6, p = 1, time = 5, R = 2, seed = 1)
	sim_l = simulate_lowrank_dbn(n = 6, p = 1, time = 5, r = 2, seed = 1)
	for (sim in list(sim_s, sim_d, sim_h, sim_l)) {
		expect_equal(dim(sim$Y), c(6, 6, 1, 5))
	}
})

test_that("simulators respect seeds", {
	a = simulate_dynamic_dbn(n = 5, time = 4, seed = 42)
	b = simulate_dynamic_dbn(n = 5, time = 4, seed = 42)
	c = simulate_dynamic_dbn(n = 5, time = 4, seed = 99)
	expect_equal(a$Y, b$Y)
	expect_false(identical(a$Y, c$Y))
})

test_that("rank likelihood preserves ordering", {
	sim = simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 1)
	# Y is ordinal-coded; latent Z should preserve rank order within each slice
	y1 = as.numeric(sim$Y[, , 1, 1])
	z1 = as.numeric(sim$Z[, , 1, 1])
	keep = !is.na(y1) & !is.na(z1)
	if (sum(keep) > 1) {
		expect_gt(cor(y1[keep], z1[keep], method = "spearman"), 0.5)
	}
})
