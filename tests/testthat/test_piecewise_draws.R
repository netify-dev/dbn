# piecewise sampler stores 4d block-operator draws
test_that("piecewise stores A_blocks and B_blocks as 4d arrays", {
	skip_on_cran()
	sim = simulate_piecewise_dbn(n = 6, time = 12, blocks = c(6, 12), p = 1, seed = 42)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Y,
			model = "piecewise",
			family = "ordinal",
			blocks = c(6, 12),
			nscan = 100,
			burn = 50,
			odens = 2,
			verbose = FALSE
		)
	))
	expect_true(is.array(fit$draws$A_blocks))
	expect_true(is.array(fit$draws$B_blocks))
	expect_equal(length(dim(fit$draws$A_blocks)), 4)
	expect_equal(dim(fit$draws$A_blocks)[2], 2)
	expect_equal(dim(fit$draws$A_blocks)[3], 6)
	expect_equal(dim(fit$draws$A_blocks)[4], 6)
})

# piecewise as_draws exposes per-block operator variables
test_that("piecewise as_draws emits A_block[k,i,j] variables", {
	skip_on_cran()
	sim = simulate_piecewise_dbn(n = 5, time = 10, blocks = c(5, 10), p = 1, seed = 1)
	mc = suppressWarnings(suppressMessages(
		dbn_multichain(
			sim$Y,
			chains = 2,
			seeds = c(1, 2),
			model = "piecewise",
			family = "ordinal",
			blocks = c(5, 10),
			nscan = 80,
			burn = 40,
			odens = 2,
			verbose = FALSE
		)
	))
	da = posterior::as_draws(mc, include_operators = TRUE)
	vars = posterior::variables(da)
	ab_vars = grep("^A_block", vars, value = TRUE)
	# 2 blocks x 5 x 5 = 50 A_block entries
	expect_equal(length(ab_vars), 50)
})

# piecewise block 1 is retained (not dropped as an anchor)
test_that("piecewise block 1 is not dropped from as_draws", {
	skip_on_cran()
	sim = simulate_piecewise_dbn(n = 5, time = 10, blocks = c(5, 10), p = 1, seed = 1)
	mc = suppressWarnings(suppressMessages(
		dbn_multichain(
			sim$Y,
			chains = 2,
			seeds = c(1, 2),
			model = "piecewise",
			family = "ordinal",
			blocks = c(5, 10),
			nscan = 80,
			burn = 40,
			odens = 2,
			verbose = FALSE
		)
	))
	da = posterior::as_draws(mc, include_operators = TRUE)
	vars = posterior::variables(da)
	block1 = grep("^A_block\\[1,", vars, value = TRUE)
	expect_gt(length(block1), 0)
})

# operator_draws helper returns block scale by default for piecewise
test_that("operator_draws returns per-block draws for piecewise", {
	skip_on_cran()
	sim = simulate_piecewise_dbn(n = 5, time = 10, blocks = c(5, 10), p = 1, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Y,
			model = "piecewise",
			family = "ordinal",
			blocks = c(5, 10),
			nscan = 80,
			burn = 40,
			odens = 2,
			verbose = FALSE
		)
	))
	od = operator_draws(fit, scale = "block")
	expect_equal(length(dim(od$A)), 4)
	expect_equal(dim(od$A)[2], 2)
})

# operator_draws scale time expands blocks across T
test_that("operator_draws scale='time' expands blocks across T", {
	skip_on_cran()
	sim = simulate_piecewise_dbn(n = 5, time = 10, blocks = c(5, 10), p = 1, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Y,
			model = "piecewise",
			family = "ordinal",
			blocks = c(5, 10),
			nscan = 80,
			burn = 40,
			odens = 2,
			verbose = FALSE
		)
	))
	od = operator_draws(fit, scale = "time")
	expect_equal(length(dim(od$A)), 4)
	expect_equal(dim(od$A)[4], 10)
	# operators at t=1 and t=5 should be identical (both block 1)
	expect_equal(od$A[, , , 1], od$A[, , , 5])
	# operators at t=5 and t=6 should differ (different blocks)
	expect_false(all(od$A[, , , 5] == od$A[, , , 6]))
})

# operator_draws on dynamic uses time scale only
test_that("operator_draws scale='block' rejects dynamic fit", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 8, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Y,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	expect_error(
		operator_draws(fit, scale = "block"),
		regexp = "piecewise"
	)
	od = operator_draws(fit, scale = "time")
	expect_equal(length(dim(od$A)), 4)
})
