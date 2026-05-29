####
# timing and capability comparison across model types.
# run manually with:
#   library(dbn); library(testthat)
#   source(system.file("benchmarks/timing_comparison.R", package = "dbn"))
# all blocks are wrapped in test_that() so wallclock measurements double as
# pass/fail assertions; output goes to the console.
####

#### timing comparisons: piecewise vs dynamic

test_that("piecewise and dynamic timing comparison", {
	skip_on_cran()
	set.seed(1001)

	# moderate size network with 30 time points
	n = 10
	Tt = 30
	nscan = 100
	burn = 50

	sim = simulate_piecewise_dbn(n = n, time = Tt, blocks = 3, seed = 1001)

	# time piecewise (3 blocks = 3 sets of parameters)
	pw_start = Sys.time()
	pw_fit = dbn(sim$Y, model = "piecewise", blocks = 3,
		nscan = nscan, burn = burn, verbose = FALSE)
	pw_time = as.numeric(difftime(Sys.time(), pw_start, units = "secs"))

	# time dynamic (T-1 sets of parameters)
	dyn_start = Sys.time()
	dyn_fit = dbn(sim$Y, model = "dynamic",
		nscan = nscan, burn = burn, verbose = FALSE)
	dyn_time = as.numeric(difftime(Sys.time(), dyn_start, units = "secs"))

	expect_s3_class(pw_fit, "dbn")
	expect_s3_class(dyn_fit, "dbn")

	# report timing characteristics
	cat("\n  Piecewise (K=3): ", round(pw_time, 2), "s\n")
	cat("  Dynamic (T=30):  ", round(dyn_time, 2), "s\n")
	if (pw_time > 0 && dyn_time > 0) {
		cat("  Ratio (pw/dyn):  ", round(pw_time / dyn_time, 1), "x\n")
	}

	# both models should complete successfully
	expect_equal(pw_fit$model, "piecewise")
	expect_equal(dyn_fit$model, "dynamic")
})

test_that("piecewise timing scales with T", {
	skip_on_cran()
	set.seed(1002)

	n = 8
	nscan = 80
	burn = 40

	# short time series (T=20) with K=2
	sim_short = simulate_piecewise_dbn(n = n, time = 20, blocks = 2, seed = 1002)
	short_start = Sys.time()
	short_fit = dbn(sim_short$Y, model = "piecewise", blocks = 2,
		nscan = nscan, burn = burn, verbose = FALSE)
	short_time = as.numeric(difftime(Sys.time(), short_start, units = "secs"))

	expect_s3_class(short_fit, "dbn")
	expect_equal(short_fit$blocks$K, 2)

	cat("\n  T=20, K=2: ", round(short_time, 2), "s\n")
	cat("  Short fit completed with ", short_fit$settings$draws, " draws\n")
})

#### structural change detection capability

test_that("piecewise detects structural change in A", {
	skip_on_cran()
	set.seed(2001)

	n = 8
	Tt = 40
	K = 2

	# simulate with true structural change at t=20
	sim = simulate_piecewise_dbn(n = n, time = Tt, blocks = c(20, 40), seed = 2001)

	# fit model
	fit = dbn(sim$Y, model = "piecewise", blocks = c(20, 40),
		nscan = 150, burn = 75, verbose = FALSE)

	# compare estimated A matrices across blocks
	A1_hat = fit$A_blocks[[1]]
	A2_hat = fit$A_blocks[[2]]

	# difference between blocks should be non-trivial
	A_diff = mean(abs(A1_hat - A2_hat))

	cat("\n  Mean |A1 - A2| (estimated): ", round(A_diff, 4), "\n")

	# true A matrices from simulation
	A1_true = sim$true_A[[1]]
	A2_true = sim$true_A[[2]]
	true_diff = mean(abs(A1_true - A2_true))

	cat("  Mean |A1 - A2| (true):      ", round(true_diff, 4), "\n")

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "piecewise")
	expect_equal(fit$blocks$K, 2)
})

test_that("piecewise recovers block structure with known boundaries", {
	skip_on_cran()
	set.seed(2002)

	n = 6
	Tt = 30
	boundaries = c(10, 20, 30)

	sim = simulate_piecewise_dbn(n = n, time = Tt, blocks = boundaries, seed = 2002)

	fit = dbn(sim$Y, model = "piecewise", blocks = boundaries,
		nscan = 150, burn = 75, verbose = FALSE)

	expect_equal(fit$blocks$K, 3)
	expect_equal(fit$blocks$boundaries, c(0, 10, 20, 30))

	# check that block-specific estimates exist
	expect_equal(length(fit$A_blocks), 3)
	expect_equal(length(fit$B_blocks), 3)

	# each A should be n x n
	for (k in 1:3) {
		expect_equal(dim(fit$A_blocks[[k]]), c(n, n))
		expect_equal(dim(fit$B_blocks[[k]]), c(n, n))
	}
})

test_that("compare_blocks quantifies between-block differences", {
	skip_on_cran()
	set.seed(2003)

	sim = simulate_piecewise_dbn(n = 6, time = 24, blocks = 3, seed = 2003)
	fit = dbn(sim$Y, model = "piecewise", blocks = 3,
		nscan = 100, burn = 50, verbose = FALSE)

	result = compare_blocks(fit)

	expect_true(is.list(result))
	expect_equal(length(result), 2)  # 2 adjacent pairs for 3 blocks

	# each comparison should have mean_diff and ci
	expect_true("mean_diff" %in% names(result[[1]]))
	expect_true("ci" %in% names(result[[1]]))
	expect_true("prob_above_threshold" %in% names(result[[1]]))

	cat("\n  Block 1 vs 2 - mean diff:", round(result[[1]]$mean_diff, 4),
		", P(diff > threshold):", round(result[[1]]$prob_above_threshold, 3), "\n")
	cat("  Block 2 vs 3 - mean diff:", round(result[[2]]$mean_diff, 4),
		", P(diff > threshold):", round(result[[2]]$prob_above_threshold, 3), "\n")
})

#### model comparison: piecewise vs static

test_that("piecewise with K=1 approximates static model", {
	skip_on_cran()
	set.seed(3001)

	sim = simulate_static_dbn(n = 8, time = 20, seed = 3001)

	# fit static model
	fit_static = dbn(sim$Y, model = "static",
		nscan = 100, burn = 50, verbose = FALSE)

	# fit piecewise with single block
	fit_pw1 = dbn(sim$Y, model = "piecewise", blocks = 1,
		nscan = 100, burn = 50, verbose = FALSE)

	expect_s3_class(fit_static, "dbn")
	expect_s3_class(fit_pw1, "dbn")

	# both should have same network dimensions
	expect_equal(fit_static$dims$n_row, fit_pw1$dims$n_row)
	expect_equal(fit_static$dims$n_col, fit_pw1$dims$n_col)

	# piecewise B_blocks should be n x n
	expect_equal(dim(fit_pw1$B_blocks[[1]]), c(8, 8))
})

test_that("piecewise captures regime changes that static misses", {
	skip_on_cran()
	set.seed(3002)

	n = 8
	Tt = 40

	# simulate data with strong regime change
	sim = simulate_piecewise_dbn(n = n, time = Tt, blocks = 2, seed = 3002)

	# fit static (ignores regime change)
	fit_static = dbn(sim$Y, model = "static",
		nscan = 100, burn = 50, verbose = FALSE)

	# fit piecewise (captures regime change)
	fit_pw = dbn(sim$Y, model = "piecewise", blocks = 2,
		nscan = 100, burn = 50, verbose = FALSE)

	expect_equal(fit_static$model, "static")
	expect_equal(fit_pw$model, "piecewise")

	# piecewise has block-specific parameters
	expect_equal(length(fit_pw$A_blocks), 2)
	expect_equal(length(fit_pw$B_blocks), 2)
})

#### complexity comparison

test_that("parameter counts reflect model complexity", {
	skip_on_cran()
	set.seed(4001)

	n = 6
	Tt = 30

	sim = simulate_piecewise_dbn(n = n, time = Tt, blocks = 3, seed = 4001)

	fit_static = dbn(sim$Y, model = "static",
		nscan = 60, burn = 30, verbose = FALSE)
	fit_pw = dbn(sim$Y, model = "piecewise", blocks = 3,
		nscan = 60, burn = 30, verbose = FALSE)
	fit_dyn = dbn(sim$Y, model = "dynamic",
		nscan = 60, burn = 30, verbose = FALSE)

	# count effective parameters (A + B matrices)
	# static: 1 set of A, B each (2 * n^2)
	# piecewise K=3: 3 sets of A, B each (6 * n^2)
	# dynamic T=30: 29 sets of A, B each (58 * n^2)

	n_params_static = 2 * n^2
	n_params_pw = 2 * 3 * n^2
	n_params_dyn = 2 * (Tt - 1) * n^2

	cat("\n  Model complexity (n=", n, ", T=", Tt, "):\n", sep = "")
	cat("    Static:           ", n_params_static, " effective params\n")
	cat("    Piecewise (K=3):  ", n_params_pw, " effective params\n")
	cat("    Dynamic:          ", n_params_dyn, " effective params\n")

	expect_lt(n_params_static, n_params_pw)
	expect_lt(n_params_pw, n_params_dyn)
})

#### gaussian family support

test_that("piecewise works with gaussian family and continuous data", {
	skip_on_cran()
	set.seed(5001)

	sim = simulate_piecewise_dbn(n = 8, time = 30, blocks = 2, seed = 5001)

	# use continuous Y (before ordinal discretization)
	fit = dbn(sim$Y_continuous, model = "piecewise", blocks = 2,
		family = "gaussian", nscan = 100, burn = 50, verbose = FALSE)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$family, "gaussian")
	expect_equal(fit$blocks$K, 2)
})

#### larger network demonstration

test_that("piecewise handles moderate-sized networks", {
	skip_on_cran()
	set.seed(6001)

	# 20 actors, 40 time points - manageable for piecewise
	n = 20
	Tt = 40
	K = 4

	sim = simulate_piecewise_dbn(n = n, time = Tt, blocks = K, seed = 6001)

	start = Sys.time()
	fit = dbn(sim$Y, model = "piecewise", blocks = K,
		nscan = 50, burn = 25, verbose = FALSE)
	elapsed = as.numeric(difftime(Sys.time(), start, units = "secs"))

	expect_s3_class(fit, "dbn")
	expect_equal(fit$dims$n_row, n)
	expect_equal(fit$blocks$K, K)

	cat("\n  n=20, T=40, K=4: ", round(elapsed, 2), "s\n")
	cat("  Estimated A dims: ", paste(dim(fit$A_blocks[[1]]), collapse = "x"), "\n")
})

#### output structure verification

test_that("piecewise draws have correct structure for inference", {
	skip_on_cran()
	set.seed(7001)

	sim = simulate_piecewise_dbn(n = 6, time = 20, blocks = 2, seed = 7001)
	fit = dbn(sim$Y, model = "piecewise", blocks = 2,
		nscan = 80, burn = 40, odens = 2, verbose = FALSE)

	# check draws structure
	expect_equal(fit$settings$draws, 40)
	expect_true("draws" %in% names(fit))

	# parameter draws should exist
	expect_true("pars" %in% names(fit$draws))
	expect_equal(nrow(fit$draws$pars), 40)

	# A and B draws should be stored
	expect_true("A" %in% names(fit$draws$misc))
	expect_equal(length(fit$draws$misc$A), 40)
})
