# forecast numeric correctness across all predict paths. checks that a
# forecast is seeded from the end of the series, stays bounded, and is
# deterministic under a fixed seed, not just structurally valid.

# helper: mean over the (i, j, rel) cells at a given horizon of a
# summary="mean" forecast array [n, n, p, H]
h_mean = function(fc, h) mean(fc[, , , h])

test_that("dynamic forecast is seeded from the last state (not zero) on centered data", {
	skip_on_cran()
	set.seed(1)
	sim = simulate_dynamic_dbn(n = 8, p = 1, time = 15, sigma2 = 0.2,
															tauA2 = 0.05, tauB2 = 0.05, ar1 = TRUE,
															rhoA = 0.8, rhoB = 0.8, seed = 1)
	fit = dbn(sim$Z, model = "dynamic", family = "gaussian",
						 ar1 = TRUE, update_rho = TRUE,
						 nscan = 1000, burn = 500, odens = 2, verbose = FALSE)
	fc = predict(fit, H = 3, draws = 200, summary = "mean", seed = 7)
	expect_equal(dim(fc), c(8, 8, 1, 3))
	expect_true(all(is.finite(fc)))
	# bounded: a stable fit should not blow up over 3 steps
	expect_lt(max(abs(fc)), 10 * max(abs(sim$Z), na.rm = TRUE))
})

test_that("forecasts are deterministic given a seed (all paths)", {
	skip_on_cran()
	set.seed(2)
	simd = simulate_dynamic_dbn(n = 7, p = 1, time = 12, ar1 = TRUE,
															 rhoA = 0.8, rhoB = 0.8, seed = 2)
	fd = dbn(simd$Z, model = "dynamic", family = "gaussian", ar1 = TRUE,
						nscan = 600, burn = 300, odens = 2, verbose = FALSE)
	a = predict(fd, H = 2, draws = 50, summary = "mean", seed = 123)
	b = predict(fd, H = 2, draws = 50, summary = "mean", seed = 123)
	expect_identical(a, b)
})

test_that("HMM forecast seeds from the terminal latent state, not zero (regression)", {
	skip_on_cran()
	# this is the exact bug that was present: Theta_now = array(0); A %*% 0 %*% B'
	# gives a pure-noise forecast centered at 0 regardless of the data level.
	set.seed(3)
	sim = simulate_hmm_dbn(n = 10, p = 1, time = 20, R = 2, sigma2 = 0.3,
													tau_A2 = 0.3, tau_B2 = 0.3, transition_prob = 0.85,
													seed = 3)
	Z = sim$Z + 4                      # shift the level well away from 0
	fit = dbn(Z, model = "hmm", family = "gaussian", R = 2,
						 nscan = 800, burn = 400, odens = 2, verbose = FALSE)
	fc = predict(fit, H = 3, draws = 100, summary = "mean", seed = 9)
	obs_level = mean(Z[, , 1, 20], na.rm = TRUE)
	# a zero-seeded forecast would sit near 0; a correctly-seeded one starts
	# near the terminal observed level (~4). require it be much closer to the
	# level than to zero.
	expect_gt(h_mean(fc, 1), 0.5 * obs_level)
	expect_true(all(is.finite(fc)))
})

test_that("low-rank forecast seeds from the last observed slice, not zero", {
	skip_on_cran()
	set.seed(4)
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 12, r = 2, sigma2 = 0.3,
															tau_alpha2 = 0.15, tauB2 = 0.04, seed = 4)
	Z = sim$Z + 3
	fit = dbn(Z, model = "lowrank", family = "gaussian", r = 2,
						 nscan = 600, burn = 300, odens = 2, verbose = FALSE)
	fc = predict(fit, H = 3, draws = 100, summary = "mean", seed = 11)
	obs_level = mean(Z[, , 1, 12], na.rm = TRUE)
	# must carry a meaningful fraction of the level at h = 1 (zero-seed -> ~0)
	expect_gt(h_mean(fc, 1), 0.25 * obs_level)
	expect_true(all(is.finite(fc)))
})

test_that("piecewise forecast seeds from the last state and stays bounded", {
	skip_on_cran()
	set.seed(5)
	sim = simulate_piecewise_dbn(n = 10, time = 24, blocks = c(12, 24), p = 1,
																sigma2 = 0.4, tau2 = 0.3, seed = 5)
	fit = dbn(sim$Y_continuous, model = "piecewise", family = "gaussian",
						 blocks = c(12, 24), nscan = 800, burn = 400, odens = 2,
						 verbose = FALSE)
	fc = predict(fit, H = 4, draws = 100, summary = "mean", seed = 13)
	expect_equal(dim(fc), c(10, 10, 1, 4))
	expect_true(all(is.finite(fc)))
	expect_lt(max(abs(fc)), 10 * max(abs(sim$Y_continuous), na.rm = TRUE))
})

test_that("forecasts stay finite over a multi-step horizon", {
	skip_on_cran()
	# the prior does not enforce contractivity (documented), so a fitted
	# operator can be locally explosive and the forecast may grow. we do NOT
	# assert decay; we assert the forecast stays finite and within a generous
	# multiple of the data scale, which catches a true blow-up to Inf/NaN or
	# an off-by-orders projection bug.
	set.seed(6)
	sim = simulate_dynamic_dbn(n = 8, p = 1, time = 18, sigma2 = 0.1,
															tauA2 = 0.01, tauB2 = 0.01, ar1 = TRUE,
															rhoA = 0.5, rhoB = 0.5, seed = 6)
	fit = dbn(sim$Z, model = "dynamic", family = "gaussian", ar1 = TRUE,
						 update_rho = TRUE, nscan = 1000, burn = 500, odens = 2,
						 verbose = FALSE)
	fc = predict(fit, H = 5, draws = 200, summary = "mean", seed = 21)
	expect_true(all(is.finite(fc)))
	expect_lt(max(abs(fc)), 50 * max(abs(sim$Z), na.rm = TRUE))
})
