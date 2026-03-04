####
# model fitting with simulated data
####

test_that("static model runs with simulated static data", {
	data <- simulate_static_dbn(n = 10, p = 2, time = 20, seed = 68861)

	results <- dbn(data$Y, model = "static", nscan = 30, burn = 10, odens = 1, verbose = FALSE)

	expect_s3_class(results, "dbn")
	expect_equal(results$model, "static")
	expect_true(is.list(results$B))
})

test_that("dynamic model runs with simulated dynamic data", {
	data <- simulate_dynamic_dbn(n = 8, p = 1, time = 15,
															 ar1 = FALSE, seed = 68862)

	results <- dbn(data$Y, model = "dynamic", nscan = 30, burn = 10, odens = 1, verbose = FALSE)

	expect_equal(results$model, "dynamic")
	expect_true("A" %in% names(results))
	expect_true("B" %in% names(results))
})

test_that("low-rank model runs with simulated low-rank data", {
	data <- simulate_lowrank_dbn(n = 10, p = 1, time = 15, r = 2, seed = 68863)

	results <- dbn(data$Y, model = "lowrank", r = 2, nscan = 20, burn = 10, odens = 1, verbose = FALSE)

	expect_equal(results$model, "lowrank")
	expect_true("U" %in% names(results))
	expect_true("alpha" %in% names(results))
})

test_that("HMM model runs with simulated HMM data", {
	data <- simulate_hmm_dbn(n = 10, p = 1, time = 15, R = 2,
													 transition_prob = 0.9, seed = 68864)

	results <- dbn(data$Y, model = "hmm", R = 2, nscan = 20, burn = 10, odens = 1, verbose = FALSE)

	expect_equal(results$model, "hmm")
	expect_true("S" %in% names(results))
	expect_true("Pi" %in% names(results))
})

test_that("simulate functions produce data with expected properties", {
	# Static data should have temporal autocorrelation
	data_static <- simulate_static_dbn(n = 10, p = 1, time = 10, seed = 68865)
	cor_t1_t2 <- cor(c(data_static$Z[,,1,1]), c(data_static$Z[,,1,2]))
	expect_true(cor_t1_t2 > 0.3)

	# HMM data should show regime persistence
	data_hmm <- simulate_hmm_dbn(n = 10, p = 1, time = 10, R = 3,
															 transition_prob = 0.8, seed = 68868)
	n_switches <- sum(diff(data_hmm$S) != 0)
	expect_true(n_switches < 8)
})
