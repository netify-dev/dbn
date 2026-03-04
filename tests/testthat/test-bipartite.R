####
# bipartite network tests
####

test_that("simulate_static_dbn produces bipartite data", {
	sim <- simulate_static_dbn(n = 6, n_col = 10, p = 1, time = 5, seed = 42)
	expect_equal(dim(sim$Y)[1], 6)
	expect_equal(dim(sim$Y)[2], 10)
	expect_equal(dim(sim$A), c(6, 6))  # sender dynamics

	expect_equal(dim(sim$B), c(10, 10))  # receiver dynamics
	expect_equal(dim(sim$M), c(6, 10, 1))
	expect_true(sim$is_bipartite)
	# bipartite should have no NA diagonal masking
	expect_false(any(is.na(diag(sim$Y[1:6, 1:6, 1, 1]))))
})

test_that("simulate_dynamic_dbn produces bipartite data", {
	sim <- simulate_dynamic_dbn(n = 5, n_col = 8, p = 1, time = 5, seed = 43)
	expect_equal(dim(sim$Y), c(5, 8, 1, 5))
	expect_equal(dim(sim$A), c(5, 5, 5))  # sender: n_row x n_row x T
	expect_equal(dim(sim$B), c(8, 8, 5))  # receiver: n_col x n_col x T
	expect_true(sim$is_bipartite)
})

test_that("simulate_lowrank_dbn produces bipartite data", {
	sim <- simulate_lowrank_dbn(n = 6, n_col = 9, p = 1, time = 5, r = 2, seed = 44)
	expect_equal(dim(sim$Y), c(6, 9, 1, 5))
	expect_equal(dim(sim$U), c(6, 2))  # sender loadings
	expect_true(sim$is_bipartite)
})

test_that("simulate_hmm_dbn produces bipartite data", {
	sim <- simulate_hmm_dbn(n = 6, n_col = 9, p = 1, time = 5, R = 2, seed = 45)
	expect_equal(dim(sim$Y), c(6, 9, 1, 5))
	expect_equal(dim(sim$A_list[[1]]), c(6, 6))
	expect_equal(dim(sim$B_list[[1]]), c(9, 9))
	expect_true(sim$is_bipartite)
})

test_that("simulate_test_data produces bipartite data", {
	Y <- simulate_test_data(n = 5, n_col = 8, p = 1, time = 3, seed = 46)
	expect_equal(dim(Y), c(5, 8, 1, 3))
})

test_that("static model fits bipartite data", {
	sim <- simulate_static_dbn(n = 6, n_col = 10, p = 1, time = 5, seed = 50)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$dims$n_row, 6)
	expect_equal(fit$dims$n_col, 10)
	expect_true(fit$dims$is_bipartite)
})

test_that("dynamic model fits bipartite data", {
	sim <- simulate_dynamic_dbn(n = 5, n_col = 8, p = 1, time = 5, seed = 51)
	fit <- dbn(sim$Y, model = "dynamic", family = "ordinal",
						 nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$dims$n_row, 5)
	expect_equal(fit$dims$n_col, 8)
	expect_true(fit$dims$is_bipartite)
	# A should be n_row x n_row
	expect_equal(dim(fit$A[[1]])[1], 5)
	expect_equal(dim(fit$A[[1]])[2], 5)
	# B should be n_col x n_col
	expect_equal(dim(fit$B[[1]])[1], 8)
	expect_equal(dim(fit$B[[1]])[2], 8)
})

test_that("hmm model rejects bipartite data with informative error", {
	sim <- simulate_hmm_dbn(n = 6, n_col = 9, p = 1, time = 5, R = 2, seed = 52)
	expect_error(
		dbn(sim$Y, model = "hmm", family = "ordinal", R = 2,
			nscan = 50, burn = 20, verbose = FALSE),
		"bipartite"
	)
})

test_that("lowrank model fits bipartite data", {
	skip("Lowrank bipartite requires n_row/n_col refactor in 3 C++ files (lowrank.cpp, lowrank_parallel.cpp, update_ab.cpp)")
	sim <- simulate_lowrank_dbn(n = 6, n_col = 9, p = 1, time = 5, r = 2, seed = 53)
	fit <- dbn(sim$Y, model = "lowrank", family = "ordinal", r = 2,
						 nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$dims$n_row, 6)
	expect_equal(fit$dims$n_col, 9)
	expect_true(fit$dims$is_bipartite)
})

test_that("unipartite backward compatibility preserved", {
	sim <- simulate_static_dbn(n = 8, p = 1, time = 5, seed = 60)
	expect_false(sim$is_bipartite)
	expect_equal(dim(sim$Y)[1], dim(sim$Y)[2])

	# self-loops should be NA
	expect_true(all(is.na(diag(sim$Y[,,1,1]))))

	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)
	expect_false(fit$dims$is_bipartite)
	expect_equal(fit$dims$n_row, fit$dims$n_col)
})
