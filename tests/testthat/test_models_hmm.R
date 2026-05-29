#### HMM model: essentials only

test_that("simulate_hmm_dbn produces valid output", {
	sim = simulate_hmm_dbn(n = 6, p = 1, time = 12, R = 2, seed = 1)
	expect_equal(dim(sim$Y), c(6, 6, 1, 12))
	expect_length(sim$A_list, 2L)
	expect_length(sim$B_list, 2L)
})

test_that("dbn_hmm fits and exposes valid regime probabilities", {
	skip_on_cran()
	sim = simulate_hmm_dbn(n = 6, p = 1, time = 12, R = 2, seed = 1)
	fit = dbn(sim$Y, model = "hmm", R = 2, nscan = 50, burn = 20,
		verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "hmm")
	rp = regime_probs(fit)
	expect_true(all(rp >= 0 & rp <= 1, na.rm = TRUE))
})

test_that("HMM rejects T < 2", {
	Y = array(rnorm(36), dim = c(6, 6, 1, 1))
	diag(Y[, , 1, 1]) = NA
	expect_error(
		dbn(Y, model = "hmm", R = 2, nscan = 20, burn = 5, verbose = FALSE),
		regexp = "at least 2 time points|T"
	)
})
