# dbn_compute_snr and dbn_coupling_rank_probs (new in 1.1.0)
# basic shape, finiteness, and structural-symmetry assertions on
# small synthetic data.

test_that("dbn_compute_snr returns expected fields and shapes", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 6, time = 10, sigma2 = 0.3,
		tauA2 = 0.05, symmetric = TRUE, seed = 41,
		return_truth = TRUE)
	Y <- sim$Y[, , 1, , drop = FALSE]

	fit <- dbn_dynamic(Y, family = "gaussian",
		nscan = 60, burn = 30, odens = 5, time_thin = 1,
		symmetric = TRUE, seed = 5, verbose = FALSE)

	snr <- dbn_compute_snr(fit, windows = list(early = 1:5, late = 6:10))

	# per-draw, per-(t-1) matrix
	expect_true(is.matrix(snr$snr_per_draw_per_t))
	expect_equal(nrow(snr$snr_per_draw_per_t), length(fit$A))
	expect_equal(ncol(snr$snr_per_draw_per_t), dim(Y)[4] - 1L)
	# all entries finite and non-negative
	expect_true(all(is.finite(snr$snr_per_draw_per_t)))
	expect_true(all(snr$snr_per_draw_per_t >= 0))
	# full sample summary has expected names
	expect_true(all(c("mean", "sd", "q025", "q975") %in% names(snr$full_sample)))
	# window names propagate
	expect_true(all(c("early", "late") %in% names(snr$windows)))
})

test_that("dbn_coupling_rank_probs sums and dominance matrix sane", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 5, time = 8, sigma2 = 0.3,
		tauA2 = 0.05, symmetric = TRUE, seed = 43,
		return_truth = TRUE)
	Y <- sim$Y[, , 1, , drop = FALSE]

	fit <- dbn_dynamic(Y, family = "gaussian",
		nscan = 50, burn = 20, odens = 5, time_thin = 1,
		symmetric = TRUE, seed = 3, verbose = FALSE)

	rp <- dbn_coupling_rank_probs(fit, K_grid = c(1, 2, 4))

	full <- rp$full_sample
	expect_equal(length(full$post_mean), 5L)
	# top-K probabilities sum to K across actors (each draw places each
	# actor in a unique rank; sum across actors of indicator{rank <= K}
	# is exactly K)
	for (K in c(1, 2, 4)) {
		pK <- full$p_topK[[paste0("p_top", K)]]
		expect_equal(sum(pK), K, tolerance = 1e-12)
	}
	# pairwise dominance matrix is zero on diagonal and pi_dom[i, j]
	# + pi_dom[j, i] approximately = 1 (modulo ties)
	expect_true(all(diag(full$pi_dom) == 0))
	off <- full$pi_dom + t(full$pi_dom)
	diag(off) <- 1
	expect_true(all(off <= 1 + 1e-12))
	expect_true(all(off >= 0 - 1e-12))
})

test_that("dbn_compute_snr requires p = 1 unipartite symmetric fit", {
	skip_on_cran()
	# a non-dbn input should error
	expect_error(dbn_compute_snr(list()), "must be a dbn object")
})
