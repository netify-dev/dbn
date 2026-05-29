# label-switching correction for HMM fits
test_that("hmm_relabel returns a relabeled Pi summary", {
	skip_on_cran()
	sim = simulate_hmm_dbn(n = 5, time = 12, R = 2, p = 1, seed = 1)
	fit = suppressMessages(suppressWarnings(
		dbn(sim$Y, model = "hmm", family = "gaussian", R = 2L,
			nscan = 80, burn = 30, odens = 2, verbose = FALSE)
	))
	res = hmm_relabel(fit, ref = "first", verbose = FALSE)
	expect_s3_class(res, "dbn_hmm_relabel")
	expect_equal(dim(res$Pi_mean), c(2L, 2L))
	expect_equal(dim(res$Pi_lo), c(2L, 2L))
	expect_equal(dim(res$Pi_hi), c(2L, 2L))
	# rows of Pi should sum to ~1 after relabeling (since each individual
	# draw's Pi sums to 1 by construction, and the mean preserves this)
	expect_true(all(abs(rowSums(res$Pi_mean) - 1) < 1e-10))
	expect_true(all(res$Pi_lo <= res$Pi_hi))
	# perms applied are valid permutations of 1..R
	expect_equal(dim(res$perms), c(length(fit$Pi), 2L))
	expect_true(all(apply(res$perms, 1, function(p) all(sort(p) == seq_len(2)))))
})

test_that("hmm_relabel refuses R > 5 with directive message", {
	# build a stub fit with R = 6 to test the guard
	stub <- list(model = "hmm", R = 6L, Pi = list(diag(6) / 1))
	class(stub) <- c("dbn", "list")
	expect_error(hmm_relabel(stub), "R <= 5")
})

test_that(".all_permutations enumerates all R! permutations", {
	p2 = dbn:::.all_permutations(2L)
	expect_equal(nrow(p2), 2L)
	expect_equal(ncol(p2), 2L)
	p3 = dbn:::.all_permutations(3L)
	expect_equal(nrow(p3), 6L)
	expect_true(all(sapply(seq_len(nrow(p3)),
		function(k) all(sort(p3[k, ]) == 1:3))))
})

test_that("hmm_relabel with mode reference picks high-trace draw", {
	skip_on_cran()
	sim = simulate_hmm_dbn(n = 4, time = 10, R = 2, p = 1, seed = 5)
	fit = suppressMessages(suppressWarnings(
		dbn(sim$Y, model = "hmm", family = "gaussian", R = 2L,
			nscan = 60, burn = 20, odens = 2, verbose = FALSE)
	))
	res_mode = hmm_relabel(fit, ref = "mode", verbose = FALSE)
	# mode reference should pick a draw with trace(Pi) >= median trace
	traces = sapply(fit$Pi, function(P) sum(diag(P)))
	expect_gte(traces[res_mode$ref_idx], stats::median(traces))
})
