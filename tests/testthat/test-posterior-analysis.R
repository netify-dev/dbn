####
# posterior extraction and summary functions
####

test_that("theta_slice extracts correct values from HMM model", {
	set.seed(101)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 101)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	# check theta draws exist
	expect_true(!is.null(fit$draws$theta))
	n_draws <- length(fit$draws$theta)

	# extract single dyad at single time
	slice1 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = 1)
	expect_equal(length(slice1), n_draws)

	# extract single dyad at multiple times
	slice2 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = c(1, 3, 5))
	expect_equal(nrow(slice2), n_draws)
	expect_equal(ncol(slice2), 3)

	# subset draws
	slice3 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = 1, draws = 1:5)
	expect_equal(length(slice3), 5)
})

test_that("theta_summary computes summaries from HMM model", {
	set.seed(102)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 102)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	# mean summary for a specific cell
	s1 <- theta_summary(fit, fun = mean, i = 1, j = 2, rel = 1, time = 1)
	expect_true(!is.null(s1))
	expect_true("value" %in% names(s1))
	expect_equal(nrow(s1), 1)

	# mean summary for a time slice (all dyads at time 1)
	s2 <- theta_summary(fit, fun = mean, time = 1)
	expect_true(nrow(s2) > 0)
})

test_that("theta_slice returns NULL for static model (no theta draws)", {
	set.seed(103)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 103)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)

	# static model does not store theta draws
	expect_warning(theta_slice(fit, i = 1, j = 2, rel = 1, time = 1),
								 "does not contain theta draws")
})

test_that("param_summary computes parameter quantiles for static model", {
	set.seed(104)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 104)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)

	ps <- param_summary(fit)
	expect_true(!is.null(ps))
	expect_true("parameter" %in% names(ps))
	expect_true("mean" %in% names(ps))
	expect_true("sd" %in% names(ps))
	expect_true(nrow(ps) >= 3) # at least s2, t2, g2
})

test_that("param_summary works for dynamic model", {
	set.seed(105)
	sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 8, seed = 105)
	fit <- dbn(sim$Y, model = "dynamic", family = "ordinal",
						 nscan = 30, burn = 10, odens = 1, verbose = FALSE)

	ps <- param_summary(fit)
	expect_true(!is.null(ps))
	expect_true(nrow(ps) >= 2) # at least sigma2, tau_A2
})

test_that("latent_summary extracts M summaries", {
	set.seed(106)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 106)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)

	ls <- latent_summary(fit, fun = mean)
	expect_true(!is.null(ls))
	expect_true(nrow(ls) > 0)
	expect_true("value" %in% names(ls))
})

test_that("latent_summary with relation filter", {
	set.seed(107)
	sim <- simulate_static_dbn(n = 6, p = 2, time = 5, seed = 107)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)

	ls_all <- latent_summary(fit, fun = mean)
	ls_r1 <- latent_summary(fit, fun = mean, rel = 1)

	expect_true(nrow(ls_all) > nrow(ls_r1))
	if ("rel" %in% names(ls_r1)) {
		expect_true(all(ls_r1$rel == 1))
	}
})

test_that("regime_probs works for HMM models", {
	set.seed(108)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 108)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	probs <- regime_probs(fit)
	expect_true(!is.null(probs))
	expect_equal(ncol(probs), 2)
	expect_equal(nrow(probs), 8)
	expect_true(all(abs(rowSums(probs) - 1) < 1e-10))
})

test_that("regime_probs returns NULL for non-HMM", {
	set.seed(109)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 109)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)

	expect_null(regime_probs(fit))
})

test_that("derive_draws computes derived quantities", {
	set.seed(110)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 110)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)

	# derive mean of M from each draw
	m_means <- derive_draws(fit, function(draw) {
		mean(draw$M[,,1], na.rm = TRUE)
	}, name = "M_mean")

	expect_true(!is.null(m_means))
	expect_true(is.numeric(m_means))
})

test_that("plot_theta produces ggplot object for HMM model", {
	set.seed(111)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 111)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	p <- plot_theta(fit, time = 1, rel = 1)
	expect_s3_class(p, "ggplot")
})

test_that("plot_trace produces ggplot object", {
	set.seed(112)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 112)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)

	p <- plot_trace(fit)
	expect_s3_class(p, "ggplot")
})

test_that("theta_credible returns correct structure", {
	set.seed(113)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 113)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	tc <- theta_credible(fit, i = 1:3, j = 1:3, time = 1:4)
	expect_true(!is.null(tc))
	expect_equal(nrow(tc), 3 * 3 * 4) # 3 senders x 3 receivers x 4 times
	expect_true(all(c("mean", "lower", "median", "upper") %in% names(tc)))
	expect_true(all(tc$lower <= tc$median))
	expect_true(all(tc$median <= tc$upper))
})

test_that("theta_credible returns NULL for static model", {
	set.seed(114)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 114)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)

	expect_warning(theta_credible(fit), "does not contain theta draws")
})

test_that("network_summary computes all stat types", {
	set.seed(115)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 115)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	for (stat in c("mean", "density", "strength")) {
		ns <- network_summary(fit, stat = stat)
		expect_true(!is.null(ns))
		expect_equal(nrow(ns), 8) # 8 time points
		expect_true(all(c("time", "mean", "lower", "upper") %in% names(ns)))
		expect_true(all(ns$lower <= ns$upper))
	}
})

test_that("edge_prob returns correct matrix", {
	set.seed(116)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 116)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	ep <- edge_prob(fit, rel = 1, time = 8)
	expect_equal(dim(ep), c(6, 6))
	expect_true(all(ep >= 0 & ep <= 1, na.rm = TRUE))
	expect_true(all(is.na(diag(ep)))) # diagonal should be NA for unipartite
})

test_that("estimate_memory returns numeric value", {
	mem <- estimate_memory(n_row = 20, n_col = 20, p = 1, Tt = 30,
												 nscan = 1000, burn = 500, odens = 1,
												 quiet = TRUE)
	expect_true(is.numeric(mem))
	expect_true(mem > 0)
})

test_that("estimate_memory scales with network size", {
	mem_small <- estimate_memory(n_row = 10, quiet = TRUE)
	mem_large <- estimate_memory(n_row = 50, quiet = TRUE)
	expect_true(mem_large > mem_small)
})

test_that("plot_ppc_ecdf works for static model", {
	set.seed(117)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 117)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)
	ppd <- posterior_predict_dbn(fit, ndraws = 5, seed = 42)

	p <- plot_ppc_ecdf(fit, ppd, Y_obs = sim$Y)
	expect_s3_class(p, "ggplot")
})

test_that("plot_ppc_density works for static model", {
	set.seed(118)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 118)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
						 nscan = 30, burn = 10, verbose = FALSE)
	ppd <- posterior_predict_dbn(fit, ndraws = 5, seed = 42)

	p <- plot_ppc_density(fit, ppd, Y_obs = sim$Y)
	expect_s3_class(p, "ggplot")
})

test_that("plot_regime_probs works for HMM model", {
	set.seed(119)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 119)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	p <- plot_regime_probs(fit)
	expect_s3_class(p, "ggplot")
})

test_that("predict works for HMM model", {
	set.seed(120)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 120)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
						 R = 2, nscan = 30, burn = 10, verbose = FALSE)

	pred <- predict(fit, H = 3, draws = 10, summary = "mean")
	expect_equal(dim(pred)[4], 3) # 3 forecast steps
})

test_that("predict works for lowrank model", {
	set.seed(121)
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 8, r = 3, seed = 121)
	fit <- dbn(sim$Y, model = "lowrank", family = "ordinal",
						 r = 3, nscan = 30, burn = 10, verbose = FALSE)

	pred <- predict(fit, H = 3, draws = 10, summary = "mean")
	expect_equal(dim(pred)[4], 3) # 3 forecast steps
})
