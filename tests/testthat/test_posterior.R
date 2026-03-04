#### s3 dispatch and plotting

test_that("plot.dbn dispatches correctly", {
	# create mock static results
	static_results <- list(
		model = "static",
		params = matrix(rnorm(30), ncol = 3),
		B = list(array(rnorm(64), dim = c(4, 4, 10)))
	)
	colnames(static_results$params) <- c("s2", "t2", "g2")
	class(static_results) <- "dbn"

	expect_silent(p <- plot(static_results))
	expect_true(inherits(p, c("gtable", "gTree", "grob")))

	# create mock dynamic results
	n <- 10
	k <- 3
	Tt <- 20
	nscan <- 10
	dynamic_results <- list(
		model = "dynamic",
		sigma2 = rnorm(nscan, 1, 0.1),
		tauA2 = rnorm(nscan, 0.5, 0.05),
		tauB2 = rnorm(nscan, 0.5, 0.05),
		g2 = rnorm(nscan, 0.3, 0.03),
		A = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
		B = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
		dims = list(m = n, n_row = n, n_col = n, Tt = Tt, p = 1, is_bipartite = FALSE)
	)
	class(dynamic_results) <- "dbn"

	expect_silent(plots <- plot(dynamic_results))
	expect_type(plots, "list")
})

test_that("tidy_dbn extracts posterior means correctly", {
	n <- 5
	k <- 3
	Tt <- 10
	nscan <- 20

	results <- list(
		model = "dynamic",
		A = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
		B = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
		M = array(rnorm(n*n*3), dim = c(n, n, 3)),
		dims = list(m = n, n_row = n, n_col = n, p = 3, Tt = Tt, is_bipartite = FALSE)
	)
	class(results) <- "dbn"

	tidy <- tidy_dbn(results)

	expect_type(tidy, "list")
	expect_equal(dim(tidy$A), c(n, k, Tt))
	expect_equal(dim(tidy$B), c(n, k, Tt))

	static_results <- list(
		model = "static",
		M = array(rnorm(n*n*3), dim = c(n, n, 3)),
		B = list(array(rnorm(n*n*10), dim = c(n, n, 10))),
		dims = list(m = n, n_row = n, n_col = n, p = 3, is_bipartite = FALSE)
	)
	class(static_results) <- "dbn"

	tidy_static <- tidy_dbn(static_results)
	expect_equal(dim(tidy_static$Theta), c(n, n, 3))
})

test_that("dyad_path handles multiple relations correctly", {
	n <- 10
	p <- 3
	Tt <- 20
	k <- 4
	nscan <- 50

	results <- list(
		model = "dynamic",
		A = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
		B = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
		dims = list(m = n, n_row = n, n_col = n, p = p, Tt = Tt, is_bipartite = FALSE)
	)
	class(results) <- "dbn"

	expect_silent(p1 <- dyad_path(results, i = 1, j = 5, rel = 1))
	expect_s3_class(p1, "ggplot")

	expect_silent(p2 <- dyad_path(results, i = 1, j = 5, rel = c(1, 2), facet = TRUE))
	expect_s3_class(p2, "ggplot")

	expect_silent(p3 <- dyad_path(results, i = 1, j = 5, rel = NULL))
	expect_s3_class(p3, "ggplot")
})

test_that("ppc_ecdf guards against out-of-range values", {
	n <- 8
	Tt <- 5
	nscan <- 50

	Y <- array(sample(1:5, n*n*Tt, replace = TRUE), dim = c(n, n, Tt))
	diag(Y[,,1]) <- NA

	results <- list(
		model = "static",
		R = Y,
		M = array(rnorm(n*n, mean = 3, sd = 0.5), dim = c(n, n)),
		B = list(
			array(rnorm(n*n*nscan), dim = c(n, n, nscan)),
			array(rnorm(n*n*nscan), dim = c(n, n, nscan))
		),
		dims = list(m = n, n_row = n, n_col = n, p = 1, Tt = Tt, is_bipartite = FALSE)
	)
	class(results) <- "dbn"

	expect_silent(p <- ppc_ecdf(results, n_rep = 20))
	expect_s3_class(p, "ggplot")
})

test_that("group influence functions work correctly", {
	n <- 20
	Tt <- 15
	nscan <- 50

	results <- list(
		model = "dynamic",
		A = lapply(1:nscan, function(i) array(rnorm(n*n*Tt), dim = c(n, n, Tt))),
		B = lapply(1:nscan, function(i) array(rnorm(n*n*Tt), dim = c(n, n, Tt))),
		dims = list(m = n, n_row = n, n_col = n, p = 1, Tt = Tt, is_bipartite = FALSE)
	)
	class(results) <- "dbn"

	group <- c(1, 3, 5, 7, 9)

	for (measure in c("rowsum", "rowmean", "l2")) {
		p <- plot_group_influence(results, group = group,
			type = "sender", measure = measure)
		expect_s3_class(p, "ggplot")
	}

	influence_data <- get_group_influence(results, group = group,
		type = "target", measure = "rowsum")
	expect_type(influence_data, "list")
	expect_true("time" %in% names(influence_data))
	expect_true("mean" %in% names(influence_data))
	expect_true("q0.025" %in% names(influence_data))
	expect_true("q0.5" %in% names(influence_data))
	expect_true("q0.975" %in% names(influence_data))

	groups <- list(
		"Group A" = 1:3,
		"Group B" = 4:6,
		"Group C" = 7:10
	)

	p_comp <- compare_group_influence(results, groups = groups,
		type = "sender", measure = "l2")
	expect_s3_class(p_comp, "ggplot")
})

test_that("check_convergence handles missing coda gracefully", {
	results <- list(
		model = "static",
		params = matrix(rnorm(300), ncol = 3)
	)
	colnames(results$params) <- c("s2", "t2", "g2")
	class(results) <- "dbn"

	output <- capture.output(check_convergence(results))
	expect_true(length(output) > 0)
})

test_that("role_trajectory works for both sender and receiver", {
	n <- 15
	k <- 5
	Tt <- 20
	nscan <- 50

	results <- list(
		model = "dynamic",
		A = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
		B = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
		dims = list(m = n, n_row = n, n_col = n, p = 1, Tt = Tt, is_bipartite = FALSE)
	)
	class(results) <- "dbn"

	expect_silent(role_trajectory(results, mat = "A", comp = 1))
	expect_silent(role_trajectory(results, mat = "B", comp = 1))
})

test_that("predict_ordinal converts correctly", {
	n <- 8
	p <- 2
	Tt <- 10

	Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
	for (t in 1:Tt) {
		for (r in 1:p) {
			diag(Y[,,r,t]) <- NA
		}
	}

	n_samples <- 10
	results <- list(
		model = "static",
		Y = Y,
		R = Y,
		M = array(3, dim = c(n, n, p)),
		B = list(
			array(rnorm(n*n*n_samples), dim = c(n, n, n_samples)),
			array(rnorm(n*n*n_samples), dim = c(n, n, n_samples)),
			array(rnorm(p*p*n_samples), dim = c(p, p, n_samples))
		),
		params = matrix(runif(n_samples * 3, 0.1, 1), ncol = 3),
		dims = list(m = n, n_row = n, n_col = n, p = p, Tt = Tt, is_bipartite = FALSE)
	)
	colnames(results$params) <- c("s2", "t2", "g2")
	class(results) <- "dbn"

	ordinal_pred <- predict_ordinal(results, draws = 20)

	vals <- unique(c(ordinal_pred))
	vals <- vals[!is.na(vals)]
	expect_true(all(vals %in% 1:5))
})

#### posterior extraction and summaries

test_that("theta_slice extracts correct values from HMM model", {
	set.seed(101)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 101)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
		R = 2, nscan = 30, burn = 10, verbose = FALSE)

	expect_true(!is.null(fit$draws$theta))
	n_draws <- length(fit$draws$theta)

	slice1 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = 1)
	expect_equal(length(slice1), n_draws)

	slice2 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = c(1, 3, 5))
	expect_equal(nrow(slice2), n_draws)
	expect_equal(ncol(slice2), 3)

	slice3 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = 1, draws = 1:5)
	expect_equal(length(slice3), 5)
})

test_that("theta_summary computes summaries from HMM model", {
	set.seed(102)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 102)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
		R = 2, nscan = 30, burn = 10, verbose = FALSE)

	s1 <- theta_summary(fit, fun = mean, i = 1, j = 2, rel = 1, time = 1)
	expect_true(!is.null(s1))
	expect_true("value" %in% names(s1))
	expect_equal(nrow(s1), 1)

	s2 <- theta_summary(fit, fun = mean, time = 1)
	expect_true(nrow(s2) > 0)
})

test_that("theta_slice returns NULL for static model (no theta draws)", {
	set.seed(103)
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 103)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 30, burn = 10, verbose = FALSE)

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
	expect_true(nrow(ps) >= 3)
})

test_that("param_summary works for dynamic model", {
	set.seed(105)
	sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 8, seed = 105)
	fit <- dbn(sim$Y, model = "dynamic", family = "ordinal",
		nscan = 30, burn = 10, odens = 1, verbose = FALSE)

	ps <- param_summary(fit)
	expect_true(!is.null(ps))
	expect_true(nrow(ps) >= 2)
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

	m_means <- derive_draws(fit, function(draw) {
		mean(draw$M[,,1], na.rm = TRUE)
	}, name = "M_mean")

	expect_true(!is.null(m_means))
	expect_true(is.numeric(m_means))
})

test_that("theta_credible returns correct structure", {
	set.seed(113)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 113)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
		R = 2, nscan = 30, burn = 10, verbose = FALSE)

	tc <- theta_credible(fit, i = 1:3, j = 1:3, time = 1:4)
	expect_true(!is.null(tc))
	expect_equal(nrow(tc), 3 * 3 * 4)
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
		expect_equal(nrow(ns), 8)
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
	expect_true(all(is.na(diag(ep))))
})

test_that("estimate_memory scales with network size", {
	mem_small <- estimate_memory(n_row = 10, quiet = TRUE)
	mem_large <- estimate_memory(n_row = 50, quiet = TRUE)
	expect_true(mem_large > mem_small)
})

test_that("predict works for HMM model", {
	set.seed(120)
	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 120)
	fit <- dbn(sim$Y, model = "hmm", family = "ordinal",
		R = 2, nscan = 30, burn = 10, verbose = FALSE)

	pred <- predict(fit, H = 3, draws = 10, summary = "mean")
	expect_equal(dim(pred)[4], 3)
})

test_that("predict works for lowrank model", {
	set.seed(121)
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 8, r = 3, seed = 121)
	fit <- dbn(sim$Y, model = "lowrank", family = "ordinal",
		r = 3, nscan = 30, burn = 10, verbose = FALSE)

	pred <- predict(fit, H = 3, draws = 10, summary = "mean")
	expect_equal(dim(pred)[4], 3)
})

#### visualization

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
