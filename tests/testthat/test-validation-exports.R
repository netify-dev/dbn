####
# Phase 3: Exported function test coverage
# Every exported function should have at least a smoke test.
# Also tests S3 dispatch for all model types.
####

# shared fixtures — small fits for each model type
setup_fits <- function() {
	sim_s <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3001)
	fit_s <- dbn(sim_s$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	sim_d <- simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 3002)
	fit_d <- dbn(sim_d$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	sim_l <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3003)
	fit_l <- dbn(sim_l$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	sim_h <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 3004)
	fit_h <- dbn(sim_h$Y, model = "hmm", R = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	list(
		static = list(sim = sim_s, fit = fit_s),
		dynamic = list(sim = sim_d, fit = fit_d),
		lowrank = list(sim = sim_l, fit = fit_l),
		hmm = list(sim = sim_h, fit = fit_h)
	)
}

####
# S3 method dispatch: print, summary, plot for all 4 model types
####

test_that("print.dbn works for all model types", {
	skip_on_cran()
	fits <- setup_fits()

	for (nm in names(fits)) {
		expect_output(print(fits[[nm]]$fit), regexp = NULL,
			label = sprintf("print works for %s", nm))
	}
})

test_that("summary works for all model types", {
	skip_on_cran()
	fits <- setup_fits()

	for (nm in names(fits)) {
		res <- tryCatch(
			{ capture.output(summary_dbn(fits[[nm]]$fit)); TRUE },
			error = function(e) e$message
		)
		expect_true(isTRUE(res),
			label = sprintf("summary works for %s (error: %s)", nm,
				if (isTRUE(res)) "none" else res))
	}
})

test_that("plot works for all model types (null device)", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")
	fits <- setup_fits()

	pdf(NULL) # null device
	on.exit(dev.off(), add = TRUE)

	for (nm in names(fits)) {
		res <- tryCatch(
			{ plot_dbn(fits[[nm]]$fit); TRUE },
			error = function(e) e$message
		)
		expect_true(isTRUE(res),
			label = sprintf("plot works for %s (error: %s)", nm,
				if (isTRUE(res)) "none" else res))
	}
})

####
# predict dispatch: PPD and forecast for applicable models
####

test_that("predict returns PPD for all model types", {
	skip_on_cran()
	fits <- setup_fits()

	for (nm in names(fits)) {
		ppd <- tryCatch(
			posterior_predict_dbn(fits[[nm]]$fit, ndraws = 5, seed = 42),
			error = function(e) NULL
		)
		expect_true(!is.null(ppd),
			label = sprintf("PPD succeeds for %s", nm))
	}
})

####
# previously untested exported functions
####

test_that("dbn_lowrank_accurate runs and returns valid output", {
	skip_on_cran()
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3010)

	fit <- dbn_lowrank_accurate(sim$Y, r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "lowrank")
	expect_true(!is.null(fit$U), label = "U stored")
	expect_true(!is.null(fit$alpha), label = "alpha stored")
})

test_that("tidy_dbn_lowrank returns data.frame with expected columns", {
	skip_on_cran()
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3011)
	fit <- dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	tidy_df <- tidy_dbn_lowrank(fit)
	expect_true(is.data.frame(tidy_df), label = "returns data.frame")
	expect_true(nrow(tidy_df) > 0, label = "non-empty data.frame")
	# should have columns for actor, factor, time, value, etc.
	expect_true(ncol(tidy_df) >= 3, label = "has at least 3 columns")
})

test_that("debug_irf runs without error", {
	skip_on_cran()
	# debug_irf calls impulse_response_dynamic which requires dynamic model
	sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 3012)
	fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	expect_no_error(
		capture.output(debug_irf(fit, draw_idx = 1, shock_i = 1, shock_j = 2))
	)
})

test_that("compare_dbn works with two fits", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3013)

	fit1 <- dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)
	fit2 <- dbn(sim$Y, model = "static", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	result <- compare_dbn(fit1, fit2)
	expect_true(!is.null(result), label = "compare_dbn returns non-NULL")
})

####
# network statistic functions
####

test_that("stat_density, stat_in_degree, stat_out_degree work", {
	n <- 6
	Theta <- matrix(rnorm(n * n), n, n)
	diag(Theta) <- NA

	d <- stat_density(Theta)
	expect_true(is.numeric(d))
	expect_length(d, 1)

	in_deg <- stat_in_degree(Theta)
	expect_true(is.numeric(in_deg))

	out_deg <- stat_out_degree(Theta)
	expect_true(is.numeric(out_deg))
})

test_that("stat_reciprocity and stat_transitivity work", {
	n <- 6
	Theta <- matrix(rnorm(n * n), n, n)
	diag(Theta) <- NA

	r <- stat_reciprocity(Theta)
	expect_true(is.numeric(r))
	expect_length(r, 1)

	tr <- stat_transitivity(Theta)
	expect_true(is.numeric(tr))
	expect_length(tr, 1)
})

####
# IRF functions
####

test_that("build_shock produces valid shock matrices", {
	n <- 8

	S_density <- build_shock(n, type = "density")
	expect_true(is.matrix(S_density))
	expect_equal(dim(S_density), c(n, n))

	S_edge <- build_shock(n, type = "unit_edge", i = 1, j = 2)
	expect_equal(S_edge[1, 2], 1)
	expect_equal(sum(S_edge != 0, na.rm = TRUE), 1)

	S_node <- build_shock(n, type = "node_out", i = 1)
	expect_true(all(S_node[1, ] == 1 | is.na(S_node[1, ])))
})

test_that("compute_irf returns list with expected structure", {
	skip_on_cran()
	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3014)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	S <- build_shock(6, type = "density")
	irf <- compute_irf(fit, shock = S, H = 3, n_draws = 3,
		stat_fun = dbn::stat_density)
	expect_true(is.list(irf), label = "compute_irf returns list")
})

####
# posterior extraction functions
####

test_that("param_summary returns data.frame for all model types", {
	skip_on_cran()
	fits <- setup_fits()

	for (nm in names(fits)) {
		ps <- param_summary(fits[[nm]]$fit)
		expect_true(is.data.frame(ps),
			label = sprintf("param_summary returns df for %s", nm))
		expect_true(nrow(ps) > 0,
			label = sprintf("param_summary non-empty for %s", nm))
	}
})

test_that("theta_summary and theta_credible work", {
	skip_on_cran()
	# use lowrank which stores draws$theta in list format
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3015)
	fit <- dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ts <- theta_summary(fit, i = 1, j = 2, rel = 1, time = 1)
	expect_true(!is.null(ts), label = "theta_summary returns non-NULL")

	tc <- theta_credible(fit, i = 1, j = 2, rel = 1, time = 1)
	expect_true(is.data.frame(tc), label = "theta_credible returns df")
})

test_that("edge_prob returns matrix", {
	skip_on_cran()
	# use lowrank which stores draws$theta in list format
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3016)
	fit <- dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ep <- edge_prob(fit)
	expect_true(is.matrix(ep) || is.array(ep), label = "edge_prob returns matrix")
})

test_that("derive_draws works", {
	skip_on_cran()
	# use lowrank which stores draws$theta in list format
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3017)
	fit <- dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	dd <- derive_draws(fit, fun = function(draw) mean(draw$theta, na.rm = TRUE))
	expect_true(!is.null(dd), label = "derive_draws returns non-NULL")
})

test_that("network_summary works", {
	skip_on_cran()
	# use lowrank which stores draws$theta in list format
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3018)
	fit <- dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ns <- network_summary(fit)
	expect_true(!is.null(ns), label = "network_summary returns non-NULL")
})

####
# simulation functions
####

test_that("simulate_test_data produces valid output", {
	td <- simulate_test_data(n = 5, time = 3)
	# simulate_test_data returns an array directly
	expect_true(is.array(td))
	expect_equal(dim(td)[1], 5)
	expect_equal(dim(td)[4], 3)
})

####
# thread management
####

test_that("get_dbn_threads returns numeric", {
	nt <- get_dbn_threads()
	expect_true(is.numeric(nt))
	expect_true(nt >= 1)
})

test_that("set_dbn_threads and get_dbn_threads round-trip", {
	old <- get_dbn_threads()
	set_dbn_threads(1)
	expect_equal(get_dbn_threads(), 1)
	set_dbn_threads(old)
})

####
# estimate_memory
####

test_that("estimate_memory returns sensible estimate", {
	est <- estimate_memory(n_row = 20, Tt = 30, p = 1, nscan = 5000, odens = 5)
	expect_true(is.numeric(est) || is.list(est))
})

####
# visualization functions (null device)
####

test_that("plot_trace produces output", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 3020)
	fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	p_trace <- plot_trace(fit)
	expect_true(!is.null(p_trace), label = "plot_trace returns non-NULL")
})

test_that("plot_theta produces output with lowrank model", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	# use lowrank which stores draws$theta
	sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3025)
	fit <- dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	p_theta <- plot_theta(fit, time = 1, rel = 1)
	expect_true(!is.null(p_theta), label = "plot_theta returns non-NULL")
})

test_that("plot_ppc_ecdf and plot_ppc_density produce output", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3021)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ppd <- tryCatch(
		posterior_predict_dbn(fit, ndraws = 5, seed = 42),
		error = function(e) NULL
	)

	if (!is.null(ppd)) {
		p_ecdf <- plot_ppc_ecdf(fit, ppd = ppd, Y_obs = sim$Y)
		expect_true(!is.null(p_ecdf), label = "plot_ppc_ecdf returns non-NULL")

		p_dens <- plot_ppc_density(fit, ppd = ppd, Y_obs = sim$Y)
		expect_true(!is.null(p_dens), label = "plot_ppc_density returns non-NULL")
	}
})

test_that("regime_probs and plot_regime_probs work for HMM", {
	skip_on_cran()

	sim <- simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 3022)
	fit <- dbn(sim$Y, model = "hmm", R = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	rp <- regime_probs(fit)
	expect_true(is.matrix(rp) || is.data.frame(rp), label = "regime_probs returns matrix/df")

	if (requireNamespace("ggplot2", quietly = TRUE)) {
		p <- plot_regime_probs(fit)
		expect_true(!is.null(p), label = "plot_regime_probs returns non-NULL")
	}
})

test_that("dyad_path and net_snapshot produce output", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 3023)
	fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	p_dyad <- dyad_path(fit, i = 1, j = 2, rel = 1)
	expect_true(!is.null(p_dyad), label = "dyad_path returns non-NULL")

	p_net <- net_snapshot(fit, t = 1, rel = 1)
	expect_true(!is.null(p_net), label = "net_snapshot returns non-NULL")
})

test_that("latent_summary works", {
	skip_on_cran()

	sim <- simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3024)
	fit <- dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ls <- latent_summary(fit)
	expect_true(!is.null(ls), label = "latent_summary returns non-NULL")
})
