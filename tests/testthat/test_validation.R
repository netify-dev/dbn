####
#### convergence diagnostics (skip_on_cran + DBN_VALIDATION)
####

skip_validation = function() {
	skip_on_cran()
	skip_if(
		!identical(Sys.getenv("DBN_VALIDATION"), "true"),
		"Set DBN_VALIDATION=true to run extended validation tests"
	)
}

test_that("validation: check_convergence produces diagnostic output", {
	skip_validation()

	sim = simulate_static_dbn(
		n = 10, p = 1, time = 10,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 9003
	)

	fit = dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 2000, burn = 1000, odens = 2,
		verbose = FALSE
	)

	diag_result = check_convergence(fit)
	expect_true(!is.null(diag_result), label = "check_convergence returns non-NULL")

	if (is.list(diag_result)) {
		expect_true(length(diag_result) > 0,
			label = "check_convergence returns non-empty result")
	}
})

####
#### PPD calibration (skip_on_cran + DBN_VALIDATION)
####

compute_ppd_coverage = function(fit, obs, ndraws = 100, prob = 0.90) {
	ppd = tryCatch(
		posterior_predict_dbn(fit, ndraws = ndraws, seed = 42),
		error = function(e) NULL
	)
	if (is.null(ppd) || is.null(ppd$y_rep) || !is.array(ppd$y_rep)) {
		return(list(coverage = NA, ppd = ppd))
	}

	yrep = ppd$y_rep
	nd = length(dim(yrep))
	alpha = (1 - prob) / 2
	q_lo = apply(yrep, 1:(nd - 1), quantile, probs = alpha, na.rm = TRUE)
	q_hi = apply(yrep, 1:(nd - 1), quantile, probs = 1 - alpha, na.rm = TRUE)

	if (length(dim(q_lo)) >= 3 && dim(q_lo)[3] >= 1) {
		in_interval = (obs >= q_lo[, , 1, ]) & (obs <= q_hi[, , 1, ])
	} else {
		in_interval = (obs >= q_lo) & (obs <= q_hi)
	}
	coverage = mean(in_interval, na.rm = TRUE)
	list(coverage = coverage, ppd = ppd)
}

test_that("validation: static ordinal PPD produces valid ordinal output", {
	skip_validation()

	sim = simulate_static_dbn(
		n = 10, p = 1, time = 15,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 6003
	)

	fit = dbn(sim$Y,
		model = "static", family = "ordinal",
		nscan = 3000, burn = 1000, odens = 2,
		verbose = FALSE
	)

	ppd = tryCatch(
		posterior_predict_dbn(fit, ndraws = 50, seed = 42),
		error = function(e) NULL
	)

	expect_true(!is.null(ppd), label = "PPD generation succeeds for ordinal")

	if (!is.null(ppd) && !is.null(ppd$y_rep)) {
		yrep = ppd$y_rep
		obs_mean = mean(sim$Y, na.rm = TRUE)
		ppd_means = apply(yrep, length(dim(yrep)), function(x) mean(x, na.rm = TRUE))
		ppd_mean_avg = mean(ppd_means)

		expect_gt(ppd_mean_avg, obs_mean / 3,
			label = "PPD mean not drastically below observed")
		expect_lt(ppd_mean_avg, obs_mean * 3,
			label = "PPD mean not drastically above observed")
	}
})

test_that("validation: static gaussian PPD p-value is non-extreme", {
	skip_validation()

	sim = simulate_static_dbn(
		n = 12, p = 1, time = 15,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 6004
	)

	fit = dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	ppd = tryCatch(
		posterior_predict_dbn(fit, ndraws = 200, seed = 42),
		error = function(e) NULL
	)

	expect_true(!is.null(ppd), label = "PPD generation succeeds")

	if (!is.null(ppd) && !is.null(ppd$y_rep)) {
		yrep = ppd$y_rep
		obs_stat = mean(sim$Z, na.rm = TRUE)
		nd = length(dim(yrep))
		ppd_stats = apply(yrep, nd, function(x) mean(x, na.rm = TRUE))
		pp_pval = mean(ppd_stats >= obs_stat)

		expect_gt(pp_pval, 0.01,
			label = sprintf("PPD p-value %.3f > 0.01 (not extreme low)", pp_pval))
		expect_lt(pp_pval, 0.99,
			label = sprintf("PPD p-value %.3f < 0.99 (not extreme high)", pp_pval))
	}
})

####
#### parameter recovery (skip_on_cran + DBN_VALIDATION)
####

extract_M_mean = function(fit) {
	if (fit$model == "static") {
		M_draws = fit$draws$misc$M
		if (!is.null(M_draws) && length(M_draws) > 0) {
			Reduce("+", M_draws) / length(M_draws)
		} else {
			NULL
		}
	} else if (fit$model == "dynamic") {
		if (!is.null(fit$M) && is.array(fit$M)) {
			nd = length(dim(fit$M))
			apply(fit$M, 1:(nd - 1), mean)
		} else if (!is.null(fit$draws$misc$M)) {
			Reduce("+", fit$draws$misc$M) / length(fit$draws$misc$M)
		} else {
			NULL
		}
	} else if (fit$model %in% c("lowrank", "hmm")) {
		M_draws = fit$draws$misc$M
		if (!is.null(M_draws) && length(M_draws) > 0) {
			Reduce("+", M_draws) / length(M_draws)
		} else if (!is.null(fit$M) && is.list(fit$M)) {
			Reduce("+", fit$M) / length(fit$M)
		} else {
			NULL
		}
	} else {
		NULL
	}
}

covers_truth = function(draws, truth, prob = 0.90) {
	q = quantile(draws, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2),
		na.rm = TRUE)
	truth >= q[1] && truth <= q[2]
}

safe_cor = function(x, y) {
	valid = complete.cases(x, y)
	if (sum(valid) < 5) return(NA_real_)
	cor(x[valid], y[valid])
}

test_that("validation: static gaussian recovers M and sigma2", {
	skip_validation()

	n = 15
	Tt = 20
	true_sigma2 = 0.3
	true_tau2 = 0.05

	sim = simulate_static_dbn(
		n = n, p = 1, time = Tt,
		sigma2 = true_sigma2, tau2 = true_tau2,
		seed = 7001
	)

	fit = dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	M_hat = extract_M_mean(fit)
	expect_false(is.null(M_hat), label = "M posterior mean extracted")
	r_M = safe_cor(as.vector(sim$M[, , 1]), as.vector(M_hat[, , 1]))
	expect_gt(r_M, 0.7,
		label = sprintf("M recovery r=%.3f > 0.7 (static gaussian)", r_M))

	sigma2_draws = fit$draws$pars$s2
	expect_true(
		covers_truth(sigma2_draws, true_sigma2, 0.90),
		label = sprintf(
			"sigma2 90%% CI [%.3f, %.3f] covers truth %.3f",
			quantile(sigma2_draws, 0.05),
			quantile(sigma2_draws, 0.95),
			true_sigma2
		)
	)

	sigma2_hat = mean(sigma2_draws)
	expect_lt(abs(sigma2_hat - true_sigma2) / true_sigma2, 1.0,
		label = "sigma2 relative error < 100%")
})

test_that("validation: static ordinal recovers M direction", {
	skip_validation()

	n = 15
	Tt = 20

	sim = simulate_static_dbn(
		n = n, p = 1, time = Tt,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 7002
	)

	fit = dbn(sim$Y,
		model = "static", family = "ordinal",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	M_hat = extract_M_mean(fit)
	expect_false(is.null(M_hat), label = "M posterior mean extracted")
	r_M = safe_cor(as.vector(sim$M[, , 1]), as.vector(M_hat[, , 1]))
	expect_gt(r_M, 0.5,
		label = sprintf("M recovery r=%.3f > 0.5 (static ordinal)", r_M))

	sigma2_draws = fit$draws$pars$s2
	expect_true(all(sigma2_draws > 0), label = "sigma2 all positive")
	expect_lt(mean(sigma2_draws), 10, label = "sigma2 mean < 10")
})

test_that("validation: static binary recovers M sign pattern", {
	skip_validation()

	n = 15
	Tt = 30

	set.seed(7003)
	M_true = matrix(rnorm(n * n, 0, 0.5), n, n)
	diag(M_true) = NA
	Z = array(NA, c(n, n, 1, Tt))
	Y = array(NA, c(n, n, 1, Tt))
	for (t in 1:Tt) {
		Z[, , 1, t] = M_true + matrix(rnorm(n * n), n, n)
		Y[, , 1, t] = ifelse(Z[, , 1, t] > 0, 1, 0)
		diag(Y[, , 1, t]) = NA
	}

	fit = dbn(Y,
		model = "static", family = "binary",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	M_hat = extract_M_mean(fit)
	expect_false(is.null(M_hat), label = "M posterior mean extracted")
	M_hat_mat = if (length(dim(M_hat)) == 3) M_hat[, , 1] else M_hat
	r_M = safe_cor(as.vector(M_true), as.vector(M_hat_mat))
	expect_gt(r_M, 0.4,
		label = sprintf("M recovery r=%.3f > 0.4 (static binary)", r_M))

	valid = !is.na(as.vector(M_true)) & !is.na(as.vector(M_hat_mat))
	sign_match = mean(sign(as.vector(M_true)[valid]) ==
		sign(as.vector(M_hat_mat)[valid]))
	expect_gt(sign_match, 0.65,
		label = sprintf("M sign recovery %.1f%% > 65%%", sign_match * 100))
})

test_that("validation: dynamic gaussian recovers Theta and variance params", {
	skip_validation()

	n = 12
	Tt = 20
	true_sigma2 = 0.3
	true_tauA2 = 0.05
	true_tauB2 = 0.05

	sim = simulate_dynamic_dbn(
		n = n, p = 1, time = Tt,
		sigma2 = true_sigma2, tauA2 = true_tauA2, tauB2 = true_tauB2,
		seed = 7004
	)

	fit = dbn(sim$Z,
		model = "dynamic", family = "gaussian",
		nscan = 8000, burn = 3000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "dynamic")

	if (!is.null(fit$Theta) && is.array(fit$Theta)) {
		theta_mean = apply(fit$Theta, 1:4, mean)
		t_last = dim(theta_mean)[4]
		t_true_last = dim(sim$Theta)[4]
		th_hat = as.vector(theta_mean[, , 1, t_last])
		th_true = as.vector(sim$Theta[, , 1, t_true_last])
		r_theta = safe_cor(th_hat, th_true)
		expect_gt(r_theta, 0.4,
			label = sprintf("Theta(T) recovery r=%.3f > 0.4 (dynamic gaussian)", r_theta))
	}

	expect_true(
		covers_truth(fit$sigma2, true_sigma2, 0.90),
		label = sprintf(
			"sigma2 90%% CI [%.3f, %.3f] covers truth %.3f",
			quantile(fit$sigma2, 0.05),
			quantile(fit$sigma2, 0.95),
			true_sigma2
		)
	)

	expect_true(all(is.finite(fit$tau_A2)), label = "tau_A2 finite")
	expect_true(all(fit$tau_A2 > 0), label = "tau_A2 positive")

	expect_true(all(is.finite(fit$tau_B2)), label = "tau_B2 finite")
	expect_true(all(fit$tau_B2 > 0), label = "tau_B2 positive")

	M_hat = extract_M_mean(fit)
	if (!is.null(M_hat)) {
		M_hat_mat = if (length(dim(M_hat)) >= 3) M_hat[, , 1] else M_hat
		M_true_mat = sim$M[, , 1]
		r_M = safe_cor(as.vector(M_true_mat), as.vector(M_hat_mat))
		expect_gt(r_M, 0.5,
			label = sprintf("M recovery r=%.3f > 0.5 (dynamic gaussian)", r_M))
	}
})

test_that("validation: dynamic ordinal completes and produces valid output", {
	skip_validation()

	sim = simulate_dynamic_dbn(
		n = 12, p = 1, time = 20,
		sigma2 = 0.3, tauA2 = 0.05, tauB2 = 0.05,
		seed = 7005
	)

	fit = dbn(sim$Y,
		model = "dynamic", family = "ordinal",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	expect_true(all(is.finite(fit$sigma2)), label = "sigma2 all finite")
	expect_true(all(fit$sigma2 > 0), label = "sigma2 all positive")
	expect_gt(mean(fit$sigma2), 0.01, label = "sigma2 mean > 0.01")

	expect_true(!is.null(fit$A), label = "A draws stored")
	expect_true(length(fit$A) > 0, label = "A draws non-empty")
})

test_that("validation: lowrank recovers non-degenerate alpha trajectories", {
	skip_validation()

	sim = simulate_lowrank_dbn(
		n = 20, p = 1, time = 15,
		r = 2, sigma2 = 0.2,
		seed = 7006
	)

	fit = dbn(sim$Y,
		model = "lowrank", r = 2,
		family = "ordinal",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	alpha_draws = fit$alpha
	expect_true(!is.null(alpha_draws), label = "alpha draws stored")
	expect_true(length(alpha_draws) > 0, label = "alpha draws non-empty")

	alpha_last = alpha_draws[[length(alpha_draws)]]
	for (r_idx in seq_len(nrow(alpha_last))) {
		alpha_sd = sd(alpha_last[r_idx, ])
		expect_gt(alpha_sd, 1e-6,
			label = sprintf("alpha[%d] has temporal variation (sd=%.4f)", r_idx, alpha_sd))
	}

	sigma2_draws = fit$draws$pars$sigma2_proc
	expect_true(all(is.finite(sigma2_draws)), label = "sigma2 finite")
	expect_true(all(sigma2_draws > 0), label = "sigma2 positive")

	U_last = fit$U[[length(fit$U)]]
	UtU = crossprod(U_last)
	expect_equal(UtU, diag(ncol(U_last)), tolerance = 0.1,
		label = "U is approximately orthonormal")
})

test_that("validation: hmm recovers regime assignments", {
	skip_validation()

	sim = simulate_hmm_dbn(
		n = 10, p = 1, time = 30,
		R = 2, transition_prob = 0.9,
		sigma2 = 0.3,
		seed = 7007
	)

	fit = dbn(sim$Y,
		model = "hmm", R = 2,
		family = "ordinal",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	rp = regime_probs(fit)
	expect_true(!is.null(rp), label = "regime_probs returns non-NULL")
	expect_equal(nrow(rp), sim$dims$Tt, label = "regime probs has correct time dim")

	modal_regime = apply(rp, 1, which.max)
	true_S = sim$S

	agreement_direct = mean(modal_regime == true_S, na.rm = TRUE)
	agreement_flipped = mean(modal_regime == (3 - true_S), na.rm = TRUE)
	best_agreement = max(agreement_direct, agreement_flipped)

	expect_gt(best_agreement, 0.55,
		label = sprintf(
			"Regime recovery %.1f%% > 55%% (best of direct/flipped)",
			best_agreement * 100
		))

	Pi_draws = fit$Pi
	if (!is.null(Pi_draws) && length(Pi_draws) > 0) {
		Pi_mean = Reduce("+", Pi_draws) / length(Pi_draws)
		for (r_idx in 1:2) {
			expect_gt(Pi_mean[r_idx, r_idx], 0.5,
				label = sprintf("Pi[%d,%d]=%.2f > 0.5 (sticky regimes)", r_idx, r_idx, Pi_mean[r_idx, r_idx]))
		}
	}

	sigma2_draws = fit$draws$pars$sigma2_proc
	expect_true(all(is.finite(sigma2_draws)), label = "sigma2 finite")
})

test_that("validation: static gaussian sigma2 coverage across seeds", {
	skip_validation()

	true_sigma2 = 0.3
	seeds = c(8001, 8002, 8003)
	covers = logical(length(seeds))

	for (i in seq_along(seeds)) {
		sim = simulate_static_dbn(
			n = 15, p = 1, time = 20,
			sigma2 = true_sigma2, tau2 = 0.05,
			seed = seeds[i]
		)

		fit = dbn(sim$Z,
			model = "static", family = "gaussian",
			nscan = 5000, burn = 2000, odens = 2,
			verbose = FALSE
		)

		covers[i] = covers_truth(fit$draws$pars$s2, true_sigma2, 0.90)
	}

	expect_gte(sum(covers), 2,
		label = sprintf(
			"sigma2 covered in %d/3 seeds (expect >= 2)",
			sum(covers)
		))
})

####
#### exported function coverage (skip_on_cran)
####

setup_fits = function() {
	sim_s = simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3001)
	fit_s = dbn(sim_s$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	sim_d = simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 3002)
	fit_d = dbn(sim_d$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	sim_l = simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3003)
	fit_l = dbn(sim_l$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	sim_h = simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 3004)
	fit_h = dbn(sim_h$Y, model = "hmm", R = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	list(
		static = list(sim = sim_s, fit = fit_s),
		dynamic = list(sim = sim_d, fit = fit_d),
		lowrank = list(sim = sim_l, fit = fit_l),
		hmm = list(sim = sim_h, fit = fit_h)
	)
}

test_that("print.dbn works for all model types", {
	skip_on_cran()
	fits = setup_fits()

	for (nm in names(fits)) {
		expect_output(print(fits[[nm]]$fit), regexp = NULL,
			label = sprintf("print works for %s", nm))
	}
})

test_that("summary works for all model types", {
	skip_on_cran()
	fits = setup_fits()

	for (nm in names(fits)) {
		res = tryCatch(
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
	fits = setup_fits()

	pdf(NULL)
	on.exit(dev.off(), add = TRUE)

	for (nm in names(fits)) {
		res = tryCatch(
			{ plot_dbn(fits[[nm]]$fit); TRUE },
			error = function(e) e$message
		)
		expect_true(isTRUE(res),
			label = sprintf("plot works for %s (error: %s)", nm,
				if (isTRUE(res)) "none" else res))
	}
})

test_that("predict returns PPD for all model types", {
	skip_on_cran()
	fits = setup_fits()

	for (nm in names(fits)) {
		ppd = tryCatch(
			posterior_predict_dbn(fits[[nm]]$fit, ndraws = 5, seed = 42),
			error = function(e) NULL
		)
		expect_true(!is.null(ppd),
			label = sprintf("PPD succeeds for %s", nm))
	}
})

test_that("dbn_lowrank_accurate runs and returns valid output", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3010)

	fit = dbn_lowrank_accurate(sim$Y, r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "lowrank")
	expect_true(!is.null(fit$U), label = "U stored")
	expect_true(!is.null(fit$alpha), label = "alpha stored")
})

test_that("tidy_dbn_lowrank returns data.frame with expected columns", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3011)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	tidy_df = tidy_dbn_lowrank(fit)
	expect_true(is.data.frame(tidy_df), label = "returns data.frame")
	expect_true(nrow(tidy_df) > 0, label = "non-empty data.frame")
	expect_true(ncol(tidy_df) >= 3, label = "has at least 3 columns")
})

test_that("debug_irf runs without error", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 3012)
	fit = dbn(sim$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	expect_no_error(
		capture.output(debug_irf(fit, draw_idx = 1, shock_i = 1, shock_j = 2))
	)
})

test_that("compare_dbn works with two fits", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim = simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3013)

	fit1 = dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)
	fit2 = dbn(sim$Y, model = "static", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	result = compare_dbn(fit1, fit2)
	expect_true(!is.null(result), label = "compare_dbn returns non-NULL")
})

test_that("stat_density, stat_in_degree, stat_out_degree work", {
	n = 6
	Theta = matrix(rnorm(n * n), n, n)
	diag(Theta) = NA

	d = stat_density(Theta)
	expect_true(is.numeric(d))
	expect_length(d, 1)

	in_deg = stat_in_degree(Theta)
	expect_true(is.numeric(in_deg))

	out_deg = stat_out_degree(Theta)
	expect_true(is.numeric(out_deg))
})

test_that("stat_reciprocity and stat_transitivity work", {
	n = 6
	Theta = matrix(rnorm(n * n), n, n)
	diag(Theta) = NA

	r = stat_reciprocity(Theta)
	expect_true(is.numeric(r))
	expect_length(r, 1)

	tr = stat_transitivity(Theta)
	expect_true(is.numeric(tr))
	expect_length(tr, 1)
})

test_that("build_shock produces valid shock matrices", {
	n = 8

	S_density = build_shock(n, type = "density")
	expect_true(is.matrix(S_density))
	expect_equal(dim(S_density), c(n, n))

	S_edge = build_shock(n, type = "unit_edge", i = 1, j = 2)
	expect_equal(S_edge[1, 2], 1)
	expect_equal(sum(S_edge != 0, na.rm = TRUE), 1)

	S_node = build_shock(n, type = "node_out", i = 1)
	expect_true(all(S_node[1, ] == 1 | is.na(S_node[1, ])))
})

test_that("compute_irf returns list with expected structure", {
	skip_on_cran()
	sim = simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3014)
	fit = dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	S = build_shock(6, type = "density")
	irf = compute_irf(fit, shock = S, H = 3, n_draws = 3,
		stat_fun = dbn::stat_density)
	expect_true(is.list(irf), label = "compute_irf returns list")
})

test_that("param_summary returns data.frame for all model types", {
	skip_on_cran()
	fits = setup_fits()

	for (nm in names(fits)) {
		ps = param_summary(fits[[nm]]$fit)
		expect_true(is.data.frame(ps),
			label = sprintf("param_summary returns df for %s", nm))
		expect_true(nrow(ps) > 0,
			label = sprintf("param_summary non-empty for %s", nm))
	}
})

test_that("theta_summary and theta_credible work", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3015)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ts = theta_summary(fit, i = 1, j = 2, rel = 1, time = 1)
	expect_true(!is.null(ts), label = "theta_summary returns non-NULL")

	tc = theta_credible(fit, i = 1, j = 2, rel = 1, time = 1)
	expect_true(is.data.frame(tc), label = "theta_credible returns df")
})

test_that("edge_prob returns matrix", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3016)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ep = edge_prob(fit)
	expect_true(is.matrix(ep) || is.array(ep), label = "edge_prob returns matrix")
})

test_that("derive_draws works", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3017)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	dd = derive_draws(fit, fun = function(draw) mean(draw$theta, na.rm = TRUE))
	expect_true(!is.null(dd), label = "derive_draws returns non-NULL")
})

test_that("network_summary works", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3018)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ns = network_summary(fit)
	expect_true(!is.null(ns), label = "network_summary returns non-NULL")
})

test_that("simulate_test_data produces valid output", {
	td = simulate_test_data(n = 5, time = 3)
	expect_true(is.array(td))
	expect_equal(dim(td)[1], 5)
	expect_equal(dim(td)[4], 3)
})

test_that("get_dbn_threads returns numeric", {
	nt = get_dbn_threads()
	expect_true(is.numeric(nt))
	expect_true(nt >= 1)
})

test_that("set_dbn_threads and get_dbn_threads round-trip", {
	old = get_dbn_threads()
	set_dbn_threads(1)
	expect_equal(get_dbn_threads(), 1)
	set_dbn_threads(old)
})

test_that("estimate_memory returns sensible estimate", {
	est = estimate_memory(n_row = 20, Tt = 30, p = 1, nscan = 5000, odens = 5)
	expect_true(is.numeric(est) || is.list(est))
})

test_that("plot_trace produces output", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim = simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 3020)
	fit = dbn(sim$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	p_trace = plot_trace(fit)
	expect_true(!is.null(p_trace), label = "plot_trace returns non-NULL")
})

test_that("plot_theta produces output with lowrank model", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 3025)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	p_theta = plot_theta(fit, time = 1, rel = 1)
	expect_true(!is.null(p_theta), label = "plot_theta returns non-NULL")
})

test_that("plot_ppc_ecdf and plot_ppc_density produce output", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim = simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3021)
	fit = dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ppd = tryCatch(
		posterior_predict_dbn(fit, ndraws = 5, seed = 42),
		error = function(e) NULL
	)

	if (!is.null(ppd)) {
		p_ecdf = plot_ppc_ecdf(fit, ppd = ppd, Y_obs = sim$Y)
		expect_true(!is.null(p_ecdf), label = "plot_ppc_ecdf returns non-NULL")

		p_dens = plot_ppc_density(fit, ppd = ppd, Y_obs = sim$Y)
		expect_true(!is.null(p_dens), label = "plot_ppc_density returns non-NULL")
	}
})

test_that("regime_probs and plot_regime_probs work for HMM", {
	skip_on_cran()

	sim = simulate_hmm_dbn(n = 6, p = 1, time = 8, R = 2, seed = 3022)
	fit = dbn(sim$Y, model = "hmm", R = 2,
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	rp = regime_probs(fit)
	expect_true(is.matrix(rp) || is.data.frame(rp), label = "regime_probs returns matrix/df")

	if (requireNamespace("ggplot2", quietly = TRUE)) {
		p = plot_regime_probs(fit)
		expect_true(!is.null(p), label = "plot_regime_probs returns non-NULL")
	}
})

test_that("dyad_path and net_snapshot produce output", {
	skip_on_cran()
	skip_if_not_installed("ggplot2")

	sim = simulate_dynamic_dbn(n = 6, p = 1, time = 5, seed = 3023)
	fit = dbn(sim$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	p_dyad = dyad_path(fit, i = 1, j = 2, rel = 1)
	expect_true(!is.null(p_dyad), label = "dyad_path returns non-NULL")

	p_net = net_snapshot(fit, t = 1, rel = 1)
	expect_true(!is.null(p_net), label = "net_snapshot returns non-NULL")
})

test_that("latent_summary works", {
	skip_on_cran()

	sim = simulate_static_dbn(n = 6, p = 1, time = 5, seed = 3024)
	fit = dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, odens = 1, verbose = FALSE)

	ls = latent_summary(fit)
	expect_true(!is.null(ls), label = "latent_summary returns non-NULL")
})
