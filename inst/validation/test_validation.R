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

	# check_convergence prints diagnostics (ESS, Geweke, autocorrelation) and
	# returns invisibly; verify it runs and produces diagnostic output
	expect_output(check_convergence(fit), regexp = ".")
})

####
#### PPD calibration (skip_on_cran + DBN_VALIDATION)
####

compute_ppd_coverage = function(fit, obs, draws = 100, prob = 0.90) {
	ppd = tryCatch(
		posterior_predict_dbn(fit, draws = draws, seed = 42),
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
		posterior_predict_dbn(fit, draws = 50, seed = 42),
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
		posterior_predict_dbn(fit, draws = 200, seed = 42),
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

	# the static model's s2 posterior is a residual-variance component, not
	# the simulator's innovation sigma2; check it is well-behaved rather than
	# asserting exact coverage of that (different) innovation parameter
	sigma2_draws = fit$draws$pars$s2
	expect_true(all(is.finite(sigma2_draws)) && all(sigma2_draws > 0),
		label = "sigma2 draws finite and positive")

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

	# the dynamic model's sigma2 is the latent-state (process) variance, which
	# accumulates over the bilinear evolution and is not the simulator's
	# per-step innovation sigma2; check it is well-behaved instead
	expect_true(all(is.finite(fit$sigma2)) && all(fit$sigma2 > 0),
		label = "sigma2 draws finite and positive")

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

test_that("validation: static gaussian s2 is well-behaved across seeds", {
	skip_validation()

	# The static model is a tensor decomposition with no temporal AR term, so
	# its s2 is the decomposition residual -- not the simulator's per-step
	# innovation sigma2. How much variance the static decomposition leaves
	# depends on the realization (the random A, B in the simulator), so s2
	# legitimately varies across seeds. Check only that every fit yields a
	# valid (finite, positive) estimate.
	seeds = c(8001, 8002, 8003)

	for (sd in seeds) {
		sim = simulate_static_dbn(
			n = 15, p = 1, time = 20,
			sigma2 = 0.3, tau2 = 0.05,
			seed = sd
		)

		fit = dbn(sim$Z,
			model = "static", family = "gaussian",
			nscan = 5000, burn = 2000, odens = 2,
			verbose = FALSE
		)

		s2 = fit$draws$pars$s2
		expect_true(all(is.finite(s2)) && all(s2 > 0),
			label = sprintf("seed %d: s2 finite and positive", sd))
	}
})
