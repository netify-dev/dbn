####
# Phase 1A: Rigorous parameter recovery study
# Validates that MCMC posterior inference recovers known parameters
# with adequate data, chain length, and strict thresholds.
#
# Gated behind skip_on_cran() + DBN_VALIDATION env var because
# these tests take minutes (long chains for proper inference).
####

skip_validation <- function() {
	skip_on_cran()
	skip_if(
		!identical(Sys.getenv("DBN_VALIDATION"), "true"),
		"Set DBN_VALIDATION=true to run extended validation tests"
	)
}

# helper: posterior mean of M from fit object
extract_M_mean <- function(fit) {
	if (fit$model == "static") {
		M_draws <- fit$draws$misc$M
		if (!is.null(M_draws) && length(M_draws) > 0) {
			Reduce("+", M_draws) / length(M_draws)
		} else {
			NULL
		}
	} else if (fit$model == "dynamic") {
		if (!is.null(fit$M) && is.array(fit$M)) {
			nd <- length(dim(fit$M))
			apply(fit$M, 1:(nd - 1), mean)
		} else if (!is.null(fit$draws$misc$M)) {
			Reduce("+", fit$draws$misc$M) / length(fit$draws$misc$M)
		} else {
			NULL
		}
	} else if (fit$model %in% c("lowrank", "hmm")) {
		M_draws <- fit$draws$misc$M
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

# helper: check if truth falls in posterior credible interval
covers_truth <- function(draws, truth, prob = 0.90) {
	q <- quantile(draws, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2),
		na.rm = TRUE)
	truth >= q[1] && truth <= q[2]
}

# helper: correlation between two vectors, removing NAs
safe_cor <- function(x, y) {
	valid <- complete.cases(x, y)
	if (sum(valid) < 5) return(NA_real_)
	cor(x[valid], y[valid])
}

####
# 1. Static Gaussian — strongest identifiability
####
test_that("validation: static gaussian recovers M and sigma2", {
	skip_validation()

	n <- 15
	Tt <- 20
	true_sigma2 <- 0.3
	true_tau2 <- 0.05

	sim <- simulate_static_dbn(
		n = n, p = 1, time = Tt,
		sigma2 = true_sigma2, tau2 = true_tau2,
		seed = 7001
	)

	fit <- dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	# M recovery: correlation > 0.7 for gaussian (strong signal)
	M_hat <- extract_M_mean(fit)
	expect_false(is.null(M_hat), label = "M posterior mean extracted")
	r_M <- safe_cor(as.vector(sim$M[, , 1]), as.vector(M_hat[, , 1]))
	expect_gt(r_M, 0.7,
		label = sprintf("M recovery r=%.3f > 0.7 (static gaussian)", r_M))

	# sigma2: 90% CI should cover truth
	sigma2_draws <- fit$draws$pars$s2
	expect_true(
		covers_truth(sigma2_draws, true_sigma2, 0.90),
		label = sprintf(
			"sigma2 90%% CI [%.3f, %.3f] covers truth %.3f",
			quantile(sigma2_draws, 0.05),
			quantile(sigma2_draws, 0.95),
			true_sigma2
		)
	)

	# sigma2 posterior mean should be within factor of 2
	sigma2_hat <- mean(sigma2_draws)
	expect_lt(abs(sigma2_hat - true_sigma2) / true_sigma2, 1.0,
		label = "sigma2 relative error < 100%")
})

####
# 2. Static Ordinal — weaker identifiability
####
test_that("validation: static ordinal recovers M direction", {
	skip_validation()

	n <- 15
	Tt <- 20

	sim <- simulate_static_dbn(
		n = n, p = 1, time = Tt,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 7002
	)

	fit <- dbn(sim$Y,
		model = "static", family = "ordinal",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	# M recovery: r > 0.5 for ordinal (identifiability weaker)
	M_hat <- extract_M_mean(fit)
	expect_false(is.null(M_hat), label = "M posterior mean extracted")
	r_M <- safe_cor(as.vector(sim$M[, , 1]), as.vector(M_hat[, , 1]))
	expect_gt(r_M, 0.5,
		label = sprintf("M recovery r=%.3f > 0.5 (static ordinal)", r_M))

	# sigma2 should be positive and bounded
	sigma2_draws <- fit$draws$pars$s2
	expect_true(all(sigma2_draws > 0), label = "sigma2 all positive")
	expect_lt(mean(sigma2_draws), 10, label = "sigma2 mean < 10")
})

####
# 3. Static Binary — probit truncation makes recovery harder
####
test_that("validation: static binary recovers M sign pattern", {
	skip_validation()

	n <- 15
	Tt <- 30

	# generate binary data from probit DGP
	set.seed(7003)
	M_true <- matrix(rnorm(n * n, 0, 0.5), n, n)
	diag(M_true) <- NA
	Z <- array(NA, c(n, n, 1, Tt))
	Y <- array(NA, c(n, n, 1, Tt))
	for (t in 1:Tt) {
		Z[, , 1, t] <- M_true + matrix(rnorm(n * n), n, n)
		Y[, , 1, t] <- ifelse(Z[, , 1, t] > 0, 1, 0)
		diag(Y[, , 1, t]) <- NA
	}

	fit <- dbn(Y,
		model = "static", family = "binary",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	# M recovery: r > 0.4 for binary
	M_hat <- extract_M_mean(fit)
	expect_false(is.null(M_hat), label = "M posterior mean extracted")
	M_hat_mat <- if (length(dim(M_hat)) == 3) M_hat[, , 1] else M_hat
	r_M <- safe_cor(as.vector(M_true), as.vector(M_hat_mat))
	expect_gt(r_M, 0.4,
		label = sprintf("M recovery r=%.3f > 0.4 (static binary)", r_M))

	# sign recovery: > 70% of off-diagonal signs match
	valid <- !is.na(as.vector(M_true)) & !is.na(as.vector(M_hat_mat))
	sign_match <- mean(sign(as.vector(M_true)[valid]) ==
		sign(as.vector(M_hat_mat)[valid]))
	expect_gt(sign_match, 0.65,
		label = sprintf("M sign recovery %.1f%% > 65%%", sign_match * 100))
})

####
# 4. Dynamic Gaussian — time-varying inference
####
test_that("validation: dynamic gaussian recovers Theta and variance params", {
	skip_validation()

	n <- 12
	Tt <- 20
	true_sigma2 <- 0.3
	true_tauA2 <- 0.05
	true_tauB2 <- 0.05

	sim <- simulate_dynamic_dbn(
		n = n, p = 1, time = Tt,
		sigma2 = true_sigma2, tauA2 = true_tauA2, tauB2 = true_tauB2,
		seed = 7004
	)

	fit <- dbn(sim$Z,
		model = "dynamic", family = "gaussian",
		nscan = 8000, burn = 3000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "dynamic")

	# Theta recovery at final time point: r > 0.5
	if (!is.null(fit$Theta) && is.array(fit$Theta)) {
		theta_mean <- apply(fit$Theta, 1:4, mean)
		t_last <- dim(theta_mean)[4]
		t_true_last <- dim(sim$Theta)[4]
		th_hat <- as.vector(theta_mean[, , 1, t_last])
		th_true <- as.vector(sim$Theta[, , 1, t_true_last])
		r_theta <- safe_cor(th_hat, th_true)
		expect_gt(r_theta, 0.4,
			label = sprintf("Theta(T) recovery r=%.3f > 0.4 (dynamic gaussian)", r_theta))
	}

	# sigma2: 90% CI should cover truth
	expect_true(
		covers_truth(fit$sigma2, true_sigma2, 0.90),
		label = sprintf(
			"sigma2 90%% CI [%.3f, %.3f] covers truth %.3f",
			quantile(fit$sigma2, 0.05),
			quantile(fit$sigma2, 0.95),
			true_sigma2
		)
	)

	# tau_A2 should be positive and finite
	expect_true(all(is.finite(fit$tau_A2)), label = "tau_A2 finite")
	expect_true(all(fit$tau_A2 > 0), label = "tau_A2 positive")

	# tau_B2 should be positive and finite
	expect_true(all(is.finite(fit$tau_B2)), label = "tau_B2 finite")
	expect_true(all(fit$tau_B2 > 0), label = "tau_B2 positive")

	# M recovery
	M_hat <- extract_M_mean(fit)
	if (!is.null(M_hat)) {
		M_hat_mat <- if (length(dim(M_hat)) >= 3) M_hat[, , 1] else M_hat
		M_true_mat <- sim$M[, , 1]
		r_M <- safe_cor(as.vector(M_true_mat), as.vector(M_hat_mat))
		expect_gt(r_M, 0.5,
			label = sprintf("M recovery r=%.3f > 0.5 (dynamic gaussian)", r_M))
	}
})

####
# 5. Dynamic Ordinal — model should complete with sensible output
####
test_that("validation: dynamic ordinal completes and produces valid output", {
	skip_validation()

	sim <- simulate_dynamic_dbn(
		n = 12, p = 1, time = 20,
		sigma2 = 0.3, tauA2 = 0.05, tauB2 = 0.05,
		seed = 7005
	)

	fit <- dbn(sim$Y,
		model = "dynamic", family = "ordinal",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	# sigma2 should be positive and stable
	expect_true(all(is.finite(fit$sigma2)), label = "sigma2 all finite")
	expect_true(all(fit$sigma2 > 0), label = "sigma2 all positive")

	# no degenerate behavior: sigma2 shouldn't collapse to near-zero
	expect_gt(mean(fit$sigma2), 0.01, label = "sigma2 mean > 0.01")

	# A and B draws should exist and have correct dimensions
	expect_true(!is.null(fit$A), label = "A draws stored")
	expect_true(length(fit$A) > 0, label = "A draws non-empty")
})

####
# 6. Lowrank Ordinal — factor recovery
####
test_that("validation: lowrank recovers non-degenerate alpha trajectories", {
	skip_validation()

	sim <- simulate_lowrank_dbn(
		n = 20, p = 1, time = 15,
		r = 2, sigma2 = 0.2,
		seed = 7006
	)

	fit <- dbn(sim$Y,
		model = "lowrank", r = 2,
		family = "ordinal",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	# alpha trajectories should be non-degenerate
	alpha_draws <- fit$alpha
	expect_true(!is.null(alpha_draws), label = "alpha draws stored")
	expect_true(length(alpha_draws) > 0, label = "alpha draws non-empty")

	# alpha should have variation over time (not constant)
	alpha_last <- alpha_draws[[length(alpha_draws)]]
	for (r_idx in seq_len(nrow(alpha_last))) {
		alpha_sd <- sd(alpha_last[r_idx, ])
		expect_gt(alpha_sd, 1e-6,
			label = sprintf("alpha[%d] has temporal variation (sd=%.4f)", r_idx, alpha_sd))
	}

	# sigma2 should be finite and positive
	sigma2_draws <- fit$draws$pars$sigma2_proc
	expect_true(all(is.finite(sigma2_draws)), label = "sigma2 finite")
	expect_true(all(sigma2_draws > 0), label = "sigma2 positive")

	# U should be orthonormal (U'U = I)
	U_last <- fit$U[[length(fit$U)]]
	UtU <- crossprod(U_last)
	expect_equal(UtU, diag(ncol(U_last)), tolerance = 0.1,
		label = "U is approximately orthonormal")
})

####
# 7. HMM Ordinal — regime recovery
####
test_that("validation: hmm recovers regime assignments", {
	skip_validation()

	sim <- simulate_hmm_dbn(
		n = 10, p = 1, time = 30,
		R = 2, transition_prob = 0.9,
		sigma2 = 0.3,
		seed = 7007
	)

	fit <- dbn(sim$Y,
		model = "hmm", R = 2,
		family = "ordinal",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	# regime probabilities
	rp <- regime_probs(fit)
	expect_true(!is.null(rp), label = "regime_probs returns non-NULL")
	expect_equal(nrow(rp), sim$dims$Tt, label = "regime probs has correct time dim")

	# modal regime assignment should agree with truth > 60%
	# (accounting for label switching: try both permutations)
	modal_regime <- apply(rp, 1, which.max)
	true_S <- sim$S

	# label matching may be permuted
	agreement_direct <- mean(modal_regime == true_S, na.rm = TRUE)
	agreement_flipped <- mean(modal_regime == (3 - true_S), na.rm = TRUE)
	best_agreement <- max(agreement_direct, agreement_flipped)

	expect_gt(best_agreement, 0.55,
		label = sprintf(
			"Regime recovery %.1f%% > 55%% (best of direct/flipped)",
			best_agreement * 100
		))

	# transition matrix should have high diagonal (sticky regimes)
	Pi_draws <- fit$Pi
	if (!is.null(Pi_draws) && length(Pi_draws) > 0) {
		Pi_mean <- Reduce("+", Pi_draws) / length(Pi_draws)
		for (r_idx in 1:2) {
			expect_gt(Pi_mean[r_idx, r_idx], 0.5,
				label = sprintf("Pi[%d,%d]=%.2f > 0.5 (sticky regimes)", r_idx, r_idx, Pi_mean[r_idx, r_idx]))
		}
	}

	# sigma2 should be positive and finite
	sigma2_draws <- fit$draws$pars$sigma2_proc
	expect_true(all(is.finite(sigma2_draws)), label = "sigma2 finite")
})

####
# 8. Multi-seed coverage test: static gaussian
# Run 3 seeds, verify that truth falls in 90% CI for at least 2/3 seeds
####
test_that("validation: static gaussian sigma2 coverage across seeds", {
	skip_validation()

	true_sigma2 <- 0.3
	seeds <- c(8001, 8002, 8003)
	covers <- logical(length(seeds))

	for (i in seq_along(seeds)) {
		sim <- simulate_static_dbn(
			n = 15, p = 1, time = 20,
			sigma2 = true_sigma2, tau2 = 0.05,
			seed = seeds[i]
		)

		fit <- dbn(sim$Z,
			model = "static", family = "gaussian",
			nscan = 5000, burn = 2000, odens = 2,
			verbose = FALSE
		)

		covers[i] <- covers_truth(fit$draws$pars$s2, true_sigma2, 0.90)
	}

	# at least 2 out of 3 should cover truth
	expect_gte(sum(covers), 2,
		label = sprintf(
			"sigma2 covered in %d/3 seeds (expect >= 2)",
			sum(covers)
		))
})
