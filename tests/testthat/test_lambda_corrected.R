# bias-corrected eb pilot produces a finite lambda
test_that("eb_corrected pilot returns valid diagnostics", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(
		n = 8,
		p = 1,
		time = 15,
		sigma2 = 0.3,
		tauA2 = 0.05,
		tauB2 = 0.05,
		seed = 42
	)
	fit = suppressWarnings(suppressMessages(
		dbn_als(sim$Z, family = "gaussian", lambda = "eb_corrected", verbose = FALSE)
	))
	p = fit$meta$tv$pilot
	expect_true(is.list(p))
	expect_true(all(c(
		"lambda_raw", "lambda_corrected", "lambda_A", "lambda_B",
		"sigma2", "tau_A2_raw", "tau_A2_corrected",
		"tau_B2_raw", "tau_B2_corrected",
		"tau_A2_meas_trace", "tau_B2_meas_trace",
		"correction_factor", "boundary_hit"
	) %in% names(p)))
	expect_true(is.finite(p$lambda_corrected) && p$lambda_corrected > 0)
	expect_true(p$tau_A2_corrected > 0)
	expect_true(p$tau_B2_corrected > 0)
})

# correction factor is bounded so a degenerate pilot does not blow up lambda
test_that("eb_corrected correction factor is bounded", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(
		n = 6,
		p = 1,
		time = 12,
		sigma2 = 0.3,
		tauA2 = 0.05,
		tauB2 = 0.05,
		seed = 1
	)
	fit = suppressWarnings(suppressMessages(
		dbn_als(sim$Z, family = "gaussian", lambda = "eb_corrected", verbose = FALSE)
	))
	p = fit$meta$tv$pilot
	# factor-3 safety cap means correction_factor in [1/3, 3] when tau_raw > 0
	expect_lte(p$correction_factor, 3.05)
	expect_gte(p$correction_factor, 1 / 3.05)
})

# corrected tau_A2 is closer to truth than raw
test_that("eb_corrected reduces tau_A2 overestimate", {
	skip_on_cran()
	true_tau_A2 = 0.05
	sim = simulate_dynamic_dbn(
		n = 8,
		p = 1,
		time = 15,
		sigma2 = 0.3,
		tauA2 = true_tau_A2,
		tauB2 = true_tau_A2,
		seed = 42
	)
	fit = suppressWarnings(suppressMessages(
		dbn_als(sim$Z, family = "gaussian", lambda = "eb_corrected", verbose = FALSE)
	))
	p = fit$meta$tv$pilot
	gap_raw = abs(p$tau_A2_raw - true_tau_A2)
	gap_corr = abs(p$tau_A2_corrected - true_tau_A2)
	expect_lt(gap_corr, gap_raw)
})

# eb_raw is no longer accepted
test_that("eb_raw lambda is rejected", {
	expect_error(
		dbn_als(
			array(rnorm(5 * 5 * 1 * 6), dim = c(5, 5, 1, 6)),
			family = "gaussian",
			lambda = "eb_raw",
			verbose = FALSE
		),
		regexp = "not recognized|Allowed"
	)
})

# cv_local picks one of three candidates around the corrected pilot
test_that("cv_local refines around the corrected pilot", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(
		n = 6,
		p = 1,
		time = 12,
		sigma2 = 0.3,
		tauA2 = 0.05,
		tauB2 = 0.05,
		seed = 7
	)
	fit = suppressWarnings(suppressMessages(
		dbn_als(sim$Z, family = "gaussian", lambda = "cv_local", verbose = FALSE)
	))
	expect_true(is.finite(fit$meta$tv$lambda_used))
	expect_true(fit$meta$tv$lambda_used > 0)
})
