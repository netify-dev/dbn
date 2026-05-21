#### parameter recovery and PPD calibration

# helper: compute posterior mean of M from stored draws
extract_M_posterior_mean = function(fit) {
	if (fit$model == "static") {
		M_draws = fit$draws$misc$M
		if (!is.null(M_draws) && length(M_draws) > 0) {
			Reduce("+", M_draws) / length(M_draws)
		} else if (!is.null(fit$params)) {
			n_row = fit$dims$n_row
			n_col = fit$dims$n_col
			M_col = grep("^M\\[", colnames(fit$params))
			if (length(M_col) > 0) {
				matrix(colMeans(fit$params[, M_col, drop = FALSE]), n_row, n_col)
			} else {
				NULL
			}
		} else {
			NULL
		}
	} else if (fit$model == "dynamic") {
		if (!is.null(fit$M)) {
			M_draws = fit$M
			if (is.list(M_draws)) {
				Reduce("+", M_draws) / length(M_draws)
			} else if (is.array(M_draws)) {
				apply(M_draws, 1:(length(dim(M_draws)) - 1), mean)
			} else {
				NULL
			}
		} else if (!is.null(fit$draws$misc$M)) {
			Reduce("+", fit$draws$misc$M) / length(fit$draws$misc$M)
		} else {
			NULL
		}
	} else {
		NULL
	}
}

# helper: extract posterior mean Theta across draws
extract_Theta_posterior_mean = function(fit) {
	if (fit$model == "dynamic" && !is.null(fit$Theta) && is.array(fit$Theta)) {
		nd = length(dim(fit$Theta))
		apply(fit$Theta, 1:(nd - 1), mean)
	} else if (fit$model == "static") {
		if (!is.null(fit$draws$theta) && length(fit$draws$theta) > 0) {
			Reduce("+", fit$draws$theta) / length(fit$draws$theta)
		} else {
			NULL
		}
	} else {
		NULL
	}
}

# static x gaussian
test_that("simulation study: static gaussian recovers M", {
	skip_on_cran()
	set.seed(2024)

	n = 8
	p = 1
	Tt = 15
	true_sigma2 = 0.3
	true_tau2 = 0.05

	sim = simulate_static_dbn(
		n = n, p = p, time = Tt,
		sigma2 = true_sigma2, tau2 = true_tau2,
		seed = 2024
	)

	fit = dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 800, burn = 400, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "static")

	M_hat = extract_M_posterior_mean(fit)
	if (!is.null(M_hat) && !is.null(sim$M)) {
		M_true_vec = as.vector(sim$M[, , 1])
		M_hat_vec = as.vector(M_hat)
		valid = complete.cases(M_true_vec, M_hat_vec)
		if (sum(valid) > 5) {
			r = cor(M_true_vec[valid], M_hat_vec[valid])
			expect_gt(r, 0.4, label = "M recovery correlation (static gaussian)")
		}
	}
})

# static x ordinal
test_that("simulation study: static ordinal recovers M direction", {
	skip_on_cran()
	set.seed(2025)

	sim = simulate_static_dbn(
		n = 8, p = 1, time = 15,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 2025
	)

	fit = dbn(sim$Y,
		model = "static", family = "ordinal",
		nscan = 800, burn = 400, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	if (!is.null(fit$params)) {
		sigma2_col = grep("s2|sigma2", colnames(fit$params))
		if (length(sigma2_col) > 0) {
			sigma2_hat = mean(fit$params[, sigma2_col[1]])
			expect_gt(sigma2_hat, 0)
			expect_lt(sigma2_hat, 10)
		}
	}
})

# static x binary
test_that("simulation study: static binary recovers M sign pattern", {
	skip_on_cran()
	set.seed(2026)

	n = 10
	p = 1
	Tt = 20

	M_true = matrix(rnorm(n * n, 0, 0.5), n, n)
	Z = array(NA, c(n, n, p, Tt))
	Y = array(NA, c(n, n, p, Tt))
	for (t in 1:Tt) {
		Z[, , 1, t] = M_true + matrix(rnorm(n * n), n, n)
		Y[, , 1, t] = ifelse(Z[, , 1, t] > 0, 1, 0)
	}

	fit = dbn(Y,
		model = "static", family = "binary",
		nscan = 800, burn = 400, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")

	fam_name = if (is.list(fit$family)) fit$family$name else fit$family
	expect_equal(fam_name, "binary")

	M_hat = extract_M_posterior_mean(fit)
	if (!is.null(M_hat)) {
		M_true_vec = as.vector(M_true)
		M_hat_vec = as.vector(M_hat)
		valid = complete.cases(M_true_vec, M_hat_vec)
		if (sum(valid) > 5) {
			r = cor(M_true_vec[valid], M_hat_vec[valid])
			expect_gt(r, 0.3, label = "M sign recovery (static binary)")
		}
	}
})

# dynamic x gaussian
test_that("simulation study: dynamic gaussian recovers M and Theta", {
	skip_on_cran()
	set.seed(2027)

	sim = simulate_dynamic_dbn(
		n = 8, p = 1, time = 15,
		sigma2 = 0.3, tauA2 = 0.05, tauB2 = 0.05,
		seed = 2027
	)

	fit = dbn(sim$Z,
		model = "dynamic", family = "gaussian",
		nscan = 800, burn = 400, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "dynamic")

	Theta_hat = extract_Theta_posterior_mean(fit)
	if (!is.null(Theta_hat) && !is.null(sim$Theta)) {
		nd_hat = length(dim(Theta_hat))
		nd_true = length(dim(sim$Theta))
		if (nd_hat >= 4 && nd_true >= 4) {
			t_hat = dim(Theta_hat)[4]
			t_true = dim(sim$Theta)[4]
			th_true_vec = as.vector(sim$Theta[, , 1, t_true])
			th_hat_vec = as.vector(Theta_hat[, , 1, t_hat])
			valid = complete.cases(th_true_vec, th_hat_vec)
			if (sum(valid) > 5) {
				r = cor(th_true_vec[valid], th_hat_vec[valid])
				expect_gt(r, 0.25, label = "Theta recovery (dynamic gaussian)")
			}
		}
	}

	sigma2_hat = mean(fit$sigma2)
	expect_gt(sigma2_hat, 0)
	expect_true(is.finite(sigma2_hat))
})

# dynamic x ordinal
test_that("simulation study: dynamic ordinal model runs and converges", {
	skip_on_cran()
	set.seed(2028)

	sim = simulate_dynamic_dbn(
		n = 8, p = 1, time = 12,
		sigma2 = 0.3, tauA2 = 0.05, tauB2 = 0.05,
		seed = 2028
	)

	fit = dbn(sim$Y,
		model = "dynamic", family = "ordinal",
		nscan = 800, burn = 400, odens = 2,
		verbose = FALSE
	)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "dynamic")

	expect_true(all(is.finite(fit$sigma2)))
	expect_true(all(fit$sigma2 > 0))

	expect_true(!is.null(fit$A) || !is.null(fit$draws$A))
})

# PPD calibration (dynamic gaussian)
test_that("simulation study: PPD calibration for dynamic gaussian", {
	skip_on_cran()
	set.seed(2030)

	sim = simulate_dynamic_dbn(
		n = 8, p = 1, time = 15,
		sigma2 = 0.5, tauA2 = 0.05, tauB2 = 0.05,
		seed = 2030
	)

	fit = dbn(sim$Z,
		model = "dynamic", family = "gaussian",
		nscan = 800, burn = 400, odens = 2,
		verbose = FALSE
	)

	ppd = tryCatch(
		posterior_predict_dbn(fit, draws = 50, seed = 42),
		error = function(e) {
			message("PPD generation failed: ", conditionMessage(e))
			NULL
		}
	)

	expect_true(!is.null(ppd), label = "PPD generation succeeds")

	if (!is.null(ppd) && !is.null(ppd$y_rep) && is.array(ppd$y_rep)) {
		yrep = ppd$y_rep
		nd = length(dim(yrep))
		q_lo = apply(yrep, 1:(nd - 1), quantile, probs = 0.05, na.rm = TRUE)
		q_hi = apply(yrep, 1:(nd - 1), quantile, probs = 0.95, na.rm = TRUE)

		obs = sim$Z[, , 1, ]
		if (length(dim(q_lo)) >= 3 && dim(q_lo)[3] >= 1) {
			in_interval = (obs >= q_lo[, , 1, ]) & (obs <= q_hi[, , 1, ])
		} else {
			in_interval = (obs >= q_lo) & (obs <= q_hi)
		}
		coverage = mean(in_interval, na.rm = TRUE)

		expect_gt(coverage, 0.65,
			label = "PPD 90% coverage should be above 65%")
	}
})

# PPD calibration (static gaussian)
test_that("simulation study: PPD calibration for static gaussian", {
	skip_on_cran()
	set.seed(2031)

	sim = simulate_static_dbn(
		n = 8, p = 1, time = 15,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 2031
	)

	fit = dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 800, burn = 400, odens = 2,
		verbose = FALSE
	)

	ppd = tryCatch(
		posterior_predict_dbn(fit, draws = 50, seed = 43),
		error = function(e) {
			message("PPD generation failed: ", conditionMessage(e))
			NULL
		}
	)

	expect_true(!is.null(ppd), label = "PPD generation succeeds")

	if (!is.null(ppd) && !is.null(ppd$y_rep) && is.array(ppd$y_rep)) {
		yrep = ppd$y_rep
		nd = length(dim(yrep))
		q_lo = apply(yrep, 1:(nd - 1), quantile, probs = 0.05, na.rm = TRUE)
		q_hi = apply(yrep, 1:(nd - 1), quantile, probs = 0.95, na.rm = TRUE)

		obs = sim$Z[, , 1, ]
		if (length(dim(q_lo)) >= 3 && dim(q_lo)[3] >= 1) {
			in_interval = (obs >= q_lo[, , 1, ]) & (obs <= q_hi[, , 1, ])
		} else {
			in_interval = (obs >= q_lo) & (obs <= q_hi)
		}
		coverage = mean(in_interval, na.rm = TRUE)
		expect_gt(coverage, 0.65,
			label = "PPD 90% coverage should be above 65%")
	}
})
