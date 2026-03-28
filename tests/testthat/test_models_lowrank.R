####
# lowrank model
####

test_that("simulate_lowrank_dbn produces valid output structure", {
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 10, r = 2, seed = 42)
	expect_true(is.list(sim))
	expect_true("Y" %in% names(sim))
	expect_equal(dim(sim$Y)[1], 8)
	expect_equal(dim(sim$Y)[2], 8)
	expect_equal(dim(sim$Y)[3], 1)
	expect_equal(dim(sim$Y)[4], 10)
	# ordinal: positive integers
	y_vals = sim$Y[!is.na(sim$Y)]
	expect_true(all(y_vals == round(y_vals)))
	expect_true(all(y_vals > 0))
})

test_that("simulate_lowrank_dbn works with multiple relations", {
	sim = simulate_lowrank_dbn(n = 6, p = 2, time = 8, r = 2, seed = 43)
	expect_equal(dim(sim$Y)[3], 2)
})

test_that("dbn_lowrank fits ordinal data and returns expected object", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 8, r = 2, seed = 44)
	fit = dbn(sim$Y, model = "lowrank", family = "ordinal",
						 nscan = 100, burn = 50, verbose = FALSE, r = 2)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "lowrank")
	# required components
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$alpha))
	expect_true(!is.null(fit$B))
	expect_true(!is.null(fit$sigma2))
})

test_that("lowrank U loadings have correct dimensions", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 8, r = 2, seed = 45)
	fit = dbn(sim$Y, model = "lowrank", family = "ordinal",
						 nscan = 100, burn = 50, verbose = FALSE, r = 2)
	# U: n_row x r
	U_1 = fit$U[[1]]
	expect_equal(nrow(U_1), 8)
	expect_equal(ncol(U_1), 2)
	# semi-orthogonal: U'U ~ I
	UtU = t(U_1) %*% U_1
	expect_true(max(abs(UtU - diag(2))) < 0.5)
})

test_that("dbn_lowrank fits gaussian data", {
	skip_on_cran()
	set.seed(46)
	Y = array(rnorm(8 * 8 * 1 * 8), dim = c(8, 8, 1, 8))
	for (t in 1:8) diag(Y[,,1,t]) = NA
	fit = dbn(Y, model = "lowrank", family = "gaussian",
						 nscan = 100, burn = 50, verbose = FALSE, r = 2)
	expect_s3_class(fit, "dbn")
	expect_true(!is.null(fit$sigma2_obs))
})

test_that("dbn_lowrank fits binary data", {
	skip_on_cran()
	set.seed(47)
	Y = array(rbinom(8 * 8 * 1 * 8, 1, 0.3), dim = c(8, 8, 1, 8))
	for (t in 1:8) diag(Y[,,1,t]) = NA
	fit = dbn(Y, model = "lowrank", family = "binary",
						 nscan = 100, burn = 50, verbose = FALSE, r = 2)
	expect_s3_class(fit, "dbn")
})

test_that("posterior_predict works for lowrank model", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 6, p = 1, time = 6, r = 2, seed = 48)
	fit = dbn(sim$Y, model = "lowrank", family = "ordinal",
						 nscan = 80, burn = 40, verbose = FALSE, r = 2)
	ppd = posterior_predict_dbn(fit, ndraws = 5, seed = 1)
	expect_s3_class(ppd, "dbn_ppd")
	expect_equal(length(ppd), 5)
})

test_that("summary and print work for lowrank model", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 6, p = 1, time = 6, r = 2, seed = 49)
	fit = dbn(sim$Y, model = "lowrank", family = "ordinal",
						 nscan = 80, burn = 40, verbose = FALSE, r = 2)
	expect_output(print(fit))
	expect_output(summary(fit))
})

test_that("predict_lowrank returns forecasts", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 6, p = 1, time = 6, r = 2, seed = 55)
	fit = dbn(sim$Y, model = "lowrank", family = "ordinal",
						 nscan = 80, burn = 40, verbose = FALSE, r = 2)
	pred = predict(fit, H = 3, draws = 5)
	expect_true(!is.null(pred))
})

test_that("lowrank rejects r > n", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 5, time = 5, r = 2, seed = 56)
	expect_error(
		dbn(sim$Y, model = "lowrank", r = 10, nscan = 30, burn = 10, verbose = FALSE),
		"cannot exceed"
	)
})

test_that("lowrank rejects T < 2", {
	skip_on_cran()
	set.seed(50)
	Y = array(rnorm(8 * 8 * 1 * 1), dim = c(8, 8, 1, 1))
	expect_error(dbn(Y, model = "lowrank", family = "ordinal", r = 2),
							 "at least 2 time points")
})

test_that("tidy_dbn_lowrank returns correct dimensions", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 10, r = 2, seed = 6886)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
						 nscan = 100, burn = 50, verbose = FALSE)
	alpha_df = tidy_dbn_lowrank(fit)
	# one row per time point per factor
	expect_equal(nrow(alpha_df), 10 * 2)
	expect_true(all(c("time", "mean", "lo", "hi", "factor") %in% names(alpha_df)))
	# time values match
	expect_equal(alpha_df$time[1:10], 1:10)
})

test_that("lowrank sigma2 alias exists", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 6, p = 1, time = 6, r = 2, seed = 6886)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
						 nscan = 80, burn = 40, verbose = FALSE)
	expect_true(!is.null(fit$sigma2))
	expect_true(!is.null(fit$sigma2_proc))
	expect_equal(fit$sigma2, fit$sigma2_proc)
})

test_that("lowrank tau_B2 stays bounded and predictions are stable", {
	skip_on_cran()
	set.seed(6886)
	sim = simulate_lowrank_dbn(n = 10, p = 1, time = 10, r = 2, seed = 6886)
	fit = dbn(sim$Y, model = "lowrank", r = 2,
						 nscan = 200, burn = 100, verbose = FALSE)
	# tau_B2 should not hit the safe_rinv_gamma ceiling of 1e8
	expect_lt(max(fit$tau_B2), 1e6)
	# predictions should be finite and not explosive
	pred = predict(fit, H = 2, draws = 10, summary = "mean")
	expect_true(all(is.finite(pred[!is.na(pred)])))
	expect_lt(max(abs(pred), na.rm = TRUE), 1000)
})
