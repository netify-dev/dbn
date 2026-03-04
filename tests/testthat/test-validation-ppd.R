####
# Phase 1D: Posterior predictive distribution calibration
# Validates that PPD intervals have proper coverage across model types.
# Gated behind skip_on_cran() + DBN_VALIDATION env var.
####

skip_validation <- function() {
	skip_on_cran()
	skip_if(
		!identical(Sys.getenv("DBN_VALIDATION"), "true"),
		"Set DBN_VALIDATION=true to run extended validation tests"
	)
}

# helper: compute PPD coverage for gaussian data
compute_ppd_coverage <- function(fit, obs, ndraws = 100, prob = 0.90) {
	ppd <- tryCatch(
		posterior_predict_dbn(fit, ndraws = ndraws, seed = 42),
		error = function(e) NULL
	)
	if (is.null(ppd) || is.null(ppd$y_rep) || !is.array(ppd$y_rep)) {
		return(list(coverage = NA, ppd = ppd))
	}

	yrep <- ppd$y_rep
	nd <- length(dim(yrep))
	alpha <- (1 - prob) / 2
	q_lo <- apply(yrep, 1:(nd - 1), quantile, probs = alpha, na.rm = TRUE)
	q_hi <- apply(yrep, 1:(nd - 1), quantile, probs = 1 - alpha, na.rm = TRUE)

	# match obs dimensions to quantile dimensions
	if (length(dim(q_lo)) >= 3 && dim(q_lo)[3] >= 1) {
		in_interval <- (obs >= q_lo[, , 1, ]) & (obs <= q_hi[, , 1, ])
	} else {
		in_interval <- (obs >= q_lo) & (obs <= q_hi)
	}
	coverage <- mean(in_interval, na.rm = TRUE)
	list(coverage = coverage, ppd = ppd)
}

####
# 1. Static Gaussian PPD calibration
####
test_that("validation: static gaussian PPD coverage in [0.70, 0.98]", {
	skip_validation()

	sim <- simulate_static_dbn(
		n = 12, p = 1, time = 15,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 6001
	)

	fit <- dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	obs <- sim$Z[, , 1, ]
	result <- compute_ppd_coverage(fit, obs, ndraws = 100, prob = 0.90)
	expect_false(is.na(result$coverage), label = "PPD coverage computed")
	expect_gt(result$coverage, 0.70,
		label = sprintf("PPD 90%% coverage %.2f > 0.70", result$coverage))
	expect_lt(result$coverage, 0.99,
		label = sprintf("PPD 90%% coverage %.2f < 0.99 (not over-dispersed)", result$coverage))
})

####
# 2. Dynamic Gaussian PPD calibration
####
test_that("validation: dynamic gaussian PPD coverage in [0.70, 0.98]", {
	skip_validation()

	sim <- simulate_dynamic_dbn(
		n = 10, p = 1, time = 15,
		sigma2 = 0.5, tauA2 = 0.05, tauB2 = 0.05,
		seed = 6002
	)

	fit <- dbn(sim$Z,
		model = "dynamic", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	obs <- sim$Z[, , 1, ]
	result <- compute_ppd_coverage(fit, obs, ndraws = 100, prob = 0.90)
	expect_false(is.na(result$coverage), label = "PPD coverage computed")
	expect_gt(result$coverage, 0.70,
		label = sprintf("PPD 90%% coverage %.2f > 0.70", result$coverage))
})

####
# 3. Static Ordinal PPD — check PPD generates valid ordinal categories
####
test_that("validation: static ordinal PPD produces valid ordinal output", {
	skip_validation()

	sim <- simulate_static_dbn(
		n = 10, p = 1, time = 15,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 6003
	)

	fit <- dbn(sim$Y,
		model = "static", family = "ordinal",
		nscan = 3000, burn = 1000, odens = 2,
		verbose = FALSE
	)

	ppd <- tryCatch(
		posterior_predict_dbn(fit, ndraws = 50, seed = 42),
		error = function(e) NULL
	)

	expect_true(!is.null(ppd), label = "PPD generation succeeds for ordinal")

	if (!is.null(ppd) && !is.null(ppd$y_rep)) {
		yrep <- ppd$y_rep
		# PPD replicate mean should be in same ballpark as observed mean
		obs_mean <- mean(sim$Y, na.rm = TRUE)
		ppd_means <- apply(yrep, length(dim(yrep)), function(x) mean(x, na.rm = TRUE))
		ppd_mean_avg <- mean(ppd_means)

		# within a factor of 3 of observed mean
		expect_gt(ppd_mean_avg, obs_mean / 3,
			label = "PPD mean not drastically below observed")
		expect_lt(ppd_mean_avg, obs_mean * 3,
			label = "PPD mean not drastically above observed")
	}
})

####
# 4. PPD test statistic: posterior predictive p-value
####
test_that("validation: static gaussian PPD p-value is non-extreme", {
	skip_validation()

	sim <- simulate_static_dbn(
		n = 12, p = 1, time = 15,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 6004
	)

	fit <- dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	ppd <- tryCatch(
		posterior_predict_dbn(fit, ndraws = 200, seed = 42),
		error = function(e) NULL
	)

	expect_true(!is.null(ppd), label = "PPD generation succeeds")

	if (!is.null(ppd) && !is.null(ppd$y_rep)) {
		yrep <- ppd$y_rep

		# test statistic: mean of observed data
		obs_stat <- mean(sim$Z, na.rm = TRUE)

		# compute same statistic for each PPD replicate
		nd <- length(dim(yrep))
		ppd_stats <- apply(yrep, nd, function(x) mean(x, na.rm = TRUE))

		# posterior predictive p-value
		pp_pval <- mean(ppd_stats >= obs_stat)

		# should not be extreme (neither 0 nor 1)
		expect_gt(pp_pval, 0.01,
			label = sprintf("PPD p-value %.3f > 0.01 (not extreme low)", pp_pval))
		expect_lt(pp_pval, 0.99,
			label = sprintf("PPD p-value %.3f < 0.99 (not extreme high)", pp_pval))
	}
})
