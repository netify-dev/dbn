#### forecast CI surface (round-15 H3 / H4 + round-16 C4 / H3 / H4)

.mk_dyn_fit = function() {
	sim = simulate_dynamic_dbn(n = 6, time = 8, p = 1L, seed = 1)
	suppressMessages(suppressWarnings(
		dbn(sim$Z, model = "dynamic", family = "gaussian",
			nscan = 30, burn = 10, verbose = FALSE)
	))
}

test_that("predict(level = N) returns dbn_forecast_ci with the right probs", {
	skip_on_cran()
	fit = .mk_dyn_fit()
	fc80 = predict(fit, H = 2, draws = 50, level = 80)
	expect_s3_class(fc80, "dbn_forecast_ci")
	expect_equal(attr(fc80, "probs"), c(0.10, 0.90))
	expect_true(!is.null(dimnames(fc80$mean)))
	expect_true(!is.null(dimnames(fc80$lower)))
	expect_true(!is.null(dimnames(fc80$upper)))
})

test_that("predict(level = ...) rejects probability-scale input and other bad values", {
	skip_on_cran()
	fit = .mk_dyn_fit()
	# probability-scale footgun
	expect_error(predict(fit, H = 2, draws = 50, level = 0.95),
		regexp = "percent|probability")
	# out-of-range
	expect_error(predict(fit, H = 2, draws = 50, level = 200),
		regexp = "level")
	# vector level: directive error
	expect_error(predict(fit, H = 2, draws = 50, level = c(80, 95)),
		regexp = "single|level")
})

test_that("predict(summary_probs = ...) honors custom probs and validates", {
	skip_on_cran()
	fit = .mk_dyn_fit()
	fc = predict(fit, H = 2, draws = 50, summary = "ci",
		summary_probs = c(0.05, 0.95))
	expect_equal(attr(fc, "probs"), c(0.05, 0.95))
	# malformed summary_probs
	expect_error(
		predict(fit, H = 2, draws = 50, summary = "ci",
			summary_probs = c(0.95, 0.05)),
		regexp = "summary_probs|lo|hi"
	)
})

test_that("predict(level=, summary_probs=) both supplied is a directive error", {
	skip_on_cran()
	fit = .mk_dyn_fit()
	expect_error(
		predict(fit, H = 2, draws = 50, level = 80,
			summary_probs = c(0.05, 0.95)),
		regexp = "both|level|summary_probs"
	)
})

test_that("draws, nsim, and seed are all validated", {
	fit = .mk_dyn_fit()
	expect_error(predict(fit, H = 2, draws = 0), regexp = "draws|positive")
	expect_error(predict(fit, H = 2, draws = NA), regexp = "draws|positive")
	expect_error(simulate(fit, nsim = 0), regexp = "nsim|positive")
	expect_error(simulate(fit, nsim = 2, seed = c(1, 2)),
		regexp = "seed|single")
	expect_error(posterior_predict_dbn(fit, draws = 0),
		regexp = "draws|positive")
})

test_that("dbn_forecast_ci summary returns a per-horizon data frame", {
	skip_on_cran()
	fit = .mk_dyn_fit()
	fc = predict(fit, H = 3, draws = 50, level = 90)
	sm = summary(fc)
	expect_s3_class(sm, "data.frame")
	expect_true(all(c("h", "mean", "lower", "upper") %in% names(sm)))
	expect_equal(nrow(sm), 3L)
})

test_that("plot / summary methods reject non-dbn_rank_probs input", {
	expect_error(plot.dbn_rank_probs(matrix(0, 3, 3)),
		regexp = "dbn_rank_probs")
	expect_error(summary.dbn_rank_probs(list()),
		regexp = "dbn_rank_probs")
})
