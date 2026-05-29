#### cross-sectional (T=1) handling

.mk_T1 = function(n = 8) {
	Y = array(sample(1:5, n * n, replace = TRUE), dim = c(n, n, 1, 1))
	diag(Y[, , 1, 1]) = NA
	Y
}

test_that("dbn() auto-detects T=1 and dispatches to static", {
	Y = .mk_T1()
	expect_message(
		fit <- dbn(Y, family = "ordinal",
			nscan = 30, burn = 10, verbose = FALSE),
		regexp = "static"
	)
	expect_equal(fit$model, "static")
})

test_that("time-varying models reject T=1", {
	Y = .mk_T1()
	for (mdl in c("dynamic", "hmm", "lowrank")) {
		expect_error(
			dbn(Y, model = mdl, nscan = 20, burn = 5, verbose = FALSE),
			regexp = "2 time points|T",
			info = mdl
		)
	}
})

test_that("static model fits T=1 and T>1 data uniformly", {
	skip_on_cran()
	# T=1
	Y1 = .mk_T1()
	fit1 = suppressMessages(dbn(Y1, model = "static", family = "ordinal",
		nscan = 30, burn = 10, verbose = FALSE))
	expect_s3_class(fit1, "dbn")
	# T>1
	Y5 = array(sample(1:5, 8 * 8 * 5, replace = TRUE), dim = c(8, 8, 1, 5))
	for (t in 1:5) diag(Y5[, , 1, t]) = NA
	fit5 = dbn(Y5, model = "static", family = "ordinal",
		nscan = 30, burn = 10, verbose = FALSE)
	expect_s3_class(fit5, "dbn")
})
