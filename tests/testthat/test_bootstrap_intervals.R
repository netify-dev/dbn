# bootstrap CI: percentile vs basic (reverse-percentile) intervals
# verify both methods are stored on dbn_boot and accessible via dbn_boot_intervals
test_that("dbn_als bootstrap stores both percentile and basic intervals", {
	skip_on_cran()
	set.seed(31L)
	sim = simulate_static_dbn(n = 5, time = 8, p = 1, seed = 17)
	fit = suppressMessages(suppressWarnings(
		dbn_als(sim$Z, family = "gaussian", bootstrap = 30L, verbose = FALSE)
	))
	bt = fit$bootstrap
	expect_s3_class(bt, "dbn_boot")
	# both interval methods should be stored
	expect_true(!is.null(bt$ci_A_lo))
	expect_true(!is.null(bt$basic_A_lo))
	expect_true(!is.null(bt$basic_A_hi))
	# basic interval invariants:
	# basic_lo = 2*theta_hat - percentile_hi
	# basic_hi = 2*theta_hat - percentile_lo
	expect_equal(as.numeric(bt$basic_A_lo),
		as.numeric(2 * bt$point_est_A - bt$ci_A_hi), tolerance = 1e-10)
	expect_equal(as.numeric(bt$basic_A_hi),
		as.numeric(2 * bt$point_est_A - bt$ci_A_lo), tolerance = 1e-10)
	expect_equal(as.numeric(bt$basic_B_lo),
		as.numeric(2 * bt$point_est_B - bt$ci_B_hi), tolerance = 1e-10)
	expect_equal(as.numeric(bt$basic_M_lo),
		as.numeric(2 * bt$point_est_M - bt$ci_M_hi), tolerance = 1e-10)
})

test_that("dbn_boot_intervals returns both methods by default", {
	skip_on_cran()
	set.seed(7L)
	sim = simulate_static_dbn(n = 4, time = 6, p = 1, seed = 3)
	fit = suppressMessages(suppressWarnings(
		dbn_als(sim$Z, family = "gaussian", bootstrap = 25L, verbose = FALSE)
	))
	# can call on the fit OR on the bootstrap component directly
	res = dbn_boot_intervals(fit, parm = "A")
	expect_true("percentile" %in% names(res))
	expect_true("basic" %in% names(res))
	expect_equal(dim(res$percentile$lo), c(4L, 4L))
	expect_equal(dim(res$basic$lo), c(4L, 4L))
	# method = "basic" only returns basic
	res_basic = dbn_boot_intervals(fit, parm = "B", method = "basic")
	expect_true("basic" %in% names(res_basic))
	expect_false("percentile" %in% names(res_basic))
})

test_that("dbn_boot_intervals errors when no bootstrap is stored", {
	# fit without bootstrap should have NULL $bootstrap
	stub <- list(bootstrap = NULL)
	class(stub) <- c("dbn", "list")
	expect_error(dbn_boot_intervals(stub), "dbn_boot")
})

test_that("basic and percentile have similar widths for symmetric bootstrap distributions", {
	# fabricate a dbn_boot with a symmetric centered distribution
	skip_on_cran()
	set.seed(99L)
	sim = simulate_static_dbn(n = 4, time = 8, p = 1, seed = 13)
	fit = suppressMessages(suppressWarnings(
		dbn_als(sim$Z, family = "gaussian", bootstrap = 50L, verbose = FALSE)
	))
	bt = fit$bootstrap
	# widths
	w_pct <- bt$ci_A_hi - bt$ci_A_lo
	w_basic <- bt$basic_A_hi - bt$basic_A_lo
	# the two widths are always equal: basic widens = (2 theta - q_lo) - (2 theta - q_hi) = q_hi - q_lo
	expect_equal(as.numeric(w_pct), as.numeric(w_basic), tolerance = 1e-10)
})
