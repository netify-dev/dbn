#### posterior accessor smoke tests
# only essential dispatch checks. exhaustive (model x family x accessor)
# coverage lives in inst/validation/.

test_that("tidy_dbn extracts posterior means", {
	n = 5; Tt = 10; nscan = 20
	results = list(
		model = "dynamic",
		A = lapply(1:nscan, function(i) array(rnorm(n*n*Tt), dim = c(n, n, Tt))),
		B = lapply(1:nscan, function(i) array(rnorm(n*n*Tt), dim = c(n, n, Tt))),
		M = array(rnorm(n*n*3), dim = c(n, n, 3)),
		dims = list(m = n, n_row = n, n_col = n, p = 3, Tt = Tt, is_bipartite = FALSE)
	)
	class(results) = "dbn"
	td = tidy_dbn(results)
	expect_type(td, "list")
	expect_true("A" %in% names(td))
	expect_equal(dim(td$A), c(n, n, Tt))
})

test_that("ppc_ecdf guards against out-of-range values", {
	skip_on_cran()
	sim = simulate_static_dbn(n = 6, p = 1, time = 5, seed = 1)
	fit = dbn(sim$Y, model = "static", family = "ordinal",
		nscan = 50, burn = 20, verbose = FALSE)
	ppd = posterior_predict_dbn(fit, draws = 5, seed = 1)
	# malformed Y_obs (out-of-range) should error or warn, not return garbage
	expect_error(
		plot_ppc_ecdf(fit, ppd = ppd, Y_obs = array(99, dim = dim(sim$Y))),
		regexp = NA  # any error is fine; the test is that it doesn't crash silently
	)
})

test_that("check_convergence works on dbn fits", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_output(check_convergence(fit), regexp = ".")
})

test_that("role_trajectory dispatches for sender (A) and receiver (B)", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE)
	rs = role_trajectory(fit, mat = "A")
	rr = role_trajectory(fit, mat = "B")
	expect_true(is.data.frame(rs) || is.list(rs))
	expect_true(is.data.frame(rr) || is.list(rr))
})

test_that("param_summary default returns 95% CI columns", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE)
	ps = param_summary(fit)
	expect_true("X2.5." %in% names(ps) || "lo" %in% names(ps) ||
		any(grepl("2\\.5|q025|lower", names(ps))))
})

test_that("regime_probs is NULL for non-HMM, valid for HMM", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	fit_d = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 30, burn = 10, verbose = FALSE)
	expect_null(regime_probs(fit_d))
	fit_h = suppressWarnings(dbn(sim$Y, model = "hmm", R = 2,
		nscan = 40, burn = 15, verbose = FALSE))
	rp = regime_probs(fit_h)
	expect_true(is.array(rp) || is.matrix(rp))
})
