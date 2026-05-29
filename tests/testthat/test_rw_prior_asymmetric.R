#### random-walk prior on the asymmetric (directed) dynamic path

test_that("RW prior produces finite, scale-stable draws on directed data", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian", symmetric = FALSE,
		nscan = 60, burn = 20, verbose = FALSE)
	expect_true(all(vapply(fit$A, function(a) all(is.finite(a)), logical(1))))
	expect_true(all(vapply(fit$B, function(b) all(is.finite(b)), logical(1))))
})

test_that("iid prior alternative is still selectable via options", {
	skip_on_cran()
	prev = getOption("dbn.prior_kind", NULL)
	on.exit(options(dbn.prior_kind = prev), add = TRUE)
	options(dbn.prior_kind = "iid")
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 30, burn = 10, verbose = FALSE)
	expect_s3_class(fit, "dbn")
})
