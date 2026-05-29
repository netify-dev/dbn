#### posterior_predict smoke tests
# exhaustive (model x family) coverage lives in inst/validation/.

test_that("posterior_predict_dbn produces draws for an MCMC fit", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 6, time = 8, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE)
	ppd = posterior_predict_dbn(fit, draws = 5, seed = 1)
	expect_s3_class(ppd, "dbn_ppd")
	expect_length(ppd, 5)
	expect_true(all(vapply(ppd, function(x) length(dim(x)) == 4L, logical(1))))
})

test_that("posterior_predict_dbn works for an ALS+bootstrap fit", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 6, time = 8, seed = 1)
	fit = suppressWarnings(dbn_als(sim$Z, family = "gaussian",
		bootstrap = 10, lambda = 1, verbose = FALSE))
	ppd = posterior_predict_dbn(fit, draws = 5, seed = 1)
	expect_s3_class(ppd, "dbn_ppd")
})
