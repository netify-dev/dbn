# shrink_rho adds a contractivity ridge that keeps the operator gain bounded.
# these tests check it is validated, recorded, reduces the gain, and that the
# off path (shrink_rho = NULL) is unchanged.

sr = function(M) max(abs(eigen(M, only.values = TRUE)$values))

make_clean = function(n = 10, Tt = 30, sigma2 = 0.3, seed = 7) {
	set.seed(seed)
	A0 = matrix(rnorm(n * n), n, n); A0 = 0.5 * A0 / sr(A0)
	B0 = matrix(rnorm(n * n), n, n); B0 = 0.5 * B0 / sr(B0)
	Theta = array(0, c(n, n, 1, Tt))
	Theta[, , 1, 1] = matrix(rnorm(n * n, 0, sqrt(sigma2)), n, n)
	for (t in 2:Tt)
		Theta[, , 1, t] = A0 %*% Theta[, , 1, t - 1] %*% t(B0) +
			matrix(rnorm(n * n, 0, sqrt(sigma2)), n, n)
	# self-ties are not modelled; blank the diagonal so dbn() does not warn
	for (t in seq_len(Tt)) diag(Theta[, , 1, t]) = NA
	Theta
}

test_that("shrink_rho rejects values outside (0, 1]", {
	Y = make_clean(n = 6, Tt = 10)
	expect_error(
		dbn(Y, model = "dynamic", family = "gaussian", nscan = 20, burn = 10,
				verbose = FALSE, shrink_rho = 0),
		"shrink_rho"
	)
	expect_error(
		dbn(Y, model = "dynamic", family = "gaussian", nscan = 20, burn = 10,
				verbose = FALSE, shrink_rho = 1.5),
		"shrink_rho"
	)
	expect_error(
		dbn(Y, model = "dynamic", family = "gaussian", nscan = 20, burn = 10,
				verbose = FALSE, shrink_rho = c(0.5, 0.6)),
		"shrink_rho"
	)
})

test_that("shrink_rho defaults to 0.9, is recorded, and NULL disables it", {
	Y = make_clean(n = 6, Tt = 10)
	fit_default = suppressWarnings(dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 60, burn = 30, odens = 2, verbose = FALSE))
	expect_equal(fit_default$settings$shrink_rho, 0.9)

	fit_custom = suppressWarnings(dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 60, burn = 30, odens = 2, verbose = FALSE, shrink_rho = 0.7))
	expect_equal(fit_custom$settings$shrink_rho, 0.7)

	fit_off = suppressWarnings(dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 60, burn = 30, odens = 2, verbose = FALSE, shrink_rho = NULL))
	expect_null(fit_off$settings$shrink_rho)
})

test_that("the default contractivity prior reduces gain vs. shrink_rho = NULL", {
	Y = make_clean(n = 10, Tt = 30)
	args = list(Y, model = "dynamic", family = "gaussian",
		nscan = 800, burn = 400, odens = 2, verbose = FALSE)
	fit0 = suppressWarnings(do.call(dbn, c(args, list(shrink_rho = NULL)))) # off
	fitS = suppressWarnings(do.call(dbn, args)) # default 0.9

	g0 = dbn_operator(fit0)$stability$gain_mean
	gS = dbn_operator(fitS)$stability$gain_mean

	# the default ridge must pull the mean gain down substantially
	expect_lt(mean(gS, na.rm = TRUE), mean(g0, na.rm = TRUE))
	# and into the contractive region on this gain-0.25 DGP
	expect_lt(mean(gS, na.rm = TRUE), 0.8)
})

test_that("default shrink_rho matches explicit 0.9; NULL (off) differs", {
	Y = make_clean(n = 6, Tt = 10)
	Abar = function(f) Reduce(`+`, f$A) / length(f$A)
	fit_default = suppressWarnings(dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 80, burn = 40, odens = 2, verbose = FALSE, seed = 99))
	fit_explicit = suppressWarnings(dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 80, burn = 40, odens = 2, verbose = FALSE, seed = 99, shrink_rho = 0.9))
	fit_off = suppressWarnings(dbn(Y, model = "dynamic", family = "gaussian",
		nscan = 80, burn = 40, odens = 2, verbose = FALSE, seed = 99, shrink_rho = NULL))
	# not passing shrink_rho is identical to passing the default value
	expect_equal(Abar(fit_default), Abar(fit_explicit))
	# disabling the prior produces a genuinely different operator
	expect_false(isTRUE(all.equal(Abar(fit_default), Abar(fit_off))))
})
