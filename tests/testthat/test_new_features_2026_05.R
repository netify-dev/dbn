# tests for the round-19 additions: smart init, custom priors,
# long-run IRF, Cholesky-ordered FEVD, non-dynamic covariate dispatch

test_that("smart init produces a usable warm-start", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 6, p = 1, time = 12, seed = 17)
	init <- dbn:::.dbn_smart_init(sim$Z, family = "gaussian",
		model = "dynamic", symmetric = FALSE,
		Tt = 12L, n_row = 6L, n_col = 6L, p = 1L, verbose = FALSE)
	expect_true(init$source %in% c("als_warm_start", "identity_default"))
	expect_equal(dim(init$A_init), c(6L, 6L, 12L))
	expect_equal(init$A_init[, , 1], diag(6L))   # t=1 pinned at identity
	expect_true(init$sigma2_init > 0)
})

test_that("dbn_static accepts a_sig / b_sig / a_g / b_g priors", {
	skip_on_cran()
	sim <- simulate_static_dbn(n = 5, time = 6, p = 1, seed = 3)
	fit_loose <- suppressMessages(suppressWarnings(
		dbn_static(sim$Y, family = "gaussian", nscan = 60, burn = 20,
			a_sig = 1, b_sig = 1, verbose = FALSE)
	))
	fit_tight <- suppressMessages(suppressWarnings(
		dbn_static(sim$Y, family = "gaussian", nscan = 60, burn = 20,
			a_sig = 20, b_sig = 0.05, verbose = FALSE)
	))
	# tighter prior shrinks toward smaller s2; fit$params is a matrix
	# with column "s2"
	mean_loose <- mean(fit_loose$params[, "s2"])
	mean_tight <- mean(fit_tight$params[, "s2"])
	expect_true(mean_tight < mean_loose)
})

test_that("dbn_dynamic exposes a_sig / b_sig / a_sig_obs / b_sig_obs", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 5, p = 1, time = 8, seed = 1)
	# call with custom obs prior to confirm signature accepts the args
	fit <- suppressMessages(suppressWarnings(
		dbn(sim$Z, model = "dynamic", family = "gaussian",
			nscan = 80, burn = 20, odens = 2,
			a_sig = 2, b_sig = 1, a_sig_obs = 5, b_sig_obs = 0.5,
			verbose = FALSE)
	))
	expect_s3_class(fit, "dbn")
})

test_that("irf_longrun computes (I - B kron A)^{-1} vec(S) when contractive", {
	skip_on_cran()
	# build a tiny known fit by hand: static, A = 0.5 I, B = 0.5 I, so the
	# Kronecker gain is 0.5 * 0.5 = 0.25 (contractive)
	n <- 4L
	fit <- list(
		model = "static",
		dims = list(n_row = n, n_col = n, p = 1L, Tt = 4L),
		A = NULL,
		B = list(array(0.5 * diag(n), c(n, n, 5L)),
			array(0.5 * diag(n), c(n, n, 5L))),
		sigma2 = rep(1, 5)
	)
	class(fit) <- "dbn"
	S <- matrix(0, n, n); S[1, 2] <- 1
	res <- irf_longrun(fit, S, n_draws = 3)
	# closed form: (I - 0.25 I)^{-1} vec(S) = (1 / 0.75) vec(S)
	expect_true(all(abs(res$Delta_inf[1, 2] - 1 / 0.75) < 1e-6))
	expect_equal(res$n_contractive, 3L)
	expect_equal(res$gain_mean, 0.25, tolerance = 1e-10)
})

test_that("irf_longrun returns NA when every draw is non-contractive", {
	skip_on_cran()
	n <- 4L
	# A = 2 I, B = I -> gain = 2 (non-contractive)
	fit <- list(
		model = "static",
		dims = list(n_row = n, n_col = n, p = 1L, Tt = 4L),
		B = list(array(2 * diag(n), c(n, n, 3L)),
			array(diag(n), c(n, n, 3L))),
		sigma2 = rep(1, 3)
	)
	class(fit) <- "dbn"
	S <- matrix(0, n, n); S[1, 2] <- 1
	res <- suppressMessages(suppressWarnings(irf_longrun(fit, S, n_draws = 3)))
	expect_equal(res$n_contractive, 0L)
})

test_that("fevd accepts Sigma_row / Sigma_col for structural decomposition", {
	skip_on_cran()
	sim <- simulate_dynamic_dbn(n = 4, p = 1, time = 8, seed = 5)
	fit <- suppressMessages(suppressWarnings(
		dbn(sim$Z, model = "dynamic", family = "gaussian",
			nscan = 60, burn = 20, odens = 2, verbose = FALSE)
	))
	# pass an explicit Sigma_row (diagonal, but not identity) and Sigma_col
	Sr <- diag(c(1, 2, 1, 1))
	Sc <- diag(c(1, 1, 2, 1))
	fv <- fevd(fit, H = 3, shock = "actor", Sigma_row = Sr, Sigma_col = Sc,
		n_draws = 10)
	expect_s3_class(fv, "dbn_fevd")
	expect_true(fv$metadata$ordering_required)
	expect_false(fv$metadata$isotropic_innovation)
})

test_that("dbn() accepts covariates on static / piecewise via shift-trick", {
	skip_on_cran()
	set.seed(13)
	n <- 5; Tt <- 8
	X <- matrix(stats::rnorm(n * n), n, n)
	Y <- array(0, c(n, n, 1, Tt))
	for (t in seq_len(Tt))
		Y[, , 1, t] <- 0.5 * X + matrix(stats::rnorm(n * n, 0, 0.5), n, n)
	for (t in seq_len(Tt)) diag(Y[, , 1, t]) <- NA
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	# static covariate fit
	fit_static <- suppressMessages(suppressWarnings(
		dbn(Y, model = "static", family = "gaussian", covariates = cv,
			nscan = 60, burn = 20, verbose = FALSE)
	))
	expect_s3_class(fit_static, "dbn_covariates_fit")
	expect_true(!is.null(fit_static$beta))
})
