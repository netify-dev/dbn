#### numerical stability

test_that("safe_rinv_gamma handles extreme values", {
	set.seed(123)
	x <- safe_rinv_gamma(shape = 1e-20, rate = 1e-20)
	expect_true(is.finite(x))
	expect_true(x > 0)
	expect_true(x >= 1e-8)
	expect_true(x <= 1e8)

	y <- safe_rinv_gamma(shape = 1e20, rate = 1e20)
	expect_true(is.finite(y))
	expect_true(y > 0)

	z <- safe_rinv_gamma(shape = 1e-10, rate = 1e10)
	expect_true(is.finite(z))
	expect_true(z > 0)

	set.seed(456)
	shape <- 5
	rate <- 10
	safe_sample <- safe_rinv_gamma(shape, rate)

	set.seed(456)
	gamma_sample <- rgamma(1, shape = shape, rate = 1)
	direct_sample <- rate / gamma_sample

	expect_equal(safe_sample, direct_sample)
})

test_that("sparse matrix helpers work correctly", {
	A <- matrix(rnorm(100), 10, 10)
	B <- matrix(rnorm(100), 10, 10)

	result1 <- sparse_mult(A, B)
	expect_equal(result1, A %*% B)

	if (requireNamespace("Matrix", quietly = TRUE)) {
		A_sparse <- Matrix::Matrix(A, sparse = TRUE)
		result2 <- sparse_mult(A_sparse, B)
		expect_equal(as.matrix(result2), A %*% B)
	}

	diag_small <- sparse_diag(5)
	expect_equal(diag_small, diag(5))

	diag_large <- sparse_diag(150)
	if (requireNamespace("Matrix", quietly = TRUE)) {
		expect_true(inherits(diag_large, "ddiMatrix") || inherits(diag_large, "diagonalMatrix"))
		expect_equal(dim(diag_large), c(150, 150))
	}

	vals <- c(1, 2, 3, 4, 5)
	diag_vals <- sparse_diag(vals)
	expect_equal(diag_vals, diag(vals))
})

test_that("models use safe_rinv_gamma instead of 1/rgamma", {
	set.seed(789)
	m <- 5
	p <- 1
	n <- 5
	Y <- array(rnorm(m*m*p*n), dim = c(m, m, p, n))

	suppressWarnings(suppressMessages({
		result <- dbn_static(Y, nscan = 10, burn = 0, odens = 1, verbose = FALSE)
	}))

	expect_true(all(is.finite(result$params[, "s2"])))
	expect_true(all(result$params[, "s2"] > 0))
	expect_true(all(is.finite(result$params[, "t2"])))
	expect_true(all(result$params[, "t2"] > 0))
	expect_true(all(is.finite(result$params[, "g2"])))
	expect_true(all(result$params[, "g2"] > 0))

	suppressWarnings(suppressMessages({
		result_dyn <- dbn_dynamic(Y, nscan = 10, burn = 0, odens = 1, verbose = FALSE)
	}))

	expect_true(all(is.finite(result_dyn$sigma2)))
	expect_true(all(result_dyn$sigma2 > 0))
	expect_true(all(is.finite(result_dyn$tau_A2)))
	expect_true(all(result_dyn$tau_A2 > 0))
	expect_true(all(is.finite(result_dyn$tau_B2)))
	expect_true(all(result_dyn$tau_B2 > 0))
	expect_true(all(is.finite(result_dyn$g2)))
	expect_true(all(result_dyn$g2 > 0))
})

test_that("C++ stability functions handle edge cases", {
	# ensure_positive_definite with negative eigenvalues
	M <- matrix(c(1, 2, 2, 1), 2, 2)
	M_pd <- ensure_positive_definite(M, 0.1)
	expect_true(all(eigen(M_pd)$values > 0))

	# ensure_positive_definite with near-symmetric input
	A_asym <- matrix(c(1, 0.5, 0.5000001, 1), 2, 2)
	expect_silent({
		A_pd <- ensure_positive_definite(A_asym, min_eigenvalue = 0.01)
	})
	expect_equal(A_pd, t(A_pd), tolerance = 1e-10)
	expect_true(all(eigen(A_pd)$values > 0))

	# safe_cholesky
	A <- matrix(c(1, 0.99, 0.99, 1), 2, 2)
	L <- matrix(0, 2, 2)
	success <- safe_cholesky(L, A, 1e-6)
	expect_true(success)
	expect_equal(L %*% t(L), A, tolerance = 1e-10)

	# stabilize_spectral_radius
	A_unstable <- matrix(c(1.5, 0.5, 0.5, 1.5), 2, 2)
	A_stable <- stabilize_spectral_radius(A_unstable, 0.99)
	sr <- max(abs(eigen(A_stable)$values))
	expect_true(sr <= 0.99)
})

test_that("edge case data handling", {
	# dynamic needs >= 2 time points
	Y_min <- array(c(1, 2, 2, 3), dim = c(2, 2, 1, 1))
	expect_error({
		dbn(Y_min, model = "dynamic", nscan = 10, burn = 0, odens = 1, verbose = FALSE)
	}, "at least 2 time points")

	# missing data patterns
	Y_miss <- array(rnorm(5*5*2*10), dim = c(5, 5, 2, 10))
	Y_miss[,,1,1:5] <- NA
	expect_no_error(suppressWarnings({
		res <- dbn(Y_miss, model = "static", nscan = 100, burn = 20, odens = 1, verbose = FALSE)
	}))

	# single node
	expect_error({
		Y_single <- array(1, dim = c(1, 1, 2, 10))
		dbn(Y_single, model = "static")
	}, "at least 2 nodes")
})

#### stress tests

test_that("near-constant gaussian data runs without crash", {
	skip_on_cran()
	set.seed(301)
	n <- 8
	Tt <- 5
	Y <- array(5.0 + rnorm(n * n * 1 * Tt, sd = 0.001), dim = c(n, n, 1, Tt))
	for (t in 1:Tt) diag(Y[,,1,t]) <- NA

	fit <- dbn(Y, model = "static", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_true(all(is.finite(fit$draws$pars$s2)),
		label = "sigma2 finite for near-constant data")
})

test_that("highly imbalanced binary data runs without crash", {
	skip_on_cran()
	set.seed(302)
	n <- 10
	Tt <- 10
	Y <- array(0, dim = c(n, n, 1, Tt))
	for (t in 1:Tt) {
		Y[,,1,t] <- matrix(rbinom(n * n, 1, 0.02), n, n)
		diag(Y[,,1,t]) <- NA
	}

	fit <- dbn(Y, model = "static", family = "binary",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_true(all(is.finite(fit$draws$pars$s2)),
		label = "sigma2 finite for sparse binary")
})

test_that("very rectangular bipartite runs without crash", {
	skip_on_cran()
	set.seed(303)
	n_row <- 4
	n_col <- 20
	Tt <- 5
	Y <- array(rnorm(n_row * n_col * 1 * Tt), dim = c(n_row, n_col, 1, Tt))

	suppressWarnings({
		fit <- dbn(Y, model = "dynamic", family = "gaussian",
			nscan = 50, burn = 20, verbose = FALSE)
	})
	expect_s3_class(fit, "dbn")
	expect_true(all(is.finite(fit$sigma2)),
		label = "sigma2 finite for very rectangular bipartite")
})

test_that("long time series with small n runs without crash", {
	skip_on_cran()
	set.seed(304)
	sim <- simulate_dynamic_dbn(n = 5, p = 1, time = 50,
		sigma2 = 0.3, tauA2 = 0.05, tauB2 = 0.05, seed = 304)

	fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
		nscan = 50, burn = 20, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_true(all(is.finite(fit$sigma2)),
		label = "sigma2 finite for long time series")
})
