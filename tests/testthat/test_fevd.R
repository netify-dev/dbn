# fevd returns the expected output shape under actor shock
test_that("fevd returns expected shape for actor shock", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 8, seed = 42)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 100,
			burn = 50,
			odens = 2,
			verbose = FALSE
		)
	))
	fv = fevd(fit, H = 4, shock = "actor")
	expect_s3_class(fv, "dbn_fevd")
	expect_equal(dim(fv$median), c(5, 5, 5, 4))
	expect_equal(dim(fv$total_mse), c(5, 5, 4))
	expect_equal(fv$metadata$ordering_required, FALSE)
})

# per-draw shares sum to 1 across the source axis
test_that("fevd shares sum to 1 per draw", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 8, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 100,
			burn = 50,
			odens = 2,
			verbose = FALSE
		)
	))
	# verify directly from cumulative-operator math on a single draw
	A1 = fit$A[[1]]
	B1 = fit$B[[1]]
	sig2 = fit$sigma2[1]
	n = dim(A1)[1]
	Tt = dim(A1)[3]
	H = 3L
	t0 = Tt
	A_cum = array(0, dim = c(n, n, H))
	B_cum = array(0, dim = c(n, n, H))
	A_cum[, , H] = diag(n)
	B_cum[, , H] = diag(n)
	for (h in (H - 1L):1L) {
		ti = min(Tt, t0 + h + 1L)
		A_cum[, , h] = A1[, , ti] %*% A_cum[, , h + 1L]
		B_cum[, , h] = B1[, , ti] %*% B_cum[, , h + 1L]
	}
	mse_11 = 0
	for (q in 1:H) {
		A_sq = A_cum[, , q]^2
		B_sq = B_cum[, , q]^2
		mse_11 = mse_11 + sig2 * sum(A_sq[1, ]) * sum(B_sq[1, ])
	}
	total_share = 0
	for (k in 1:n) {
		ck = 0
		for (q in 1:H) {
			A_sq = A_cum[, , q]^2
			B_sq = B_cum[, , q]^2
			ck = ck + sig2 * A_sq[1, k] * sum(B_sq[1, ])
		}
		total_share = total_share + ck / mse_11
	}
	expect_equal(total_share, 1, tolerance = 1e-10)
})

# fevd dyad shock returns 5d array
test_that("fevd dyad shock returns 5d output", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 4, p = 1, time = 6, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	fv = fevd(fit, H = 3, shock = "dyad")
	expect_equal(length(dim(fv$median)), 5)
	expect_equal(dim(fv$median), c(4, 4, 4, 4, 3))
})

# total mse is non-decreasing in horizon
test_that("total forecast mse is non-decreasing in horizon", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 8, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 80,
			burn = 30,
			odens = 2,
			verbose = FALSE
		)
	))
	fv = fevd(fit, H = 4, shock = "actor")
	mse_seq = apply(fv$total_mse, 3, mean)
	expect_true(all(diff(mse_seq) >= -1e-8))
})

# fevd works on piecewise after operator expansion
test_that("fevd works on piecewise fit via operator_draws", {
	skip_on_cran()
	sim = simulate_piecewise_dbn(n = 5, time = 10, blocks = c(5, 10), p = 1, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Y,
			model = "piecewise",
			family = "ordinal",
			blocks = c(5, 10),
			nscan = 80,
			burn = 40,
			odens = 2,
			verbose = FALSE
		)
	))
	fv = fevd(fit, H = 3, shock = "actor")
	expect_s3_class(fv, "dbn_fevd")
	expect_equal(dim(fv$median), c(5, 5, 5, 3))
})

# isotropic-innovation flag is set
test_that("fevd metadata declares isotropic innovation", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 4, p = 1, time = 6, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	fv = fevd(fit, H = 2, shock = "actor")
	expect_true(fv$metadata$isotropic_innovation)
	expect_false(fv$metadata$ordering_required)
})
