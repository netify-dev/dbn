#### IRF essentials

test_that("impulse_response_const and _dynamic produce sensible cubes", {
	n = 4; H = 3
	A = matrix(rnorm(n * n) * 0.3, n, n); diag(A) = 0.5
	B = matrix(rnorm(n * n) * 0.3, n, n); diag(B) = 0.5
	S = matrix(0, n, n); S[1, 2] = 1
	d_const = impulse_response_const(A, B, S, H)
	expect_equal(dim(d_const), c(n, n, H + 1))
	expect_equal(d_const[, , 1], S, tolerance = 1e-10)

	Aarr = array(0, dim = c(n, n, 5))
	Barr = array(0, dim = c(n, n, 5))
	for (t in 1:5) { Aarr[, , t] = A; Barr[, , t] = B }
	d_dyn = impulse_response_dynamic(Aarr, Barr, S, 0L, H)
	expect_equal(dim(d_dyn), c(n, n, H + 1))
})

test_that("build_shock validates input and produces correct shapes", {
	# n_row <= 0 should error
	expect_error(build_shock(0, type = "unit_edge"), regexp = "n_row|positive")
	# unknown type should error (match.arg)
	expect_error(build_shock(5, type = "bogus"))
	s_unit = build_shock(5, type = "unit_edge", i = 1, j = 2)
	expect_equal(dim(s_unit), c(5, 5))
	expect_equal(s_unit[1, 2], 1)
})

test_that("compute_irf returns a data.frame on a dynamic fit", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 8, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 30, burn = 10, verbose = FALSE)
	irf = compute_irf(fit, shock = "unit_edge", H = 3, n_draws = 4,
		shock_pars = list(i = 1, j = 2), seed = 1)
	expect_s3_class(irf, "dbn_irf")
	expect_s3_class(irf, "data.frame")
})

# static IRF should use the sender factor A = fit$B[[1]] (not identity)
# so Delta_h = A^h S (B^h)' reflects both sender and receiver dynamics
test_that("static IRF closed form uses sender factor A, not identity", {
	skip_on_cran()
	sim = simulate_static_dbn(n = 5, time = 6, p = 1, seed = 7)
	fit = suppressWarnings(suppressMessages(
		dbn(sim$Y, model = "static", family = "gaussian",
			nscan = 60, burn = 20, verbose = FALSE)
	))
	# inspect the post-hoc IRF: A = posterior mean of B[[1]], B = posterior mean of B[[2]]
	A_post <- apply(fit$B[[1]], c(1, 2), mean)
	B_post <- apply(fit$B[[2]], c(1, 2), mean)
	# pick a shock with non-zero entry; compute the closed form on the
	# posterior-mean operators and confirm it matches Delta_1 = A S B'
	S <- matrix(0, 5, 5); S[1, 2] <- 1
	expected_h1 <- A_post %*% S %*% t(B_post)
	# call the C++ helper directly on the posterior mean operators
	delta <- impulse_response_const(A_post, B_post, S, H = 3L)
	expect_equal(delta[, , 1], S, tolerance = 1e-10)
	expect_equal(delta[, , 2], expected_h1, tolerance = 1e-10)
	# the IRF should NOT equal S at h = 1 (identity contraction would force this);
	# instead it should reflect A %*% S %*% B'
	expect_false(isTRUE(all.equal(delta[, , 2], S)))
})

# the user-facing compute_irf wrapper should pick up the sender factor too
test_that("compute_irf on a static fit propagates through A", {
	skip_on_cran()
	sim = simulate_static_dbn(n = 4, time = 6, p = 1, seed = 11)
	fit = suppressWarnings(suppressMessages(
		dbn(sim$Y, model = "static", family = "gaussian",
			nscan = 50, burn = 20, verbose = FALSE)
	))
	options(dbn.suppress_static_irf_warning = TRUE)
	on.exit(options(dbn.suppress_static_irf_warning = FALSE), add = TRUE)
	irf = compute_irf(fit, shock = "unit_edge", H = 3, n_draws = 4,
		shock_pars = list(i = 1, j = 2), seed = 1)
	expect_s3_class(irf, "dbn_irf")
	# the IRF profile should not be flat (would be flat under A = B = I)
	stat_vals <- irf$mean
	expect_gt(diff(range(stat_vals)), 0)
})
