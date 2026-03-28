####
# IRF functions
####

test_that("impulse_response_const works correctly", {
		m = 5
		A = diag(0.5, m)
		B = diag(0.3, m)
		S = matrix(0, m, m)
		S[1, 2] = 1  # unit shock
		H = 3
		
		Delta = impulse_response_const(A, B, S, H)
		
		# dimensions
		expect_equal(dim(Delta), c(m, m, H + 1))

		# initial shock
		expect_equal(Delta[,,1], S)

		# first propagation: A * S * B'
		expected_1 = A %*% S %*% t(B)
		expect_equal(Delta[,,2], expected_1, tolerance = 1e-6)
		
		# response decays for stable A, B
		for (h in 2:(H+1)) {
				expect_true(max(abs(Delta[,,h])) <= max(abs(Delta[,,h-1])) + 1e-10)
		}
})

test_that("impulse_response_dynamic works correctly", {
		m = 4
		T_len = 10
		
		# time-varying A and B
		Aarray = array(0, dim = c(m, m, T_len))
		Barray = array(0, dim = c(m, m, T_len))
		
		for (t in 1:T_len) {
				Aarray[,,t] = diag(0.5 - 0.01 * t, m)  # decreasing
				Barray[,,t] = diag(0.3 + 0.01 * t, m)  # increasing
		}
		
		S = matrix(0, m, m)
		S[1, 2] = 1
		t0 = 2  # 0-based
		H = 3
		
		Delta = impulse_response_dynamic(Aarray, Barray, S, t0, H)
		
		# dimensions
		expect_equal(dim(Delta), c(m, m, H + 1))

		# initial shock
		expect_equal(Delta[,,1], S)

		# first propagation (t0 is 0-based in C++)
		expected_1 = Aarray[,,t0+1] %*% S %*% t(Barray[,,t0+1])
		
		# non-zero element close (numerical precision)
		expect_true(abs(Delta[1,2,2] - expected_1[1,2]) < 0.01)
		expect_equal(sum(abs(Delta[,,2])), sum(abs(expected_1)), tolerance = 0.01)
})

test_that("build_shock creates correct shock matrices", {
		m = 5
		
		# unit edge shock
		S_edge = build_shock(m, "unit_edge", i = 2, j = 3)
		expect_equal(sum(S_edge), 1)
		expect_equal(S_edge[2, 3], 1)
		
		# node out shock
		S_out = build_shock(m, "node_out", i = 2)
		expect_equal(sum(S_out), m)
		expect_equal(S_out[2, ], rep(1, m))
		
		# node in shock
		S_in = build_shock(m, "node_in", i = 3)
		expect_equal(sum(S_in), m)
		expect_equal(S_in[, 3], rep(1, m))
		
		# density shock
		S_dens = build_shock(m, "density")
		expect_equal(sum(S_dens), 1, tolerance = 1e-10)
		expect_equal(diag(S_dens), rep(0, m))
})

test_that("network statistics are computed correctly", {
		m = 4
		X = matrix(c(
				0, 1, 0, 1,
				1, 0, 1, 0,
				0, 1, 0, 1,
				1, 0, 1, 0
		), m, m, byrow = TRUE)
		
		# density: 8 edges / 12 possible
		expect_equal(stat_density(X), 8/12)
		
		# degrees
		expect_equal(stat_in_degree(X), c(2, 2, 2, 2))
		expect_equal(stat_out_degree(X), c(2, 2, 2, 2))
		
		# reciprocity (perfect)
		expect_equal(stat_reciprocity(X), 1)
		
		# transitivity
		expect_equal(stat_transitivity(X), 0)  # no triangles
})

test_that("compute_irf works with static model", {
		skip_if_not_installed("dbn")
		set.seed(6886)

		# minimal static model fit
		m = 5
		n_draws = 10

		fit = list(
				model = "static",
				dims = list(m = m, n_row = m, n_col = m, p = 1, Tt = 1, is_bipartite = FALSE),
				B = array(rnorm(n_draws * m * m, 0, 0.1), dim = c(n_draws, m, m)),
				M = matrix(0, m, m),
				class = c("dbn_static", "dbn")
		)
		class(fit) = c("dbn_static", "dbn")
		
		# unit edge shock
		irf = compute_irf(fit, shock = "unit_edge", H = 5, n_draws = 5,
											shock_pars = list(i = 1, j = 2))
		
		expect_s3_class(irf, "dbn_irf")
		expect_equal(nrow(irf), 6)  # H + 1
		expect_true(all(c("horizon", "mean", "median", "sd", "q025", "q975") %in% names(irf)))
		
		# horizon 0 effect bounded
		expect_true(abs(irf$mean[1]) < 1)
})

test_that("error handling works correctly", {
		m = 5
		
		# dimension mismatch
		A = diag(m)
		B = diag(m + 1)  # wrong size
		S = matrix(0, m, m)
		
		expect_error(impulse_response_const(A, B, S, 5), "dimensions")
		
		# invalid shock parameters
		expect_error(build_shock(m, "unit_edge", i = 0, j = 1), "between 1 and")
		expect_error(build_shock(m, "unit_edge", i = 1, j = m + 1), "between 1 and")
})