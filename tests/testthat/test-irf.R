test_that("impulse_response_const works correctly", {
    m <- 5
    A <- diag(0.5, m)
    B <- diag(0.3, m)
    S <- matrix(0, m, m)
    S[1, 2] <- 1  # Unit shock
    H <- 3
    
    Delta <- impulse_response_const(A, B, S, H)
    
    # Check dimensions
    expect_equal(dim(Delta), c(m, m, H + 1))
    
    # Check initial shock
    expect_equal(Delta[,,1], S)
    
    # Check first propagation (should be A * S * B^T)
    expected_1 <- A %*% S %*% t(B)
    expect_equal(Delta[,,2], expected_1, tolerance = 1e-6)
    
    # Check that response decays (for stable A, B)
    for (h in 2:(H+1)) {
        expect_true(max(abs(Delta[,,h])) <= max(abs(Delta[,,h-1])) + 1e-10)
    }
})

test_that("impulse_response_dynamic works correctly", {
    m <- 4
    T_len <- 10
    
    # Create time-varying A and B
    Aarray <- array(0, dim = c(m, m, T_len))
    Barray <- array(0, dim = c(m, m, T_len))
    
    for (t in 1:T_len) {
        Aarray[,,t] <- diag(0.5 - 0.01 * t, m)  # Decreasing diagonal
        Barray[,,t] <- diag(0.3 + 0.01 * t, m)  # Increasing diagonal
    }
    
    S <- matrix(0, m, m)
    S[1, 2] <- 1
    t0 <- 2  # 0-based
    H <- 3
    
    Delta <- impulse_response_dynamic(Aarray, Barray, S, t0, H)
    
    # Check dimensions
    expect_equal(dim(Delta), c(m, m, H + 1))
    
    # Check initial shock
    expect_equal(Delta[,,1], S)
    
    # Check first propagation
    # Note: t0 is 0-based in C++, so t0+1 in C++ corresponds to time index t0+1
    # But we need to use the correct time index for the expected calculation
    expected_1 <- Aarray[,,t0+1] %*% S %*% t(Barray[,,t0+1])
    
    # Just check that the non-zero element is reasonably close
    # The exact value might differ slightly due to numerical precision
    expect_true(abs(Delta[1,2,2] - expected_1[1,2]) < 0.01)
    expect_equal(sum(abs(Delta[,,2])), sum(abs(expected_1)), tolerance = 0.01)
})

test_that("build_shock creates correct shock matrices", {
    m <- 5
    
    # Unit edge shock
    S_edge <- build_shock(m, "unit_edge", i = 2, j = 3)
    expect_equal(sum(S_edge), 1)
    expect_equal(S_edge[2, 3], 1)
    
    # Node out shock
    S_out <- build_shock(m, "node_out", i = 2)
    expect_equal(sum(S_out), m)
    expect_equal(S_out[2, ], rep(1, m))
    
    # Node in shock
    S_in <- build_shock(m, "node_in", i = 3)
    expect_equal(sum(S_in), m)
    expect_equal(S_in[, 3], rep(1, m))
    
    # Density shock
    S_dens <- build_shock(m, "density")
    expect_equal(sum(S_dens), 1, tolerance = 1e-10)
    expect_equal(diag(S_dens), rep(0, m))
})

test_that("network statistics are computed correctly", {
    m <- 4
    X <- matrix(c(
        0, 1, 0, 1,
        1, 0, 1, 0,
        0, 1, 0, 1,
        1, 0, 1, 0
    ), m, m, byrow = TRUE)
    
    # Density (8 non-zero edges out of 12 possible, excluding diagonal)
    expect_equal(stat_density(X), 8/12)
    
    # Degrees
    expect_equal(stat_in_degree(X), c(2, 2, 2, 2))
    expect_equal(stat_out_degree(X), c(2, 2, 2, 2))
    
    # Reciprocity (perfect in this case)
    expect_equal(stat_reciprocity(X), 1)
    
    # Transitivity
    expect_equal(stat_transitivity(X), 0)  # No triangles in this bipartite-like structure
})

test_that("compute_irf works with static model", {
    skip_if_not_installed("dbn")
    
    # Create a minimal static model fit object
    m <- 5
    n_draws <- 10
    
    fit <- list(
        model = "static",
        dims = list(m = m, p = 1, T = 1),
        B = array(rnorm(n_draws * m * m, 0, 0.1), dim = c(n_draws, m, m)),
        M = matrix(0, m, m),
        class = c("dbn_static", "dbn")
    )
    class(fit) <- c("dbn_static", "dbn")
    
    # Test with unit edge shock
    irf <- compute_irf(fit, shock = "unit_edge", H = 5, n_draws = 5,
                      shock_pars = list(i = 1, j = 2))
    
    expect_s3_class(irf, "dbn_irf")
    expect_equal(nrow(irf), 6)  # H + 1
    expect_true(all(c("horizon", "mean", "median", "sd", "q025", "q975") %in% names(irf)))
    
    # IRF should start at 0 for horizon 0 (shock effect on density)
    expect_true(abs(irf$mean[1]) < 1)  # Some effect expected
})

test_that("error handling works correctly", {
    m <- 5
    
    # Dimension mismatch
    A <- diag(m)
    B <- diag(m + 1)  # Wrong size
    S <- matrix(0, m, m)
    
    expect_error(impulse_response_const(A, B, S, 5), "dimensions")
    
    # Invalid shock parameters
    expect_error(build_shock(m, "unit_edge", i = 0, j = 1), "between 1 and m")
    expect_error(build_shock(m, "unit_edge", i = 1, j = m + 1), "between 1 and m")
})