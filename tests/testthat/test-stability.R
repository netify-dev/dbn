test_that("edge case data handling", {
  
  # test with minimal data
  Y_min <- array(c(1, 2, 2, 3), dim = c(2, 2, 1, 1))
  
  expect_error({
    dbn(Y_min, model = "dynamic", nscan = 10, burn = 0, odens = 1, verbose = FALSE)
  }, "at least 2 time points")
  
  # test with missing data patterns
  Y_miss <- array(rnorm(5*5*2*10), dim = c(5, 5, 2, 10))
  Y_miss[,,1,1:5] <- NA  # first relation missing for half the time
  
  expect_no_error({
    res <- dbn(Y_miss, model = "static", nscan = 100, burn = 20, odens = 1, verbose = FALSE)
  })
  
  # test with single node
  expect_error({
    Y_single <- array(1, dim = c(1, 1, 2, 10))
    dbn(Y_single, model = "static")
  }, "zero variance")
})

test_that("C++ functions handle edge cases", {
  
  # test ensure_positive_definite
  M <- matrix(c(1, 2, 2, 1), 2, 2)  # not positive definite
  M_pd <- ensure_positive_definite(M, 0.1)
  
  expect_true(all(eigen(M_pd)$values > 0))
  
  # test safe_cholesky
  A <- matrix(c(1, 0.99, 0.99, 1), 2, 2)
  L <- matrix(0, 2, 2)
  success <- safe_cholesky(L, A, 1e-6)
  
  expect_true(success)
  expect_equal(L %*% t(L), A, tolerance = 1e-10)
  
  # test stabilize_spectral_radius
  A_unstable <- matrix(c(1.5, 0.5, 0.5, 1.5), 2, 2)
  A_stable <- stabilize_spectral_radius(A_unstable, 0.99)
  
  sr <- max(abs(eigen(A_stable)$values))
  expect_true(sr <= 0.99)
})

test_that("model outputs are consistent", {
  
  # generate stable test data
  set.seed(6886)
  Y <- array(sample(1:5, 5*5*2*10, replace = TRUE), dim = c(5, 5, 2, 10))
  
  # run models with same seed
  set.seed(6887)
  res1 <- dbn(Y, model = "static", nscan = 100, burn = 20, odens = 1, verbose = FALSE)
  
  set.seed(6887)
  res2 <- dbn(Y, model = "static", nscan = 100, burn = 20, odens = 1, verbose = FALSE)
  
  # results should be identical
  expect_equal(res1$A[[1]], res2$A[[1]])
  expect_equal(res1$sigma2, res2$sigma2)
})