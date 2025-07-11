test_that("ensure_positive_definite handles non-symmetric matrices", {
  skip_if_not(exists("ensure_positive_definite"))
  
  # create a matrix that's almost symmetric but not quite
  # (this mimics numerical precision issues)
  A <- matrix(c(1, 0.5, 0.5000001, 1), 2, 2)
  
  # this should not produce a warning
  expect_silent({
    A_pd <- ensure_positive_definite(A, min_eigenvalue = 0.01)
  })
  
  # result should be symmetric
  expect_equal(A_pd, t(A_pd), tolerance = 1e-10)
  
  # result should be positive definite
  eigenvals <- eigen(A_pd)$values
  expect_true(all(eigenvals > 0))
})

test_that("ensure_positive_definite preserves exact symmetry", {
  skip_if_not(exists("ensure_positive_definite"))
  
  # create a perfectly symmetric matrix
  A <- matrix(c(2, 1, 1, 2), 2, 2)
  
  A_pd <- ensure_positive_definite(A, min_eigenvalue = 0.01)
  
  # should remain symmetric
  expect_equal(A_pd, t(A_pd))
  
  # should be positive definite
  expect_true(all(eigen(A_pd)$values > 0))
})

test_that("ensure_positive_definite handles negative eigenvalues", {
  skip_if_not(exists("ensure_positive_definite"))
  
  # create a matrix with negative eigenvalues
  A <- matrix(c(1, 2, 2, 1), 2, 2)  # eigenvalues are 3 and -1
  
  A_pd <- ensure_positive_definite(A, min_eigenvalue = 0.1)
  
  # all eigenvalues should be at least min_eigenvalue
  eigenvals <- eigen(A_pd)$values
  expect_true(all(eigenvals >= 0.1))
  
  # should be symmetric
  expect_equal(A_pd, t(A_pd), tolerance = 1e-10)
})