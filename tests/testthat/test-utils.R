# Test suite for utility functions

test_that("mat unfold operation works correctly", {
  # create test array
  X <- array(1:24, dim = c(2, 3, 4))
  
  # test mode 1 unfold
  X1 <- mat(X, 1)
  expect_equal(dim(X1), c(2, 12))
  expect_equal(X1[1,1], 1)
  expect_equal(X1[2,1], 2)
  
  # test mode 2 unfold
  X2 <- mat(X, 2)
  expect_equal(dim(X2), c(3, 8))
  
  # test mode 3 unfold
  X3 <- mat(X, 3)
  expect_equal(dim(X3), c(4, 6))
})

test_that("amprod array-matrix product works", {
  # create test data
  A <- array(1:24, dim = c(2, 3, 4))
  M <- matrix(1:6, nrow = 3, ncol = 2)
  
  # test multiplication along mode 1
  result <- amprod(A, M, 1)
  expect_equal(dim(result), c(3, 3, 4))
  
  # test with identity matrix
  I <- diag(2)
  result_I <- amprod(A, I, 1)
  expect_equal(result_I, A)
})

test_that("tprod tensor product works", {
  # create test tensor
  A <- array(rnorm(24), dim = c(2, 3, 4))
  
  # create transformation matrices
  B <- list(
    matrix(rnorm(4), 2, 2),
    matrix(rnorm(9), 3, 3),
    matrix(rnorm(16), 4, 4)
  )
  
  # test full tensor product
  result <- tprod(A, B)
  expect_equal(dim(result), dim(A))
  
  # test partial modes
  result_partial <- tprod(A, B[1:2], modes = 1:2)
  expect_equal(dim(result_partial), dim(A))
})

test_that("rsan generates correct dimensions", {
  # test various dimensions
  d1 <- c(5, 5)
  X1 <- rsan(d1)
  expect_equal(dim(X1), d1)
  expect_true(all(!is.na(X1)))
  
  d2 <- c(3, 4, 5)
  X2 <- rsan(d2)
  expect_equal(dim(X2), d2)
  
  # check approximately normal
  d3 <- c(100, 100)
  X3 <- rsan(d3)
  expect_true(abs(mean(X3)) < 0.1)  # should be near 0
  expect_true(abs(sd(c(X3)) - 1) < 0.1)  # should be near 1
})

test_that("zscores handles 2D and 3D arrays correctly", {
  # test 2d array
  y2d <- matrix(c(1, 3, 2, 5, 4, NA), nrow = 2)
  z2d <- zscores(y2d)
  expect_equal(dim(z2d), dim(y2d))
  expect_true(is.na(z2d[2, 3]))
  expect_true(all(!is.na(z2d[!is.na(y2d)])))
  
  # test 3d array
  y3d <- array(c(1:24, NA, 26:30), dim = c(5, 3, 2))
  z3d <- zscores(y3d)
  expect_equal(dim(z3d), dim(y3d))
  expect_equal(sum(is.na(z3d)), 1)
  
  # test ties
  y_ties <- matrix(c(1, 1, 2, 2, 3, 3), nrow = 2)
  z_ties <- zscores(y_ties, ties.method = "average")
  expect_equal(z_ties[1,1], z_ties[2,1])  # same rank -> same z-score
})

# Note: rz_fc has been completely rewritten with different signature
# Old tests removed as they no longer apply to the new implementation




test_that("matrix operations preserve structure", {
  # test that unfold-refold preserves array
  X <- array(rnorm(60), dim = c(3, 4, 5))
  
  # unfold and refold mode 1
  X_mat <- mat(X, 1)
  X_back <- array(X_mat, dim = c(3, 4, 5))
  expect_equal(X_back, X)
  
  # test amprod with identity preserves structure
  I <- diag(3)
  X_same <- amprod(X, I, 1)
  expect_equal(X_same, X)
})

test_that("edge cases for array functions", {
  # empty arrays
  expect_error(mat(array(dim = c(0, 5, 3)), 1))
  
  # single element arrays
  X_single <- array(42, dim = c(1, 1, 1))
  expect_equal(mat(X_single, 1), matrix(42, 1, 1))
  
  # large dimensions
  X_large <- array(1, dim = c(2, 2, 1000))
  X_mat <- mat(X_large, 3)
  expect_equal(dim(X_mat), c(1000, 4))
})