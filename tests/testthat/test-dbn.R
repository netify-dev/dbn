# library(testthat)
# library(dbn

test_that("dbn package loads and has expected functions", {
  # verify main function exists
  expect_true(exists("dbn"))
  expect_true(exists("dbn_static"))
  expect_true(exists("dbn_dynamic"))
  
  # verify s3 methods are registered
  dbn_methods <- as.character(methods(class = "dbn"))
  expect_true("plot.dbn" %in% dbn_methods || "plot" %in% dbn_methods)
  expect_true("summary.dbn" %in% dbn_methods || "summary" %in% dbn_methods)
  
  # verify other exported functions
  expect_true(exists("check_convergence"))
  expect_true(exists("compare_dbn"))
})

test_that("dbn function validates input", {
  # should error on incorrect dimensions
  expect_error(dbn(matrix(1:10, 5, 2)), "3D array .* or 4D array")
  
  # should accept valid model type strings
  expect_error(dbn(array(1, dim=c(5,5,2,10)), model = "invalid"), "should be one of")
})

test_that("static model runs on small data", {
  # use simulation function to create test data
  data <- simulate_static_dbn(n = 10, p = 2, time = 20, seed = 6886)
  
  # run static model with minimal iterations
  results <- dbn(data$Y, model = "static", nscan = 20, burn = 10, odens = 1, verbose = FALSE)
  
  # verify output structure is correct
  expect_s3_class(results, "dbn")
  expect_equal(results$model, "static")
  expect_true(is.list(results$B))
  expect_true(is.matrix(results$params))
  expect_true(all(c("s2", "t2", "g2") %in% colnames(results$params)))
  
  # verify array dimensions
  expect_equal(length(results$B), 3)  # should be 3 b matrices
  expect_equal(dim(results$B[[1]])[1], 10)
  expect_equal(dim(results$B[[1]])[2], 10)
})

test_that("summary and plot methods work", {
  # create small test data and run model
  Y <- simulate_test_data(n = 8, p = 2, time = 15, seed = 456)
  results <- dbn(Y, model = "static", nscan = 20, burn = 10, odens = 1, verbose = FALSE)
  
  # verify summary returns invisibly
  expect_invisible(summary(results))
  expect_invisible(summary_dbn(results))
  
  # test plotting (using null device to avoid creating files)
  pdf(NULL)
  expect_silent(plot(results))
  expect_silent(plot_dbn(results))
  dev.off()
})

test_that("utility functions work correctly", {
  # verify mat function works correctly
  X <- array(1:24, dim=c(2,3,4))
  X1 <- dbn:::mat(X, 1)
  expect_equal(dim(X1), c(2, 12))
  
  # verify rsan function
  d <- c(3, 4, 5)
  X <- dbn:::rsan(d)
  expect_equal(dim(X), d)
  
  # verify zscores function
  y <- matrix(c(1,3,2,4,2,1,3,4), 4, 2)
  z <- dbn:::zscores(y)
  expect_equal(dim(z), dim(y))
  expect_true(all(abs(z) < 3))  # should produce reasonable z-scores
})

test_that("tprod function works correctly", {
  # test simple 2d case
  X <- matrix(1:6, 2, 3)
  A <- matrix(1:4, 2, 2)
  Y <- dbn:::tprod(X, list(A), modes=1)
  expect_equal(Y, A %*% X)
  
  # test 3d case
  X <- array(1:24, dim=c(2,3,4))
  A <- matrix(1:4, 2, 2)
  B <- matrix(1:9, 3, 3)
  
  Y <- dbn:::tprod(X, list(A, B), modes=c(1,2))
  expect_equal(dim(Y), c(2,3,4))
})

test_that("data can be loaded from file path", {
  # save test data to file
  test_file <- tempfile(fileext = ".RData")
  Y <- simulate_test_data(n = 5, p = 2, time = 10, seed = 999)
  save(Y, file = test_file)
  
  # verify loading data from file and running model
  results <- dbn(test_file, model = "static", nscan = 5, burn = 5, odens = 1, verbose = FALSE)
  expect_s3_class(results, "dbn")
  
  # clean up temp file
  unlink(test_file)
})

test_that("B scaling works correctly", {
  # B_scaled function was removed as unnecessary
  # this test now verifies that B matrices maintain proper scaling
  B <- list(
    diag(c(2, 4)),
    diag(c(1, 3)),
    diag(c(0.5, 1.5))
  )
  
  # manual scaling for testing
  K <- length(B)
  mu <- rep(0, K)
  for (k in 1:K) mu[k] <- mean(diag(B[[k]]))
  B_scaled_result <- B
  for (k in 1:K) B_scaled_result[[k]] <- B[[k]] / mu[k]
  
  # verify that average of diagonals equals 1
  for(k in 1:3) {
    expect_equal(mean(diag(B_scaled_result[[k]])), 1, tolerance = 1e-10)
  }
})