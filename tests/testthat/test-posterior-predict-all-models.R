# Test posterior prediction for all model types
# This ensures posterior_predict_dbn works correctly after structural changes

test_that("posterior_predict_dbn works for static model", {
  set.seed(123)
  n <- 8
  p <- 2
  Tt <- 15
  
  # Simulate data
  Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  
  # Fit static model
  fit <- dbn_static(Y, family = "ordinal", nscan = 100, burn = 50, verbose = FALSE)
  
  # Generate posterior predictions
  ppd <- posterior_predict_dbn(fit, ndraws = 20, seed = 123)
  
  # Check structure
  expect_type(ppd, "list")
  expect_equal(length(ppd), 20)
  expect_equal(dim(ppd[[1]]), dim(Y))
  
  # Check values are in valid range for ordinal
  vals <- unique(unlist(ppd))
  vals <- vals[!is.na(vals)]
  expect_true(all(vals %in% 1:5))
})

test_that("posterior_predict_dbn works for dynamic model", {
  set.seed(456)
  n <- 6
  p <- 1
  Tt <- 10
  
  # Simulate data
  Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  
  # Fit dynamic model
  fit <- dbn_dynamic(Y, family = "ordinal", nscan = 100, burn = 50, verbose = FALSE)
  
  # Generate posterior predictions
  ppd <- posterior_predict_dbn(fit, ndraws = 15, seed = 456)
  
  # Check structure
  expect_type(ppd, "list")
  expect_equal(length(ppd), 15)
  expect_equal(dim(ppd[[1]]), dim(Y))
  
  # Check values are in valid range
  vals <- unique(unlist(ppd))
  vals <- vals[!is.na(vals)]
  expect_true(all(vals %in% 1:5))
})

test_that("posterior_predict_dbn works for dynamic model with time thinning", {
  set.seed(789)
  n <- 5
  p <- 2
  Tt <- 20
  
  # Simulate data
  Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  
  # Fit dynamic model with time thinning
  fit <- dbn_dynamic(Y, family = "ordinal", nscan = 100, burn = 50, 
                     time_thin = 2, verbose = FALSE)
  
  # Check that A and B have reduced time dimension
  expect_equal(dim(fit$A[[1]])[3], 10)  # Tt/2
  
  # Generate posterior predictions - should still work
  ppd <- posterior_predict_dbn(fit, ndraws = 10, seed = 789)
  
  # Check structure
  expect_type(ppd, "list")
  expect_equal(length(ppd), 10)
  expect_equal(dim(ppd[[1]]), dim(Y))  # Should match original Y dimensions
})

test_that("posterior_predict_dbn works for gaussian family", {
  set.seed(111)
  n <- 5
  p <- 1
  Tt <- 10
  
  # Simulate gaussian data
  Y <- array(rnorm(n*n*p*Tt), dim = c(n, n, p, Tt))
  
  # Test static model
  fit_static <- dbn_static(Y, family = "gaussian", nscan = 80, burn = 40, verbose = FALSE)
  ppd_static <- posterior_predict_dbn(fit_static, ndraws = 10, seed = 111)
  
  expect_type(ppd_static, "list")
  expect_equal(length(ppd_static), 10)
  expect_equal(dim(ppd_static[[1]]), dim(Y))
  
  # Test dynamic model
  fit_dynamic <- dbn_dynamic(Y, family = "gaussian", nscan = 80, burn = 40, verbose = FALSE)
  ppd_dynamic <- posterior_predict_dbn(fit_dynamic, ndraws = 10, seed = 222)
  
  expect_type(ppd_dynamic, "list")
  expect_equal(length(ppd_dynamic), 10)
  expect_equal(dim(ppd_dynamic[[1]]), dim(Y))
})

test_that("posterior_predict_dbn works for binary family", {
  set.seed(333)
  n <- 6
  p <- 1
  Tt <- 8
  
  # Simulate binary data
  Y <- array(sample(0:1, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  
  # Test static model
  fit_static <- dbn_static(Y, family = "binary", nscan = 60, burn = 30, verbose = FALSE)
  ppd_static <- posterior_predict_dbn(fit_static, ndraws = 8, seed = 333)
  
  expect_type(ppd_static, "list")
  expect_equal(length(ppd_static), 8)
  expect_equal(dim(ppd_static[[1]]), dim(Y))
  
  # Check binary values
  vals <- unique(unlist(ppd_static))
  vals <- vals[!is.na(vals)]
  expect_true(all(vals %in% 0:1))
})

test_that("posterior_predict_dbn works for lowrank model", {
  skip_if_not_installed("Matrix")
  
  set.seed(444)
  n <- 8
  p <- 1
  Tt <- 12
  
  # Simulate data
  Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  
  # Fit lowrank model
  fit <- dbn_lowrank(Y, family = "ordinal", r = 2, nscan = 80, burn = 40, verbose = FALSE)
  
  # Generate posterior predictions
  ppd <- posterior_predict_dbn(fit, ndraws = 10, seed = 444)
  
  # Check structure
  expect_type(ppd, "list")
  expect_equal(length(ppd), 10)
  expect_equal(dim(ppd[[1]]), dim(Y))
})

test_that("posterior_predict_dbn works for hmm model", {
  set.seed(555)
  n <- 6
  p <- 1
  Tt <- 15
  
  # Simulate data
  Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  
  # Fit HMM model with 2 regimes
  fit <- dbn_hmm(Y, family = "ordinal", R = 2, nscan = 100, burn = 50, verbose = FALSE)
  
  # Generate posterior predictions
  ppd <- posterior_predict_dbn(fit, ndraws = 12, seed = 555)
  
  # Check structure
  expect_type(ppd, "list")
  expect_equal(length(ppd), 12)
  expect_equal(dim(ppd[[1]]), dim(Y))
})

test_that("posterior_predict_dbn handles edge cases", {
  set.seed(666)
  n <- 4
  p <- 1
  Tt <- 5
  
  # Small data
  Y <- array(sample(1:3, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  
  # Fit with minimal iterations
  fit <- dbn_static(Y, family = "ordinal", nscan = 20, burn = 10, verbose = FALSE)
  
  # Request more draws than available
  ppd <- posterior_predict_dbn(fit, ndraws = 50, seed = 666)
  
  # Should sample with replacement when requesting more than available
  expect_equal(length(ppd), 50)
  
  # Request specific draws
  ppd_specific <- posterior_predict_dbn(fit, draws = c(1, 5, 10), seed = 777)
  expect_equal(length(ppd_specific), 3)
})

test_that("posterior_predict_dbn works with 3D input data", {
  set.seed(888)
  n <- 5
  Tt <- 10
  
  # Create 3D array (single relation)
  Y_3d <- array(sample(1:5, n*n*Tt, replace = TRUE), dim = c(n, n, Tt))
  
  # Fit static model with 3D input
  suppressMessages({
    fit <- dbn(Y_3d, family = "ordinal", model = "static", 
               nscan = 50, burn = 25, verbose = FALSE)
  })
  
  # Generate predictions
  ppd <- posterior_predict_dbn(fit, ndraws = 10, seed = 888)
  
  # Check that predictions match the expanded 4D structure
  expect_type(ppd, "list")
  expect_equal(length(ppd), 10)
  expect_equal(dim(ppd[[1]]), c(n, n, 1, Tt))  # Should be 4D with p=1
})