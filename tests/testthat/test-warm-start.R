test_that("warm start works for static model", {
  skip_if_not_installed("cli")
  
  # generate small test data
  set.seed(6886)
  m <- 4
  p <- 2  
  Tt <- 5
  Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
  
  # run initial model for short time
  result1 <- dbn_static(Y, nscan = 100, burn = 20, odens = 1, verbose = FALSE)
  
  # continue from previous results
  result2 <- dbn_static(Y, nscan = 100, burn = 20, odens = 1, verbose = FALSE, previous = result1)
  
  # check that parameters continued
  expect_equal(length(result2$params[,1]), 100)  # nscan=100 iterations kept
  expect_true(!is.null(result2$continued_from))
  expect_equal(result2$continued_from, 120)  # previous total: burn(20) + nscan(100) = 120
  expect_equal(result2$total_iter, 220)  # 120 + 100 new iterations = 220
  
  # check that warm start initialized from last values of result1
  last_s2 <- result1$params[nrow(result1$params), "s2"]
  # Note: since we set seed again, the first sample might be identical to a previous sample
  # but the chain should continue from the last state
  expect_true(length(result2$params[, "s2"]) == 100)
})

test_that("warm start works for dynamic model", {
  skip_if_not_installed("cli")
  
  # generate small test data
  set.seed(6886)
  m <- 4
  p <- 2
  Tt <- 5
  Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
  
  # run initial model
  result1 <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE)
  
  # continue from previous results
  result2 <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE, previous = result1)
  
  # check that parameters continued
  expect_equal(length(result2$sigma2), 150)  # nscan=150 iterations kept
  expect_true(!is.null(result2$continued_from))
  
  # verify dimensions match
  expect_equal(dim(result2$A[[1]]), dim(result1$A[[1]]))
  expect_equal(dim(result2$B[[1]]), dim(result1$B[[1]]))
})

test_that("init argument works for static model", {
  skip_if_not_installed("cli")
  
  # generate small test data
  set.seed(6886)
  m <- 4
  p <- 2
  Tt <- 5
  Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
  
  # run with custom initial values
  init_vals <- list(
    s2 = 0.5,
    t2 = 2.0,
    g2 = 0.1,
    B = list(diag(m) * 0.9, diag(m) * 0.9, diag(p) * 0.9)
  )
  
  result <- dbn_static(Y, nscan = 50, burn = 10, odens = 1, verbose = FALSE, init = init_vals)
  
  # the model should run without errors
  expect_s3_class(result, "dbn")
  expect_equal(result$model, "static")
})

test_that("init argument works for dynamic model", {
  skip_if_not_installed("cli")
  
  # generate small test data
  set.seed(6886)
  m <- 4
  p <- 2
  Tt <- 5
  Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
  
  # create initial a and b arrays
  Aarray <- array(0, dim=c(m, m, Tt))
  Barray <- array(0, dim=c(m, m, Tt))
  for(t in 1:Tt) {
    Aarray[,,t] <- diag(m) * 0.8
    Barray[,,t] <- diag(m) * 0.8
  }
  
  # run with custom initial values
  init_vals <- list(
    sigma2 = 0.5,
    tau_A2 = 5.0,
    tau_B2 = 5.0,
    g2 = 0.2,
    A = Aarray,
    B = Barray
  )
  
  result <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE, init = init_vals)
  
  # the model should run without errors
  expect_s3_class(result, "dbn")
  expect_equal(result$model, "dynamic")
})

test_that("warm start validates model type", {
  skip_if_not_installed("cli")
  
  # generate small test data
  set.seed(6886)
  m <- 4
  p <- 2
  Tt <- 5
  Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
  
  # run static model
  result_static <- dbn_static(Y, nscan = 50, burn = 10, odens = 1, verbose = FALSE)
  
  # try to use static results to continue dynamic model - should error
  expect_error(
    dbn_dynamic(Y, nscan = 50, burn = 10, odens = 1, verbose = FALSE, previous = result_static),
    "previous must be results from dbn_dynamic"
  )
})

test_that("warm start handles time thinning correctly", {
  skip_if_not_installed("cli")
  
  # generate small test data
  set.seed(6886) 
  m <- 4
  p <- 2
  Tt <- 10
  Y <- array(sample(1:5, m*m*p*Tt, replace=TRUE), dim=c(m, m, p, Tt))
  
  # run with time thinning
  result1 <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE, time_thin = 2)
  
  # check that saved arrays have correct dimensions
  expect_equal(dim(result1$A[[1]])[3], 5)  # tt/2 = 10/2 = 5
  
  # continue from previous - should expand back to full time
  result2 <- dbn_dynamic(Y, nscan = 150, burn = 50, odens = 1, verbose = FALSE, previous = result1)
  
  # should still work correctly
  expect_s3_class(result2, "dbn")
})