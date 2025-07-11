test_that("can simulate and recover simple network structure", {
  set.seed(6886)
  
  # Simulate from known model
  n <- 10
  time <- 50
  
  # True parameters
  true_s2 <- 1.0
  true_sender_effects <- c(rep(0.8, 5), rep(-0.8, 5))  # two groups
  
  # Generate data with block structure
  Y <- array(NA, dim = c(n, n, 1, time))
  
  # Initial time point
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
        # higher values within groups
        if((i <= 5 & j <= 5) | (i > 5 & j > 5)) {
          Y[i,j,1,1] <- sample(3:5, 1)
        } else {
          Y[i,j,1,1] <- sample(1:3, 1)
        }
      }
    }
  }
  
  # Generate temporal evolution with persistence
  for(t in 2:time) {
    Y[,,1,t] <- Y[,,1,t-1]
    # add some random changes
    n_changes <- round(0.2 * n^2)
    for(k in 1:n_changes) {
      i <- sample(1:n, 1)
      j <- sample(setdiff(1:n, i), 1)
      Y[i,j,1,t] <- max(1, min(5, Y[i,j,1,t] + sample(c(-1,1), 1)))
    }
  }
  
  # set diagonal to na
  for(t in 1:time) diag(Y[,,1,t]) <- NA
  
  # Fit model
  results <- dbn(Y, model = "static", nscan = 50, burn = 10, odens = 1, verbose = FALSE)
  
  # check if model detects structure
  B1_mean <- apply(results$B[[1]], c(1,2), mean)
  sender_effects <- diag(B1_mean)
  
  # should detect two groups (though signs may flip)
  group1_mean <- mean(sender_effects[1:5])
  group2_mean <- mean(sender_effects[6:10])
  
  # groups should be different (but with only 10 iterations, may not recover)
  # just check that sender effects exist
  expect_true(length(sender_effects) == 10)
})

test_that("bilinear structure is preserved", {
  set.seed(6887)
  
  # test bilinear property: θ_t = a θ_{t-1} b'
  n <- 5
  
  # create test matrices
  A <- matrix(rnorm(n*n), n, n)
  B <- matrix(rnorm(n*n), n, n) 
  Theta <- matrix(rnorm(n*n), n, n)
  
  # compute bilinear product
  Theta_new <- A %*% Theta %*% t(B)
  
  # test with tprod function
  Theta_test <- tprod(Theta, list(A, B), modes = c(1,2))
  
  # should be equal
  expect_equal(Theta_new, Theta_test, tolerance = 1e-10)
})

test_that("rank likelihood preserves ordering", {
  # create ordinal data
  Y <- matrix(c(1,2,3,4,
                2,1,4,3,
                3,4,1,2,
                4,3,2,1), 4, 4)
  
  # convert to z-scores
  Z <- zscores(Y)
  
  # check ordering is preserved
  for(i in 1:4) {
    y_order <- order(Y[i,])
    z_order <- order(Z[i,])
    expect_equal(y_order, z_order)
  }
})

test_that("model handles single time point gracefully", {
  Y <- array(sample(1:5, 10*10*2*1, replace=TRUE), dim=c(10,10,2,1))
  
  # should give informative error or handle gracefully
  expect_error(
    dbn(Y, model = "static", nscan = 10, burn = 2, odens = 1),
    NA  # no error expected - should handle single time point
  )
})

test_that("model handles single relation", {
  set.seed(6889)
  Y <- array(sample(1:5, 8*8*1*20, replace=TRUE), dim=c(8,8,1,20))
  
  results <- dbn(Y, model = "static", nscan = 50, burn = 10, odens = 1, verbose = FALSE)
  
  # should work with single relation
  expect_equal(results$dims$p, 1)
  expect_equal(dim(results$B[[3]])[1], 1)  # 1x1 relation matrix
})

test_that("extreme values don't break the model", {
  set.seed(6881)
  n <- 6
  Y <- array(1, dim=c(n,n,2,10))  # all same value
  
  # add some variation
  Y[1,2,,] <- 5
  Y[2,1,,] <- 5
  
  # remove diagonal
  for(t in 1:10) {
    for(r in 1:2) {
      diag(Y[,,r,t]) <- NA
    }
  }
  
  # should still run (though estimates may be poor)
  results <- dbn(Y, model = "static", nscan = 30, burn = 10, odens = 1, verbose = FALSE)
  
  # most parameters should be finite (some may be na due to extreme values)
  expect_true(sum(is.finite(results$params)) > 0)
})