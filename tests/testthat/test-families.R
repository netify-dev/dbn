test_that("Gaussian family works correctly", {
  # simulate very small gaussian network with stable data
  set.seed(6886)
  m <- 3
  Tt <- 5
  
  # create simple, stable test data
  Y_gau <- array(0, c(m, m, 1, Tt))
  for(t in 1:Tt) {
    # create a simple pattern with small noise
    base <- matrix(0.1 * (1:m) %*% t(1:m), m, m)
    Y_gau[,,1,t] <- base + matrix(rnorm(m*m, 0, 0.1), m, m)
  }
  
  # fit low-rank dbn with fewer iterations
  fit_g <- dbn_lowrank(Y_gau, r = 1, nscan = 50, burn = 10, odens = 1,
                       family = "gaussian", verbose = FALSE)
  
  expect_s3_class(fit_g, "dbn")
  expect_true("sigma2_proc" %in% names(fit_g))
  expect_true("sigma2_obs" %in% names(fit_g))
  expect_equal(fit_g$settings$family, "gaussian")
})

# test_that("Binary family works correctly", {
#   skip_if_not_installed("truncnorm")
#   
#   # Simulate probit-binary network
#   set.seed(6885)
#   m <- 4
#   Tt <- 6
#   Theta <- array(rnorm(m*m*Tt, 0, 0.5), c(m, m, Tt))
#   eta <- Theta  # use same latent mean
#   prob <- pnorm(eta)
#   Y_bin <- array(rbinom(length(prob), 1, prob), dim = c(m, m, 1, Tt))
#   
#   # Fit HMM DBN
#   fit_b <- dbn_hmm(Y_bin, R = 2, nscan = 300, burn = 100, odens = 1,
#                    family = "binary", verbose = FALSE)
#   
#   expect_s3_class(fit_b, "dbn")
#   expect_equal(fit_b$settings$family, "binary")
#   expect_true("S" %in% names(fit_b))
# })

