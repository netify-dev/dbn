# Test suite for simulation functions

test_that("simulate_static_dbn generates valid data", {
  # Generate data
  data <- simulate_static_dbn(n = 10, p = 2, time = 10, seed = 6886)
  
  # Check structure
  expect_type(data, "list")
  expect_true(all(c("Y", "Z", "A", "B", "M", "sigma2", "tau2") %in% names(data)))
  
  # Check dimensions
  expect_equal(dim(data$Y), c(10, 10, 2, 10))
  expect_equal(dim(data$Z), c(10, 10, 2, 10))
  expect_equal(dim(data$A), c(10, 10))
  expect_equal(dim(data$B), c(10, 10))
  expect_equal(dim(data$M), c(10, 10, 2))
  
  # Check Y is ordinal
  y_vals <- unique(c(data$Y[!is.na(data$Y)]))
  expect_true(all(y_vals %in% 1:5))
  
  # Check self-loops are NA
  for(t in 1:10) {
    for(p in 1:2) {
      expect_true(all(is.na(diag(data$Y[,,p,t]))))
    }
  }
  
  # Check parameters
  expect_equal(data$sigma2, 0.5)
  expect_equal(data$tau2, 0.1)
})

test_that("simulate_dynamic_dbn generates valid data", {
  # Test without AR(1)
  data1 <- simulate_dynamic_dbn(n = 10, p = 2, time = 10, 
                               ar1 = FALSE, seed = 6886)
  
  # Check structure
  expect_equal(dim(data1$Y), c(10, 10, 2, 10))
  expect_equal(dim(data1$A), c(10, 10, 10))  # Time-varying
  expect_equal(dim(data1$B), c(10, 10, 10))
  expect_false(data1$ar1)
  expect_null(data1$rhoA)
  
  # Test with AR(1)
  data2 <- simulate_dynamic_dbn(n = 10, p = 2, time = 10,
                               ar1 = TRUE, rhoA = 0.9, rhoB = 0.8,
                               seed = 6886)
  
  expect_true(data2$ar1)
  expect_equal(data2$rhoA, 0.9)
  expect_equal(data2$rhoB, 0.8)
  
  # Check that A and B vary over time
  A_var <- var(c(data2$A))
  B_var <- var(c(data2$B))
  expect_true(A_var > 0)
  expect_true(B_var > 0)
})

test_that("simulate_lowrank_dbn generates valid data", {
  # Generate low-rank data
  data <- simulate_lowrank_dbn(n = 10, p = 2, time = 10, r = 3, seed = 6886)
  
  # Check structure
  expect_equal(dim(data$Y), c(10, 10, 2, 10))
  expect_equal(dim(data$U), c(10, 3))  # n x r
  expect_equal(dim(data$alpha), c(3, 10))  # r x time
  expect_equal(data$r, 3)
  
  # Check U is orthonormal
  UtU <- t(data$U) %*% data$U
  expect_equal(UtU, diag(3), tolerance = 1e-10)
  
  # Check A matrices are correctly constructed
  for(t in 1:5) {  # Check first few time points
    A_t <- data$U %*% diag(data$alpha[,t]) %*% t(data$U)
    expect_equal(A_t, data$A[,,t], tolerance = 1e-10)
  }
})

test_that("simulate_hmm_dbn generates valid data", {
  # Generate HMM data
  data <- simulate_hmm_dbn(n = 10, p = 2, time = 10, R = 3,
                          transition_prob = 0.8, seed = 6886)
  
  # Check structure
  expect_equal(dim(data$Y), c(10, 10, 2, 10))
  expect_length(data$S, 10)
  expect_equal(length(data$A_list), 3)
  expect_equal(length(data$B_list), 3)
  expect_equal(dim(data$Pi), c(3, 3))
  expect_equal(data$R, 3)
  
  # Check states are valid
  expect_true(all(data$S %in% 1:3))
  
  # Check transition matrix
  expect_true(all(data$Pi >= 0))
  expect_equal(rowSums(data$Pi), rep(1, 3), tolerance = 1e-10)
  expect_equal(diag(data$Pi), rep(0.8, 3))  # Diagonal should be transition_prob
  
  # Check that different regimes have different dynamics
  A1_norm <- norm(data$A_list[[1]] - diag(10), "F")
  A2_norm <- norm(data$A_list[[2]] - diag(10), "F")
  A3_norm <- norm(data$A_list[[3]] - diag(10), "F")
  expect_true(abs(A1_norm - A2_norm) > 0.01)  # Should be different
})

test_that("simulate_test_data generates simple data", {
  # Generate simple test data
  Y <- simulate_test_data(n = 10, p = 2, time = 20, seed = 333)
  
  # Check dimensions
  expect_equal(dim(Y), c(10, 10, 2, 20))
  
  # Check values are ordinal
  y_vals <- unique(c(Y[!is.na(Y)]))
  expect_true(all(y_vals %in% 1:5))
  
  # Check self-loops are NA
  for(t in 1:20) {
    for(p in 1:2) {
      expect_true(all(is.na(diag(Y[,,p,t]))))
    }
  }
  
  # Check relations have different baselines
  mean_r1 <- mean(Y[,,1,], na.rm = TRUE)
  mean_r2 <- mean(Y[,,2,], na.rm = TRUE)
  expect_true(mean_r1 > mean_r2)  # Relation 1 should have higher values
})

test_that("simulation functions respect seeds", {
  # Static model
  data1 <- simulate_static_dbn(n = 10, p = 1, time = 5, seed = 999)
  data2 <- simulate_static_dbn(n = 10, p = 1, time = 5, seed = 999)
  expect_equal(data1$Y, data2$Y)
  
  # Dynamic model
  data3 <- simulate_dynamic_dbn(n = 10, p = 1, time = 5, seed = 999)
  data4 <- simulate_dynamic_dbn(n = 10, p = 1, time = 5, seed = 999)
  expect_equal(data3$Y, data4$Y)
  
  # Low-rank model
  data5 <- simulate_lowrank_dbn(n = 10, p = 1, time = 5, r = 2, seed = 999)
  data6 <- simulate_lowrank_dbn(n = 10, p = 1, time = 5, r = 2, seed = 999)
  expect_equal(data5$Y, data6$Y)
  
  # HMM model
  data7 <- simulate_hmm_dbn(n = 10, p = 1, time = 5, R = 2, seed = 999)
  data8 <- simulate_hmm_dbn(n = 10, p = 1, time = 5, R = 2, seed = 999)
  expect_equal(data7$Y, data8$Y)
  expect_equal(data7$S, data8$S)  # States should also match
})

test_that("simulation functions handle edge cases", {
  # Single time point
  data1 <- simulate_static_dbn(n = 5, p = 1, time = 1, seed = 444)
  expect_equal(dim(data1$Y), c(5, 5, 1, 1))
  
  # Single relation
  data2 <- simulate_dynamic_dbn(n = 5, p = 1, time = 10, seed = 555)
  expect_equal(dim(data2$Y)[3], 1)
  
  # Rank 1 low-rank model
  data3 <- simulate_lowrank_dbn(n = 5, p = 1, time = 10, r = 1, seed = 666)
  expect_equal(ncol(data3$U), 1)
  expect_equal(nrow(data3$alpha), 1)
  
  # Two regimes HMM
  data4 <- simulate_hmm_dbn(n = 5, p = 1, time = 10, R = 2, seed = 777)
  expect_equal(length(data4$A_list), 2)
})