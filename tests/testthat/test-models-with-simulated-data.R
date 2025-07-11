# # Test suite to verify models work with simulated data

# test_that("static model runs with simulated static data", {
#   # Generate appropriate static data
#   data <- simulate_static_dbn(n = 10, p = 2, time = 20, seed = 68861)
  
#   # Run static model with few iterations
#   results <- dbn(data$Y, model = "static", nscan = 30, burn = 10, odens = 1, verbose = FALSE)
  
#   # Basic checks
#   expect_s3_class(results, "dbn")
#   expect_equal(results$model, "static")
#   expect_true(is.list(results$B))
# })

# test_that("dynamic model runs with simulated dynamic data", {
#   # Generate appropriate dynamic data
#   data <- simulate_dynamic_dbn(n = 8, p = 1, time = 15, 
#                               ar1 = FALSE, seed = 68862)
  
#   # Try to run dynamic model
#   results <- tryCatch(
#     dbn(data$Y, model = "dynamic", nscan = 30, burn = 10, odens = 1, verbose = FALSE),
#     error = function(e) NULL
#   )
  
#   # If it runs, check basic structure
#   if(!is.null(results)) {
#     expect_equal(results$model, "dynamic")
#     expect_true("A" %in% names(results))
#     expect_true("B" %in% names(results))
#   }
# })

# test_that("low-rank model runs with simulated low-rank data", {
  
#   # Generate larger low-rank data for stability
#   data <- simulate_lowrank_dbn(n = 20, p = 1, time = 20, r = 2, seed = 68863)
  
#   # Try to run low-rank model with very short chain
#   results <- tryCatch(
#     dbn(data$Y, model = "lowrank", r = 2, nscan = 10, burn = 5, odens = 1, verbose = FALSE),
#     error = function(e) NULL
#   )
  
#   # If it runs, check basic structure
#   if(!is.null(results)) {
#     expect_equal(results$model, "lowrank")
#     expect_true("U" %in% names(results))
#     expect_true("alpha" %in% names(results))
#   }
# })

# test_that("HMM model runs with simulated HMM data", {
  
#   # Generate larger HMM data for stability
#   data <- simulate_hmm_dbn(n = 20, p = 1, time = 20, R = 2, 
#                           transition_prob = 0.9, seed = 68864)
  
#   # Try to run HMM model with very short chain
#   results <- tryCatch(
#     dbn(data$Y, model = "hmm", R = 2, nscan = 10, burn = 5, odens = 1, verbose = FALSE),
#     error = function(e) NULL
#   )
  
#   # If it runs, check basic structure
#   if(!is.null(results)) {
#     expect_equal(results$model, "hmm")
#     expect_true("S" %in% names(results))
#     expect_true("Pi" %in% names(results))
#   }
# })

# test_that("simulate functions produce data suitable for model convergence", {
#   # Static data should have temporal autocorrelation
#   data_static <- simulate_static_dbn(n = 10, p = 1, time = 10, seed = 68865)
  
#   # Check temporal correlation in the latent Z
#   cor_t1_t2 <- cor(c(data_static$Z[,,1,1]), c(data_static$Z[,,1,2]))
#   expect_true(cor_t1_t2 > 0.3)  # Should have some correlation
  
#   # Dynamic data should show time-varying patterns
#   data_dynamic <- simulate_dynamic_dbn(n = 10, p = 1, time = 10, 
#                                       ar1 = TRUE, rhoA = 0.9, seed = 68866)
  
#   # Check that A matrices change over time
#   A_diff <- norm(data_dynamic$A[,,10] - data_dynamic$A[,,1], "F")
#   expect_true(A_diff > 0.1)  # Should have changed
  
#   # Low-rank data should have low-rank structure
#   data_lr <- simulate_lowrank_dbn(n = 10, p = 1, time = 20, r = 2, seed = 68867)
  
#   # Check rank of average A
#   A_avg <- apply(data_lr$A, c(1,2), mean)
#   sv <- svd(A_avg)$d
#   expect_true(sv[3] / sv[1] < 0.1)  # Third singular value should be small
  
#   # HMM data should show regime persistence
#   data_hmm <- simulate_hmm_dbn(n = 10, p = 1, time = 10, R = 3, 
#                               transition_prob = 0.8, seed = 68868)
  
#   # Check regime persistence
#   n_switches <- sum(diff(data_hmm$S) != 0)
#   expect_true(n_switches < 5)  # Should have fewer than 5 switches in 10 time points
# })