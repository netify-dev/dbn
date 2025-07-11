# # Test suite for new posterior analysis layer

# test_that("theta_slice extracts correct values", {
#   # create mock dynamic results with new format
#   n <- 5
#   p <- 2
#   Tt <- 10
#   n_draws <- 20
  
#   # create theta arrays
#   theta_arrays <- lapply(1:n_draws, function(s) {
#     array(rnorm(n*n*p*Tt), dim = c(n, n, p, Tt))
#   })
  
#   fit <- list(
#     draws = list(
#       theta = theta_arrays,
#       pars = data.frame(
#         sigma2_proc = runif(n_draws, 0.5, 1.5),
#         tau_A2 = runif(n_draws, 0.1, 1),
#         tau_B2 = runif(n_draws, 0.1, 1)
#       )
#     ),
#     meta = list(
#       dims = list(m = n, p = p, Tt = Tt),
#       draws = n_draws,
#       model = "dynamic"
#     ),
#     model = "dynamic"
#   )
#   class(fit) <- "dbn"
  
#   # test extraction of specific dyad
#   slice1 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = 5)
#   expect_equal(length(slice1), n_draws)
  
#   # test extraction of multiple time points
#   slice2 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = c(1, 5, 10))
#   expect_equal(dim(slice2), c(n_draws, 3))
  
#   # test extraction with draws subset
#   slice3 <- theta_slice(fit, i = 1, j = 2, rel = 1, time = 5, draws = 1:5)
#   expect_equal(length(slice3), 5)
# })

# test_that("theta_summary computes summaries correctly", {
#   # create mock fit object
#   n <- 4
#   p <- 1
#   Tt <- 8
#   n_draws <- 50
  
#   theta_arrays <- lapply(1:n_draws, function(s) {
#     array(rnorm(n*n*p*Tt, mean = s/n_draws), dim = c(n, n, p, Tt))
#   })
  
#   fit <- list(
#     draws = list(
#       theta = theta_arrays
#     ),
#     meta = list(
#       dims = list(m = n, p = p, Tt = Tt),
#       draws = n_draws,
#       model = "dynamic"
#     ),
#     model = "dynamic"
#   )
#   class(fit) <- "dbn"
  
#   # test mean summary
#   summary1 <- theta_summary(fit, fun = mean, i = 1, j = 2, rel = 1, time = 1)
#   expect_equal(nrow(summary1), 1)
#   expect_true("value" %in% names(summary1))
  
#   # test quantile summary
#   summary2 <- theta_summary(fit, fun = function(x) quantile(x, c(0.025, 0.5, 0.975)))
#   expect_true(nrow(summary2) > 0)
  
#   # test with subset
#   summary3 <- theta_summary(fit, fun = mean, time = 1:3)
#   expect_true(all(summary3$time %in% 1:3))
# })

# test_that("param_summary extracts parameters correctly", {
#   n_draws <- 100
  
#   # test with new format
#   fit_new <- list(
#     draws = list(
#       pars = data.frame(
#         sigma2_proc = runif(n_draws, 0.5, 1.5),
#         tau_A2 = runif(n_draws, 0.1, 1),
#         tau_B2 = runif(n_draws, 0.1, 1),
#         g2 = runif(n_draws, 0.5, 2)
#       )
#     ),
#     model = "dynamic"
#   )
#   class(fit_new) <- "dbn"
  
#   summary_new <- param_summary(fit_new)
#   expect_equal(nrow(summary_new), 4)
#   expect_true(all(c("parameter", "mean", "q5", "q50", "q95") %in% names(summary_new)))
  
#   # test with legacy format
#   fit_legacy <- list(
#     sigma2 = runif(n_draws, 0.5, 1.5),
#     tau_A2 = runif(n_draws, 0.1, 1),
#     tau_B2 = runif(n_draws, 0.1, 1),
#     model = "dynamic"
#   )
#   class(fit_legacy) <- "dbn"
  
#   summary_legacy <- param_summary(fit_legacy)
#   expect_true(nrow(summary_legacy) >= 3)
# })

# test_that("latent_summary extracts M arrays correctly", {
#   n <- 5
#   p <- 2
#   n_draws <- 30
  
#   # create m arrays
#   M_arrays <- lapply(1:n_draws, function(s) {
#     array(rnorm(n*n*p), dim = c(n, n, p))
#   })
  
#   fit <- list(
#     draws = list(
#       misc = list(M = M_arrays)
#     ),
#     meta = list(
#       dims = list(m = n, p = p),
#       draws = n_draws,
#       model = "dynamic"
#     ),
#     model = "dynamic"
#   )
#   class(fit) <- "dbn"
  
#   # test extraction
#   summary1 <- latent_summary(fit, fun = mean)
#   expect_equal(nrow(summary1), n * n * p)
  
#   # test specific relation
#   summary2 <- latent_summary(fit, fun = mean, rel = 1)
#   expect_equal(nrow(summary2), n * n)
#   expect_true(all(summary2$rel == 1))
# })

# test_that("regime_probs works for HMM models", {
#   Tt <- 20
#   R <- 3
#   n_draws = 50
  
#   # create mock hmm results
#   S_arrays <- lapply(1:n_draws, function(s) {
#     sample(1:R, Tt, replace = TRUE)
#   })
  
#   fit <- list(
#     draws = list(
#       misc = list(S = S_arrays)
#     ),
#     meta = list(
#       dims = list(Tt = Tt),
#       R = R,
#       model = "hmm"
#     ),
#     settings = list(R = R),
#     model = "hmm"
#   )
#   class(fit) <- "dbn"
  
#   # test extraction
#   probs <- regime_probs(fit)
#   expect_equal(dim(probs), c(Tt, R))
#   expect_true(all(abs(rowSums(probs) - 1) < 1e-10))
  
#   # test non-hmm model
#   fit_dynamic <- list(model = "dynamic")
#   class(fit_dynamic) <- "dbn"
#   expect_null(regime_probs(fit_dynamic))
# })

# test_that("posterior_predict_dbn generates correct predictions", {
#   n <- 6
#   p <- 2
#   Tt <- 10
#   n_draws_fit <- 30
  
#   # create mock ordinal data
#   Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
#   for(t in 1:Tt) {
#     for(r in 1:p) {
#       diag(Y[,,r,t]) <- NA
#     }
#   }
  
#   # create theta and z arrays for ordinal family
#   theta_arrays <- lapply(1:n_draws_fit, function(s) {
#     array(rnorm(n*n*p*Tt), dim = c(n, n, p, Tt))
#   })
  
#   z_arrays <- lapply(1:n_draws_fit, function(s) {
#     array(rnorm(n*n*p*Tt), dim = c(n, n, p, Tt))
#   })
  
#   # create m arrays
#   M_arrays <- lapply(1:n_draws_fit, function(s) {
#     array(rnorm(n*n*p, mean = 3), dim = c(n, n, p))
#   })
  
#   # For dynamic model, create Theta as 5D array (m x m x p x time x iter)
#   Theta_array <- array(NA, dim = c(n, n, p, Tt, n_draws_fit))
#   for (s in 1:n_draws_fit) {
#     Theta_array[,,,,s] <- theta_arrays[[s]]
#   }
  
#   # Create A and B arrays for dynamic model
#   A_arrays <- lapply(1:n_draws_fit, function(s) {
#     A_arr <- array(0, dim = c(n, n, Tt))
#     for (t in 1:Tt) {
#       A_arr[,,t] <- diag(n)  # Identity matrices
#     }
#     A_arr
#   })
  
#   B_arrays <- lapply(1:n_draws_fit, function(s) {
#     B_arr <- array(0, dim = c(n, n, Tt))
#     for (t in 1:Tt) {
#       B_arr[,,t] <- diag(n)  # Identity matrices
#     }
#     B_arr
#   })
  
#   # Create M as 4D array for dynamic model
#   M_4d <- array(NA, dim = c(n, n, p, n_draws_fit))
#   for (s in 1:n_draws_fit) {
#     M_4d[,,,s] <- M_arrays[[s]]
#   }
  
#   fit <- list(
#     # Theta = Theta_array,  # No longer stored to save memory
#     A = A_arrays,  # Dynamic model requires these
#     B = B_arrays,
#     M = M_4d,
#     sigma2 = rep(1, n_draws_fit),  # Dynamic model uses this to count draws
#     draws = list(
#       theta = theta_arrays,
#       z = z_arrays,
#       misc = list(M = M_arrays)
#     ),
#     meta = list(
#       dims = list(m = n, p = p, Tt = Tt),
#       draws = n_draws_fit,
#       model = "dynamic"
#     ),
#     dims = list(m = n, p = p, Tt = Tt),  # Add dims at top level too
#     family = family_ordinal(),
#     Y = Y,
#     model = "dynamic"
#   )
#   class(fit) <- "dbn"
  
#   # test prediction
#   n_pred_draws <- 20
#   yrep <- posterior_predict_dbn(fit, ndraws = n_pred_draws)
  
#   expect_equal(length(yrep), n_pred_draws)
#   # Check if first prediction exists and has correct dimensions
#   expect_false(is.null(yrep[[1]]))
#   expect_equal(dim(yrep[[1]]), dim(Y))
  
#   # check predictions are in valid range
#   vals <- unique(unlist(yrep))
#   vals <- vals[!is.na(vals)]
#   expect_true(all(vals %in% 1:5))
# })

# test_that("plot functions create valid plots", {
#   # create minimal fit object
#   n <- 5
#   p <- 1
#   Tt <- 10
#   n_draws <- 50
  
#   fit <- list(
#     draws = list(
#       theta = lapply(1:n_draws, function(s) array(rnorm(n*n*p*Tt), dim = c(n, n, p, Tt))),
#       pars = data.frame(
#         sigma2_proc = runif(n_draws, 0.5, 1.5),
#         tau_A2 = runif(n_draws, 0.1, 1)
#       ),
#       misc = list(
#         M = lapply(1:n_draws, function(s) array(rnorm(n*n*p), dim = c(n, n, p)))
#       )
#     ),
#     meta = list(
#       dims = list(m = n, p = p, Tt = Tt),
#       draws = n_draws,
#       model = "dynamic"
#     ),
#     model = "dynamic"
#   )
#   class(fit) <- "dbn"
  
#   # test trace plot
#   expect_silent(plot_trace(fit, pars = c("sigma2_proc", "tau_A2")))
  
#   # test theta plot
#   expect_silent(plot_theta(fit, time = 5, rel = 1))
  
#   # test dyad path (should create ggplot if available)
#   if(requireNamespace("ggplot2", quietly = TRUE)) {
#     p <- plot_dyad_path(fit, i = 1, j = 2, rel = 1)
#     expect_s3_class(p, "ggplot")
#   }
# })

# test_that("derive_draws computes derived quantities", {
#   n <- 4
#   n_draws <- 20
  
#   # create fit with a and b arrays
#   A_arrays <- lapply(1:n_draws, function(s) {
#     A <- diag(n) + matrix(rnorm(n*n, 0, 0.1), n, n)
#     A
#   })
  
#   B_arrays <- lapply(1:n_draws, function(s) {
#     B <- diag(n) + matrix(rnorm(n*n, 0, 0.1), n, n) 
#     B
#   })
  
#   fit <- list(
#     draws = list(
#       misc = list(
#         A = lapply(1:n_draws, function(s) array(A_arrays[[s]], dim = c(n, n, 1))),
#         B = lapply(1:n_draws, function(s) array(B_arrays[[s]], dim = c(n, n, 1)))
#       )
#     ),
#     meta = list(
#       dims = list(m = n),
#       draws = n_draws,
#       model = "dynamic"
#     ),
#     model = "dynamic"
#   )
#   class(fit) <- "dbn"
  
#   # define function to compute persistence (diagonal of a)
#   persistence_fn <- function(draw) {
#     A <- draw$A[,,1]
#     diag(A)
#   }
  
#   # compute derived quantity
#   persistence <- derive_draws(fit, persistence_fn)
  
#   expect_equal(dim(persistence), c(n_draws, n))
#   expect_true(all(is.finite(persistence)))
# })

# # test_that("extraction functions handle missing data gracefully", {
# #   # Create fit with some missing components
# #   fit <- list(
# #     draws = list(
# #       pars = data.frame(sigma2 = runif(10))
# #     ),
# #     meta = list(
# #       dims = list(m = 5),
# #       model = "static"
# #     ),
# #     model = "static"
# #   )
# #   class(fit) <- "dbn"
# #   
# #   # Should return NULL for missing theta
# #   expect_null(theta_slice(fit, i = 1, j = 2))
# #   
# #   # Should still work for parameters
# #   param_sum <- param_summary(fit)
# #   expect_equal(nrow(param_sum), 1)
# #   
# #   # Should return NULL for missing M
# #   expect_null(latent_summary(fit))
# # })