# Test suite for posterior analysis methods

test_that("plot.dbn dispatches correctly", {
  # create mock static results
  static_results <- list(
    model = "static",
    params = matrix(rnorm(30), ncol = 3),
    B = list(array(rnorm(64), dim = c(4, 4, 10)))
  )
  colnames(static_results$params) <- c("s2", "t2", "g2")
  class(static_results) <- "dbn"
  
  # test static dispatch
  expect_silent(p <- plot(static_results))
  # plot_static returns a gtable/grob, not a ggplot
  expect_true(inherits(p, c("gtable", "gTree", "grob")))
  
  # create mock dynamic results
  n <- 10
  k <- 3
  Tt <- 20
  nscan <- 10
  dynamic_results <- list(
    model = "dynamic",
    sigma2 = rnorm(nscan, 1, 0.1),
    tauA2 = rnorm(nscan, 0.5, 0.05),
    tauB2 = rnorm(nscan, 0.5, 0.05),
    g2 = rnorm(nscan, 0.3, 0.03),
    A = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
    B = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
    dims = list(m = n, Tt = Tt)
  )
  class(dynamic_results) <- "dbn"
  
  # test dynamic dispatch
  expect_silent(plots <- plot(dynamic_results))
  expect_type(plots, "list")
})

test_that("tidy_dbn extracts posterior means correctly", {
  # create mock results for dynamic model (since tidy_dbn only extracts a/b for dynamic)
  n <- 5
  k <- 3
  Tt <- 10
  nscan <- 20
  
  # dynamic model structure
  results <- list(
    model = "dynamic",
    A = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
    B = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
    M = array(rnorm(n*n*3), dim = c(n, n, 3)),
    dims = list(m = n, p = 3)
  )
  class(results) <- "dbn"
  
  # test extraction
  tidy <- tidy_dbn(results)
  
  expect_type(tidy, "list")
  expect_equal(dim(tidy$A), c(n, k, Tt))
  expect_equal(dim(tidy$B), c(n, k, Tt))
  
  # for static model, only theta (which is m) should be extracted
  static_results <- list(
    model = "static",
    M = array(rnorm(n*n*3), dim = c(n, n, 3)),
    B = list(array(rnorm(n*n*10), dim = c(n, n, 10))),
    dims = list(m = n, p = 3)
  )
  class(static_results) <- "dbn"
  
  tidy_static <- tidy_dbn(static_results)
  expect_equal(dim(tidy_static$Theta), c(n, n, 3))
})

test_that("dyad_path handles multiple relations correctly", {
  # create mock dynamic results with multiple relations
  n <- 10
  p <- 3
  Tt <- 20
  k <- 4
  nscan <- 50
  
  results <- list(
    model = "dynamic",
    A = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
    B = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
    dims = list(m = n, p = p, Tt = Tt)
  )
  class(results) <- "dbn"
  
  # test single relation
  expect_silent(p1 <- dyad_path(results, i = 1, j = 5, rel = 1))
  expect_s3_class(p1, "ggplot")
  
  # test multiple relations with faceting
  expect_silent(p2 <- dyad_path(results, i = 1, j = 5, rel = c(1, 2), facet = TRUE))
  expect_s3_class(p2, "ggplot")
  
  # test all relations
  expect_silent(p3 <- dyad_path(results, i = 1, j = 5, rel = NULL))
  expect_s3_class(p3, "ggplot")
})

test_that("ppc_ecdf guards against out-of-range values", {
  # create results with continuous predictions
  n <- 8
  Tt <- 5
  nscan <- 50
  n_iter <- nscan  # define n_iter for B array dimensions
  
  # mock data with ordinal values 1-5
  Y <- array(sample(1:5, n*n*Tt, replace = TRUE), dim = c(n, n, Tt))
  diag(Y[,,1]) <- NA
  
  results <- list(
    model = "static",
    R = Y,  # ppc_ecdf expects r, not y
    M = array(rnorm(n*n, mean = 3, sd = 0.5), dim = c(n, n)),  # baseline mean for static model
    B = list(
      array(rnorm(n*n*n_iter), dim = c(n, n, n_iter)),
      array(rnorm(n*n*n_iter), dim = c(n, n, n_iter))
    ),
    dims = list(m = n, p = 1, n = Tt)
  )
  class(results) <- "dbn"
  
  # should not error even with extreme values
  expect_silent(p <- ppc_ecdf(results, n_rep = 20))
  expect_s3_class(p, "ggplot")
})

test_that("predict.dbn works for static and dynamic models", {
  # test static model
  n <- 6
  k <- 3
  p <- 2
  Tt <- 10
  nscan <- 50
  
  Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  for(t in 1:Tt) {
    for(r in 1:p) {
      diag(Y[,,r,t]) <- NA
    }
  }
  
  static_results <- list(
    model = "static",
    Y = Y,
    R = Y,  # add r field
    M = array(rnorm(n*n*p, mean = 0, sd = 0.5), dim = c(n, n, p)),  # baseline mean
    B = list(
      array(rnorm(n*n*nscan), dim = c(n, n, nscan)),  # b[[1]] has dims [n, n, nscan]
      array(rnorm(n*n*nscan), dim = c(n, n, nscan)),  # b[[2]] has dims [n, n, nscan]
      array(rnorm(p*p*nscan), dim = c(p, p, nscan))   # b[[3]] has dims [p, p, nscan]
    ),
    params = matrix(runif(nscan * 3, 0.1, 1), ncol = 3),
    dims = list(m = n, p = p, n = Tt)
  )
  colnames(static_results$params) <- c("s2", "t2", "g2")
  class(static_results) <- "dbn"
  
  # test prediction
  pred <- predict(static_results, S = 10)
  # check dimensions are correct
  expect_equal(length(dim(pred)), 5)
  expect_equal(prod(dim(pred)), n * n * p * Tt * 10)
  
  # test dynamic model
  dynamic_results <- list(
    model = "dynamic",
    Y = Y,
    A = lapply(1:nscan, function(i) array(rnorm(n*n*Tt), dim = c(n, n, Tt))),
    B = lapply(1:nscan, function(i) array(rnorm(n*n*Tt), dim = c(n, n, Tt))),
    Theta = lapply(1:nscan, function(i) array(rnorm(n*n*p*Tt), dim = c(n, n, p, Tt))),
    sigma2 = runif(nscan, 0.1, 1),
    dims = list(m = n, p = p, Tt = Tt)
  )
  class(dynamic_results) <- "dbn"
  
  # test forecasting
  H <- 5
  forecast <- predict(dynamic_results, H = H, S = 10)
  expect_equal(dim(forecast), c(n, n, p, H, 10))
})

test_that("net_snapshot handles sparse networks", {
  # small network - should use regular heatmap
  n_small <- 20
  Tt <- 10
  
  results_small <- list(
    model = "dynamic",  # net_snapshot requires dynamic model
    A = list(array(rnorm(n_small*3*Tt), dim = c(n_small, 3, Tt))),
    B = list(array(rnorm(n_small*3*Tt), dim = c(n_small, 3, Tt))),
    dims = list(m = n_small, p = 1, Tt = Tt)
  )
  class(results_small) <- "dbn"
  
  p_small <- net_snapshot(results_small, t = 5)
  expect_s3_class(p_small, "ggplot")
  
  # large network - should use sparse visualization
  n_large <- 200
  results_large <- list(
    model = "dynamic",  # net_snapshot requires dynamic model
    A = list(array(rnorm(n_large*5*Tt, sd = 0.1), dim = c(n_large, 5, Tt))),
    B = list(array(rnorm(n_large*5*Tt, sd = 0.1), dim = c(n_large, 5, Tt))),
    dims = list(m = n_large, p = 1, Tt = Tt)
  )
  class(results_large) <- "dbn"
  
  # suppress expected warning about no edges above threshold
  suppressWarnings({
    p_large <- net_snapshot(results_large, t = 5, eps = 0.5)
  })
  expect_s3_class(p_large, "ggplot")
})

test_that("group influence functions work correctly", {
  # create dynamic results with correct dimensions
  n <- 20
  Tt <- 15
  nscan <- 50
  
  results <- list(
    model = "dynamic",
    A = lapply(1:nscan, function(i) array(rnorm(n*n*Tt), dim = c(n, n, Tt))),
    B = lapply(1:nscan, function(i) array(rnorm(n*n*Tt), dim = c(n, n, Tt))),
    dims = list(m = n, p = 1, Tt = Tt)
  )
  class(results) <- "dbn"
  
  # test plot_group_influence
  group <- c(1, 3, 5, 7, 9)
  
  # test different measures
  for(measure in c("rowsum", "rowmean", "l2")) {
    p <- plot_group_influence(results, group = group, 
                             type = "sender", measure = measure)
    expect_s3_class(p, "ggplot")
  }
  
  # test get_group_influence
  influence_data <- get_group_influence(results, group = group, 
                                       type = "target", measure = "rowsum")
  expect_type(influence_data, "list")
  # check that it has time, mean, and quantile columns
  expect_true("time" %in% names(influence_data))
  expect_true("mean" %in% names(influence_data))
  expect_true("q0.025" %in% names(influence_data))
  expect_true("q0.5" %in% names(influence_data))
  expect_true("q0.975" %in% names(influence_data))
  
  # test compare_group_influence
  groups <- list(
    "Group A" = 1:3,
    "Group B" = 4:6,
    "Group C" = 7:10
  )
  
  p_comp <- compare_group_influence(results, groups = groups,
                                   type = "sender", measure = "l2")
  expect_s3_class(p_comp, "ggplot")
})

test_that("check_convergence handles missing coda gracefully", {
  # create results with params
  results <- list(
    model = "static",
    params = matrix(rnorm(300), ncol = 3)
  )
  colnames(results$params) <- c("s2", "t2", "g2")
  class(results) <- "dbn"
  
  # should work even without coda
  output <- capture.output(check_convergence(results))
  expect_true(length(output) > 0)
})

test_that("role_trajectory works for both sender and receiver", {
  # Create dynamic results
  n <- 15
  k <- 5
  Tt <- 20
  nscan <- 50
  
  results <- list(
    model = "dynamic",
    A = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
    B = lapply(1:nscan, function(i) array(rnorm(n*k*Tt), dim = c(n, k, Tt))),
    dims = list(m = n, p = 1, Tt = Tt)
  )
  class(results) <- "dbn"
  
  # test first component of a matrix (senders)
  # role_trajectory uses base graphics and returns null
  expect_silent(role_trajectory(results, mat = "A", comp = 1))
  
  # test first component of b matrix (receivers)
  expect_silent(role_trajectory(results, mat = "B", comp = 1))
})

test_that("predict_ordinal converts correctly", {
  # create results
  n <- 8
  p <- 2
  Tt <- 10
  
  Y <- array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
  for(t in 1:Tt) {
    for(r in 1:p) {
      diag(Y[,,r,t]) <- NA
    }
  }
  
  # create results with proper number of posterior samples
  n_samples <- 10
  results <- list(
    model = "static", 
    Y = Y,
    R = Y,  # Ordinal data
    M = array(3, dim = c(n, n, p)),  # baseline mean
    B = list(
      array(rnorm(n*n*n_samples), dim = c(n, n, n_samples)),
      array(rnorm(n*n*n_samples), dim = c(n, n, n_samples)),
      array(rnorm(p*p*n_samples), dim = c(p, p, n_samples))
    ),
    params = matrix(runif(n_samples * 3, 0.1, 1), ncol = 3),  # n_samples rows
    dims = list(m = n, p = p, n = Tt)
  )
  colnames(results$params) <- c("s2", "t2", "g2")
  class(results) <- "dbn"
  
  # test conversion
  ordinal_pred <- predict_ordinal(results, draws = 20)
  
  # check all values are in 1:5
  vals <- unique(c(ordinal_pred))
  vals <- vals[!is.na(vals)]
  expect_true(all(vals %in% 1:5))
})

# test_that("summary.dbn works for both model types", {
#   # Static model
#   static_results <- list(
#     model = "static",
#     params = matrix(rnorm(300), ncol = 3),
#     nscan = 50,
#     dims = list(m = 10, p = 2, n = 20),
#     settings = list(nscan = 50, burn = 10, odens = 1)
#   )
#   colnames(static_results$params) <- c("s2", "t2", "g2")
#   class(static_results) <- "dbn"
#   
#   # Summary should return invisibly
#   expect_invisible(summary(static_results))
#   
#   # Dynamic model
#   dynamic_results <- list(
#     model = "dynamic",
#     sigma2 = rnorm(100, 1, 0.1),
#     tauA2 = rnorm(100, 0.5, 0.05),  # Use tauA2 not tau2_A
#     tauB2 = rnorm(100, 0.5, 0.05),  # Use tauB2 not tau2_B
#     nscan = 50,
#     dims = list(m = 10, p = 2, Tt = 20),
#     settings = list(nscan = 100, burn = 20, odens = 1, ar1 = FALSE)
#   )
#   class(dynamic_results) <- "dbn"
#   
#   # Summary should return invisibly
#   expect_invisible(summary(dynamic_results))
# })