test_that("models handle missing data", {
  
  # create data with missing values
  set.seed(6886)
  Y <- array(sample(1:5, 5*5*2*10, replace = TRUE), dim = c(5, 5, 2, 10))
  
  # add 10% missing data randomly
  n_missing <- round(0.1 * length(Y))
  Y[sample(length(Y), n_missing)] <- NA
  
  # test static model with missing data
  expect_no_error({
    res_static <- dbn(Y, model = "static", nscan = 50, burn = 10, odens = 1, verbose = FALSE)
  })
  
  # test dynamic model with missing data  
  expect_no_error({
    res_dynamic <- dbn(Y, model = "dynamic", nscan = 50, burn = 10, odens = 1, verbose = FALSE)
  })
  
  # test low-rank model with smaller network
  Y_small <- array(sample(1:5, 3*3*1*10, replace = TRUE), dim = c(3, 3, 1, 10))
  Y_small[sample(length(Y_small), 5)] <- NA
  
  expect_no_error({
    res_lr <- dbn_lowrank(Y_small, r = 2, nscan = 50, burn = 10, odens = 1, verbose = FALSE)
  })
  
  # test hmm model with missing data
  expect_no_error({
    res_hmm <- dbn_hmm(Y, R = 2, nscan = 50, burn = 10, odens = 1, verbose = FALSE)
  })
})

test_that("rank likelihood handles missing values correctly", {
  
  # create rank data with na values
  R <- matrix(c(1, 2, NA, 3, 4, NA, 2, 1, 3), 3, 3)
  
  # verify zscores function handles nas
  z <- zscores(R)
  expect_equal(sum(is.na(z)), sum(is.na(R)))
  expect_true(all(!is.na(z[!is.na(R)])))
  
  # verify rank index building
  IR <- precompute_ranks(array(R, dim = c(3, 3, 1, 1)))
  expect_true("NA" %in% names(IR[[1]]))
  expect_equal(length(IR[[1]][["NA"]]), sum(is.na(R)))
  
  # verify rz_fc handles missing values
  Z <- matrix(rnorm(9), 3, 3)
  EZ <- matrix(rnorm(9), 3, 3)
  Z_new <- rz_fc(R, Z, EZ, IR[[1]])
  
  expect_equal(dim(Z_new), dim(R))
  expect_true(all(is.finite(Z_new)))
})

test_that("missing data patterns are preserved", {
  
  # create specific missing data pattern
  Y <- array(1:40, dim = c(2, 2, 2, 5))
  Y[1, 1, 1, 1:2] <- NA
  Y[2, 2, 2, 5] <- NA
  
  # verify preprocessing preserves missing pattern
  pre <- shared_preprocess(Y)
  
  # original na positions should map to finite values in z
  expect_true(all(is.finite(pre$Z)))
  
  # ir should track na positions correctly
  expect_equal(length(pre$IR[[1]][["NA"]]), 2)
  expect_equal(length(pre$IR[[2]][["NA"]]), 1)
})