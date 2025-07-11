test_that("verbose parameter controls output in dbn() wrapper", {
  # Generate small test data
  Y <- simulate_test_data(n = 6, p = 1, time = 5, seed = 123)
  
  # Test verbose = FALSE suppresses output
  expect_silent({
    results <- dbn(Y, model = "static", nscan = 10, burn = 2, odens = 1, verbose = FALSE)
  })
  expect_s3_class(results, "dbn")
  
  # Test verbose = TRUE runs without error
  results <- dbn(Y, model = "static", nscan = 10, burn = 2, odens = 1, verbose = TRUE)
  expect_s3_class(results, "dbn")
})

test_that("verbose parameter is accepted by all model functions", {
  Y <- simulate_test_data(n = 6, p = 1, time = 5, seed = 456)
  
  # Test static model accepts verbose parameter
  expect_silent({
    results <- dbn_static(Y, nscan = 10, burn = 2, odens = 1, verbose = FALSE)
  })
  expect_type(results, "list")
  
  # Just verify the functions accept the verbose parameter without error
  # (actual model runs may fail due to data issues, but parameter should be accepted)
  
  # Check function signatures include verbose
  expect_true("verbose" %in% names(formals(dbn)))
  expect_true("verbose" %in% names(formals(dbn_static)))
  expect_true("verbose" %in% names(formals(dbn_dynamic)))
  expect_true("verbose" %in% names(formals(dbn_lowrank)))
  expect_true("verbose" %in% names(formals(dbn_hmm)))
})

test_that("verbose defaults are set correctly", {
  # Check default values
  expect_equal(formals(dbn)$verbose, TRUE)
  expect_equal(formals(dbn_static)$verbose, 100)
  expect_equal(formals(dbn_dynamic)$verbose, 100)
  expect_equal(formals(dbn_lowrank)$verbose, 100)
  expect_equal(formals(dbn_hmm)$verbose, 100)
})