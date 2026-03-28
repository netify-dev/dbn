#### piecewise model tests

test_that("parse_blocks works with integer", {
	block_info <- parse_blocks(4, Tt = 40, n = 10)
	expect_equal(block_info$K, 4)
	expect_equal(length(block_info$block_indices), 4)
	expect_equal(sum(block_info$lengths), 40)
})

test_that("parse_blocks works with vector", {
	block_info <- parse_blocks(c(10, 25, 40), Tt = 40, n = 10)
	expect_equal(block_info$K, 3)
	expect_equal(block_info$boundaries, c(0, 10, 25, 40))
})

test_that("parse_blocks works with named vector", {
	block_info <- parse_blocks(c(pre = 15, crisis = 30, post = 50), Tt = 50, n = 10)
	expect_equal(block_info$K, 3)
	expect_true("pre" %in% block_info$names)
})

test_that("parse_blocks validates input", {
	expect_error(parse_blocks(0, Tt = 40, n = 10))
	expect_error(parse_blocks(50, Tt = 40, n = 10))
})

test_that("simulate_piecewise_dbn generates correct dimensions", {
	sim <- simulate_piecewise_dbn(n = 10, time = 40, blocks = 4, p = 1, seed = 123)
	expect_equal(dim(sim$Y), c(10, 10, 1, 40))
	expect_equal(sim$block_info$K, 4)
	expect_equal(length(sim$true_A), 4)
	expect_equal(length(sim$true_B), 4)
})

test_that("simulate_piecewise_dbn with boundary vector", {
	sim <- simulate_piecewise_dbn(n = 8, time = 30, blocks = c(10, 20, 30), seed = 456)
	expect_equal(dim(sim$Y), c(8, 8, 1, 30))
	expect_equal(sim$block_info$K, 3)
})

test_that("dbn piecewise model runs with integer blocks", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 8, time = 30, blocks = 3, seed = 789)
	fit <- dbn(sim$Y, model = "piecewise", blocks = 3,
			   nscan = 50, burn = 25, verbose = FALSE)

	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "piecewise")
	expect_equal(fit$blocks$K, 3)
})

test_that("dbn piecewise model runs with boundary vector", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 8, time = 40, blocks = c(15, 30, 40), seed = 111)
	fit <- dbn(sim$Y, model = "piecewise", blocks = c(15, 30, 40),
			   nscan = 50, burn = 25, verbose = FALSE)

	expect_equal(fit$blocks$boundaries, c(0, 15, 30, 40))
})

test_that("piecewise model requires blocks parameter", {
	sim <- simulate_piecewise_dbn(n = 8, time = 30, blocks = 3, seed = 222)
	expect_error(dbn(sim$Y, model = "piecewise", verbose = FALSE), "blocks")
})

test_that("piecewise model requires sufficient time points", {
	Y <- array(rnorm(64), dim = c(8, 8, 1, 3))
	expect_error(dbn(Y, model = "piecewise", blocks = 2, verbose = FALSE), "4 time points")
})

test_that("piecewise summary works", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 6, time = 24, blocks = 3, seed = 333)
	fit <- dbn(sim$Y, model = "piecewise", blocks = 3,
			   nscan = 50, burn = 25, verbose = FALSE)

	expect_output(summary(fit), "Piecewise-Static")
	expect_output(summary(fit), "Block Structure")
})

test_that("compare_blocks works", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 6, time = 24, blocks = 3, seed = 6886)
	fit <- dbn(sim$Y, model = "piecewise", blocks = 3,
			   nscan = 50, burn = 25, verbose = FALSE)

	result <- compare_blocks(fit)
	expect_true(is.list(result))
	expect_equal(length(result), 2)

	# verify s3 class and clean print output
	expect_s3_class(result, "dbn_block_comparison")
	out <- capture.output(print(result), type = "message")
	expect_false(any(grepl("diff_norms", out)))

	# underlying data still accessible
	expect_true(length(result[[1]]$diff_norms) > 0)
})

test_that("compare_blocks requires piecewise model", {
	skip_on_cran()
	Y <- array(rnorm(64 * 10), dim = c(8, 8, 1, 10))
	fit_static <- dbn(Y, model = "static", nscan = 30, burn = 15, verbose = FALSE)
	expect_error(compare_blocks(fit_static), "piecewise")
})

test_that("piecewise print shows block info", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 6, time = 24, blocks = 3, seed = 555)
	fit <- dbn(sim$Y, model = "piecewise", blocks = 3,
			   nscan = 50, burn = 25, verbose = FALSE)

	output <- capture.output(print(fit))
	expect_true(any(grepl("PIECEWISE", output)))
	expect_true(any(grepl("Blocks:", output)))
})

test_that("piecewise model with gaussian family", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 6, time = 24, blocks = 2, seed = 666)
	fit <- dbn(sim$Y_continuous, model = "piecewise", blocks = 2,
			   family = "gaussian", nscan = 50, burn = 25, verbose = FALSE)

	expect_equal(fit$family, "gaussian")
	expect_s3_class(fit, "dbn")
})

test_that("piecewise model stores block-specific A and B", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 6, time = 20, blocks = 2, seed = 777)
	fit <- dbn(sim$Y, model = "piecewise", blocks = 2,
			   nscan = 50, burn = 25, verbose = FALSE)

	expect_equal(length(fit$A_blocks), 2)
	expect_equal(length(fit$B_blocks), 2)
	expect_equal(dim(fit$A_blocks[[1]]), c(6, 6))
})

test_that("piecewise draws have correct structure", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 6, time = 20, blocks = 2, seed = 888)
	fit <- dbn(sim$Y, model = "piecewise", blocks = 2,
			   nscan = 60, burn = 30, odens = 2, verbose = FALSE)

	expect_equal(fit$settings$draws, 30)
	expect_equal(nrow(fit$draws$pars), 30)
	expect_equal(length(fit$draws$misc$A), 30)
})

test_that("piecewise stores Theta draws by default for actor position inference", {
	skip_on_cran()
	n <- 6
	Tt <- 20
	sim <- simulate_piecewise_dbn(n = n, time = Tt, blocks = 2, seed = 999)
	fit <- dbn(sim$Y, model = "piecewise", blocks = 2,
			   nscan = 60, burn = 30, odens = 2, verbose = FALSE)

	# Theta should be stored by default
	expect_true("Theta" %in% names(fit$draws$misc))
	expect_equal(length(fit$draws$misc$Theta), 30)

	# each Theta draw should have correct dimensions
	expect_equal(dim(fit$draws$misc$Theta[[1]]), c(n, n, 1, Tt))

	# can access actor positions at specific time points
	theta_draw1 <- fit$draws$misc$Theta[[1]]
	actor_positions_t5 <- theta_draw1[, , 1, 5]
	expect_equal(dim(actor_positions_t5), c(n, n))
})

test_that("piecewise store_theta=FALSE reduces memory", {
	skip_on_cran()
	sim <- simulate_piecewise_dbn(n = 8, time = 30, blocks = 3, seed = 1001)

	# with Theta storage (default)
	fit_with_theta <- dbn(sim$Y, model = "piecewise", blocks = 3,
						  nscan = 50, burn = 25, verbose = FALSE,
						  store_theta = TRUE)

	# without Theta storage
	fit_no_theta <- dbn(sim$Y, model = "piecewise", blocks = 3,
						nscan = 50, burn = 25, verbose = FALSE,
						store_theta = FALSE)

	# verify Theta is stored/not stored
	expect_true("Theta" %in% names(fit_with_theta$draws$misc))
	expect_false("Theta" %in% names(fit_no_theta$draws$misc))

	# fit without Theta should be smaller
	size_with <- object.size(fit_with_theta)
	size_without <- object.size(fit_no_theta)
	expect_lt(as.numeric(size_without), as.numeric(size_with))
})
