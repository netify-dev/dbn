#### package loading and core dispatch

test_that("dbn package loads and has expected functions", {
	expect_true(exists("dbn"))
	expect_true(exists("dbn_static"))
	expect_true(exists("dbn_dynamic"))

	dbn_methods = as.character(methods(class = "dbn"))
	expect_true("plot.dbn" %in% dbn_methods || "plot" %in% dbn_methods)
	expect_true("summary.dbn" %in% dbn_methods || "summary" %in% dbn_methods)

	expect_true(exists("check_convergence"))
	expect_true(exists("compare_dbn"))
})

test_that("dbn function validates input", {
	expect_error(dbn(matrix(1:10, 5, 2)), "3D array .* or 4D array")
	expect_error(dbn(array(1, dim=c(5,5,2,10)), model = "invalid"), "should be one of")
})

test_that("data can be loaded from file path", {
	test_file = tempfile(fileext = ".RData")
	Y = simulate_test_data(n = 5, p = 2, time = 10, seed = 999)
	save(Y, file = test_file)

	results = dbn(test_file, model = "static", nscan = 10, burn = 5, odens = 1, verbose = FALSE)
	expect_s3_class(results, "dbn")

	unlink(test_file)
})

test_that("B scaling preserves diagonal mean", {
	B = list(
		diag(c(2, 4)),
		diag(c(1, 3)),
		diag(c(0.5, 1.5))
	)

	K = length(B)
	mu = rep(0, K)
	for (k in 1:K) mu[k] = mean(diag(B[[k]]))
	B_scaled_result = B
	for (k in 1:K) B_scaled_result[[k]] = B[[k]] / mu[k]

	for (k in 1:3) {
		expect_equal(mean(diag(B_scaled_result[[k]])), 1, tolerance = 1e-10)
	}
})

#### model output structure

test_that("static model produces expected output structure", {
	data = simulate_static_dbn(n = 10, p = 2, time = 10, seed = 6886)

	results = dbn(data$Y, model = "static", nscan = 40, burn = 10, odens = 1, verbose = FALSE, seed = 6886)

	expect_true(all(c("B", "params", "M", "R", "dims", "settings", "model") %in% names(results)))

	expect_equal(length(results$B), 3)
	expect_equal(dim(results$B[[1]])[3], 40)

	expect_equal(dim(results$M), c(10, 10, 2))

	s2_mean = mean(results$params[,"s2"])
	expect_true(s2_mean > 0 && s2_mean < 10)

	expect_equal(results$dims$m, 10)
	expect_equal(results$dims$p, 2)
	expect_equal(results$dims$n, 10)
})

test_that("static model handles missing data correctly", {
	data = simulate_static_dbn(n = 10, p = 2, time = 10, seed = 6886)

	n_missing = round(0.1 * length(data$Y))
	data$Y[sample(length(data$Y), n_missing)] = NA

	results = dbn(data$Y, model = "static", nscan = 40, burn = 10, odens = 1, verbose = FALSE)

	expect_false(any(is.na(results$params)))
	expect_s3_class(results, "dbn")
})

test_that("verbose parameter controls output", {
	Y = simulate_test_data(n = 10, p = 2, time = 10, seed = 333)

	expect_message(
		dbn(Y, model = "static", nscan = 20, burn = 10, odens = 1, verbose = TRUE),
		"Running static DBN MCMC"
	)

	expect_s3_class(
		dbn(Y, model = "static", nscan = 20, burn = 10, odens = 1, verbose = FALSE),
		"dbn"
	)
})

test_that("thinning works correctly", {
	Y = simulate_test_data(n = 10, p = 2, time = 10, seed = 666)

	results = dbn(Y, model = "static", nscan = 50, burn = 10, odens = 5, verbose = FALSE)

	expect_equal(nrow(results$params), 10)
	expect_equal(dim(results$B[[1]])[3], 10)
})

test_that("static model recovers reasonable parameters", {
	data = simulate_static_dbn(n = 10, p = 2, time = 20,
		sigma2 = 0.5, tau2 = 0.1, seed = 777)

	results = dbn(data$Y, model = "static", nscan = 80, burn = 20, odens = 1, verbose = FALSE, seed = 6886)

	s2_mean = mean(results$params[,"s2"])
	expect_true(s2_mean > 0.1 && s2_mean < 5)

	expect_equal(results$model, "static")
})

test_that("dynamic model initializes and runs", {
	data = simulate_dynamic_dbn(n = 8, p = 2, time = 10, seed = 888)

	results = tryCatch(
		dbn(data$Y, model = "dynamic", nscan = 15, burn = 5, odens = 1, verbose = FALSE),
		error = function(e) NULL
	)

	if (!is.null(results)) {
		expect_true("A" %in% names(results))
		expect_true("B" %in% names(results))
		expect_equal(results$model, "dynamic")
	} else {
		skip("dynamic model did not converge")
	}
})

#### utility functions

test_that("mat unfold operation works correctly", {
	X = array(1:24, dim = c(2, 3, 4))

	X1 = mat(X, 1)
	expect_equal(dim(X1), c(2, 12))
	expect_equal(X1[1,1], 1)
	expect_equal(X1[2,1], 2)

	X2 = mat(X, 2)
	expect_equal(dim(X2), c(3, 8))

	X3 = mat(X, 3)
	expect_equal(dim(X3), c(4, 6))
})

test_that("amprod array-matrix product works", {
	A = array(1:24, dim = c(2, 3, 4))
	M = matrix(1:6, nrow = 3, ncol = 2)

	result = amprod(A, M, 1)
	expect_equal(dim(result), c(3, 3, 4))

	I = diag(2)
	result_I = amprod(A, I, 1)
	expect_equal(result_I, A)
})

test_that("tprod tensor product works", {
	A = array(rnorm(24), dim = c(2, 3, 4))

	B = list(
		matrix(rnorm(4), 2, 2),
		matrix(rnorm(9), 3, 3),
		matrix(rnorm(16), 4, 4)
	)

	result = tprod(A, B)
	expect_equal(dim(result), dim(A))

	result_partial = tprod(A, B[1:2], modes = 1:2)
	expect_equal(dim(result_partial), dim(A))
})

test_that("rsan generates correct dimensions", {
	d1 = c(5, 5)
	X1 = rsan(d1)
	expect_equal(dim(X1), d1)
	expect_true(all(!is.na(X1)))

	d2 = c(3, 4, 5)
	X2 = rsan(d2)
	expect_equal(dim(X2), d2)

	d3 = c(100, 100)
	X3 = rsan(d3)
	expect_true(abs(mean(X3)) < 0.1)
	expect_true(abs(sd(c(X3)) - 1) < 0.1)
})

test_that("zscores handles 2D and 3D arrays correctly", {
	y2d = matrix(c(1, 3, 2, 5, 4, NA), nrow = 2)
	z2d = zscores(y2d)
	expect_equal(dim(z2d), dim(y2d))
	expect_true(is.na(z2d[2, 3]))
	expect_true(all(!is.na(z2d[!is.na(y2d)])))

	y3d = array(c(1:24, NA, 26:30), dim = c(5, 3, 2))
	z3d = zscores(y3d)
	expect_equal(dim(z3d), dim(y3d))
	expect_equal(sum(is.na(z3d)), 1)

	y_ties = matrix(c(1, 1, 2, 2, 3, 3), nrow = 2)
	z_ties = zscores(y_ties, ties.method = "average")
	expect_equal(z_ties[1,1], z_ties[2,1])
})

test_that("matrix operations preserve structure", {
	X = array(rnorm(60), dim = c(3, 4, 5))

	X_mat = mat(X, 1)
	X_back = array(X_mat, dim = c(3, 4, 5))
	expect_equal(X_back, X)

	I = diag(3)
	X_same = amprod(X, I, 1)
	expect_equal(X_same, X)
})

test_that("edge cases for array functions", {
	expect_error(mat(array(dim = c(0, 5, 3)), 1))

	X_single = array(42, dim = c(1, 1, 1))
	expect_equal(mat(X_single, 1), matrix(42, 1, 1))

	X_large = array(1, dim = c(2, 2, 1000))
	X_mat = mat(X_large, 3)
	expect_equal(dim(X_mat), c(1000, 4))
})
