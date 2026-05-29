# thread setting round-trips
test_that("set_dbn_threads round-trips", {
	prev = get_dbn_threads()
	on.exit(set_dbn_threads(prev), add = TRUE)
	suppressMessages(set_dbn_threads(2L))
	expect_equal(get_dbn_threads(), 2L)
})

# keep argument drops Theta and shrinks the fit
test_that("keep argument drops Theta and shrinks the fit", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 4)
	fit_full = dbn(
		sim$Y,
		model = "dynamic",
		family = "gaussian",
		nscan = 30,
		burn = 10,
		verbose = FALSE
	)
	fit_thin = dbn(
		sim$Y,
		model = "dynamic",
		family = "gaussian",
		nscan = 30,
		burn = 10,
		verbose = FALSE,
		keep = c("A", "B", "M")
	)
	expect_null(fit_thin$Theta)
	expect_lt(
		length(serialize(fit_thin, NULL)),
		length(serialize(fit_full, NULL))
	)
})

# detects symmetric adjacency
test_that(".dbn_data_looks_undirected detects symmetric adjacency", {
	Y_sym = array(NA, dim = c(4, 4, 1, 3))
	for (t in 1:3) {
		M = matrix(rnorm(16), 4, 4)
		M = (M + t(M)) / 2
		diag(M) = NA
		Y_sym[, , 1, t] = M
	}
	expect_true(dbn:::.dbn_data_looks_undirected(Y_sym, 4, 4))
})
