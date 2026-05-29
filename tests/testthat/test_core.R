#### core dispatch and input validation

test_that("dbn rejects malformed data", {
	# a 2D matrix is neither a 3D nor a 4D array
	expect_error(dbn(matrix(1:9, 3, 3)), regexp = "3D|4D|array")
	# a character array (length > 1) is not a file path either
	expect_error(dbn(c("a", "b", "c")), regexp = "numeric array|character")
})

test_that("dbn can load Y from an .RData file path", {
	skip_on_cran()
	tmp_file = tempfile(fileext = ".RData")
	on.exit(unlink(tmp_file), add = TRUE)
	Y = simulate_static_dbn(n = 5, p = 1, time = 5, seed = 1)$Y
	save(Y, file = tmp_file)
	expect_message(
		fit <- dbn(tmp_file, model = "static",
			nscan = 30, burn = 10, verbose = FALSE),
		regexp = "Loading data"
	)
	expect_s3_class(fit, "dbn")
})

test_that("verbose = FALSE silences chatter; verbose = TRUE informs", {
	sim = simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	msgs_off = character()
	withCallingHandlers(
		suppressWarnings(dbn(sim$Y, family = "gaussian",
			nscan = 10, burn = 5, verbose = FALSE)),
		message = function(m) { msgs_off <<- c(msgs_off, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_false(any(grepl("defaulting to", msgs_off)))

	msgs_on = character()
	withCallingHandlers(
		suppressWarnings(dbn(sim$Y, family = "gaussian",
			nscan = 10, burn = 5, verbose = TRUE)),
		message = function(m) { msgs_on <<- c(msgs_on, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_true(any(grepl("defaulting to", msgs_on)))
})

test_that("static model produces a valid fit object", {
	skip_on_cran()
	sim = simulate_static_dbn(n = 5, p = 1, time = 5, seed = 1)
	fit = dbn(sim$Y, model = "static", family = "gaussian",
		nscan = 30, burn = 10, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "static")
	# core slots every fit must carry; let the implementation grow other slots
	# without flagging a stale-test failure.
	for (slot in c("model", "family", "dims", "settings", "Y", "M", "meta")) {
		expect_true(slot %in% names(fit), info = slot)
	}
})

test_that("dynamic model produces a valid fit object", {
	skip_on_cran()
	sim = simulate_dynamic_dbn(n = 5, time = 8, seed = 1)
	fit = dbn(sim$Y, model = "dynamic", family = "gaussian",
		nscan = 30, burn = 10, verbose = FALSE)
	expect_s3_class(fit, "dbn")
	expect_equal(fit$model, "dynamic")
	expect_true(is.list(fit$A) && length(fit$A) >= 2L)
})

test_that("thinning reduces stored draws", {
	skip_on_cran()
	sim = simulate_static_dbn(n = 5, p = 1, time = 5, seed = 1)
	fit_thin1 = dbn(sim$Y, model = "static", nscan = 60, burn = 20, odens = 1,
		verbose = FALSE)
	fit_thin5 = dbn(sim$Y, model = "static", nscan = 60, burn = 20, odens = 5,
		verbose = FALSE)
	n_thin1 = fit_thin1$settings$draws %||% fit_thin1$meta$draws
	n_thin5 = fit_thin5$settings$draws %||% fit_thin5$meta$draws
	expect_gt(n_thin1, n_thin5)
})
