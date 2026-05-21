####
# PPD for all model types
####

test_that("posterior_predict_dbn works for static model", {
	skip_on_cran()
	set.seed(123)
	n = 8
	p = 2
	Tt = 15
	
	# simulate data
	Y = array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
	
	# fit static
	fit = dbn_static(Y, family = "ordinal", nscan = 100, burn = 50, verbose = FALSE)
	
	# PPD
	ppd = posterior_predict_dbn(fit, draws = 20, seed = 123)
	
	# structure
	expect_type(ppd, "list")
	expect_equal(length(ppd), 20)
	expect_equal(dim(ppd[[1]]), dim(Y))
	
	# ordinal range
	vals = unique(unlist(ppd))
	vals = vals[!is.na(vals)]
	expect_true(all(vals %in% 1:5))
})

test_that("posterior_predict_dbn works for dynamic model", {
	skip_on_cran()
	set.seed(456)
	n = 6
	p = 1
	Tt = 10
	
	# simulate data
	Y = array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
	
	# fit dynamic
	fit = dbn_dynamic(Y, family = "ordinal", nscan = 100, burn = 50, verbose = FALSE)
	
	# PPD
	ppd = posterior_predict_dbn(fit, draws = 15, seed = 456)
	
	# structure
	expect_type(ppd, "list")
	expect_equal(length(ppd), 15)
	expect_equal(dim(ppd[[1]]), dim(Y))
	
	# valid range
	vals = unique(unlist(ppd))
	vals = vals[!is.na(vals)]
	expect_true(all(vals %in% 1:5))
})

test_that("posterior_predict_dbn works for dynamic model with time thinning", {
	skip_on_cran()
	set.seed(789)
	n = 5
	p = 2
	Tt = 20
	
	# simulate data
	Y = array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
	
	# fit dynamic with time thinning
	fit = dbn_dynamic(Y, family = "ordinal", nscan = 100, burn = 50, 
										 time_thin = 2, verbose = FALSE)
	
	# A and B have reduced time dim
	expect_equal(dim(fit$A[[1]])[3], 10)  # Tt / 2
	
	# PPD works with thinned draws
	ppd = posterior_predict_dbn(fit, draws = 10, seed = 789)
	
	# structure
	expect_type(ppd, "list")
	expect_equal(length(ppd), 10)
	expect_equal(dim(ppd[[1]]), dim(Y))  # matches original Y
})

test_that("posterior_predict_dbn works for gaussian family", {
	skip_on_cran()
	set.seed(111)
	n = 5
	p = 1
	Tt = 10
	
	# gaussian data
	Y = array(rnorm(n*n*p*Tt), dim = c(n, n, p, Tt))
	
	# static
	fit_static = dbn_static(Y, family = "gaussian", nscan = 80, burn = 40, verbose = FALSE)
	ppd_static = posterior_predict_dbn(fit_static, draws = 10, seed = 111)
	
	expect_type(ppd_static, "list")
	expect_equal(length(ppd_static), 10)
	expect_equal(dim(ppd_static[[1]]), dim(Y))
	
	# dynamic
	fit_dynamic = dbn_dynamic(Y, family = "gaussian", nscan = 80, burn = 40, verbose = FALSE)
	ppd_dynamic = posterior_predict_dbn(fit_dynamic, draws = 10, seed = 222)
	
	expect_type(ppd_dynamic, "list")
	expect_equal(length(ppd_dynamic), 10)
	expect_equal(dim(ppd_dynamic[[1]]), dim(Y))
})

test_that("posterior_predict_dbn works for binary family", {
	skip_on_cran()
	set.seed(333)
	n = 6
	p = 1
	Tt = 8
	
	# binary data
	Y = array(sample(0:1, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
	
	# static
	fit_static = dbn_static(Y, family = "binary", nscan = 60, burn = 30, verbose = FALSE)
	ppd_static = posterior_predict_dbn(fit_static, draws = 8, seed = 333)
	
	expect_type(ppd_static, "list")
	expect_equal(length(ppd_static), 8)
	expect_equal(dim(ppd_static[[1]]), dim(Y))
	
	# binary values
	vals = unique(unlist(ppd_static))
	vals = vals[!is.na(vals)]
	expect_true(all(vals %in% 0:1))
})

test_that("posterior_predict_dbn works for lowrank model", {
	skip_on_cran()
	skip_if_not_installed("Matrix")

	set.seed(444)
	n = 8
	p = 1
	Tt = 12
	
	# simulate data
	Y = array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
	
	# fit lowrank
	fit = dbn_lowrank(Y, family = "ordinal", r = 2, nscan = 80, burn = 40, verbose = FALSE)
	
	# PPD
	ppd = posterior_predict_dbn(fit, draws = 10, seed = 444)
	
	# structure
	expect_type(ppd, "list")
	expect_equal(length(ppd), 10)
	expect_equal(dim(ppd[[1]]), dim(Y))
})

test_that("posterior_predict_dbn works for hmm model", {
	skip_on_cran()
	set.seed(555)
	n = 6
	p = 1
	Tt = 15
	
	# simulate data
	Y = array(sample(1:5, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
	
	# fit HMM, R=2
	fit = dbn_hmm(Y, family = "ordinal", R = 2, nscan = 100, burn = 50, verbose = FALSE)
	
	# PPD
	ppd = posterior_predict_dbn(fit, draws = 12, seed = 555)
	
	# structure
	expect_type(ppd, "list")
	expect_equal(length(ppd), 12)
	expect_equal(dim(ppd[[1]]), dim(Y))
})

test_that("posterior_predict_dbn handles edge cases", {
	skip_on_cran()
	set.seed(666)
	n = 4
	p = 1
	Tt = 5
	
	# small data
	Y = array(sample(1:3, n*n*p*Tt, replace = TRUE), dim = c(n, n, p, Tt))
	
	# minimal iterations
	fit = dbn_static(Y, family = "ordinal", nscan = 20, burn = 10, verbose = FALSE)
	
	# more draws than available
	ppd = posterior_predict_dbn(fit, draws = 50, seed = 666)
	
	# samples with replacement
	expect_equal(length(ppd), 50)
	
	# specific draws
	ppd_specific = posterior_predict_dbn(fit, draw_indices = c(1, 5, 10), seed = 777)
	expect_equal(length(ppd_specific), 3)
})

test_that("posterior_predict_dbn works with 3D input data", {
	skip_on_cran()
	set.seed(888)
	n = 5
	Tt = 10
	
	# 3D array (single relation)
	Y_3d = array(sample(1:5, n*n*Tt, replace = TRUE), dim = c(n, n, Tt))
	
	# fit static with 3D input
	suppressMessages({
		fit = dbn(Y_3d, family = "ordinal", model = "static", 
							 nscan = 50, burn = 25, verbose = FALSE)
	})
	
	# predictions
	ppd = posterior_predict_dbn(fit, draws = 10, seed = 888)
	
	# matches expanded 4D structure
	expect_type(ppd, "list")
	expect_equal(length(ppd), 10)
	expect_equal(dim(ppd[[1]]), c(n, n, 1, Tt))  # 4D with p=1
})