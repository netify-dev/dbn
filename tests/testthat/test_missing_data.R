####
# missing data handling
####

test_that("models handle missing data", {
	
	# data with missing values
	set.seed(6886)
	Y <- array(sample(1:5, 5*5*2*10, replace = TRUE), dim = c(5, 5, 2, 10))
	
	# 10% missing at random
	n_missing <- round(0.1 * length(Y))
	Y[sample(length(Y), n_missing)] <- NA
	
	# static with missing
	expect_no_error({
		res_static <- dbn(Y, model = "static", nscan = 50, burn = 10, odens = 1, verbose = FALSE)
	})
	
	# dynamic with missing
	expect_no_error({
		res_dynamic <- dbn(Y, model = "dynamic", nscan = 50, burn = 10, odens = 1, verbose = FALSE)
	})
	
	# lowrank with smaller network
	Y_small <- array(sample(1:5, 3*3*1*10, replace = TRUE), dim = c(3, 3, 1, 10))
	Y_small[sample(length(Y_small), 5)] <- NA
	
	expect_no_error({
		res_lr <- dbn_lowrank(Y_small, r = 2, nscan = 50, burn = 10, odens = 1, verbose = FALSE)
	})
	
	# HMM with missing
	expect_no_error({
		res_hmm <- dbn_hmm(Y, R = 2, nscan = 50, burn = 10, odens = 1, verbose = FALSE)
	})
})

test_that("rank likelihood handles missing values correctly", {
	
	# rank data with NA
	R <- matrix(c(1, 2, NA, 3, 4, NA, 2, 1, 3), 3, 3)
	
	# zscores handles NA
	z <- zscores(R)
	expect_equal(sum(is.na(z)), sum(is.na(R)))
	expect_true(all(!is.na(z[!is.na(R)])))
	
	# rank index building
	IR <- precompute_ranks(array(R, dim = c(3, 3, 1, 1)))
	expect_true("NA" %in% names(IR[[1]]))
	expect_equal(length(IR[[1]][["NA"]]), sum(is.na(R)))
	
	# rz_fc handles NA
	Z <- matrix(rnorm(9), 3, 3)
	EZ <- matrix(rnorm(9), 3, 3)
	Z_new <- rz_fc(R, Z, EZ, IR[[1]])
	
	expect_equal(dim(Z_new), dim(R))
	expect_true(all(is.finite(Z_new)))
})

test_that("missing data patterns are preserved", {
	
	# specific missing pattern
	Y <- array(1:40, dim = c(2, 2, 2, 5))
	Y[1, 1, 1, 1:2] <- NA
	Y[2, 2, 2, 5] <- NA
	
	# preprocessing preserves NA pattern
	pre <- shared_preprocess(Y)
	
	# NA positions map to finite Z values
	expect_true(all(is.finite(pre$Z)))
	
	# IR tracks NA positions
	expect_equal(length(pre$IR[[1]][["NA"]]), 2)
	expect_equal(length(pre$IR[[2]][["NA"]]), 1)
})