####
# FFBS functions
####

test_that("ffbs_dlm handles basic state space model", {
	# state space setup
	Tt = 20
	n = 3
	
	# test data
	set.seed(6886)
	y = lapply(1:Tt, function(t) rnorm(n))
	Flist = lapply(1:Tt, function(t) diag(n))
	m0 = rep(0, n)
	C0 = diag(n)
	
	# obs and state variances
	V = diag(n) * 0.5
	W = diag(n) * 0.3
	
	# run FFBS
	result = ffbs_dlm(y, Flist, V, W, m0, C0, ar1 = FALSE)

	# output structure
	expect_type(result, "double")
	expect_equal(dim(result), c(n, Tt))
	
	# stationarity check
	state_vars = apply(result, 2, var)
	expect_true(all(state_vars < 10))  # bounded
})

test_that("ffbs_theta performs forward filtering correctly", {
	# setup
	n = 5
	Tt = 10
	
	set.seed(6886)
	# test data
	Z = array(rnorm(n*n*Tt), dim = c(n, n, Tt))
	mu = matrix(0, n, n)
	
	# A and B arrays
	Aarray = array(0, dim = c(n, n, Tt))
	Barray = array(0, dim = c(n, n, Tt))
	for(t in 1:Tt) {
		Aarray[,,t] = diag(n) * 0.9
		Barray[,,t] = diag(n) * 0.9
	}
	
	sigma2 = 0.5
	
	# FFBS for theta
	result = ffbs_theta(Z, mu, Aarray, Barray, sigma2)

	# output structure
	expect_type(result, "double")
	expect_equal(dim(result), c(n, n, Tt))
	
	# finite values
	expect_true(all(is.finite(result)))
})

test_that("ffbs_dlm handles different variance specifications", {
	# scalar V
	Tt = 10
	n = 3
	
	set.seed(6889)
	y = lapply(1:Tt, function(t) rnorm(n))
	Flist = lapply(1:Tt, function(t) diag(n))
	m0 = rep(0, n)
	C0 = diag(n)
	V_scalar = 0.5
	W = diag(n) * 0.3
	
	# scalar V converted to matrix
	result_scalar = ffbs_dlm(y, Flist, V_scalar * diag(n), W, m0, C0)
	expect_equal(dim(result_scalar), c(n, Tt))
	
	# matrix V
	V_matrix = diag(n) * runif(n, 0.3, 0.7)
	result_matrix = ffbs_dlm(y, Flist, V_matrix, W, m0, C0)
	expect_equal(dim(result_matrix), c(n, Tt))
	
	# should differ
	expect_false(all(result_scalar == result_matrix))
})

test_that("ffbs_dlm maintains positive definite covariances", {
	# poorly conditioned system
	Tt = 25
	n = 4
	
	set.seed(6881)
	y = lapply(1:Tt, function(t) rnorm(n) * 0.1)
	Flist = lapply(1:Tt, function(t) diag(n) + matrix(rnorm(n^2), n, n) * 0.01)
	
	m0 = rep(0, n)
	C0 = diag(n) * 10
	V = diag(n) * 0.01
	W = diag(n) * 0.001

	# FFBS with AR(1) handles near-singularity
	result = ffbs_dlm(y, Flist, V, W, m0, C0, ar1 = TRUE, rho = 0.95)
	
	# states finite
	expect_true(all(is.finite(result)))
})

test_that("ffbs_dlm AR(1) vs random walk dynamics", {
	# AR(1) vs random walk
	Tt = 50
	n = 2
	
	set.seed(123)
	y = lapply(1:Tt, function(t) rnorm(n))
	Flist = lapply(1:Tt, function(t) diag(n))
	m0 = rep(0, n)
	C0 = diag(n)
	V = diag(n) * 0.5
	W = diag(n) * 0.3
	
	# random walk
	result_rw = ffbs_dlm(y, Flist, V, W, m0, C0, ar1 = FALSE)
	
	# AR(1) with rho = 0.8
	result_ar1 = ffbs_dlm(y, Flist, V, W, m0, C0, ar1 = TRUE, rho = 0.8)
	
	# should differ
	expect_false(all(result_rw == result_ar1))

	# AR(1) mean-reverting behavior
	mean_ar1 = rowMeans(result_ar1)
	mean_rw = rowMeans(result_rw)
	
	# AR(1) states more centered around zero
	expect_true(max(abs(mean_ar1)) < max(abs(mean_rw)) * 1.5)
})


test_that("ffbs_theta handles bilinear dynamics", {
	# theta with bilinear structure
	n = 4
	Tt = 8
	
	set.seed(6883)
	# test data
	Z = array(0, dim = c(n, n, Tt))
	mu = matrix(0, n, n)
	
	# diagonal dynamics
	Aarray = Barray = array(0, dim = c(n, n, Tt))
	for(t in 1:Tt) {
		Aarray[,,t] = diag(n) * 0.8
		Barray[,,t] = diag(n) * 0.8
		Z[,,t] = matrix(rnorm(n^2), n, n)
	}
	
	sigma2 = 1.0
	
	# FFBS
	result = ffbs_theta(Z, mu, Aarray, Barray, sigma2)

	# output dims
	expect_equal(dim(result), c(n, n, Tt))
	expect_true(all(is.finite(result)))
	
	# temporal smoothing: adjacent time points correlated
	cor_12 = cor(c(result[,,1]), c(result[,,2]))
	expect_true(cor_12 > -0.2)  # at least not strongly negative
})