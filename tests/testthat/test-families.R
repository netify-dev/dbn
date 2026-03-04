####
# family-specific tests
####

test_that("Gaussian family works correctly", {
	# simulate very small gaussian network with stable data
	set.seed(6886)
	m <- 3
	Tt <- 5
	
	# create simple, stable test data
	Y_gau <- array(0, c(m, m, 1, Tt))
	for(t in 1:Tt) {
		# create a simple pattern with small noise
		base <- matrix(0.1 * (1:m) %*% t(1:m), m, m)
		Y_gau[,,1,t] <- base + matrix(rnorm(m*m, 0, 0.1), m, m)
	}
	
	# fit low-rank dbn with fewer iterations
	fit_g <- dbn_lowrank(Y_gau, r = 1, nscan = 50, burn = 10, odens = 1,
											 family = "gaussian", verbose = FALSE)
	
	expect_s3_class(fit_g, "dbn")
	expect_true("sigma2_proc" %in% names(fit_g))
	expect_true("sigma2_obs" %in% names(fit_g))
	expect_equal(fit_g$settings$family, "gaussian")
})

test_that("Binary family works correctly", {
	skip_if_not_installed("truncnorm")

	# Simulate probit-binary network
	set.seed(6885)
	m <- 4
	Tt <- 6
	prob <- array(0.3, dim = c(m, m, 1, Tt))
	Y_bin <- array(rbinom(length(prob), 1, prob), dim = c(m, m, 1, Tt))
	for (t in 1:Tt) diag(Y_bin[,,1,t]) <- NA

	# Fit HMM DBN with binary family
	fit_b <- dbn(Y_bin, model = "hmm", family = "binary",
							 nscan = 100, burn = 50, verbose = FALSE, R = 2)

	expect_s3_class(fit_b, "dbn")
	expect_equal(fit_b$model, "hmm")
	expect_true(!is.null(fit_b$S))
})

