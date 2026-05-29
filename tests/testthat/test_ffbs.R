#### FFBS C++ smoke tests

test_that("ffbs_dlm handles a basic state-space model", {
	Tt = 10; n = 4
	# y is a list of per-time observation vectors of length n
	y = lapply(seq_len(Tt), function(t) rnorm(n))
	Flist = replicate(Tt, diag(n), simplify = FALSE)
	V = diag(n); W = diag(n) * 0.1
	m0 = rep(0, n); C0 = diag(n) * 10
	res = ffbs_dlm(y, Flist, V, W, m0, C0)
	expect_true(is.array(res) || is.matrix(res))
	expect_true(all(is.finite(res)))
})

test_that("ffbs_theta performs forward filtering on bilinear dynamics", {
	skip_on_cran()
	n = 4; Tt = 8
	# ffbs_theta(Z, mu, Aarray, Barray, sigma2): per-time A and B as
	# n x n x Tt cubes, baseline mu as n x n, scalar sigma2.
	Z = array(rnorm(n * n * Tt), dim = c(n, n, Tt))
	A = array(0, dim = c(n, n, Tt))
	B = array(0, dim = c(n, n, Tt))
	for (t in seq_len(Tt)) {
		A[, , t] = diag(n) + rnorm(n * n) * 0.1
		B[, , t] = diag(n) + rnorm(n * n) * 0.1
	}
	mu = matrix(0, n, n)
	res = tryCatch(ffbs_theta(Z, mu, A, B, sigma2 = 1),
		error = function(e) e)
	expect_false(inherits(res, "error"))
})
