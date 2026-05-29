# the fit-time variance-saturation diagnostic must catch both full saturation
# (every draw at the 1e8 ceiling) and partial saturation (a subset clamped,
# dragging the posterior mean or upper ci into the millions).

test_that(".warn_variance_saturation flags a fully-saturated chain", {
	fake = list(model = "dynamic",
							 sigma2 = rep(1e8, 200),         # every draw at the ceiling
							 tau_A2 = rnorm(200, 0.1, 0.01),
							 tau_B2 = rnorm(200, 0.1, 0.01))
	expect_warning(hit <- dbn:::.warn_variance_saturation(fake),
								 regexp = "saturated the sampler ceiling for every")
	expect_true("sigma2" %in% hit)
})

test_that(".warn_variance_saturation flags a partially-saturated chain", {
	# 20% of draws clamped at the ceiling, rest small -> mean is in the millions
	v = c(rep(1e8, 40), rnorm(160, 0.05, 0.01))
	fake = list(model = "hmm", g2 = v,
							 sigma2 = rnorm(200, 1, 0.1),
							 tau_A2 = rnorm(200, 0.3, 0.05),
							 tau_B2 = rnorm(200, 0.3, 0.05))
	expect_warning(hit <- dbn:::.warn_variance_saturation(fake),
								 regexp = "hit the sampler ceiling on a subset")
	expect_true("g2" %in% hit)
})

test_that(".warn_variance_saturation flags a heavy right tail even below 5%", {
	# a couple of draws orders of magnitude above the median
	v = c(1e8, 5e7, rnorm(198, 0.05, 0.01))
	fake = list(model = "hmm", tau_B2 = v)
	expect_warning(dbn:::.warn_variance_saturation(fake),
								 regexp = "subset of draws")
})

test_that(".warn_variance_saturation stays silent on a healthy fit", {
	fake = list(model = "dynamic",
							 sigma2 = rnorm(200, 1, 0.1),
							 tau_A2 = rnorm(200, 0.3, 0.05),
							 tau_B2 = rnorm(200, 0.3, 0.05),
							 g2 = rnorm(200, 0.5, 0.05))
	expect_silent(hit <- dbn:::.warn_variance_saturation(fake))
	expect_length(hit, 0L)
})

test_that(".warn_variance_saturation reads static/piecewise params traces", {
	# static/piecewise store scalar draws in a $params data.frame, not at the
	# top level; the diagnostic must scan there too.
	pars = data.frame(s2 = rnorm(200, 1, 0.1),
										 t2 = c(rep(1e8, 30), rnorm(170, 0.1, 0.02)),  # partial
										 g2 = rnorm(200, 0.5, 0.05))
	fake = list(model = "piecewise", params = pars)
	expect_warning(hit <- dbn:::.warn_variance_saturation(fake),
								 regexp = "subset")
	expect_true("t2" %in% hit)
})
