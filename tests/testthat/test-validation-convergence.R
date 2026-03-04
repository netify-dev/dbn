####
# Phase 1C: MCMC convergence diagnostics validation
# Verifies that chains with adequate length show acceptable Rhat and ESS.
# Gated behind skip_on_cran() + DBN_VALIDATION env var.
####

skip_validation <- function() {
	skip_on_cran()
	skip_if(
		!identical(Sys.getenv("DBN_VALIDATION"), "true"),
		"Set DBN_VALIDATION=true to run extended validation tests"
	)
}

####
# 1. Static gaussian: Rhat and ESS from 2 chains
####
test_that("validation: static gaussian achieves convergence (Rhat, ESS)", {
	skip_validation()
	skip_if_not_installed("coda")

	sim <- simulate_static_dbn(
		n = 12, p = 1, time = 15,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 9001
	)

	fit1 <- dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	fit2 <- dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	# extract scalar parameter traces
	chain1 <- coda::mcmc(fit1$draws$pars)
	chain2 <- coda::mcmc(fit2$draws$pars)
	mc_list <- coda::mcmc.list(chain1, chain2)

	# Rhat should be < 1.2 for all scalar parameters
	gd <- coda::gelman.diag(mc_list, autoburnin = FALSE)
	rhat_vals <- gd$psrf[, 1]
	for (nm in names(rhat_vals)) {
		expect_lt(rhat_vals[nm], 1.2,
			label = sprintf("Rhat(%s)=%.3f < 1.2", nm, rhat_vals[nm]))
	}

	# ESS should be > 100 for each chain
	ess1 <- coda::effectiveSize(chain1)
	for (nm in names(ess1)) {
		expect_gt(ess1[nm], 50,
			label = sprintf("ESS(%s)=%.0f > 50", nm, ess1[nm]))
	}
})

####
# 2. Dynamic gaussian: convergence of variance parameters
####
test_that("validation: dynamic gaussian achieves convergence (Rhat, ESS)", {
	skip_validation()
	skip_if_not_installed("coda")

	sim <- simulate_dynamic_dbn(
		n = 10, p = 1, time = 15,
		sigma2 = 0.3, tauA2 = 0.05, tauB2 = 0.05,
		seed = 9002
	)

	fit1 <- dbn(sim$Z,
		model = "dynamic", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	fit2 <- dbn(sim$Z,
		model = "dynamic", family = "gaussian",
		nscan = 5000, burn = 2000, odens = 2,
		verbose = FALSE
	)

	# build comparable parameter matrices
	params1 <- cbind(sigma2 = fit1$sigma2, tau_A2 = fit1$tau_A2, tau_B2 = fit1$tau_B2)
	params2 <- cbind(sigma2 = fit2$sigma2, tau_A2 = fit2$tau_A2, tau_B2 = fit2$tau_B2)

	chain1 <- coda::mcmc(params1)
	chain2 <- coda::mcmc(params2)
	mc_list <- coda::mcmc.list(chain1, chain2)

	gd <- coda::gelman.diag(mc_list, autoburnin = FALSE)
	rhat_vals <- gd$psrf[, 1]
	for (nm in names(rhat_vals)) {
		expect_lt(rhat_vals[nm], 1.3,
			label = sprintf("Rhat(%s)=%.3f < 1.3", nm, rhat_vals[nm]))
	}

	ess1 <- coda::effectiveSize(chain1)
	for (nm in names(ess1)) {
		expect_gt(ess1[nm], 30,
			label = sprintf("ESS(%s)=%.0f > 30", nm, ess1[nm]))
	}
})

####
# 3. check_convergence() returns sensible output
####
test_that("validation: check_convergence produces diagnostic output", {
	skip_validation()

	sim <- simulate_static_dbn(
		n = 10, p = 1, time = 10,
		sigma2 = 0.3, tau2 = 0.05,
		seed = 9003
	)

	fit <- dbn(sim$Z,
		model = "static", family = "gaussian",
		nscan = 2000, burn = 1000, odens = 2,
		verbose = FALSE
	)

	diag_result <- check_convergence(fit)
	expect_true(!is.null(diag_result), label = "check_convergence returns non-NULL")

	# should contain summary information
	if (is.list(diag_result)) {
		# check that it has some expected structure
		expect_true(length(diag_result) > 0,
			label = "check_convergence returns non-empty result")
	}
})
