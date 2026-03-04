####
# End-to-end validation: Dynamic Gaussian model
# Simulates a moderately sized network and runs the full analysis pipeline.
####

library(dbn)
set_dbn_threads(min(4, parallel::detectCores(logical = FALSE)))

cat("=== Dynamic Gaussian End-to-End Validation ===\n\n")

# 1. Simulate
cat("1. Simulating trade-like network (n=30, Tt=20)...\n")
sim <- simulate_dynamic_dbn(
	n = 30, p = 1, time = 20,
	sigma2 = 0.3, tauA2 = 0.05, tauB2 = 0.05,
	seed = 42
)
cat(sprintf("   Data dimensions: %s\n", paste(dim(sim$Z), collapse = " x ")))

# 2. Memory estimation
cat("\n2. Memory estimation...\n")
est <- estimate_memory(n_row = 30, n_col = 30, p = 1, Tt = 20,
	nscan = 5000, burn = 2000, odens = 5)
print(est)

# 3. Fit
cat("\n3. Fitting dynamic gaussian model...\n")
t0 <- proc.time()
fit <- dbn(sim$Z,
	model = "dynamic", family = "gaussian",
	nscan = 5000, burn = 2000, odens = 5,
	verbose = 500L
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("   Fit completed in %.1f seconds\n", elapsed))

# 4. Summary
cat("\n4. Model summary:\n")
summary(fit)

# 5. Convergence diagnostics
cat("\n5. Convergence diagnostics:\n")
ps <- param_summary(fit)
print(ps)

# 6. Posterior predictive check
cat("\n6. Posterior predictive distribution...\n")
ppd <- posterior_predict_dbn(fit, ndraws = 100, seed = 42)

obs <- sim$Z[, , 1, ]
yrep <- ppd$y_rep
nd <- length(dim(yrep))
q_lo <- apply(yrep, 1:(nd - 1), quantile, probs = 0.05, na.rm = TRUE)
q_hi <- apply(yrep, 1:(nd - 1), quantile, probs = 0.95, na.rm = TRUE)
if (length(dim(q_lo)) >= 3 && dim(q_lo)[3] >= 1) {
	in_interval <- (obs >= q_lo[, , 1, ]) & (obs <= q_hi[, , 1, ])
} else {
	in_interval <- (obs >= q_lo) & (obs <= q_hi)
}
coverage <- mean(in_interval, na.rm = TRUE)
cat(sprintf("   90%% PPD coverage: %.1f%%\n", coverage * 100))

# 7. Impulse response
cat("\n7. Impulse response analysis...\n")
S <- build_shock(30, type = "density")
irf <- compute_irf(fit, shock = S, H = 10, n_draws = 50, stat_fun = dbn::stat_density)
cat("   IRF computed successfully\n")

# 8. Forecasting
cat("\n8. Forecasting 5 steps ahead...\n")
pred <- predict(fit, H = 5, draws = 50, summary = "mean")
cat(sprintf("   Forecast dimensions: %s\n", paste(dim(pred), collapse = " x ")))

cat("\n=== Validation complete ===\n")
