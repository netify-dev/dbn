####
# End-to-end validation: Low-rank model scalability
# Tests with a larger network to demonstrate scalability.
####

library(dbn)
set_dbn_threads(min(4, parallel::detectCores(logical = FALSE)))

cat("=== Low-rank Scalability Validation ===\n\n")

# 1. Simulate large sparse network
cat("1. Simulating large network (n=50, Tt=10, r=3)...\n")
sim <- simulate_lowrank_dbn(
	n = 50, p = 1, time = 10,
	r = 3,
	sigma2 = 0.2,
	tau_alpha2 = 0.05,
	tauB2 = 0.03,
	seed = 42
)
cat(sprintf("   Data dimensions: %s\n", paste(dim(sim$Y), collapse = " x ")))

# 2. Memory estimation
cat("\n2. Memory estimation...\n")
est <- estimate_memory(n_row = 50, n_col = 50, p = 1, Tt = 10,
	nscan = 3000, burn = 1000, odens = 3,
	family = "ordinal")
print(est)

# 3. Fit
cat("\n3. Fitting lowrank model (r=3)...\n")
t0 <- proc.time()
fit <- dbn(sim$Y,
	model = "lowrank", family = "ordinal",
	r = 3,
	nscan = 3000, burn = 1000, odens = 3,
	verbose = 500L
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("   Fit completed in %.1f seconds\n", elapsed))

# 4. Summary
cat("\n4. Model summary:\n")
summary(fit)

# 5. Diagnostics
cat("\n5. Parameter summary:\n")
ps <- param_summary(fit)
print(ps)

# 6. Factor structure
cat("\n6. Factor trajectories (tidy format):\n")
alpha_df <- tidy_dbn_lowrank(fit)
cat(sprintf("   %d rows, %d columns\n", nrow(alpha_df), ncol(alpha_df)))
print(head(alpha_df))

# 7. U orthonormality check
cat("\n7. U orthonormality check:\n")
U_last <- fit$U[[length(fit$U)]]
UtU <- crossprod(U_last)
offdiag_max <- max(abs(UtU - diag(ncol(U_last))))
cat(sprintf("   Max off-diagonal deviation from I: %.6f\n", offdiag_max))

# 8. Forecasting
cat("\n8. Forecasting 3 steps ahead...\n")
pred <- predict(fit, H = 3, draws = 30, summary = "mean")
cat(sprintf("   Forecast dimensions: %s\n", paste(dim(pred), collapse = " x ")))

cat("\n=== Scalability validation complete ===\n")
