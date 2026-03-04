####
# End-to-end validation: HMM vs Dynamic model comparison
# Demonstrates model comparison workflow with ordinal data.
####

library(dbn)
set_dbn_threads(min(4, parallel::detectCores(logical = FALSE)))

cat("=== HMM vs Dynamic Model Comparison ===\n\n")

# 1. Simulate regime-switching data
cat("1. Simulating regime-switching network (n=10, Tt=50, R=2)...\n")
sim <- simulate_hmm_dbn(
	n = 10, p = 1, time = 50,
	R = 2, transition_prob = 0.9,
	sigma2 = 0.3,
	seed = 42
)
cat(sprintf("   Regime distribution: %s\n",
	paste(names(table(sim$S)), table(sim$S), sep = "=", collapse = ", ")))

# 2. Fit dynamic model
cat("\n2. Fitting dynamic ordinal model...\n")
t0 <- proc.time()
fit_dyn <- dbn(sim$Y,
	model = "dynamic", family = "ordinal",
	nscan = 3000, burn = 1000, odens = 3,
	verbose = FALSE
)
cat(sprintf("   Dynamic fit: %.1fs\n", (proc.time() - t0)[3]))

# 3. Fit HMM model
cat("\n3. Fitting HMM ordinal model (R=2)...\n")
t0 <- proc.time()
fit_hmm <- dbn(sim$Y,
	model = "hmm", family = "ordinal",
	R = 2,
	nscan = 3000, burn = 1000, odens = 3,
	verbose = FALSE
)
cat(sprintf("   HMM fit: %.1fs\n", (proc.time() - t0)[3]))

# 4. Compare
cat("\n4. Model summaries:\n")
cat("\n--- Dynamic ---\n")
ps_dyn <- param_summary(fit_dyn)
print(ps_dyn)

cat("\n--- HMM ---\n")
ps_hmm <- param_summary(fit_hmm)
print(ps_hmm)

# 5. Regime recovery
cat("\n5. Regime recovery:\n")
rp <- regime_probs(fit_hmm)
modal <- apply(rp, 1, which.max)
agree_direct <- mean(modal == sim$S)
agree_flip <- mean(modal == (3 - sim$S))
best <- max(agree_direct, agree_flip)
cat(sprintf("   Best regime agreement: %.1f%%\n", best * 100))

# 6. Dyad tracking
cat("\n6. Dyad trajectory comparison (actors 1 -> 3):\n")
tc_dyn <- theta_credible(fit_dyn, i = 1, j = 3, rel = 1)
cat(sprintf("   Dynamic: %d time points with credible intervals\n", nrow(tc_dyn)))

cat("\n=== Comparison complete ===\n")
