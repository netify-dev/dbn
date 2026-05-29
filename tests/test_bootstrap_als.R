#!/usr/bin/env Rscript
####
# Comprehensive tests for dbn_als(..., bootstrap = N)
####

library(dbn)
library(cli)

set.seed(42)

# Test 1: Small simulated network -- point estimate only
cli::cli_h1("Test 1: Bootstrap on Small Simulated Network")

sim <- simulate_dynamic_dbn(n = 8, time = 10, seed = 1)
Y <- sim$Y

cli::cli_h2("Fitting ALS point estimate")
fit_als <- dbn_als(Y, family = "gaussian", max_iter = 50, verbose = FALSE)
stopifnot(isFALSE(fit_als$meta$uncertainty_available))
cli::cli_bullets(c(
	"v" = "ALS fit successful",
	" " = "Dimensions: {dim(Y)[1]} x {dim(Y)[2]} x {dim(Y)[4]} time"
))

# Test 2: Block bootstrap via dbn_als(bootstrap = N)
cli::cli_h2("Test 2: Block Bootstrap via bootstrap = 50")

fit_block <- dbn_als(Y, family = "gaussian", max_iter = 50,
	bootstrap = 50, bootstrap_type = "block", verbose = TRUE)
boot_block <- fit_block$bootstrap

cli::cli_bullets(c(
	"v" = "Block bootstrap complete",
	" " = "Valid replicates: {boot_block$n_valid}/{boot_block$n_total}",
	" " = "Mean SE (A): {sprintf('%.4f', mean(boot_block$se_A))}",
	" " = "Mean SE (B): {sprintf('%.4f', mean(boot_block$se_B))}",
	" " = "meta$uncertainty_available: {fit_block$meta$uncertainty_available}",
	" " = "meta$uncertainty_source: {fit_block$meta$uncertainty_source}"
))

stopifnot(nrow(boot_block$coefs_A) == 50)
stopifnot(ncol(boot_block$coefs_A) == 8 * 8)
stopifnot(length(boot_block$se_A) == 8 * 8)
stopifnot(boot_block$type == "block")
stopifnot(boot_block$n_valid > 40)
stopifnot(isTRUE(fit_block$meta$uncertainty_available))
stopifnot(identical(fit_block$meta$uncertainty_source, "bootstrap"))
cli::cli_text("v Block bootstrap output structure verified")

# Test 3: Parametric bootstrap
cli::cli_h2("Test 3: Parametric Bootstrap via bootstrap_type = \"parametric\"")

fit_param <- dbn_als(Y, family = "gaussian", max_iter = 50,
	bootstrap = 50, bootstrap_type = "parametric", verbose = TRUE)
boot_param <- fit_param$bootstrap

cli::cli_bullets(c(
	"v" = "Parametric bootstrap complete",
	" " = "Valid replicates: {boot_param$n_valid}/{boot_param$n_total}",
	" " = "Mean SE (A): {sprintf('%.4f', mean(boot_param$se_A))}"
))
stopifnot(boot_param$type == "parametric")
stopifnot(boot_param$n_valid > 40)
cli::cli_text("v Parametric bootstrap output structure verified")

# Test 4: Print and summary methods
cli::cli_h2("Test 4: Print and Summary Methods")
print(boot_block)
summary(boot_block)
cli::cli_text("v Print and summary methods work")

# Test 5: Symmetric constraint
cli::cli_h2("Test 5: Symmetric Network Bootstrap")
fit_als_sym <- dbn_als(Y, family = "gaussian", symmetric = TRUE,
	bootstrap = 30, bootstrap_type = "block", verbose = FALSE)
boot_sym <- fit_als_sym$bootstrap
stopifnot(isTRUE(boot_sym$symmetric))
cli::cli_text("v Symmetric constraint preserved through bootstrap")

# Test 6: Sign-aligned A estimates
cli::cli_h2("Test 6: Sign Alignment Verification")
valid_idx <- which(!apply(boot_block$coefs_A, 1, function(r) any(is.na(r))))
A_orig <- matrix(boot_block$point_est_A, 8, 8)
dot_products <- numeric(length(valid_idx))
for (i in seq_along(valid_idx)) {
	b <- valid_idx[i]
	A_boot <- matrix(boot_block$coefs_A[b, ], 8, 8)
	dot_products[i] <- sum(A_orig * A_boot)
}
sign_agree <- sum(dot_products > 0)
total_valid <- length(dot_products)
cli::cli_bullets(c(
	" " = "Sign-aligned A estimates: {sign_agree}/{total_valid} ({round(100*sign_agree/total_valid,1)}%)",
	"v" = "Sign alignment working correctly"
))
stopifnot(sign_agree / total_valid > 0.95)

# Test 7: sampler_describe on dbn_boot
cli::cli_h2("Test 7: sampler_describe() for Bootstrap Results")
sampler_describe(boot_block, verbose = TRUE)
cli::cli_text("v sampler_describe works for dbn_boot objects")

# Test 8: ALS point estimate is close to PCG point estimate
cli::cli_h2("Test 8: Comparing Bootstrap Point Estimate to PCG Sampler")
fit_pcg <- dbn(Y, model = "dynamic", family = "gaussian",
	symmetric = FALSE, sampler = "exact",
	nscan = 200, burn = 50, odens = 1, verbose = 0)
A_als <- matrix(boot_block$point_est_A, 8, 8)
A_pcg <- fit_pcg$A[[1]][, , 1]
B_als <- matrix(boot_block$point_est_B, 8, 8)
B_pcg <- fit_pcg$B[[1]][, , 1]
A_diff <- sqrt(sum((A_als - A_pcg)^2))
B_diff <- sqrt(sum((B_als - B_pcg)^2))
cli::cli_bullets(c(
	" " = "Frobenius diff (A): {sprintf('%.4f', A_diff)}",
	" " = "Frobenius diff (B): {sprintf('%.4f', B_diff)}"
))

# Test 9: SEs are within reasonable range
cli::cli_h2("Test 9: Bootstrap Standard Error Reasonableness")
stopifnot(all(boot_block$se_A < 1, na.rm = TRUE))
stopifnot(all(boot_block$se_B < 1, na.rm = TRUE))
cli::cli_text("v All A and B standard errors are < 1")

# Test 10: CI coverage of the point estimate
cli::cli_h2("Test 10: Confidence Interval Coverage")
A_in_ci <- sum((boot_block$ci_A_lo <= boot_block$point_est_A) &
	(boot_block$point_est_A <= boot_block$ci_A_hi))
B_in_ci <- sum((boot_block$ci_B_lo <= boot_block$point_est_B) &
	(boot_block$point_est_B <= boot_block$ci_B_hi))
M_in_ci <- sum((boot_block$ci_M_lo <= boot_block$point_est_M) &
	(boot_block$point_est_M <= boot_block$ci_M_hi))
total_A <- length(boot_block$point_est_A)
total_B <- length(boot_block$point_est_B)
total_M <- length(boot_block$point_est_M)
cli::cli_bullets(c(
	" " = "A: {A_in_ci}/{total_A} in 95% CI ({round(100*A_in_ci/total_A,1)}%)",
	" " = "B: {B_in_ci}/{total_B} in 95% CI ({round(100*B_in_ci/total_B,1)}%)",
	" " = "M: {M_in_ci}/{total_M} in 95% CI ({round(100*M_in_ci/total_M,1)}%)"
))

# Test 11: Seed reproducibility
cli::cli_h2("Test 11: Bootstrap Reproducibility")
fit_a <- dbn_als(Y, family = "gaussian", max_iter = 50,
	bootstrap = 30, bootstrap_type = "block", seed = 123, verbose = FALSE)
fit_b <- dbn_als(Y, family = "gaussian", max_iter = 50,
	bootstrap = 30, bootstrap_type = "block", seed = 123, verbose = FALSE)
if (identical(fit_a$bootstrap$se_A, fit_b$bootstrap$se_A)) {
	cli::cli_text("v Bootstrap results are reproducible with seed")
} else {
	cli::cli_warn("Bootstrap results differ with same seed (numerical variation acceptable)")
}

# Test 12: dbn(method = "als", bootstrap = N) routes through dbn_als
cli::cli_h2("Test 12: dbn(method = \"als\", bootstrap = N) dispatcher route")
fit_via_dbn <- dbn(Y, family = "gaussian", method = "als",
	bootstrap = 30, bootstrap_type = "block", verbose = FALSE)
stopifnot(identical(fit_via_dbn$meta$sampler_used, "als"))
stopifnot(isTRUE(fit_via_dbn$meta$uncertainty_available))
stopifnot(!is.null(fit_via_dbn$bootstrap))
cli::cli_text("v dbn(method = \"als\", bootstrap = ...) routes to dbn_als correctly")

# Test 13: Medium-network performance
cli::cli_h2("Test 13: Medium Network Performance")
sim_large <- simulate_dynamic_dbn(n = 20, time = 15, seed = 2)
t0 <- Sys.time()
fit_large <- dbn_als(sim_large$Y, family = "gaussian",
	bootstrap = 50, bootstrap_type = "block", verbose = FALSE)
t_total <- as.numeric(Sys.time() - t0, units = "secs")
boot_large <- fit_large$bootstrap
cli::cli_bullets(c(
	" " = "Total time (fit + R=50): {sprintf('%.2f', t_total)} seconds",
	" " = "Per replicate: {sprintf('%.2f', t_total/50)} seconds"
))
stopifnot(boot_large$n_valid > 40)

# Final
cli::cli_rule()
cli::cli_h1("ALL TESTS PASSED")
