# Bootstrap Implementation: Known Issues and Limitations

## Summary

The `dbn_explore_bootstrap()` implementation is **functionally correct and working**, but users should be aware of a fundamental limitation: **ALS can converge to different local optima on different bootstrap samples**, leading to inflated variance estimates.

## Issue: Multimodality in Bootstrap Replicates

### What We Found

When testing bootstrap replicates, we discovered that the magnitude of coefficient estimates can vary dramatically across samples:

```
Simple network (5x5):
  A norm range: [0.33, 62.66]
  Coefficient of variation: 1.83 (extremely high!)

Better network (6x6):
  A norm range: [1.04, 7.08]  
  Coefficient of variation: 0.58 (still concerning)
```

### Root Cause

The ALS algorithm has multiple local optima, especially when applied to resampled data. With different random initializations, the algorithm can converge to solutions with very different magnitudes (scales) that all fit the data reasonably well.

This is **NOT a bug in the bootstrap code** - it's a fundamental property of the ALS optimization landscape.

### Impact on Inference

- Bootstrap standard errors may be **inflated** due to averaging over different modes
- Bootstrap confidence intervals may be **too wide**
- The bootstrap mean may **differ substantially from the point estimate**
- Posterior uncertainty is **not properly quantified**

## Diagnostic Signs

The implementation includes automatic detection of replicates in different modes:

```r
boot <- dbn_explore_bootstrap(fit_als, R = 50, verbose = TRUE)
# Warning: Bootstrap: 5 replicates with large coefficient norms (possible local optima)
```

High variability in Frobenius norms across replicates is a warning sign:
- CV(norms) > 0.5: concerning
- CV(norms) > 1.0: serious multimodality issue

## Recommendations

### 1. For Point Estimation Only
If you only care about point estimates (not uncertainty):
- `dbn_explore()` is **fine** - it provides fast, reasonable point estimates
- Bootstrap adds little value in this case

### 2. For Uncertainty Quantification
**Use MCMC instead:**

```r
# Instead of:
fit_als <- dbn_explore(Y, ...)
boot <- dbn_explore_bootstrap(fit_als, R = 100)

# Use:
fit_mcmc <- dbn(Y, model = "dynamic", sampler = "exact", ...)
sampler_describe(fit_mcmc)  # Credible intervals available
```

The MCMC sampler (PCG for symmetric models, FFBS for asymmetric) properly accounts for posterior uncertainty.

### 3. If You Must Use Bootstrap
- Use **block bootstrap only** (parametric bootstrap is worse due to scale ambiguity)
- Increase `R` substantially (500-1000) to average over noise
- Set `verbose = TRUE` to see multimodality warnings
- **Interpret results with extreme caution**
- Validate against MCMC results

## When Bootstrap Works Well

Bootstrap for ALS tends to work better when:
- **Network is small** (n < 10) and **well-identified**
- **Time series is long** (T > 20) relative to network size
- **Bootstrap samples are less divergent** (less resampling variation)
- **You accept inflated uncertainty bounds** as conservative estimates

Bootstrap tends to **fail** when:
- Network is moderately large (n > 10)
- Time series is short relative to network size
- Resampling creates substantially different dynamics
- Exact uncertainty quantification is required

## Technical Details

### Why ALS is Sensitive to Resampling

The ALS algorithm solves the optimization problem:
```
minimize ||Z_t - A * Z_{t-1} * B'||^2
```

For resampled data:
1. Different time periods have different "energy levels"
2. The resampling can change this energy distribution significantly
3. Multiple scale-equivalent solutions exist: (A, B) ≈ (cA, c^-1B)
4. Different initializations find different modes
5. Bootstrap averaging over modes → inflated variance

### Alignment Attempts

The implementation tries to detect and align divergent replicates using:
- Correlation-based flip detection (sign alignment)
- Frobenius norm outlier detection (warning system)

However, these are **band-aids**, not solutions. Alignment can fail when correlations are moderate (0.4-0.7), and we can't reliably identify the "true" mode.

## Alternative Approaches

### 1. Multiple Random Starts (Not Implemented)
Run ALS from multiple random initializations and keep the best:
```r
# Pseudo-code (not available in dbn)
best_fit <- NULL
best_objective <- Inf
for (start in 1:10) {
  fit <- als(data, random_init = start)
  if (fit$objective < best_objective) {
    best_fit <- fit
    best_objective <- fit$objective
  }
}
```

Benefit: More stable convergence
Drawback: Even more computation

### 2. Procrustes Alignment (Not Implemented)
Align all bootstrap replicates to the original estimate using orthogonal transformations:
```r
# Pseudo-code
for (rep in 1:R) {
  # Find orthogonal P minimizing ||A_orig - A_boot %*% P||_F
  P <- optimal_rotation(A_orig, A_boot)
  A_boot_aligned <- A_boot %*% P
}
```

Benefit: Removes rotation ambiguity
Drawback: Still doesn't address scale ambiguity

### 3. Constrained ALS (Not Implemented)
Fix one element (e.g., A[1,1] = 1) to remove scale ambiguity:
```r
# Pseudo-code
als_constrained <- function(data, constraint_value = 1.0) {
  # Enforce A[1,1] = constraint_value throughout optimization
}
```

Benefit: Removes scale ambiguity entirely
Drawback: Changes the optimization landscape significantly

## References

This issue is related to:
- **Local optima in matrix factorization**: [See bilinear model literature]
- **Bootstrap properties of dependent data**: [See block bootstrap literature]
- **Identifiability in latent variable models**: [See structural equation modeling literature]

## Conclusion

The bootstrap implementation is **correct and useful for diagnostics**, but **NOT recommended for obtaining posterior uncertainty intervals**. Users seeking proper uncertainty quantification should use MCMC methods instead.

The warning system (`verbose = TRUE`) will alert users to detected multimodality. When warnings appear frequently, it's a sign to switch to MCMC.

---

**Date**: May 20, 2026
**Status**: Documented and operational with clear warnings
