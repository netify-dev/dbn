# Code Review Prompt: Critical Fixes to Dynamic Bilinear Network Package

## Context

You are reviewing four critical fixes to an R package for Bayesian inference on dynamic bilinear network models. The package implements:
- **Core model**: Θ_t = A_t Θ_{t-1} B_t^T + M + ε_t where A_t, B_t are time-varying influence matrices
- **Samplers**: Preconditioned conjugate gradient (PCG), forward-filtering backward-sampling (FFBS), alternating least squares (ALS)
- **Data types**: Gaussian (continuous), Ordinal (rank-based), Binary (probit)
- **Network types**: Unipartite (n × n square) and bipartite (n_row × n_col rectangular)

The package uses:
- **Frontend**: Pure R with S3 method dispatch, vectorized operations
- **Computational backend**: Rcpp + RcppArmadillo for likelihood, sampling, matrix operations
- **Statistical approach**: Bayesian MCMC with data-driven priors

## The Four Fixes

### FIX #1: Singular Matrix Handling in Ordinal ALS

**Problem**: Ordinal data requires latent score augmentation via shared_preprocess(). The resulting Z matrices are transformed continuous values from the rank likelihood. When fitting alternating least squares (ALS) model Θ = AΘ_{t-1}B^T, the row-wise OLS solves:
```
a_i^* = argmin_a ||Z_i,: - a^T (B ⊙ Z_{t-1})^T||^2
```
produce rank-deficient design matrices X_design = [Z_{t-1} ⊗ B]^T with condition numbers 50-200+.

**Current fix**: Added try-catch in ALS loops (dbn_explore, als_refit_warmstart):
```r
tryCatch({
  solve(crossprod(X_design), crossprod(X_design, y_rhs))
}, error = function(e) {
  MASS::ginv(XtX) %*% Xty  # Moore-Penrose pseudoinverse
})
```

**Test results**:
- Gaussian: condition κ(A) ≈ 14.4 (normal range)
- Ordinal integer: κ(A) ≈ 108.9 (ill-conditioned, solve fails)
- Binary: All-zero A matrix from pseudoinverse (extreme rank deficiency)
- All results numerically stable (no NaN/Inf propagation)

**Questions for reviewer**:
1. **Regularization vs pseudoinverse tradeoff**: The pseudoinverse gives the minimum-norm least-squares solution, but doesn't penalize magnitude. Would Ridge regularization (λI term) be theoretically preferable for this problem? How would you choose λ?
2. **Warm-start initialization bias**: When pseudoinverse is used, does the minimum-norm solution bias subsequent bootstrap replicates? Should we add perturbation or use different initialization for replicates 2+?
3. **Diagnostic reporting**: Should we flag when pseudoinverse is used (high-condition matrix) and report it in sampler_describe() output for transparency?
4. **Ordinal preprocessing**: Is the rank transformation in shared_preprocess() the bottleneck? Would standardizing Z scores before ALS (zscore normalization) reduce conditioning?

---

### FIX #2: S3 Method Dispatch for sampler_describe

**Problem**: In development environment (load_dev.R), R functions sourced into isolated namespace with parent=getNamespace("dbn"). S3 method dispatch failed because:
- Generic `sampler_describe()` defined in sourced env
- Methods `sampler_describe.dbn()`, `sampler_describe.dbn_boot()` also in sourced env
- UseMethod() search couldn't resolve method table in locked namespace

**Current fix** (load_dev.R):
```r
sampler_describe_wrapper <- function(object, verbose = TRUE) {
  if (inherits(object, "dbn_boot")) {
    return(sampler_describe.dbn_boot(object, verbose = verbose))
  } else if (inherits(object, "dbn")) {
    return(sampler_describe.dbn(object, verbose = verbose))
  }
}
assign("sampler_describe", sampler_describe_wrapper, envir = .GlobalEnv)
```

**Test results**:
- ALS fit (class dbn, sampler_used="als"): Output correct ✓
- Block bootstrap (class dbn_boot, type="block"): Output correct ✓
- Parametric bootstrap (class dbn_boot, type="parametric"): Output correct ✓
- Multiple calls stable (no dispatch failure on retry) ✓

**Questions for reviewer**:
1. **S3 vs NextMethod()**: Should we be using NextMethod() to chain methods (e.g., sampler_describe.dbn_boot calls sampler_describe.default then adds bootstrap-specific info)? Current approach duplicates baseline info.
2. **Scope of wrapper**: Is this wrapper-based solution appropriate for development only, or should it be the production approach? How does package-side S3 registration (Roxygen @export, NAMESPACE) differ?
3. **Method robustness**: What if a new sampler type is added (e.g., exact PCG for asymmetric)? Should we refactor dispatch to use a lookup table (list of methods by class) instead of explicit if-else?
4. **Namespace initialization order**: The wrapper assigns at load time. What if dbn_explore or dbn_explore_bootstrap haven't been defined yet? Should dispatch be lazy (check at call time)?

---

### FIX #3: Package Initialization & C++ Loading

**Problem**: Compiled C++ functions (build_rank_indices, theta_ffbs_*, sample_z_*, etc.) weren't accessible at package load due to failures in .onLoad() hook:
```r
library.dynam("dbn", pkgname, libname)  # Failed silently
```

**Current fix** (zzz.R):
```r
.onLoad <- function(libname, pkgname) {
  tryCatch({
    library.dynam("dbn", pkgname, libname)
  }, error = function(e) {
    so_paths <- c(
      file.path(libname, pkgname, "libs", "dbn.so"),
      file.path(libname, pkgname, "libs", paste0("dbn", .Platform$dynlib.ext)),
      "src/dbn.so",
      "libs/dbn.so"
    )
    for (so_path in so_paths) {
      if (file.exists(so_path)) {
        tryCatch(dyn.load(so_path), error = function(e2) {})
        break
      }
    }
  })
}
```

**Test results**:
- Ordinal rank likelihood (build_rank_indices): ✓ C++ called, 20 draws
- Dynamic FFBS (theta_ffbs_*): ✓ C++ called, Theta shape 6×6×2×8×20
- All families (Gaussian, Ordinal, Binary): ✓ C++ executed
- Numerical validity: ✓ 100% finite results

**Questions for reviewer**:
1. **Error suppression**: Wrapping dyn.load() in error handler silently swallows failures. Should we at least warn user if load fails? What's the failure mode if C++ functions become unavailable mid-session?
2. **Platform-specific paths**: The multi-path fallback works on Linux but relies on file.exists(). What about Windows (.dll) or macOS (.dylib)? Should this use .Platform$dynlib.ext consistently?
3. **Session-level state**: If dyn.load fails once, is there a way to retry? Or is the .onLoad failure permanent for that R session?
4. **Symbol visibility**: Does dyn.load() make all Rcpp-exported symbols visible, or just entry points? How does this interact with .Call() in wrapper functions?

---

### FIX #4: Bipartite Network Support

**Problem**: Package claimed to support bipartite networks (n_row ≠ n_col) but implementation assumptions weren't verified. Concerns:
- A matrix should always be n_row × n_row (sender influence)
- B matrix should always be n_col × n_col (receiver influence)
- M matrix should be n_row × n_col (dyadic baseline)
- ALS, MCMC, and bootstrap all need to handle this correctly

**Current fix**: Verification that existing code properly handles rectangular dims (no code changes needed, but comprehensive testing added).

**Test results**:
- Small bipartite (4×6): ✓ Static, A shape 4×4×5, B shape 6×6×5
- Medium bipartite (8×12): ✓ Static, M shape 8×12×1
- Extreme asymmetry (2×10): ✓ 5:1 ratio handled
- Bipartite + ordinal: ✓ Rank likelihood computed correctly
- Bipartite + dynamic: ✓ Theta shape 4×6×2×5×20 (correct)
- Bipartite + ALS: ✓ Converges, warm-start works
- Bipartite + bootstrap: ✓ SE_A length 16 (correct for 4×4)

**Questions for reviewer**:
1. **Theoretical model validity**: For bipartite networks, do A and B have the same stationarity/identifiability constraints as unipartite? Should we add checks to prevent symmetric=TRUE with n_row ≠ n_col?
2. **Dimension asymmetry implications**: In extreme asymmetry (2 senders, 10 receivers), does the lower-rank A (2×2) matrix create identifiability issues? Should priors adapt to n_row vs n_col?
3. **Baseline parameterization**: M matrix is n_row × n_col, but posterior draws are typically n_row × n_col × draws. Is sparse storage (Θ - A Θ_{t-1} B^T mean) more efficient for large bipartite networks?
4. **Missing diagonal handling**: Unipartite requires diag(Y) = NA. For bipartite, are off-diagonal entries always observed? Should we enforce this in data validation?

---

## Integration & Cross-Fix Concerns

### Singular matrices + Bipartite
When ordinal bipartite network has extreme asymmetry AND sparse events (few unique values), pseudoinverse solution for A (small 2×2) matrix becomes near-zero. Does bootstrap properly propagate uncertainty, or does pseudoinverse mask true posterior variability?

### S3 dispatch + sampler types
Current sampler_describe wrapper checks only for dbn, dbn_boot. If asymmetric PCG (Phase 2) adds dbn_pcg_asym class, will dispatch still work? Should it?

### C++ loading + ordinal + bipartite
Rank likelihood (build_rank_indices) is called for n_row × n_col × T array. Does Rcpp handle this dimension correctly, or are there off-by-one errors in vectorization?

---

## Reviewer Instructions

As a Bayesian statistical computing expert with specialization in optimization of statistical routines, please provide:

1. **Conceptual soundness** (Bayesian perspective):
   - Are pseudoinverse solutions theoretically justified for bilinear ALS models?
   - Does the prior structure in MCMC complement these fixes?
   - Are there identifiability issues introduced by these fixes?

2. **Numerical stability** (optimization perspective):
   - Condition number analysis: when should we use regularization vs pseudoinverse vs iterative refinement?
   - How do these fixes affect convergence rates of ALS and MCMC?
   - Any risk of numerical instability in downstream computations (e.g., IRF, rank probs)?

3. **Implementation quality**:
   - Are the tryCatch blocks sufficient, or should there be fallback strategies?
   - Does error handling mask important failures?
   - Are there memory/efficiency concerns?

4. **Diagnostic transparency**:
   - What should users know about when these fixes are engaged?
   - Should we add verbose flags or diagnostic outputs?
   - Any risk of silent degradation of model fit quality?

5. **Testing adequacy**:
   - Are 28 scenarios sufficient? Any edge cases missed?
   - Should we add stress tests (very large networks, extreme sparsity)?
   - How would you validate that results remain statistically consistent?

6. **Recommended improvements** (actionable, specific):
   - For each fix, suggest the single most important enhancement
   - Prioritize by impact and implementation cost
   - Suggest metrics to monitor in production

Please structure your response by fix number, and flag any concerns that cross fixes or affect the overall statistical validity of the package.
