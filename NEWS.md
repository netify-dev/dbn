# dbn 0.2.0

## Major features

* Four model types: static, dynamic, low-rank, and HMM (regime-switching).
* Three outcome families: ordinal (rank likelihood), Gaussian, and binary (probit).
* Full bipartite network support for static and dynamic models (`n_row != n_col`).
* OpenMP-parallelized A/B row/column updates with thread-safe RNG.
* Impulse response analysis (`compute_irf()`, `build_shock()`).
* Posterior predictive distribution generation and PPC plots.
* Memory-aware dynamic model: auto time-thinning, optional Z storage, `estimate_memory()`.
* Warm start support via the `previous` argument for continuing MCMC chains.

## Infrastructure

* C++ backend via Rcpp/RcppArmadillo with safe matrix operations (regularized inversions, Cholesky fallbacks).
* Comprehensive test suite (1200+ tests).
* Five vignettes covering all model types plus impulse response analysis.
* pkgdown documentation site.

## Known limitations

* HMM and low-rank models are still under construction.
* Dynamic binary models may encounter numerical singularities with small networks (n < 15).
