# Changelog

## dbn 1.0.0

### Models

- Static, dynamic, and piecewise-static bilinear network models.
- Low-rank and HMM (regime-switching) variants available but still under
  development.
- Three outcome families: ordinal (rank likelihood), Gaussian, binary
  (probit).
- Bipartite network support for all models.

### Features

- Impulse response analysis for counterfactual shock propagation.
- Block comparison for piecewise models via
  [`compare_blocks()`](https://netify-dev.github.io/dbn/reference/compare_blocks.md).
- Posterior predictive checks and visualization.
- OpenMP parallelization for row/column FFBS updates.
- Memory-aware fitting with
  [`estimate_memory()`](https://netify-dev.github.io/dbn/reference/estimate_memory.md)
  and `store_theta` options.
- Warm start support for continuing MCMC chains.

### Infrastructure

- C++ backend with numerical stability safeguards.
- Comprehensive test suite.
- Eight vignettes covering methodology and all model types.
- pkgdown documentation site.
