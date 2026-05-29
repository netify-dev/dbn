# dbn: Dynamic Bilinear Network Models

An R package for analyzing temporal relational data using Dynamic Bilinear Network (DBN) models. Supports ordinal, Gaussian, and binary outcomes across static, dynamic, piecewise, low-rank, and HMM model variants for unipartite, bipartite, and symmetric (undirected) networks. Dyadic, sender, and receiver covariates and random actor effects are handled by a block-Gibbs sampler.

### Scope of the package

The DBN model is **a latent-state model of the network's own dynamics**. The
default predictor of `Theta_t` is the lagged latent state `Theta_{t-1}` through
the bilinear operator `A_t Theta_{t-1} B_t^T`. Two things to know:

- **Exogenous covariates are supported.** Dyadic, sender, and receiver
  covariates and random actor effects ride on top of the bilinear dynamics
  via a block-Gibbs sampler. Build a covariate object with
  `dbn_covariates()` and pass it as `covariates = ...` to `dbn()`. Priors,
  time-varying coefficients, and actor effects are all configurable; see
  the "Covariates and actor effects" section below.
- **MCMC is single-chain conjugate Gibbs by default.** A single `dbn()`
  call runs one chain; for multi-chain inference (and rank-normalised Rhat
  / ESS) use `dbn_multichain(Y, chains = K, seeds = ...)`, then
  `check_convergence(mc)` for cross-chain diagnostics or
  `as_draws(mc)` for the `posterior` / `bayesplot` workflow. The
  per-chain output still reports Geweke z-scores and effective sample
  size.

## Installation

### Prerequisites

This package uses OpenMP for parallel computing. Setup by platform:

**macOS:**

1. Install Xcode Command Line Tools:

```bash
xcode-select --install
```

2. Install gfortran (includes OpenMP support):

- For R 4.5.0+: Download [gfortran-14.2-universal.pkg](https://github.com/R-macos/gfortran-for-macOS/releases)
- For R 4.4.x and earlier: Download [gfortran-12.2-universal.pkg](https://mac.r-project.org/tools/gfortran-12.2-universal.pkg)
- Double-click the downloaded .pkg file to install

**Linux:**

- OpenMP is typically pre-installed with GCC
- If compilation fails, install:

```bash
# Ubuntu/Debian
sudo apt-get install libomp-dev

# RHEL/CentOS/Fedora
sudo yum install libomp-devel
```

**Windows:**

- Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version
- OpenMP is included automatically

### Install Package

```r
# Recommended: install with the vignettes (the main learning resource)
devtools::install_github("netify-dev/dbn", build_vignettes = TRUE)

# Faster install without vignettes
devtools::install_github("netify-dev/dbn")
```

The vignettes are the primary tutorial material. Install with
`build_vignettes = TRUE`, then list them with `browseVignettes("dbn")` or open
one with `vignette("static_dbn")`. If you installed without that flag,
`browseVignettes("dbn")` returns nothing; reinstall with the flag, or read
the `.Rmd` sources in the package's `vignettes/` directory.

### Troubleshooting

If installation fails, verify your setup:

```r
# Test if compilers work
install.packages("minpack.lm", type = "source")

# Check your R version
R.version.string
```

**Common macOS Issues:**

- If you see "gfortran not found", add to your `~/.zshrc`:

```bash
export PATH="/opt/gfortran/bin:$PATH"
```

Then restart your terminal and R.

- For persistent OpenMP errors on macOS, you may need to create `~/.R/Makevars` with appropriate compiler flags. Here an example:

```bash
# gfortran 14.2 configuration for R 4.5.1 on ARM64 Mac
FC = /opt/gfortran/bin/gfortran
F77 = /opt/gfortran/bin/gfortran

# Use the actual path where libgfortran.dylib is located
FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0 -L/opt/gfortran/lib -lgfortran -lquadmath -lm

# Compiler flags
FCFLAGS = -mmacosx-version-min=11.0
FFLAGS = -mmacosx-version-min=11.0

# brew
CPPFLAGS = -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib
```

## Quick Start

```r
library(dbn)

# Optional: set the number of threads for the C++ backend
set_dbn_threads(min(4, parallel::detectCores(logical = FALSE)))

# Simulate a 10-period dynamic network on 15 actors
# (sim$Z is the continuous series, sim$Y the ordinal one)
sim <- simulate_dynamic_dbn(n = 15, p = 1, time = 10, seed = 42)

# Fit the dynamic model (FFBS-based Gibbs sampler)
fit <- dbn(
  sim$Z,
  model  = "dynamic",
  family = "gaussian",
  nscan  = 2000,
  burn   = 1000,
  odens  = 2
)

# Diagnostics
summary(fit)
check_convergence(fit)
plot_trace(fit)

# Posterior predictive checks
ppd <- posterior_predict_dbn(fit, draws = 100)
plot_ppc_ecdf(fit, ppd, Y_obs = sim$Z)

# H-step-ahead forecast with credible intervals
fc <- predict(fit, H = 3, draws = 200, level = 90)
summary(fc)             # per-horizon mean / lower / upper
plot(fc)                # ggplot facets per (sender, receiver) cell

# Impulse response: shock actor 1's tie to actor 2 by +1, propagate forward
shock <- build_shock(m = 15, type = "unit_edge", i = 1, j = 2)
irf   <- compute_irf(fit, shock = shock, H = 4, n_draws = 100,
                     stat_fun = dbn::stat_density)
print(irf)
```

### Fast point estimate (no MCMC)

For exploratory work, `dbn_als()` returns a point estimate of the same
operator structure in seconds rather than minutes. Pair it with
`bootstrap = N` for entry-wise standard errors. Downstream accessors
(`compute_irf()`, `predict()`, `dbn_compute_snr()`, ...) work uniformly
on MCMC and ALS+bootstrap fits.

```r
als_fit <- dbn_als(sim$Z, family = "gaussian", bootstrap = 200)
coef(als_fit, "A")       # point estimate of the influence operator
confint(als_fit)         # bootstrap CIs
compute_irf(als_fit, shock = shock, H = 4, stat_fun = dbn::stat_density)
```

### Covariates and actor effects

Dyadic, sender, and receiver covariates are wrapped into a
`dbn_covariates` object and handed to `dbn()`. The block-Gibbs sampler
keeps the bilinear dynamics on the residual state and updates the
covariate coefficients in a conjugate Gaussian block. Random actor
effects (sender, receiver, or both) and time-varying coefficients are
configurable, and the prior on the coefficients is exposed via
`prior_beta_scale` and `prior_kind`.

```r
covars <- dbn_covariates(
  dyad = list(distance = D_array),    # [n, n, T] or [n, n]
  row  = list(gdp = gdp_matrix)       # [n, T] or [n]; senders
)

fit_cov <- dbn(
  Y, model = "dynamic", family = "gaussian",
  covariates    = covars,
  actor_effects = "both",             # random sender + receiver effects
  prior_beta_scale = 2.5, prior_kind = "rw",
  nscan = 3000, burn = 1500
)

tidy(fit_cov)              # coefficient summaries
predict(fit_cov, H = 2)    # forecasts that include the covariate path
```

## Models

All models share the bilinear AR form, parameterised in deviation from the baseline `M`:

`Theta_t = M + A_t * (Theta_{t-1} - M) * B_t' + noise`

(equivalently, the centred state `Theta_t - M` follows an AR(1) with operators
`A_t` (left) and `B_t` (right)). The simulator and the FFBS sampler both use
this deviation form; `M` is the long-run baseline, and the time-varying
operators `A_t`, `B_t` govern how perturbations propagate forward.

### Static DBN (`model = "static"`)

- Fixed sender/receiver effects (A, B constant)
- Good for smaller datasets or cross-sectional data (T=1)

### Dynamic DBN (`model = "dynamic"`)

- Time-varying A_t and B_t matrices estimated via FFBS (default) or
  exact PCG on the symmetric path. Choose via `sampler = "auto"` /
  `"approx"` / `"exact"`.
- Optional AR(1) persistence with `ar1 = TRUE` and `update_rho = TRUE`
  to estimate `rhoA`, `rhoB` from the data
- Auto-selects time-thinning for large networks to manage memory

### Low-Rank DBN (`model = "lowrank"`)

- `A_t = U diag(alpha_t) U'` with U on the Stiefel manifold
- Better scalability for networks with many nodes (50+)

### HMM DBN (`model = "hmm"`)

- Regime-switching: A_t and B_t selected from R discrete regimes
- Regime transitions governed by a Markov chain with estimated transition matrix
- Use `plot_regime_probs()` to visualize regime assignments

### Piecewise DBN (`model = "piecewise"`)

- Operator constant within known time blocks, free to jump at block boundaries
- For data with a known structural break (a regime change at a known date)
- See `vignette("piecewise_dbn")`

### Symmetric Networks (`symmetric = TRUE`)

For undirected networks where sender and receiver dynamics are identical
(`B = A`), pass `symmetric = TRUE`. Supported on the **static, dynamic,
piecewise, and HMM** model paths; only the low-rank model rejects it with
an informative error (the Tucker factorisation does not yet implement the
constraint). All paths require a square network (`n_row == n_col`).

The constraint is two-fold on every supported path:

- `B_t = A_t` is enforced at every sampling iteration, so the latent state
  evolves as `Theta_t = A_t %*% Theta_{t-1} %*% t(A_t) + noise`.
- `A_t = t(A_t)` is enforced at storage: the package averages `A_t` with
  its transpose at each iteration, so the stored operator is itself
  symmetric.

The second constraint makes the diagonal-penalty constraint-counting
argument rigorous (Proposition 5 of the companion methods paper).
`fit$symmetric` is set to `TRUE` whenever this path is taken (it mirrors
`fit$dims$is_symmetric`).

## Outcome Families

- `family = "ordinal"`: Ordered categorical data via rank likelihood
- `family = "gaussian"`: Continuous data with estimated observation variance
- `family = "binary"`: Binary (0/1) data via probit link with data augmentation

## Main Functions

### Model Fitting

- `dbn()`: Main wrapper (dispatches to static/dynamic/piecewise/lowrank/hmm).
  Pass `method = "als"` to route to the fast point-estimate path. Pass
  `sampler = "auto" | "approx" | "exact"` to choose the underlying
  sampler explicitly; `"auto"` (default) picks the best for the model.
- `dbn_static()`, `dbn_dynamic()`, `dbn_piecewise()`, `dbn_lowrank()`, `dbn_hmm()`: Direct MCMC model functions
- `dbn_multichain()`: Run `K` independent chains with seed control;
  returns a `dbn_multichain` object that `check_convergence()` and
  `as_draws()` operate on for cross-chain diagnostics.
- `dbn_als()`: Fast alternating-least-squares estimator. Supports a static
  point estimate (`lambda = "static"`, default; the SIR estimator of Minhas
  & Hoff 2025) and an RW(1)-smoothed time-varying estimator with
  `lambda = "cv"` (rolling-origin CV, 1-SE rule), `lambda = "eb"`
  (empirical-Bayes pilot), `lambda = "path"` (full λ-path), or any positive
  numeric. Pair with `bootstrap = N` to attach bootstrap CIs in the same call.
- `dbn_covariates()`: Build a covariates object (`dyad`, `row`, `col`,
  or `actor` slots) to pass into `dbn()` as `covariates = ...`. Random
  actor effects are toggled separately via `dbn(..., actor_effects = ...)`.

### Simulation

- `simulate_static_dbn()`: Simulate from static model
- `simulate_dynamic_dbn()`: Simulate from dynamic model
- `simulate_lowrank_dbn()`: Simulate from low-rank model
- `simulate_hmm_dbn()`: Simulate from HMM model
- `simulate_piecewise_dbn()`: Simulate from piecewise model

### Actor-Level Diagnostics

- `dbn_coupling_rank_probs()`: Posterior rank probabilities for actor coupling
  (top-K and pairwise dominance), for symmetric dynamic fits
- `dbn_leverage()`: Actor leverage, measuring how widely a shock involving the actor
  propagates over a finite horizon (distinct from coupling)
- `dbn_operator()`: Posterior-mean time-varying operator `W_t = A_t B_t'` plus
  a per-period stability (operator-gain) summary
- `dbn_compute_snr()`: Posterior signal-to-residual audit for a symmetric fit
- `dbn_identifiability_diagnostics()`: Jacobian-rank identification checks

### Posterior Analysis

- `posterior_predict_dbn()`: Posterior predictive samples
- `as_draws.dbn()` / `as_draws.dbn_multichain()`: Convert a fit into a
  `posterior::draws_array` for `posterior`, `tidybayes`, `bayesplot`
- `tidy(fit)`, `glance(fit)`, `augment(fit)`: `broom`-style accessors
- `coef()`, `confint()`, `fitted()`, `residuals()`, `vcov()` where defined
- `param_summary()`: Summarize scalar parameter traces
- `theta_summary()`: Summarize latent Theta arrays
- `theta_credible()`: Dyad-level credible intervals over time
- `latent_summary()`: Summarize baseline mean M
- `theta_slice()`: Extract specific Theta posterior draws
- `derive_draws()`: Compute derived quantities from posterior
- `network_summary()`: Network-level statistics (mean, density, strength) with CIs
- `edge_prob()`: Posterior probability of positive edges
- `regime_probs()`: Extract HMM regime probabilities
- `gof()`: Per-time posterior-predictive interval coverage of network
  statistics (density, reciprocity, transitivity, mean in/out degree)
- `fevd()`: Forecast-error variance decomposition over a chosen horizon

### Visualization

- `plot()` / `summary()` / `print()`: S3 methods for dbn objects
- `plot_trace()`: Parameter trace plots with running means
- `plot_theta()`: Heatmap of posterior mean Theta
- `plot_ppc_ecdf()` / `plot_ppc_density()`: Posterior predictive checks
- `plot_regime_probs()`: HMM regime probability areas
- `plot_group_influence()`: Group-level influence trajectories
- `dyad_path()`: Individual dyad trajectory plots
- `net_snapshot()`: Network visualization at a time point

### Forecasting and IRF

- `predict(fit, H = h)`: H-step-ahead forecast; returns a `dbn_forecast`
  array. Pass `level = N` (percent) or `summary = "ci"` for credible
  intervals; the return is then a `dbn_forecast_ci` list with `$mean`,
  `$lower`, `$upper`. Built-in `print()`, `summary()`, `plot()`, and
  `autoplot()` (via ggplot2) methods.
- `forecast(fit, h = h)`: forecast-package convention; identical
  return shape, accepts `H` as an alias.
- `forecast_gain()`: In-sample one-step predictive gain per actor over
  a chosen window
- `forecast_with_covariates()`: Roll-forward forecasts for a
  covariate-equipped fit
- `compute_irf()`: Impulse response functions; works on dynamic and
  piecewise fits (piecewise is auto-expanded to per-time operators).
- `irf_longrun()`: Closed-form long-run cumulative IRF for static,
  dynamic, piecewise, and HMM fits.
- `build_shock()`: Construct shock matrices for IRF analysis
  (unit-edge, node-out, full-row, custom).
- `debug_irf()`: Per-draw IRF diagnostic; auto-expands piecewise fits.

### Utilities

- `estimate_memory()`: RAM usage estimation before fitting
- `check_convergence()`: MCMC convergence diagnostics (single chain or
  multichain; reports Rhat / bulk-ESS / tail-ESS when given a
  `dbn_multichain` and `posterior` is installed)
- `compare_dbn()` / `compare_blocks()` / `compare_samplers()`: pairwise
  model and sampler comparisons
- `sampler_describe()`: One-line description of how a fit was produced
- `set_dbn_threads()` / `get_dbn_threads()`: OpenMP thread control
- `role_trajectory()`: Track sender/receiver role changes
- `as_dbn_array()` / `as_dbn_array_edgelist()`: Coerce `network`,
  `networkDynamic`, igraph panels, or tidy edgelists into the 4D array
  layout `dbn()` expects

## Data Format

Input data should be a 4-dimensional array:

- Dimension 1: Sender actors (n_row)
- Dimension 2: Receiver actors (n_col, can differ from n_row for bipartite networks)
- Dimension 3: Relation types
- Dimension 4: Time points

3D arrays (actors x actors x time) are auto-promoted to 4D with a single relation.

For cross-sectional data (T=1), `dbn()` automatically selects the static model.

## Vignettes

The vignettes are built only when the package is installed with
`build_vignettes = TRUE` (see Install Package above); afterwards browse them
with:

```r
browseVignettes("dbn")
```

- **Methodology and estimands**: `vignette("methodology")`
- **Getting started: static model**: `vignette("static_dbn")`
- **Dynamic bilinear network models**: `vignette("dynamic_dbn")` (includes bipartite example)
- **Piecewise (structural-break) models**: `vignette("piecewise_dbn")`
- **Regime-switching (HMM) models**: `vignette("hmm_dbn")`
- **Impulse response analysis**: `vignette("impulse_response")`
- **Applied impulse-response walk-through**: `vignette("applied_ir")`

## Author

- [Tosin Salau](https://polisci.msu.edu/people/directory/salau-tosin.html) <salaubol@msu.edu>
- [Shahryar Minhas](https://s7minhas.com) <minhassh@msu.edu>


## License

MIT License - see [LICENSE](LICENSE) file for details.
