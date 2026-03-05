# dbn: Dynamic Bilinear Network Models

An R package for analyzing temporal relational data using Dynamic
Bilinear Network (DBN) models. Supports ordinal, Gaussian, and binary
outcomes across static, dynamic, low-rank, and HMM model variants for
unipartite, bipartite, and symmetric (undirected) networks.

## Installation

### Prerequisites

This package uses OpenMP for parallel computing. Setup by platform:

**macOS:**

1.  Install Xcode Command Line Tools:

``` bash
xcode-select --install
```

2.  Install gfortran (includes OpenMP support):

- For R 4.5.0+: Download
  [gfortran-14.2-universal.pkg](https://github.com/R-macos/gfortran-for-macOS/releases)
- For R 4.4.x and earlier: Download
  [gfortran-12.2-universal.pkg](https://mac.r-project.org/tools/gfortran-12.2-universal.pkg)
- Double-click the downloaded .pkg file to install

**Linux:**

- OpenMP is typically pre-installed with GCC
- If compilation fails, install:

``` bash
# Ubuntu/Debian
sudo apt-get install libomp-dev

# RHEL/CentOS/Fedora
sudo yum install libomp-devel
```

**Windows:**

- Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  matching your R version
- OpenMP is included automatically

### Install Package

``` r
# Install from GitHub
devtools::install_github("netify-dev/dbn")
```

### Troubleshooting

If installation fails, verify your setup:

``` r
# Test if compilers work
install.packages("minpack.lm", type = "source")

# Check your R version
R.version.string
```

**Common macOS Issues:**

- If you see “gfortran not found”, add to your `~/.zshrc`:

``` bash
export PATH="/opt/gfortran/bin:$PATH"
```

Then restart your terminal and R.

- For persistent OpenMP errors on macOS, you may need to create
  `~/.R/Makevars` with appropriate compiler flags. Here an example:

``` bash
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

``` r
library(dbn)

# Set number of threads for parallel computation (optional)
set_dbn_threads(min(4, parallel::detectCores(logical = FALSE)))

# Simulate ordinal network data
sim <- simulate_static_dbn(n = 15, p = 1, time = 10, seed = 42)

# Fit static model
fit <- dbn(
  sim$Y,
  model = "static",
  family = "ordinal",
  nscan = 1000,
  burn = 500,
  odens = 1
)

# Diagnostics
summary(fit)
plot(fit)
check_convergence(fit)
plot_trace(fit)

# Posterior predictive checks
ppd <- posterior_predict_dbn(fit, ndraws = 100)
plot_ppc_ecdf(fit, ppd, Y_obs = sim$Y)
```

## Models

All models share the bilinear form
`Theta_t = A_t * Theta_{t-1} * B_t' + M + noise`.

### Static DBN (`model = "static"`)

- Fixed sender/receiver effects (A, B constant)
- Good for smaller datasets or cross-sectional data (T=1)

### Dynamic DBN (`model = "dynamic"`)

- Time-varying A_t and B_t matrices estimated via FFBS
- Auto-selects time-thinning for large networks to manage memory
- Optional AR(1) dynamics for A and B

### Low-Rank DBN (`model = "lowrank"`)

- `A_t = U diag(alpha_t) U'` with U on the Stiefel manifold
- Better scalability for networks with many nodes (50+)

### HMM DBN (`model = "hmm"`)

- Regime-switching: A_t and B_t selected from R discrete regimes
- Regime transitions governed by a Markov chain with estimated
  transition matrix
- Use
  [`plot_regime_probs()`](https://netify-dev.github.io/dbn/reference/plot_regime_probs.md)
  to visualize regime assignments

### Symmetric Networks (`symmetric = TRUE`)

All models except low-rank support a `symmetric = TRUE` option that
constrains B = A, appropriate for undirected networks where sender and
receiver dynamics are identical. Requires square networks (n_row ==
n_col).

## Outcome Families

- `family = "ordinal"`: Ordered categorical data via rank likelihood
- `family = "gaussian"`: Continuous data with estimated observation
  variance
- `family = "binary"`: Binary (0/1) data via probit link with data
  augmentation

## Main Functions

### Model Fitting

- [`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md): Main
  wrapper for model estimation (dispatches to
  static/dynamic/lowrank/hmm)
- [`dbn_static()`](https://netify-dev.github.io/dbn/reference/dbn_static.md),
  [`dbn_dynamic()`](https://netify-dev.github.io/dbn/reference/dbn_dynamic.md),
  [`dbn_lowrank()`](https://netify-dev.github.io/dbn/reference/dbn_lowrank.md),
  [`dbn_hmm()`](https://netify-dev.github.io/dbn/reference/dbn_hmm.md):
  Direct model functions

### Simulation

- [`simulate_static_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_static_dbn.md):
  Simulate from static model
- [`simulate_dynamic_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_dynamic_dbn.md):
  Simulate from dynamic model
- [`simulate_lowrank_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_lowrank_dbn.md):
  Simulate from low-rank model
- [`simulate_hmm_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_hmm_dbn.md):
  Simulate from HMM model

### Posterior Analysis

- [`posterior_predict_dbn()`](https://netify-dev.github.io/dbn/reference/posterior_predict_dbn.md):
  Posterior predictive samples
- [`param_summary()`](https://netify-dev.github.io/dbn/reference/param_summary.md):
  Summarize scalar parameter traces
- [`theta_summary()`](https://netify-dev.github.io/dbn/reference/theta_summary.md):
  Summarize latent Theta arrays
- [`theta_credible()`](https://netify-dev.github.io/dbn/reference/theta_credible.md):
  Dyad-level credible intervals over time
- [`latent_summary()`](https://netify-dev.github.io/dbn/reference/latent_summary.md):
  Summarize baseline mean M
- [`theta_slice()`](https://netify-dev.github.io/dbn/reference/theta_slice.md):
  Extract specific Theta posterior draws
- [`derive_draws()`](https://netify-dev.github.io/dbn/reference/derive_draws.md):
  Compute derived quantities from posterior
- [`network_summary()`](https://netify-dev.github.io/dbn/reference/network_summary.md):
  Network-level statistics (mean, density, strength) with CIs
- [`edge_prob()`](https://netify-dev.github.io/dbn/reference/edge_prob.md):
  Posterior probability of positive edges
- [`regime_probs()`](https://netify-dev.github.io/dbn/reference/regime_probs.md):
  Extract HMM regime probabilities

### Visualization

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) /
  [`summary()`](https://rdrr.io/r/base/summary.html) /
  [`print()`](https://rdrr.io/r/base/print.html): S3 methods for dbn
  objects
- [`plot_trace()`](https://netify-dev.github.io/dbn/reference/plot_trace.md):
  Parameter trace plots with running means
- [`plot_theta()`](https://netify-dev.github.io/dbn/reference/plot_theta.md):
  Heatmap of posterior mean Theta
- [`plot_ppc_ecdf()`](https://netify-dev.github.io/dbn/reference/plot_ppc_ecdf.md)
  /
  [`plot_ppc_density()`](https://netify-dev.github.io/dbn/reference/plot_ppc_density.md):
  Posterior predictive checks
- [`plot_regime_probs()`](https://netify-dev.github.io/dbn/reference/plot_regime_probs.md):
  HMM regime probability areas
- [`plot_group_influence()`](https://netify-dev.github.io/dbn/reference/plot_group_influence.md):
  Group-level influence trajectories
- [`dyad_path()`](https://netify-dev.github.io/dbn/reference/dyad_path.md):
  Individual dyad trajectory plots
- [`net_snapshot()`](https://netify-dev.github.io/dbn/reference/net_snapshot.md):
  Network visualization at a time point

### Forecasting and IRF

- [`predict()`](https://rdrr.io/r/stats/predict.html): Multi-step-ahead
  forecasting
- [`compute_irf()`](https://netify-dev.github.io/dbn/reference/compute_irf.md):
  Impulse response functions
- [`build_shock()`](https://netify-dev.github.io/dbn/reference/build_shock.md):
  Construct shock matrices for IRF analysis

### Utilities

- [`estimate_memory()`](https://netify-dev.github.io/dbn/reference/estimate_memory.md):
  RAM usage estimation before fitting
- [`check_convergence()`](https://netify-dev.github.io/dbn/reference/check_convergence.md):
  MCMC convergence diagnostics
- [`set_dbn_threads()`](https://netify-dev.github.io/dbn/reference/set_dbn_threads.md)
  /
  [`get_dbn_threads()`](https://netify-dev.github.io/dbn/reference/get_dbn_threads.md):
  OpenMP thread control
- [`role_trajectory()`](https://netify-dev.github.io/dbn/reference/role_trajectory.md):
  Track sender/receiver role changes

## Data Format

Input data should be a 4-dimensional array:

- Dimension 1: Sender actors (n_row)
- Dimension 2: Receiver actors (n_col, can differ from n_row for
  bipartite networks)
- Dimension 3: Relation types
- Dimension 4: Time points

3D arrays (actors x actors x time) are auto-promoted to 4D with a single
relation.

For cross-sectional data (T=1),
[`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md)
automatically selects the static model.

## Vignettes

After installation, browse the vignettes for detailed tutorials:

``` r
browseVignettes("dbn")
```

- **Getting started: static model** –
  [`vignette("static_dbn")`](https://netify-dev.github.io/dbn/articles/static_dbn.md)
- **Dynamic bilinear network models** –
  [`vignette("dynamic_dbn")`](https://netify-dev.github.io/dbn/articles/dynamic_dbn.md)
  (includes bipartite example)
- **Regime-switching (HMM) models** –
  [`vignette("hmm_dbn")`](https://netify-dev.github.io/dbn/articles/hmm_dbn.md)
- **Low-rank models** –
  [`vignette("lowrank_dbn")`](https://netify-dev.github.io/dbn/articles/lowrank_dbn.md)
- **Impulse response analysis** –
  [`vignette("impulse_response")`](https://netify-dev.github.io/dbn/articles/impulse_response.md)

## Author

- [Tosin
  Salau](https://polisci.msu.edu/people/directory/salau-tosin.html)
  <salaubol@msu.edu>
- [Shahryar Minhas](https://s7minhas.com) <sminhas@msu.edu>

## License

MIT License - see [LICENSE](https://netify-dev.github.io/dbn/LICENSE)
file for details.
