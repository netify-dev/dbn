# Dynamic Bilinear Network Analysis

Main entry point for fitting Dynamic Bilinear Network (DBN) models.
These models estimate how past network interactions predict future
interactions, recovering time-varying influence structures from temporal
relational data.

The core model is: \\\Theta_t = A_t \Theta\_{t-1} B_t' + M +
\varepsilon_t\\, where \\A_t\\ captures sender influence, \\B_t\\
captures receiver influence, and \\M\\ captures stable dyad-specific
tendencies.

## Usage

``` r
dbn(
  data,
  family = c("ordinal", "gaussian", "binary"),
  model = c("static", "dynamic", "lowrank", "hmm", "piecewise"),
  nscan = 10000,
  burn = 1000,
  odens = 1,
  verbose = TRUE,
  symmetric = FALSE,
  ...
)
```

## Arguments

- data:

  Numeric array of network data, or a file path to an `.RData` file that
  contains an object named `Y`. The array should be 3-dimensional
  `[actors, actors, time]` for a single relation type, or 4-dimensional
  `[actors, actors, relations, time]` for multiple relation types.
  Diagonal entries (self-ties) should be `NA` for unipartite networks.
  For bipartite networks, pass a rectangular array where the first
  dimension (senders) differs from the second (receivers).

- family:

  Character string specifying the outcome distribution. See **Details**
  for guidance on choosing:

  - `"ordinal"`: For ordinal/ranked data (positive integers)

  - `"gaussian"`: For continuous data (any real numbers)

  - `"binary"`: For binary data (0/1 or logical)

- model:

  Character string specifying the model type. See **Details** for
  guidance on choosing:

  - `"static"`: Fixed sender/receiver effects across time

  - `"dynamic"`: Time-varying sender/receiver effects

  - `"lowrank"`: Low-rank factorization of sender effects (large
    networks)

  - `"hmm"`: Regime-switching with data-driven regime discovery

  - `"piecewise"`: Block-constant influence with known break points

- nscan:

  Number of posterior samples to draw after burn-in

- burn:

  Number of initial MCMC samples to discard (warm-up period)

- odens:

  Thinning interval: save every odens-th sample (reduces autocorrelation
  and memory)

- verbose:

  Logical or numeric. If TRUE, show progress. If numeric, print detailed
  info every n iterations (default: TRUE)

- symmetric:

  Logical. If TRUE, enforce B = A (symmetric/undirected network).
  Requires square network (n_row == n_col). Not supported for lowrank
  models. Default: FALSE.

- ...:

  Additional model-specific parameters:

  `r`

  :   Rank for lowrank model (default: 2)

  `R`

  :   Number of regimes for HMM model (default: 3)

  `blocks`

  :   Block specification for piecewise model: integer (number of equal
      blocks), numeric vector (block boundaries), named vector (labeled
      boundaries), or "auto" for automatic selection

  `ar1`

  :   Use AR(1) dynamics for dynamic model (default: FALSE)

  `update_rho`

  :   Update AR coefficient in dynamic model (default: FALSE)

  `seed`

  :   Random seed for reproducibility (default: 6886)

  `previous`

  :   Previous fit object to continue MCMC from

  `init`

  :   List of initial values for parameters

  `time_thin`

  :   Time thinning factor for dynamic/lowrank/HMM (default: auto for
      dynamic, 1 for others)

  `store_z`

  :   Store Z draws for dynamic model (default: auto based on memory)

  `store_theta`

  :   Store full Theta trajectory draws for piecewise model (default:
      TRUE). **Critical for large networks:** Set to FALSE for networks
      with 100+ actors to avoid memory issues. Theta storage scales as
      O(n^2 \* T \* draws) – a 200-actor network with 50 time points and
      500 draws requires ~40 GB. With `store_theta = FALSE`, you retain
      posterior draws for A, B, M and variance parameters,
      [`compare_blocks()`](https://netify-dev.github.io/dbn/reference/compare_blocks.md)
      functionality, and convergence diagnostics. You lose full
      posterior uncertainty on individual Theta entries and
      [`posterior_predict_dbn()`](https://netify-dev.github.io/dbn/reference/posterior_predict_dbn.md)
      with uncertainty propagation.

## Value

A list of class `"dbn"` with model-specific contents. Common elements:

- model:

  Character string indicating which model was fit

- family:

  Character string indicating the outcome family

- dims:

  List of data dimensions (n_row, n_col, p, Tt)

- settings:

  List of MCMC settings used

- Y:

  Original data array

- M:

  Posterior draws for baseline mean M

- Theta:

  Posterior draws for latent network state

Model-specific elements include:

- A:

  Posterior draws for sender influence matrices

- B:

  Posterior draws for receiver influence matrices

- sigma2, tau_A2, tau_B2, g2:

  Posterior draws for variance parameters

- rhoA, rhoB:

  AR(1) persistence parameters (dynamic model with ar1=TRUE)

- A_blocks:

  List of regime-specific posterior mean A matrices (piecewise)

- time_kept:

  Which time indices are stored (dynamic/lowrank/HMM)

Use [`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`param_summary()`](https://netify-dev.github.io/dbn/reference/param_summary.md),
and
[`check_convergence()`](https://netify-dev.github.io/dbn/reference/check_convergence.md)
to inspect results. See model-specific vignettes for full workflows.

## Details

**Choosing a family:**

- `"gaussian"`: Use when your relational data are continuous
  measurements (e.g., trade volumes, similarity scores). This is the
  simplest family and converges fastest.

- `"ordinal"`: Use when your data are ordered categories or counts
  (e.g., conflict severity 1-5, event counts) and you trust the ordering
  but not the exact values. Uses a rank likelihood.

- `"binary"`: Use for presence/absence data (0/1), e.g., whether a tie
  exists. Uses a probit link with data augmentation.

**Choosing a model:**

- `"static"`: Simplest. Influence structure is fixed over time. Good
  starting point and for short time series.

- `"dynamic"`: Influence structure changes over time. Use when you
  expect shifting alliances, evolving trade patterns, etc.

- `"piecewise"`: Influence is constant within known regimes but differs
  across them. Use when you know when structural breaks occurred (e.g.,
  before/after a crisis).

- `"hmm"`: Like piecewise but discovers regimes from data. Use when
  breaks are unknown.

- `"lowrank"`: Like dynamic but with dimensionality reduction for large
  networks (50+ actors).

**MCMC settings:** The sampler draws `nscan` posterior samples after
discarding the first `burn` as warm-up. Setting `odens > 1` thins the
output by saving every k-th sample. For initial exploration,
`nscan = 5000, burn = 2000, odens = 5` is a reasonable starting point.
For publication, use longer chains (`nscan = 10000+`) and verify
convergence with
[`check_convergence()`](https://netify-dev.github.io/dbn/reference/check_convergence.md).

## Examples

``` r
# \donttest{
sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)

# static model with gaussian family
fit <- dbn(sim$Z, model = "static", family = "gaussian",
    nscan = 200, burn = 100, verbose = FALSE)

# dynamic model
fit_dyn <- dbn(sim$Z, model = "dynamic", family = "gaussian",
    nscan = 200, burn = 100, verbose = FALSE)

# piecewise model with 2 blocks
fit_pw <- dbn(sim$Y, model = "piecewise", blocks = 2,
    nscan = 200, burn = 100, verbose = FALSE)
# }
```
