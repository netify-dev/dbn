# Dynamic Bilinear Network Analysis

Main wrapper function for Dynamic Bilinear Network (DBN) analysis. This
is the primary interface for fitting DBN models to network data. DBN
models capture complex dependencies in network data through bilinear
interactions between latent sender and receiver effects.

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

  Data array or path to .RData file containing Y array. Array should be
  3-dimensional (actors x actors x time) for single relation or
  4-dimensional (actors x actors x relations x time) for multiple
  relations

- family:

  Character string specifying the data family/distribution:

  - "ordinal": For ordinal/ranked data (e.g., ratings 1-5). Data should
    be positive integers.

  - "gaussian": For continuous data. Data can be any real numbers.

  - "binary": For binary data. Data should be 0/1 or logical values.

- model:

  Character string specifying model type:

  - "static": Fixed sender/receiver effects across time

  - "dynamic": Time-varying sender/receiver effects

  - "lowrank": Low-rank factorization of sender effects

  - "hmm": Regime-switching model with hidden Markov states

  - "piecewise": Block-constant influence matrices for structural change

- nscan:

  Number of iterations of the Markov chain (beyond burn-in)

- burn:

  Burn-in for the Markov chain

- odens:

  Output density for the Markov chain (save every odens-th iteration)

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

A list of class "dbn" containing:

- B:

  List of posterior samples for B matrices (static model)

- A:

  List of posterior samples for time-varying A matrices (dynamic model)

- params:

  Matrix of parameter traces (static model)

- sigma2:

  Vector of sigma^2 samples (dynamic model)

- model:

  Character string indicating which model was run

- dims:

  List containing data dimensions

- settings:

  List of model settings used

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
data(example_data)

# Run static model with default settings
results <- dbn(example_data, model = "static")

# Run dynamic model with custom MCMC settings
results <- dbn(example_data,
    model = "dynamic",
    nscan = 5000, burn = 1000, odens = 10
)

# Run HMM model with 3 regimes
results <- dbn(example_data, model = "hmm", R = 3)

# Run low-rank model with rank 2
results <- dbn(example_data, model = "lowrank", r = 2)

# Run quietly without progress output
results <- dbn(example_data, model = "static", verbose = FALSE)

# Run with detailed output every 100 iterations
results <- dbn(example_data, model = "dynamic", verbose = 100)

# Run piecewise model with 4 blocks
results <- dbn(example_data, model = "piecewise", blocks = 4)

# Run piecewise model with specific block boundaries
results <- dbn(example_data, model = "piecewise",
    blocks = c(pre = 25, crisis = 50, post = 100))

# Run piecewise model with automatic block selection
results <- dbn(example_data, model = "piecewise", blocks = "auto")
} # }
```
