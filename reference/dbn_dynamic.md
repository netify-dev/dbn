# Dynamic DBN MCMC

Fits DBN model with time-varying sender/receiver effects

## Usage

``` r
dbn_dynamic(
  Y,
  family = c("ordinal", "gaussian", "binary"),
  nscan = 10000,
  burn = 1000,
  odens = 1,
  ar1 = FALSE,
  update_rho = FALSE,
  seed = 6886,
  verbose = TRUE,
  time_thin = NULL,
  store_z = NULL,
  previous = NULL,
  init = NULL,
  symmetric = FALSE
)
```

## Arguments

- Y:

  Data array (nodes x nodes x relations x time)

- family:

  Character string specifying the data family/distribution:

  - "ordinal": Ordinal data (ordered categories). Data should be
    positive integers.

  - "gaussian": Continuous data with Gaussian errors. Data can be any
    real numbers.

  - "binary": Binary (0/1) data. Data should be 0/1 or logical values.

- nscan:

  Number of iterations of the Markov chain (beyond burn-in)

- burn:

  Burn-in for the Markov chain

- odens:

  Output density for the Markov chain

- ar1:

  Use AR(1) dynamics instead of random walk (default: FALSE)

- update_rho:

  Update AR coefficients (default: FALSE)

- seed:

  Random seed

- verbose:

  Logical or numeric. TRUE prints every 100 iterations, numeric prints
  every n iterations, FALSE suppresses output.

- time_thin:

  Save every nth time point to reduce memory (default: NULL = auto)

- store_z:

  Whether to store Z draws (default: NULL = auto based on memory)

- previous:

  Previous dbn_dynamic results to continue from (optional)

- init:

  List of initial values: A, B, sigma2, tau_A2, tau_B2, g2, rho_A,
  rho_B, Theta, M, Z (optional)

- symmetric:

  Logical. If TRUE, enforce B = A after each MCMC update. Default:
  FALSE.

## Value

List containing MCMC results

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md) for the main
dispatcher,
[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md)
for posterior summaries

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn_dynamic(sim$Y, nscan = 200, burn = 100, verbose = FALSE)
# }
```
