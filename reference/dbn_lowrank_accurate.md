# DBN Low-rank Model (Accurate)

Fits DBN with low-rank sender effects using batch FFBS

## Usage

``` r
dbn_lowrank_accurate(
  Y,
  family = c("ordinal", "gaussian", "binary"),
  r = 2,
  nscan = 10000,
  burn = 1000,
  odens = 1,
  ar1_alpha = TRUE,
  update_rho_alpha = FALSE,
  ar1_B = FALSE,
  update_rho_B = FALSE,
  seed = 6886,
  verbose = TRUE,
  time_thin = 1,
  previous = NULL,
  init = NULL,
  symmetric = FALSE,
  ...
)
```

## Arguments

- Y:

  Data array (nodes x nodes x relations x time)

- family:

  Data family

- r:

  Rank for low-rank factorization

- nscan:

  Number of iterations of the Markov chain (beyond burn-in)

- burn:

  Burn-in for the Markov chain

- odens:

  Output density for the Markov chain

- ar1_alpha:

  Use AR(1) for alpha dynamics

- update_rho_alpha:

  Update AR coefficient for alpha

- ar1_B:

  Use AR(1) for B dynamics

- update_rho_B:

  Update AR coefficient for B

- seed:

  Random seed

- verbose:

  Logical or numeric. TRUE prints every 100 iterations, numeric prints
  every n iterations, FALSE suppresses output.

- time_thin:

  Save every nth time point

- previous:

  Previous results to continue from

- init:

  Initial values

- symmetric:

  Logical. Not supported for low-rank models (will error). Default:
  FALSE.

- ...:

  Additional arguments (currently unused)

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
sim <- simulate_lowrank_dbn(n = 8, time = 5, r = 2, seed = 1)
fit <- dbn_lowrank(sim$Y, r = 2, nscan = 200, burn = 100, verbose = FALSE)
# }
```
