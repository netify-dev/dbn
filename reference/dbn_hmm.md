# DBN HMM Model

Regime-switching DBN for large-scale networks (100+ actors, 50+ time
steps)

## Usage

``` r
dbn_hmm(
  Y,
  family = c("ordinal", "gaussian", "binary"),
  R = 3,
  nscan = 10000,
  burn = 1000,
  odens = 1,
  delta = rep(1, R),
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

  Data family (ordinal, gaussian, or binary)

- R:

  Number of regimes

- nscan:

  Number of iterations of the Markov chain (beyond burn-in)

- burn:

  Burn-in for the Markov chain

- odens:

  Output density for the Markov chain

- delta:

  Dirichlet prior for transition matrix

- seed:

  Random seed

- verbose:

  Print progress

- time_thin:

  Save every nth time point, when 1, save all time points

- previous:

  Previous dbn_hmm results to continue from (optional)

- init:

  List of initial values: S, A_list, B_list, Pi, sigma2_proc, tau_A2,
  tau_B2, g2, pi0 (optional)

- symmetric:

  Logical. If TRUE, enforce B = A for each regime. Default: FALSE.

- ...:

  Additional arguments

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
sim <- simulate_hmm_dbn(n = 8, time = 10, R = 2, seed = 1)
fit <- dbn_hmm(sim$Y, R = 2, nscan = 200, burn = 100, verbose = FALSE)
# }
```
