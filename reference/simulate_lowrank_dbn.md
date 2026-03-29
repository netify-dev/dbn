# Simulate from Low-Rank DBN Model

Generates data from a low-rank DBN with factored A = U diag(alpha) U'

## Usage

``` r
simulate_lowrank_dbn(
  n = 30,
  n_col = n,
  p = 2,
  time = 50,
  r = 3,
  sigma2 = 0.5,
  tau_alpha2 = 0.1,
  tauB2 = 0.05,
  ar1_alpha = TRUE,
  rho_alpha = 0.9,
  seed = NULL,
  return_truth = TRUE
)
```

## Arguments

- n:

  Number of sender actors

- n_col:

  Number of receiver actors (default: n)

- p:

  Number of relation types

- time:

  Number of time points

- r:

  Rank of the factorization

- sigma2:

  Process noise variance

- tau_alpha2:

  Variance for alpha factor innovations

- tauB2:

  Variance for B innovations

- ar1_alpha:

  Use AR(1) for alpha dynamics (default TRUE)

- rho_alpha:

  AR(1) persistence for alpha (default 0.9)

- seed:

  Random seed for reproducibility

- return_truth:

  If TRUE (default), include true parameters in output

## Value

A list containing:

- Y:

  Observed ordinal data array `[n, n, p, time]`

- Z:

  Continuous latent values (use with `family = "gaussian"`)

- Theta:

  True latent network state at each time point

- U:

  True orthogonal factor matrix `[n, r]`

- alpha:

  True factor trajectories `[r, time]`

- A:

  True time-varying sender influence `[n, n, time]` (reconstructed from
  U and alpha)

- B:

  True time-varying receiver influence `[n_col, n_col, time]`

- M:

  True baseline mean array `[n, n_col, p]`

- sigma2, tau_alpha2, tauB2, r:

  True parameter values used in simulation

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md) for model
fitting,
[`simulate_test_data`](https://netify-dev.github.io/dbn/reference/simulate_test_data.md)
for quick test data

## Examples

``` r
sim <- simulate_lowrank_dbn(n = 8, p = 1, time = 5, r = 2, seed = 6886)
dim(sim$Y)
#> [1] 8 8 1 5
dim(sim$alpha)
#> [1] 2 5
```
