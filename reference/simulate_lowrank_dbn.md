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

  Rank

- sigma2:

  Innovation variance

- tau_alpha2:

  Variance for alpha innovations

- tauB2:

  Variance for B innovations

- ar1_alpha:

  Use AR(1) for alpha dynamics

- rho_alpha:

  AR coefficient for alpha

- seed:

  Random seed

- return_truth:

  Return true latent factors and parameters

## Value

List containing simulated data and true parameters

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md) for model
fitting,
[`simulate_test_data`](https://netify-dev.github.io/dbn/reference/simulate_test_data.md)
for quick test data

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_lowrank_dbn(n = 25, p = 1, time = 12, r = 1,
    sigma2 = 0.2, tau_alpha2 = 0.02, seed = 42)
fit <- dbn_lowrank(sim$Y, r = 1, n_iter = 600, burn = 200, thin = 2)
} # }
```
