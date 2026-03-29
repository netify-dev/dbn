# Simulate from HMM DBN Model

Generates data from a regime-switching HMM DBN

## Usage

``` r
simulate_hmm_dbn(
  n = 30,
  n_col = n,
  p = 2,
  time = 50,
  R = 3,
  sigma2 = 0.5,
  tau_A2 = 0.2,
  tau_B2 = 0.2,
  transition_prob = 0.8,
  stickiness = NULL,
  seed = NULL,
  return_truth = TRUE,
  symmetric = FALSE
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

- R:

  Number of regimes

- sigma2:

  Innovation variance

- tau_A2:

  Prior variance for regime A matrices

- tau_B2:

  Prior variance for regime B matrices

- transition_prob:

  Diagonal transition probability

- stickiness:

  Deprecated alias for transition_prob

- seed:

  Random seed

- return_truth:

  Return true latent states and parameters

- symmetric:

  Logical. If TRUE, set B = A for each regime. Default: FALSE.

## Value

List containing simulated data and true parameters

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md) for model
fitting,
[`simulate_test_data`](https://netify-dev.github.io/dbn/reference/simulate_test_data.md)
for quick test data

## Examples

``` r
sim <- simulate_hmm_dbn(n = 8, p = 1, time = 10, R = 2, seed = 6886)
dim(sim$Y)
#> [1]  8  8  1 10
table(sim$S)
#> 
#> 1 2 
#> 4 6 
```
