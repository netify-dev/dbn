# Simulate from Dynamic DBN Model

Generates data from a dynamic DBN with time-varying A and B

## Usage

``` r
simulate_dynamic_dbn(
  n = 30,
  n_col = n,
  p = 2,
  time = 50,
  sigma2 = 0.5,
  tauA2 = 0.05,
  tauB2 = 0.05,
  ar1 = FALSE,
  rhoA = 0.9,
  rhoB = 0.9,
  K = 5,
  return_truth = TRUE,
  seed = NULL,
  symmetric = FALSE
)
```

## Arguments

- n:

  Number of row actors / senders

- n_col:

  Number of column actors / receivers (default: n)

- p:

  Number of relation types

- time:

  Number of time points

- sigma2:

  Innovation variance

- tauA2:

  Variance for A innovations

- tauB2:

  Variance for B innovations

- ar1:

  Use AR(1) dynamics instead of random walk

- rhoA:

  AR coefficient for A

- rhoB:

  AR coefficient for B

- K:

  Number of ordinal categories

- return_truth:

  Return true parameters in a truth sub-list

- seed:

  Random seed

- symmetric:

  Logical. If TRUE, set B = A at each time point. Default: FALSE.

## Value

List containing simulated data and true parameters

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md) for model
fitting,
[`simulate_test_data`](https://netify-dev.github.io/dbn/reference/simulate_test_data.md)
for quick test data

## Examples

``` r
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 42)
str(sim$Y)
#>  int [1:6, 1:6, 1:2, 1:5] NA 2 2 3 2 5 5 NA 2 5 ...
```
