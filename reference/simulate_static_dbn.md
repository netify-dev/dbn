# Simulate from Static DBN Model

Generates data from a static DBN with fixed A and B matrices

## Usage

``` r
simulate_static_dbn(
  n = 30,
  n_col = n,
  p = 2,
  time = 50,
  sigma2 = 0.5,
  tau2 = 0.1,
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

- tau2:

  Variance for A/B deviations from identity

- K:

  Number of ordinal categories

- return_truth:

  Return true parameters in a truth sub-list

- seed:

  Random seed

- symmetric:

  Logical. If TRUE, set B = A for symmetric/undirected networks.
  Default: FALSE.

## Value

List containing simulated data and true parameters

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md) for model
fitting,
[`simulate_test_data`](https://netify-dev.github.io/dbn/reference/simulate_test_data.md)
for quick test data

## Examples

``` r
sim <- simulate_static_dbn(n = 8, time = 5, seed = 42)
str(sim$Y)
#>  int [1:8, 1:8, 1:2, 1:5] NA 4 3 2 4 5 5 1 3 NA ...
```
