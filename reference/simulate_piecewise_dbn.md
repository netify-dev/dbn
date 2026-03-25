# Simulate Data from Piecewise DBN

generates synthetic data from piecewise-static model

## Usage

``` r
simulate_piecewise_dbn(
  n = 20,
  time = 50,
  blocks = 4,
  p = 1,
  sigma2 = 1,
  tau2 = 0.5,
  seed = NULL
)
```

## Arguments

- n:

  number of actors

- time:

  number of time points

- blocks:

  block specification

- p:

  number of relations

- sigma2:

  observation variance

- tau2:

  A/B prior variance

- seed:

  random seed

## Value

list with Y, true_A, true_B, true_M, block_info

## Examples

``` r
# \donttest{
sim <- simulate_piecewise_dbn(n = 10, time = 40, blocks = c(15, 30, 40))
dim(sim$Y)
#> [1] 10 10  1 40
# }
```
