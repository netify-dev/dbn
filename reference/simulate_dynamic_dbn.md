# Simulate from Dynamic DBN Model

Generates synthetic network data from a dynamic DBN with time-varying A
and B influence matrices. Each period's influence structure evolves from
the previous period's via a random walk (default) or AR(1) process.

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

  Number of actors (senders)

- n_col:

  Number of receivers (default: same as `n`)

- p:

  Number of relation types (default: 2)

- time:

  Number of time periods to simulate

- sigma2:

  Process noise variance

- tauA2:

  Innovation variance for A (how fast sender influence changes). Larger
  values produce more volatile influence dynamics.

- tauB2:

  Innovation variance for B (how fast receiver influence changes)

- ar1:

  If TRUE, use AR(1) dynamics (smooth evolution). If FALSE (default),
  use random walk (less smooth).

- rhoA:

  AR(1) persistence parameter for A (only used if `ar1 = TRUE`). Values
  near 1 produce slowly changing influence.

- rhoB:

  AR(1) persistence parameter for B

- K:

  Number of ordinal categories for observed data (default: 5)

- return_truth:

  If TRUE (default), include true parameters in output

- seed:

  Random seed for reproducibility

- symmetric:

  If TRUE, set B = A at each time point.

## Value

A list containing:

- Y:

  Observed ordinal data array `[n_row, n_col, p, time]`

- Z:

  Continuous latent values (use with `family = "gaussian"`)

- Theta:

  True latent network state at each time point

- A:

  True time-varying sender influence `[n_row, n_row, time]`

- B:

  True time-varying receiver influence `[n_col, n_col, time]`

- M:

  True baseline mean array

- sigma2, tauA2, tauB2, rhoA, rhoB:

  True parameter values

## See also

[`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md) for model
fitting,
[`simulate_static_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_static_dbn.md)
for fixed-influence version

## Examples

``` r
sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 42)
dim(sim$Y)    # observed ordinal data
#> [1]  6  6  2 10
dim(sim$A)    # true A matrices [n, n, time]
#> [1]  6  6 10
```
