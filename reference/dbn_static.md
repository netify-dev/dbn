# Static DBN MCMC

Fits DBN model with fixed sender/receiver effects

## Usage

``` r
dbn_static(
  Y,
  family = c("ordinal", "gaussian", "binary"),
  nscan = 10000,
  burn = 1000,
  odens = 1,
  seed = 6886,
  verbose = TRUE,
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

- seed:

  Random seed for reproducibility

- verbose:

  Logical or numeric. TRUE prints every 100 iterations, numeric prints
  every n iterations, FALSE suppresses output.

- previous:

  Previous dbn_static results to continue from (optional)

- init:

  List of initial values for parameters: B, s2, t2, g2, M, Z (optional)

- symmetric:

  Logical. If TRUE, store symmetric flag in output dims. Default: FALSE.

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
sim <- simulate_static_dbn(n = 8, time = 5, seed = 1)
fit <- dbn_static(sim$Y, nscan = 200, burn = 100, verbose = FALSE)
# }
```
