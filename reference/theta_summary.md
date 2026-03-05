# Summarize Theta over posterior draws

Compute summary statistics for Theta parameters

## Usage

``` r
theta_summary(
  fit,
  fun = mean,
  draws = NULL,
  i = NULL,
  j = NULL,
  rel = NULL,
  time = NULL,
  chunk = 20
)
```

## Arguments

- fit:

  A dbn model fit object

- fun:

  Summary function (default: mean)

- draws:

  Draw indices (default: all)

- i:

  Row indices

- j:

  Column indices

- rel:

  Relation indices

- time:

  Time indices

- chunk:

  Chunk size for memory-efficient processing

## Value

Data frame with summarized values

## See also

[`theta_slice`](https://netify-dev.github.io/dbn/reference/theta_slice.md),
[`theta_credible`](https://netify-dev.github.io/dbn/reference/theta_credible.md),
[`plot_theta`](https://netify-dev.github.io/dbn/reference/plot_theta.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
ts <- theta_summary(fit)
# }
```
