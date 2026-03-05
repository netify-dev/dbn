# Summarize scalar parameters

Compute quantiles for scalar parameter traces

## Usage

``` r
param_summary(fit, probs = c(0.05, 0.5, 0.95))
```

## Arguments

- fit:

  A dbn model fit object

- probs:

  Probability levels for quantiles

## Value

Data frame with parameter summaries

## See also

[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md),
[`plot_trace`](https://netify-dev.github.io/dbn/reference/plot_trace.md),
[`derive_draws`](https://netify-dev.github.io/dbn/reference/derive_draws.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
ps <- param_summary(fit)
# }
```
