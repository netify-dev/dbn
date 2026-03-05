# Plot Parameter Trace Plots

Trace plots with running mean and posterior mean overlay

## Usage

``` r
plot_trace(fit, pars = NULL, ncol = 2)
```

## Arguments

- fit:

  A dbn model fit object

- pars:

  Character vector of parameter names to plot

- ncol:

  Number of columns for multi-panel plot

## Value

A ggplot object or NULL

## See also

[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md),
[`plot_theta`](https://netify-dev.github.io/dbn/reference/plot_theta.md),
[`plot_regime_probs`](https://netify-dev.github.io/dbn/reference/plot_regime_probs.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
if (requireNamespace("ggplot2", quietly = TRUE)) plot_trace(fit)

# }
```
