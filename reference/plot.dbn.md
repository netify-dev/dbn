# Plot Diagnostics for a Fitted DBN Model

Creates model-specific diagnostic plots for a fitted DBN object. The
type of plot depends on the model variant:

- **Static**: trace plots, posterior histograms, and a network summary
  of the estimated B matrix.

- **Dynamic**: trace plots for variance parameters and, if available,
  time-varying A/B summaries.

- **Lowrank**: trace plots, estimated factor trajectories
  (\\\alpha_t\\), and the posterior mean node-loading matrix U.

- **HMM**: regime probabilities over time, the estimated transition
  matrix, and MCMC trace plots.

- **Piecewise**: trace plots for each regime block.

## Usage

``` r
# S3 method for class 'dbn'
plot(x, ...)
```

## Arguments

- x:

  A fitted `dbn` object returned by
  [`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md).

- ...:

  Additional arguments passed to model-specific plot functions (e.g.,
  `alpha` for edge significance in the static model).

## Value

A ggplot2 object (or arranged multi-panel plot) is printed and returned
invisibly.

## See also

[`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`plot_trace()`](https://netify-dev.github.io/dbn/reference/plot_trace.md),
[`check_convergence()`](https://netify-dev.github.io/dbn/reference/check_convergence.md),
[`summary_dbn()`](https://netify-dev.github.io/dbn/reference/summary_dbn.md)

## Examples

``` r
# \donttest{
sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
plot(fit)

# }
```
