# Plot Regime Probabilities

Stacked area plot of posterior regime probabilities over time (HMM only)

## Usage

``` r
plot_regime_probs(fit)
```

## Arguments

- fit:

  A dbn_hmm model fit object

## Value

A ggplot object or NULL

## See also

[`regime_probs`](https://netify-dev.github.io/dbn/reference/regime_probs.md),
[`plot_trace`](https://netify-dev.github.io/dbn/reference/plot_trace.md),
[`dbn_hmm`](https://netify-dev.github.io/dbn/reference/dbn_hmm.md)

## Examples

``` r
# \donttest{
sim <- simulate_hmm_dbn(n = 6, time = 10, R = 2, seed = 1)
fit <- dbn(sim$Y, model = "hmm", R = 2, nscan = 200, burn = 100, verbose = FALSE)
if (requireNamespace("ggplot2", quietly = TRUE)) plot_regime_probs(fit)

# }
```
