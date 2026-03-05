# Posterior predictive density plot

Plot density of observed vs replicated data

## Usage

``` r
plot_ppc_density(fit, ppd = NULL, rel = 1, time = NULL, Y_obs = NULL)
```

## Arguments

- fit:

  A dbn model fit object

- ppd:

  Posterior predictive samples

- rel:

  Relation index

- time:

  Time indices

- Y_obs:

  Observed data array (required)

## Value

A ggplot object or NULL

## See also

[`posterior_predict_dbn`](https://netify-dev.github.io/dbn/reference/posterior_predict_dbn.md),
[`plot_ppc_ecdf`](https://netify-dev.github.io/dbn/reference/plot_ppc_ecdf.md),
[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_static_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
ppd <- posterior_predict_dbn(fit, ndraws = 5)
if (requireNamespace("ggplot2", quietly = TRUE)) plot_ppc_density(fit, ppd, Y_obs = sim$Y)

# }
```
