# Generate posterior predictive samples

Generate new observations from the posterior predictive distribution

## Usage

``` r
posterior_predict_dbn(fit, ndraws = 100, seed = NULL, draws = NULL)
```

## Arguments

- fit:

  A dbn model fit object

- ndraws:

  Number of posterior draws to use (default: 100)

- seed:

  Random seed for reproducibility

- draws:

  Specific draw indices to use (overrides ndraws)

## Value

List of predicted observations with class "dbn_ppd"

## See also

[`plot_ppc_ecdf`](https://netify-dev.github.io/dbn/reference/plot_ppc_ecdf.md),
[`plot_ppc_density`](https://netify-dev.github.io/dbn/reference/plot_ppc_density.md),
[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
ppd <- posterior_predict_dbn(fit, ndraws = 5)
# }
```
