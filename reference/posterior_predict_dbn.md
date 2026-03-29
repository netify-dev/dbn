# Generate Posterior Predictive Samples

Simulates new datasets from the fitted model. For each posterior draw of
the model parameters, generates a complete replicated dataset. These
replications can be compared to the observed data using
[`plot_ppc_ecdf()`](https://netify-dev.github.io/dbn/reference/plot_ppc_ecdf.md)
or
[`plot_ppc_density()`](https://netify-dev.github.io/dbn/reference/plot_ppc_density.md)
to check whether the model captures the key features of the data (a
"posterior predictive check").

## Usage

``` r
posterior_predict_dbn(fit, ndraws = 100, seed = NULL, draws = NULL)
```

## Arguments

- fit:

  A dbn model fit object (output from
  [`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md))

- ndraws:

  Number of replicated datasets to generate (default: 100)

- seed:

  Random seed for reproducibility

- draws:

  Specific posterior draw indices to use (overrides `ndraws`)

## Value

A list of `ndraws` replicated data arrays, each with the same dimensions
as the original data. Has class `"dbn_ppd"`.

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
