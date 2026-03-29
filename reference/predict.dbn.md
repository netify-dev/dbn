# Predict from a Fitted DBN Model

Generates predictions from a fitted DBN model. With no arguments (or
PPD-only arguments), returns posterior predictive replications of the
observed data for model checking. With forecasting arguments (`H`, `S`,
`summary`), propagates the estimated dynamics forward to produce
multi-step-ahead forecasts.

## Usage

``` r
# S3 method for class 'dbn'
predict(object, ...)
```

## Arguments

- object:

  A fitted `dbn` object returned by
  [`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md).

- ...:

  Model-specific arguments:

  H

  :   Forecast horizon (number of steps ahead). Triggers forecasting
      mode.

  S, draws

  :   Number of posterior draws to use for forecasting.

  summary

  :   If `"mean"`, return the pointwise mean across draws. If `NULL`
      (default), return the full array of simulated trajectories.

  ndraws

  :   Number of draws for posterior predictive checks (default: 100).

  seed

  :   Random seed for reproducibility.

## Value

For forecasting (`H` specified): an array of dimensions
`[n_row, n_col, p, H]` (if `summary = "mean"`) or
`[n_row, n_col, p, H, draws]`. For posterior predictive checks: a list
of class `"dbn_ppd"` containing replicated datasets.

## Details

**Posterior predictive distribution (default):** Simulates new datasets
from the fitted model to compare against the observed data. Use
[`plot_ppc_ecdf()`](https://netify-dev.github.io/dbn/reference/plot_ppc_ecdf.md)
or
[`plot_ppc_density()`](https://netify-dev.github.io/dbn/reference/plot_ppc_density.md)
for visual checks.

**Forecasting:** Pass `H` (horizon), `S` or `draws` (number of
simulations), and optionally `summary = "mean"` to average across draws.
The forecast propagates \\A_t\\, \\B_t\\, and process noise forward from
the final observed time point.

## See also

[`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`posterior_predict_dbn()`](https://netify-dev.github.io/dbn/reference/posterior_predict_dbn.md),
[`plot_ppc_ecdf()`](https://netify-dev.github.io/dbn/reference/plot_ppc_ecdf.md),
[`plot_ppc_density()`](https://netify-dev.github.io/dbn/reference/plot_ppc_density.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 6886)
fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
           nscan = 200, burn = 100, verbose = FALSE)

# forecasting: 5 steps ahead
fc <- predict(fit, H = 5, S = 50, summary = "mean")
dim(fc)
#> [1] 6 6 2 5

# posterior predictive check
ppd <- predict(fit, ndraws = 10)
# }
```
