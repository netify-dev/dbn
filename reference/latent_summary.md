# Summarize latent means (M arrays)

Compute summaries for latent mean arrays M

## Usage

``` r
latent_summary(fit, fun = mean, draws = NULL, rel = NULL, chunk = 20)
```

## Arguments

- fit:

  A dbn model fit object

- fun:

  Summary function

- draws:

  Draw indices

- rel:

  Relation indices (optional)

- chunk:

  Chunk size for processing

## Value

Data frame with M summaries

## See also

[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md),
[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md),
[`derive_draws`](https://netify-dev.github.io/dbn/reference/derive_draws.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
ls <- latent_summary(fit)
# }
```
