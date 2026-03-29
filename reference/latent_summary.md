# Summarize Baseline Mean M

Computes posterior summaries of the baseline mean matrix M, which
captures persistent dyad-specific tendencies (e.g., stable alliances or
rivalries) that do not change over time. Returns a data frame with one
row per dyad, suitable for plotting or comparison against known ground
truth in simulation studies.

## Usage

``` r
latent_summary(fit, fun = mean, draws = NULL, rel = NULL, chunk = 20)
```

## Arguments

- fit:

  A fitted `dbn` object returned by
  [`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md).

- fun:

  Summary function applied across posterior draws (default:
  [mean](https://rdrr.io/r/base/mean.html)). Use
  [median](https://rdrr.io/r/stats/median.html) for a robust
  alternative, or a custom function.

- draws:

  Integer vector of posterior draw indices to use. If `NULL` (default),
  all available draws are used.

- rel:

  Integer vector of relation indices to summarize. If `NULL` (default),
  all relations are included.

- chunk:

  Integer controlling memory-efficient processing. Draws are processed
  in blocks of this size. Only relevant for very large fits.

## Value

Data frame with columns `i` (sender), `j` (receiver), `rel`, and `value`
(the summary statistic).

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
