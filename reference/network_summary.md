# Network-level posterior summary

Compute network-level statistics at each time point across posterior
draws (e.g., mean edge weight, density of positive edges)

## Usage

``` r
network_summary(
  fit,
  stat = c("mean", "density", "strength"),
  rel = 1,
  time = NULL
)
```

## Arguments

- fit:

  A dbn model fit object

- stat:

  Statistic to compute: "mean" (average theta), "density" (fraction of
  positive values), or "strength" (mean absolute theta)

- rel:

  Relation index (default: 1)

- time:

  Time indices (default: all)

## Value

Data frame with columns: time, mean, lower, upper

## See also

[`edge_prob`](https://netify-dev.github.io/dbn/reference/edge_prob.md),
[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md),
[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
ns <- network_summary(fit)
# }
```
