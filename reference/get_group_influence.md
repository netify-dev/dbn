# Extract Group Influence Trajectories

Computes group influence trajectories with posterior quantiles

## Usage

``` r
get_group_influence(
  fit,
  group,
  type = c("sender", "target"),
  measure = c("rowsum", "rowmean", "l2"),
  fun = c("mean", "sum"),
  probs = c(0.025, 0.5, 0.975)
)
```

## Arguments

- fit:

  A "dbn" object from dbn_dynamic()

- group:

  Integer vector of actor indices

- type:

  "sender" or "target"

- measure:

  Per-actor metric: "rowsum", "rowmean", "l2"

- fun:

  Aggregation across actors: "mean" or "sum"

- probs:

  Quantile probabilities to compute

## Value

Data frame with time, posterior quantiles, and mean

## See also

[`plot_group_influence`](https://netify-dev.github.io/dbn/reference/plot_group_influence.md),
[`compare_group_influence`](https://netify-dev.github.io/dbn/reference/compare_group_influence.md),
[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 10, time = 10, seed = 6886)
fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
    nscan = 200, burn = 100, verbose = FALSE)
inf_data <- get_group_influence(fit, group = c(1, 3, 5), type = "sender")
# }
```
