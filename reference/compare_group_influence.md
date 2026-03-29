# Compare Group Influences

Compares influence trajectories of multiple groups

## Usage

``` r
compare_group_influence(
  fit,
  groups,
  group_names = NULL,
  type = c("sender", "target"),
  measure = c("rowsum", "rowmean", "l2"),
  fun = c("mean", "sum"),
  cred = 0.95
)
```

## Arguments

- fit:

  A "dbn" object from dbn_dynamic()

- groups:

  List of integer vectors, each defining a group

- group_names:

  Optional character vector of group names

- type:

  "sender" or "target"

- measure:

  Per-actor metric: "rowsum", "rowmean", "l2"

- fun:

  Aggregation: "mean" or "sum"

- cred:

  Credible band level

## Value

A ggplot2 object

## See also

[`plot_group_influence`](https://netify-dev.github.io/dbn/reference/plot_group_influence.md),
[`get_group_influence`](https://netify-dev.github.io/dbn/reference/get_group_influence.md),
[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 10, time = 10, seed = 6886)
fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
    nscan = 200, burn = 100, verbose = FALSE)
compare_group_influence(fit,
    groups = list(c(1, 3, 5), c(2, 4, 6)),
    group_names = c("Group A", "Group B"))
#> Error in if (type == "sender") "Sender" else "Target": the condition has length > 1
# }
```
