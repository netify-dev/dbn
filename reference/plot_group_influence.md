# Plot Group Influence Profile

Plots posterior group influence over time for dynamic models

## Usage

``` r
plot_group_influence(
  fit,
  group,
  type = c("sender", "target"),
  fun = c("mean", "sum"),
  measure = c("rowsum", "rowmean", "l2"),
  cred = 0.95
)
```

## Arguments

- fit:

  A "dbn" object from dbn_dynamic()

- group:

  Integer vector of actor indices

- type:

  "sender" (rows of A_t) or "target" (columns of B_t)

- fun:

  Aggregation across actors: "mean" or "sum"

- measure:

  Per-actor metric: "rowsum" (default), "rowmean", "l2"

- cred:

  Credible band level (0.95 gives 95% bands)

## Value

A ggplot2 object

## See also

[`get_group_influence`](https://netify-dev.github.io/dbn/reference/get_group_influence.md),
[`compare_group_influence`](https://netify-dev.github.io/dbn/reference/compare_group_influence.md),
[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot sender influence for actors 1, 3, 5
plot_group_influence(fit, group = c(1, 3, 5), type = "sender")

# Plot target influence using L2 norm
plot_group_influence(fit,
    group = c(1, 3, 5), type = "target",
    fun = "sum", measure = "l2", cred = 0.8
)
} # }
```
