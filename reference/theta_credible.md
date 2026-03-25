# Posterior credible intervals for Theta

Compute pointwise credible intervals for each dyad over time

## Usage

``` r
theta_credible(
  fit,
  probs = c(0.05, 0.5, 0.95),
  i = NULL,
  j = NULL,
  rel = 1,
  time = NULL
)
```

## Arguments

- fit:

  A dbn model fit object

- probs:

  Probability levels (default: 90 percent interval + median)

- i:

  Row indices (sender nodes, default: all)

- j:

  Column indices (receiver nodes, default: all)

- rel:

  Relation index (default: 1)

- time:

  Time indices (default: all)

## Value

Data frame with columns: i, j, rel, time, mean, lower, median, upper

## See also

[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md),
[`theta_slice`](https://netify-dev.github.io/dbn/reference/theta_slice.md),
[`edge_prob`](https://netify-dev.github.io/dbn/reference/edge_prob.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
tc <- theta_credible(fit)
# }
```
