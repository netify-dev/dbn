# Plot Dyad Trajectory

Plots posterior mean and 95% bands for a single dyad through time

## Usage

``` r
dyad_path(fit, i, j, rel = NULL, facet = TRUE, cred = c(0.025, 0.975))
```

## Arguments

- fit:

  Dynamic dbn object

- i:

  Actor i index

- j:

  Actor j index

- rel:

  Relation indices (default: NULL = all relations)

- facet:

  Whether to facet by relation (default: TRUE)

- cred:

  Credible interval quantiles (default: c(0.025, 0.975))

## Value

A ggplot2 object

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`net_snapshot`](https://netify-dev.github.io/dbn/reference/net_snapshot.md),
[`role_trajectory`](https://netify-dev.github.io/dbn/reference/role_trajectory.md),
[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
dyad_path(fit, i = 1, j = 2)

# }
```
