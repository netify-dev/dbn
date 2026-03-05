# Network Snapshot

Heat map of Theta at given time (posterior mean)

## Usage

``` r
net_snapshot(
  fit,
  t,
  rel = 1,
  sparse = NULL,
  eps = 1e-04,
  show_significant = FALSE,
  cred_level = 0.025
)
```

## Arguments

- fit:

  Dynamic dbn object

- t:

  Time point

- rel:

  Relation index (default: 1)

- sparse:

  Auto-switch to sparse visualization for large networks

- eps:

  Threshold for sparse plotting

- show_significant:

  Logical, whether to show only significant effects (default: FALSE)

- cred_level:

  Credible level for significance (default corresponds to 95% CI)

## Value

ggplot2 object or base R plot

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`dyad_path`](https://netify-dev.github.io/dbn/reference/dyad_path.md),
[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md),
[`role_trajectory`](https://netify-dev.github.io/dbn/reference/role_trajectory.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
net_snapshot(fit, t = 5)

# }
```
