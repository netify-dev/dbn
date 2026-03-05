# Posterior edge probability

For each dyad, compute the posterior probability that the latent theta
is positive (or exceeds a threshold)

## Usage

``` r
edge_prob(fit, threshold = 0, rel = 1, time = NULL)
```

## Arguments

- fit:

  A dbn model fit object

- threshold:

  Threshold value (default: 0)

- rel:

  Relation index (default: 1)

- time:

  Time index (default: last time point)

## Value

Matrix of posterior probabilities (n_row x n_col)

## See also

[`network_summary`](https://netify-dev.github.io/dbn/reference/network_summary.md),
[`theta_credible`](https://netify-dev.github.io/dbn/reference/theta_credible.md),
[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
ep <- edge_prob(fit)
#> Warning: Model fit does not contain theta draws
# }
```
