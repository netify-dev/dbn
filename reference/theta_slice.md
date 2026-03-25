# Extract Theta slices from posterior draws

Extract specific slices of Theta arrays from posterior draws

## Usage

``` r
theta_slice(fit, draws = NULL, i = NULL, j = NULL, rel = NULL, time = NULL)
```

## Arguments

- fit:

  A dbn model fit object

- draws:

  Integer vector of draw indices to extract

- i:

  Row indices (sender nodes)

- j:

  Column indices (receiver nodes)

- rel:

  Relation indices

- time:

  Time indices

## Value

List of Theta slices

## See also

[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md),
[`theta_credible`](https://netify-dev.github.io/dbn/reference/theta_credible.md),
[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
ts <- theta_slice(fit, time = 1)
# }
```
