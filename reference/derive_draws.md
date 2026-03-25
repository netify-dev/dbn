# Derive new quantities from posterior draws

Apply a transformation function to each posterior draw

## Usage

``` r
derive_draws(fit, fun, draws = NULL, chunk = 20, name = "derived")
```

## Arguments

- fit:

  A dbn model fit object

- fun:

  Function to apply to each Theta draw

- draws:

  Draw indices to process

- chunk:

  Chunk size for memory efficiency

- name:

  Name for the derived quantity

## Value

List of derived quantities with class "dbn_derived"

## See also

[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md),
[`theta_slice`](https://netify-dev.github.io/dbn/reference/theta_slice.md),
[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
dd <- derive_draws(fit, function(x) x$sigma2)
# }
```
