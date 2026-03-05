# Plot IRF

IRF plot with credible intervals

## Usage

``` r
# S3 method for class 'dbn_irf'
plot(x, ci_level = 0.95, title = NULL, ...)
```

## Arguments

- x:

  A dbn_irf object from compute_irf

- ci_level:

  Credible interval level (0.90 or 0.95)

- title:

  Plot title

- ...:

  Ignored

## Value

A ggplot2 object

## See also

[`compute_irf`](https://netify-dev.github.io/dbn/reference/compute_irf.md),
[`build_shock`](https://netify-dev.github.io/dbn/reference/build_shock.md),
[`print.dbn_irf`](https://netify-dev.github.io/dbn/reference/print.dbn_irf.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
S <- build_shock(m = 6, type = "unit_edge", i = 1, j = 2)
irf <- compute_irf(fit, shock = S, H = 5)
plot(irf)

# }
```
