# Compare Multiple DBN Models

Creates comparative plots for multiple DBN results using ggplot2

## Usage

``` r
compare_dbn(...)
```

## Arguments

- ...:

  Multiple dbn objects to compare

## Value

A ggplot2 object or list of plots

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`check_convergence`](https://netify-dev.github.io/dbn/reference/check_convergence.md),
[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_static_dbn(n = 6, time = 10, seed = 1)
fit1 <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
fit2 <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
compare_dbn(fit1, fit2)

# }
```
