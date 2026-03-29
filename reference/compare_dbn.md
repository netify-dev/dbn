# Compare Multiple DBN Models

Creates side-by-side trace plots of scalar variance parameters from two
or more fitted DBN models. Useful for comparing convergence behavior
across different model specifications (e.g., static vs. dynamic,
different ranks, different families).

## Usage

``` r
compare_dbn(...)
```

## Arguments

- ...:

  Two or more fitted `dbn` objects to compare. Objects are labeled
  "Model 1", "Model 2", etc. in the plot legend.

## Value

A ggplot2 object showing overlaid trace plots, faceted by parameter.

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
