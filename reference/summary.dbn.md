# Summarize a Fitted DBN Model

Prints a structured summary of a fitted DBN model including data
dimensions, posterior means and 95\\ parameters, and model-specific
details (e.g., transition matrix for HMM, block structure for piecewise,
DIC for Gaussian models).

## Usage

``` r
# S3 method for class 'dbn'
summary(object, ...)
```

## Arguments

- object:

  A fitted `dbn` object returned by
  [`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `object`. The summary is printed to the console.

## See also

[`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`param_summary()`](https://netify-dev.github.io/dbn/reference/param_summary.md),
[`check_convergence()`](https://netify-dev.github.io/dbn/reference/check_convergence.md),
[`plot_dbn()`](https://netify-dev.github.io/dbn/reference/plot_dbn.md)

## Examples

``` r
# \donttest{
sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
summary(fit)
#> 
#> ── Static Bilinear Network Model ───────────────────────────────────────────────
#> 
#> ── Dimensions 
#>   Nodes: 8
#>   Relations: 2
#>   Time points: 10
#> 
#> ── Parameter estimates (mean [95% CI]) 
#> s2: 1 [1, 1]
#> t2: 0.183 [0.098, 0.315]
#> g2: 53.044 [19.036, 100.586]
#> 
#> ── Settings 
#> family: ordinal
#> draws: 200
#> settings: 200, 100, and 1
# }
```
