# Check MCMC Convergence

Provides convergence diagnostics for MCMC chains

## Usage

``` r
check_convergence(results)
```

## Arguments

- results:

  Output from dbn()

## Value

Invisible NULL (diagnostics are printed and plotted)

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`compare_dbn`](https://netify-dev.github.io/dbn/reference/compare_dbn.md),
[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_static_dbn(n = 6, time = 10, seed = 1)
fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
check_convergence(fit)
#> 
#> ── Effective Sample Sizes 
#>         s2         t2         g2 
#>   0.000000 114.447593   9.374666 
#> 
#> ── Geweke Diagnostic 
#> 
#> Fraction in 1st window = 0.1
#> Fraction in 2nd window = 0.5 
#> 
#>     s2     t2     g2 
#>    NaN  1.263 -2.793 
#> 
# }
```
