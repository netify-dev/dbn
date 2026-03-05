# Debug IRF

Diagnostic output for troubleshooting IRF computation

## Usage

``` r
debug_irf(fit, draw_idx = 1, shock_i = 8, shock_j = 12)
```

## Arguments

- fit:

  A dbn model fit object

- draw_idx:

  Posterior draw index to debug

- shock_i:

  Source node

- shock_j:

  Target node

## See also

[`compute_irf`](https://netify-dev.github.io/dbn/reference/compute_irf.md),
[`build_shock`](https://netify-dev.github.io/dbn/reference/build_shock.md),
[`plot.dbn_irf`](https://netify-dev.github.io/dbn/reference/plot.dbn_irf.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 15, time = 50, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
debug_irf(fit, draw_idx = 1, shock_i = 1, shock_j = 2)
#> === IRF Debug Information ===
#> 
#> 1. Model Structure:
#>    Model type: dynamic 
#>    n_row = 15 
#>    n_col = 15 
#>    p = 2 
#>    T = 50 
#> 
#> 2. A and B matrices:
#>    A is list: TRUE 
#>    B is list: TRUE 
#>    Length of A: 200 
#>    Dim of A[[1]]: 15 15 25 
#>    Length of B: 200 
#>    Dim of B[[1]]: 15 15 25 
#> 
#> 3. M matrix:
#>    M dimensions: 15 15 2 200 
#> 
#> 4. Attempting to extract matrices for draw 1 :
#>    A_array dim: 15 15 25 
#>    B_array dim: 15 15 25 
#>    M for draw 1 extracted, dim: 15 15 
#> 
#> 5. Testing C++ impulse_response_dynamic:
#>    ERROR: t0 must be between 0 and T-1 
#> 
#> === End Debug ===
# }
```
