# Compute Network-Level Impulse Response Functions

Computes IRFs for network-level statistics given a shock.

**Important:** `ggplot2` exports its own `stat_density` function, which
masks
[`dbn::stat_density`](https://netify-dev.github.io/dbn/reference/stat_density.md)
when both packages are loaded. If you use `ggplot2`, pass the statistic
explicitly: `stat_fun = dbn::stat_density`.

## Usage

``` r
compute_irf(
  fit,
  shock,
  H = 20,
  t0 = 1,
  stat_fun = stat_density,
  n_draws = NULL,
  shock_pars = list(),
  ...
)
```

## Arguments

- fit:

  A dbn model fit object

- shock:

  Shock matrix or shock type (see
  [`build_shock`](https://netify-dev.github.io/dbn/reference/build_shock.md))

- H:

  Number of horizons to compute

- t0:

  Shock time for dynamic models

- stat_fun:

  Network statistic function (default:
  [`stat_density`](https://netify-dev.github.io/dbn/reference/stat_density.md)).
  Must accept a matrix and return a scalar. When `ggplot2` is loaded,
  use
  [`dbn::stat_density`](https://netify-dev.github.io/dbn/reference/stat_density.md)
  explicitly to avoid name masking.

- n_draws:

  Number of posterior draws to use

- shock_pars:

  List of parameters for
  [`build_shock`](https://netify-dev.github.io/dbn/reference/build_shock.md)
  if `shock` is a character string

- ...:

  Additional arguments

## Value

Data frame with IRF results including posterior summaries

## See also

[`build_shock`](https://netify-dev.github.io/dbn/reference/build_shock.md),
[`plot.dbn_irf`](https://netify-dev.github.io/dbn/reference/plot.dbn_irf.md),
[`print.dbn_irf`](https://netify-dev.github.io/dbn/reference/print.dbn_irf.md),
[`debug_irf`](https://netify-dev.github.io/dbn/reference/debug_irf.md),
[`stat_density`](https://netify-dev.github.io/dbn/reference/stat_density.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
S <- build_shock(m = 6, type = "unit_edge", i = 1, j = 2)
irf <- compute_irf(fit, shock = S, H = 5, stat_fun = dbn::stat_density)
print(irf)
#> Network Impulse Response Function
#> Model: dynamic 
#> Statistic: custom_function 
#> Shock time: 1 
#> 
#> Summary:
#>   horizon   mean median    sd   q025  q975    q10   q90
#> 1       0  0.033  0.033 0.000  0.033 0.033  0.033 0.033
#> 2       1  0.002  0.001 0.036 -0.065 0.096 -0.042 0.043
#> 3       2 -0.003  0.000 0.056 -0.096 0.096 -0.055 0.041
#> 4       3 -0.005 -0.003 0.057 -0.106 0.088 -0.057 0.051
#> 5       4 -0.014 -0.001 0.061 -0.192 0.075 -0.086 0.038
#> 6       5  0.011  0.000 0.097 -0.128 0.237 -0.050 0.064
# }
```
