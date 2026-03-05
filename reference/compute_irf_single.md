# Compute IRF for a Single Posterior Draw

Compute IRF for a Single Posterior Draw

## Usage

``` r
compute_irf_single(fit, draw_idx, shock, H, t0 = 1, stat_fun = stat_density)
```

## Arguments

- fit:

  dbn model fit object

- draw_idx:

  Index of posterior draw

- shock:

  Shock matrix

- H:

  Number of horizons

- t0:

  Shock time (for dynamic models, 1-based)

- stat_fun:

  Network statistic function

## Value

Vector of IRF values at each horizon
