# Safe Inverse-Gamma Sampling

Draws from inverse-Gamma on log-scale to avoid under/overflow

## Usage

``` r
safe_rinv_gamma(shape, rate, floor = 1e-08, ceiling = 1e+08)
```

## Arguments

- shape:

  Shape parameter

- rate:

  Rate parameter

- floor:

  Minimum return value

- ceiling:

  Maximum return value

## Value

Sample from inverse-gamma distribution
