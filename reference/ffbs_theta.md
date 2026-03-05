# FFBS for Theta Matrix

FFBS for the latent Theta process in bilinear model

## Usage

``` r
ffbs_theta(Z, mu, Aarray, Barray, sigma2)
```

## Arguments

- Z:

  Observed data (m x m x Tt)

- mu:

  Baseline mean (m x m)

- Aarray:

  Time-varying A matrices (m x m x Tt)

- Barray:

  Time-varying B matrices (m x m x Tt)

- sigma2:

  Innovation variance

## Value

Array of sampled Theta matrices (m x m x Tt)
