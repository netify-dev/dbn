# Compute impulse response for time-varying A,B matrices

Compute impulse response for time-varying A,B matrices

## Usage

``` r
impulse_response_dynamic(Aarray, Barray, S, t0, H)
```

## Arguments

- Aarray:

  Cube of A matrices over time (m x m x T)

- Barray:

  Cube of B matrices over time (m x m x T)

- S:

  Shock matrix (m x m)

- t0:

  Time index of shock (0-based)

- H:

  Number of horizons to compute

## Value

Cube of impulse responses (m x m x H+1)
