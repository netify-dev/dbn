# Compute impulse response for constant A,B matrices

Compute impulse response for constant A,B matrices

## Usage

``` r
impulse_response_const(A, B, S, H)
```

## Arguments

- A:

  Transition matrix A (m x m)

- B:

  Transition matrix B (m x m)

- S:

  Shock matrix (m x m)

- H:

  Number of horizons to compute

## Value

Cube of impulse responses (m x m x H+1)
