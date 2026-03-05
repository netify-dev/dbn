# Forward-Filter Backward-Sample for DLM

FFBS for dynamic linear model with optional AR(1) dynamics

## Usage

``` r
ffbs_dlm(y, Flist, V, W, m0, C0, ar1 = FALSE, rho = 0)
```

## Arguments

- y:

  List of observation vectors

- Flist:

  List of design matrices

- V:

  Observation variance matrix

- W:

  State innovation variance matrix

- m0:

  Prior mean vector

- C0:

  Prior covariance matrix

- ar1:

  Logical: use AR(1) dynamics

- rho:

  AR(1) coefficient (used if ar1=TRUE)

## Value

Matrix of sampled state vectors (columns are time points)
