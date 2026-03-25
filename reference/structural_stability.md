# Structural Stability Scan

rolling window analysis to detect structural changes

## Usage

``` r
structural_stability(
  Y,
  family = "ordinal",
  window_size = NULL,
  step_size = NULL,
  nscan = 500,
  verbose = TRUE
)
```

## Arguments

- Y:

  data array (n_row x n_col x p x Tt)

- family:

  distribution family

- window_size:

  rolling window size

- step_size:

  step between windows

- nscan:

  MCMC iterations per window

- verbose:

  show progress

## Value

list with change magnitudes and suggested boundaries
