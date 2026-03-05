# Update A and B with Pre-allocated Memory

A/B updates using batch C++ implementation

## Usage

``` r
update_AB_optimized(
  Theta_all,
  Aarray,
  Barray,
  sigma2,
  tauA2,
  tauB2,
  ar1 = FALSE,
  rhoA = 0,
  rhoB = 0
)
```

## Arguments

- Theta_all:

  Full Theta array (n_row x n_col x p x Tt)

- Aarray:

  Current A array (n_row x n_row x Tt)

- Barray:

  Current B array (n_col x n_col x Tt)

- sigma2:

  Observation variance

- tauA2:

  A innovation variance

- tauB2:

  B innovation variance

- ar1:

  Use AR(1) dynamics

- rhoA:

  AR(1) coefficient for A

- rhoB:

  AR(1) coefficient for B

## Value

List with updated A and B arrays
