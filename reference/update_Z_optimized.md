# Update Z for Ordinal or Binary Data

Dispatches Z updates using either Gaussian approximation or exact
truncated normal sampling based on data characteristics

## Usage

``` r
update_Z_optimized(R, Z, Theta, M, IR = NULL, family = "ordinal")
```

## Arguments

- R:

  Rank data array

- Z:

  Current latent values

- Theta:

  Current theta values

- M:

  Mean array

- IR:

  Rank indices (only needed for exact sampling)

- family:

  Data family

## Value

Updated Z array
