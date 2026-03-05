# Tidy Extractor for Low-Rank Factor Paths

Extracts factor trajectories in tidy data-frame format

## Usage

``` r
tidy_dbn_lowrank(fit, factors = NULL)
```

## Arguments

- fit:

  A dbn object with model="lowrank"

- factors:

  Which factors to extract (default: all)

## Value

Data frame with columns: time, mean, lo, hi, factor
