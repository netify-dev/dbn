# Automatic Block Selection

two-stage algorithm for selecting number of blocks

## Usage

``` r
select_K_auto(
  Y,
  family = "ordinal",
  K_min = 1L,
  K_max = NULL,
  nscan_stage1 = 500,
  nscan_stage2 = 2000,
  verbose = TRUE,
  ...
)
```

## Arguments

- Y:

  data array (n_row x n_col x p x Tt)

- family:

  distribution family

- K_min:

  minimum blocks to consider

- K_max:

  maximum blocks to consider

- nscan_stage1:

  MCMC iterations for stage 1

- nscan_stage2:

  MCMC iterations for stage 2

- verbose:

  show progress

- ...:

  additional arguments passed to dbn_piecewise

## Value

list with selected_K, selected_boundaries, comparison results
