# Piecewise DBN MCMC

fits DBN model with block-constant influence matrices

## Usage

``` r
dbn_piecewise(
  Y,
  family = c("ordinal", "gaussian", "binary"),
  blocks = 4,
  nscan = 10000,
  burn = 1000,
  odens = 1,
  seed = 6886,
  verbose = TRUE,
  symmetric = FALSE,
  previous = NULL,
  store_theta = TRUE
)
```

## Arguments

- Y:

  data array (n_row x n_col x p x Tt)

- family:

  distribution family: "ordinal", "gaussian", or "binary"

- blocks:

  block specification: integer, vector, or parsed block_info

- nscan:

  number of MCMC iterations after burn-in

- burn:

  burn-in iterations

- odens:

  output density (save every odens-th iteration)

- seed:

  random seed

- verbose:

  show progress

- symmetric:

  enforce B = A (symmetric networks)

- previous:

  previous fit to continue from

- store_theta:

  Logical. Store full Theta trajectory draws (default TRUE). **For large
  networks (100+ actors), set to FALSE.** Theta storage scales as O(n^2
  \* T \* draws) and becomes prohibitive for large networks (e.g., 200
  actors with 50 time points and 500 draws requires ~40 GB). With FALSE,
  you retain A, B, M draws and
  [`compare_blocks()`](https://netify-dev.github.io/dbn/reference/compare_blocks.md)
  functionality but lose
  [`theta_slice()`](https://netify-dev.github.io/dbn/reference/theta_slice.md),
  [`theta_summary()`](https://netify-dev.github.io/dbn/reference/theta_summary.md),
  and
  [`posterior_predict_dbn()`](https://netify-dev.github.io/dbn/reference/posterior_predict_dbn.md)
  with uncertainty.

## Value

list containing MCMC results

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md) for main
dispatcher
