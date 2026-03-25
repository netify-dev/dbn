# Compare Influence Matrices Across Blocks

computes posterior probability of difference between blocks

## Usage

``` r
compare_blocks(fit, blocks = NULL, parameter = "A", threshold = 0.1)
```

## Arguments

- fit:

  piecewise dbn fit object

- blocks:

  pair of block indices to compare (default: adjacent blocks)

- parameter:

  which parameter to compare: "A" or "B"

- threshold:

  minimum difference to be considered meaningful

## Value

list with comparison results

## Examples

``` r
# \donttest{
sim <- simulate_piecewise_dbn(n = 10, time = 40, blocks = c(15, 30, 40))
fit <- dbn(sim$Y, model = "piecewise", blocks = 3,
           nscan = 200, burn = 100, verbose = FALSE)
compare_blocks(fit)
#> 
#> ── Block Comparison Results (A) 
#> ℹ block_1 vs block_2: ||ΔA|| = 1.895 [1.615, 2.192]
#> →   P(||ΔA|| > 0.1) = 1
#> ℹ block_2 vs block_3: ||ΔA|| = 1.91 [1.628, 2.259]
#> →   P(||ΔA|| > 0.1) = 1
# }
```
