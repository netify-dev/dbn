# Network Statistic: In-Degree

Column sums of network matrix. Includes the diagonal entry for each
node. In typical DBN usage the diagonal is `NA` or zero, so the result
matches the conventional in-degree. If your matrix has non-zero diagonal
entries, subtract them manually.

## Usage

``` r
stat_in_degree(X)
```

## Arguments

- X:

  Network matrix

## Value

Vector of in-degrees

## See also

[`compute_irf`](https://netify-dev.github.io/dbn/reference/compute_irf.md),
[`build_shock`](https://netify-dev.github.io/dbn/reference/build_shock.md),
[`stat_density`](https://netify-dev.github.io/dbn/reference/stat_density.md),
[`stat_out_degree`](https://netify-dev.github.io/dbn/reference/stat_out_degree.md)

## Examples

``` r
X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
stat_in_degree(X)
#> [1] 1 2 0
```
