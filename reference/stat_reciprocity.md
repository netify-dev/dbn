# Network Statistic: Reciprocity

Correlation between X\[i,j\] and X\[j,i\]. Only defined for square
(unipartite) networks. For bipartite networks, use
[`stat_density`](https://netify-dev.github.io/dbn/reference/stat_density.md),
[`stat_in_degree`](https://netify-dev.github.io/dbn/reference/stat_in_degree.md),
or
[`stat_out_degree`](https://netify-dev.github.io/dbn/reference/stat_out_degree.md)
instead.

## Usage

``` r
stat_reciprocity(X)
```

## Arguments

- X:

  Square network matrix

## Value

Scalar reciprocity value

## See also

[`compute_irf`](https://netify-dev.github.io/dbn/reference/compute_irf.md),
[`build_shock`](https://netify-dev.github.io/dbn/reference/build_shock.md),
[`stat_density`](https://netify-dev.github.io/dbn/reference/stat_density.md),
[`stat_transitivity`](https://netify-dev.github.io/dbn/reference/stat_transitivity.md)

## Examples

``` r
X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
stat_reciprocity(X)
#> [1] 0.5
```
