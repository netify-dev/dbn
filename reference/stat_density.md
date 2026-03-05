# Network Statistic: Density

Mean of off-diagonal elements.

**Note:** `ggplot2` also exports a function called `stat_density`. If
both packages are loaded, use `dbn::stat_density` to refer to this
function unambiguously (e.g., when passing to
[`compute_irf`](https://netify-dev.github.io/dbn/reference/compute_irf.md)).

## Usage

``` r
stat_density(X)
```

## Arguments

- X:

  Network matrix

## Value

Scalar density value

## See also

[`compute_irf`](https://netify-dev.github.io/dbn/reference/compute_irf.md),
[`build_shock`](https://netify-dev.github.io/dbn/reference/build_shock.md),
[`stat_in_degree`](https://netify-dev.github.io/dbn/reference/stat_in_degree.md),
[`stat_out_degree`](https://netify-dev.github.io/dbn/reference/stat_out_degree.md),
[`stat_reciprocity`](https://netify-dev.github.io/dbn/reference/stat_reciprocity.md),
[`stat_transitivity`](https://netify-dev.github.io/dbn/reference/stat_transitivity.md)

## Examples

``` r
X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
dbn::stat_density(X)
#> [1] 0.5
```
