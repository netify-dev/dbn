# Build Shock Matrix for IRF Analysis

Creates shock matrices for different types of network interventions

## Usage

``` r
build_shock(
  m,
  type = c("unit_edge", "node_out", "node_in", "density"),
  i = 1,
  j = 2,
  magnitude = 1,
  n_col = m
)
```

## Arguments

- m:

  Number of sender nodes

- type:

  Type of shock: "unit_edge", "node_out", "node_in", or "density"

- i:

  Source node index

- j:

  Target node index (for unit_edge)

- magnitude:

  Shock magnitude

- n_col:

  Number of receiver nodes (default: m)

## Value

n_row x n_col shock matrix

## See also

[`compute_irf`](https://netify-dev.github.io/dbn/reference/compute_irf.md),
[`plot.dbn_irf`](https://netify-dev.github.io/dbn/reference/plot.dbn_irf.md)

## Examples

``` r
# Unit edge shock: activate edge from node 1 to node 2
S <- build_shock(m = 5, type = "unit_edge", i = 1, j = 2)
str(S)
#>  num [1:5, 1:5] 0 0 0 0 0 1 0 0 0 0 ...

# Node-level shock: activate all outgoing edges from node 1
S_out <- build_shock(m = 5, type = "node_out", i = 1)

# Density shock: small shock spread across the network
S_dens <- build_shock(m = 5, type = "density", magnitude = 0.1)

# Bipartite shock: 4 senders, 6 receivers
S_bip <- build_shock(m = 4, n_col = 6, type = "unit_edge", i = 1, j = 3)
```
