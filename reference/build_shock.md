# Build Shock Matrix for IRF Analysis

Creates structured perturbation matrices for impulse response analysis.
The shock matrix S is added to the latent network at a single time
point, and
[`compute_irf()`](https://netify-dev.github.io/dbn/reference/compute_irf.md)
tracks how the perturbation propagates.

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

  Number of actors (senders). This should match the fitted model's
  network size.

- type:

  Type of shock:

  - `"unit_edge"`: perturb a single directed edge (i -\> j)

  - `"node_out"`: perturb all outgoing edges from actor i

  - `"node_in"`: perturb all incoming edges to actor i

  - `"density"`: uniform perturbation to all edges

- i:

  Source actor index (sender)

- j:

  Target actor index (receiver), used only for `"unit_edge"`

- magnitude:

  Size of the perturbation. Scale relative to your data: e.g., if data
  ranges from -1 to 1, a magnitude of 0.5 is a large shock.

- n_col:

  Number of receiver nodes (default: same as `m`)

## Value

An `m x n_col` matrix of zeros with the specified entries set to
`magnitude`. Pass this to
[`compute_irf()`](https://netify-dev.github.io/dbn/reference/compute_irf.md)
as the `shock` argument.

## Details

**When to use each shock type:**

- `"unit_edge"`: Shock a single bilateral tie (e.g., "what if USA-Russia
  relations improve?"). Use when you care about a specific dyad.

- `"node_out"`: Shock all outgoing ties from one actor (e.g., "what if
  China engages all partners at once?"). Use when you care about one
  actor's overall engagement.

- `"node_in"`: Shock all incoming ties to one actor.

- `"density"`: Apply a uniform shock to all edges. Use to study how a
  system-wide shift propagates.

For **symmetric (undirected) networks**, symmetrize the shock after
building it: `S = S + t(S)`.

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
