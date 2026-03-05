# dbn: Dynamic Bilinear Network Models for Ordinal Relational Data

Implements Dynamic Bilinear Network (DBN) models for analyzing temporal
relational data. Supports ordinal (rank likelihood), continuous
(Gaussian), and binary (probit) outcome types across static, dynamic,
low-rank, and HMM model variants for both unipartite and bipartite
directed networks. Includes functions for model estimation via MCMC,
visualization, impulse response analysis, and diagnostics.

## Known limitations

- Bipartite HMM/lowrank:

  The HMM and low-rank models currently support only unipartite (square)
  networks. Bipartite data will produce an informative error. Use
  `model = "dynamic"` for bipartite networks.

- Dynamic binary with small networks:

  The dynamic model with `family = "binary"` may encounter numerical
  singularities when the network has fewer than ~15 nodes. A warning is
  issued at model entry. Consider `model = "static"` or a larger
  network.

- HMM label switching:

  Regime numbering (1, 2, ..., R) in the HMM model is arbitrary and may
  differ across MCMC chains. Compare regimes by their estimated A/B
  matrices, not by label.

- Lowrank Stiefel identifiability:

  The orthonormal factor matrix U in the low-rank model is identified
  only up to orthogonal rotation. Factor loadings \\\alpha_t\\ and U
  should be interpreted together.

## See also

Useful links:

- <https://netify-dev.github.io/dbn/>

- <https://github.com/netify-dev/dbn>

- Report bugs at <https://github.com/netify-dev/dbn/issues>

## Author

**Maintainer**: Shahryar Minhas <minhassh@msu.edu>

Authors:

- Tosin Salau <salaubol@msu.edu>
