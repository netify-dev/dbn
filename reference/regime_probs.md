# Extract regime probabilities for HMM models

Compute posterior probabilities of regime assignments

## Usage

``` r
regime_probs(fit)
```

## Arguments

- fit:

  A dbn_hmm model fit object

## Value

Matrix of regime probabilities (T x R) or NULL

## See also

[`plot_regime_probs`](https://netify-dev.github.io/dbn/reference/plot_regime_probs.md),
[`param_summary`](https://netify-dev.github.io/dbn/reference/param_summary.md),
[`dbn_hmm`](https://netify-dev.github.io/dbn/reference/dbn_hmm.md)

## Examples

``` r
# \donttest{
sim <- simulate_hmm_dbn(n = 6, time = 10, R = 2, seed = 1)
fit <- dbn(sim$Y, model = "hmm", R = 2, nscan = 200, burn = 100, verbose = FALSE)
rp <- regime_probs(fit)
# }
```
