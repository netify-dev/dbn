# Summarize Scalar Parameters

Returns posterior mean, standard deviation, and quantiles for the scalar
variance parameters estimated by the model. These typically include:

- `sigma2` or `s2`: process noise variance

- `tau_A2` / `tau_B2` or `t2`: innovation variance for A/B

- `g2`: latent variance

- `rhoA` / `rhoB`: AR(1) persistence (dynamic model with `ar1 = TRUE`)

- `sigma2_obs`: observation variance (gaussian family)

## Usage

``` r
param_summary(fit, probs = c(0.025, 0.5, 0.975))
```

## Arguments

- fit:

  A dbn model fit object (output from
  [`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md))

- probs:

  Quantile probabilities (default: 5th, 50th, 95th percentiles)

## Value

Data frame with columns: `parameter`, `mean`, `sd`, and one column per
requested quantile

## See also

[`theta_summary`](https://netify-dev.github.io/dbn/reference/theta_summary.md),
[`plot_trace`](https://netify-dev.github.io/dbn/reference/plot_trace.md),
[`derive_draws`](https://netify-dev.github.io/dbn/reference/derive_draws.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
ps <- param_summary(fit)
# }
```
