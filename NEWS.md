# dbn 1.2.0 (2026-05-28)

* Dynamic models now apply a light contractivity prior to the lag operator by default (`shrink_rho = 0.9`), keeping impulse responses and forecasts stable out of the box. Lower it to regularise more, or set `shrink_rho = NULL` to turn it off.
* `dbn()` now exposes the underlying sampler via `sampler = "auto" | "approx" | "exact"`. The default chooses based on the model.
* AR(1) persistence on the lag operators: pass `ar1 = TRUE` to give `A_t` and `B_t` mean-reverting dynamics, and `update_rho = TRUE` to estimate the persistence parameters.
* Covariates: dyadic, sender, and receiver covariates handled via a block-Gibbs sampler, with configurable priors (`prior_beta_scale`, `prior_kind`), random actor effects, and time-varying coefficients.
* Smart initialisation (`init = "smart"`, the default) runs a brief ALS warm-start before MCMC.
* `posterior` / `tidybayes` / `bayesplot` interop via `as_draws.dbn()` and `as_draws.dbn_multichain()`, with optional per-cell operator entries.
* ALS+bootstrap fits expand into per-replicate operator draws so `compute_irf()`, `predict()`, `dbn_compute_snr()`, and the rest of the accessor surface work uniformly across MCMC and ALS fits.
* Vignettes refreshed: methodology opens with a worked example, the dynamic vignette includes an AR(1) recovery demo, and the applied / IRF / piecewise vignettes build noticeably faster.
