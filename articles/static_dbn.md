# Getting started: static bilinear network model

## 1 Simulate a small static network

The static model estimates time-invariant sender and receiver effects.
We start by simulating ordinal relational data.

``` r
set.seed(1)
sim = simulate_static_dbn(n = 15,            # network size: 15 actors
                          p   =  1,            # single relation type
                          time  =  5)          # 5 replicate time periods
Y = sim$Y
dim(Y)   # [n_row, n_col, p, time]
#> [1] 15 15  1  5
```

## 2 Fit the model

We pass the simulated array to
[`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md) and specify
the model type and outcome family. The MCMC sampler draws from the
posterior over latent positions and variance parameters.

``` r
fit_static = dbn(sim$Y,
                  model = 'static',
                  family = 'ordinal',   # rank-likelihood for ordinal data
                  nscan = 600,          # MCMC iterations after burn-in
                  burn = 200,           # burn-in period
                  odens = 1             # save every iteration
                  )
```

## 3 Convergence diagnostics

Always check convergence before interpreting results. You want to see
stable, well-mixed chains with no visible drift or trend. If a trace
plot shows the chain wandering or stuck in one region, increase `nscan`
and `burn`.

``` r
check_convergence(fit_static)
#>         s2         t2         g2 
#>   0.000000 600.000000   1.615836
#> 
#> Fraction in 1st window = 0.1
#> Fraction in 2nd window = 0.5 
#> 
#>    s2    t2    g2 
#>   NaN  1.76 -5.79
plot_trace(fit_static, pars = c("s2", "t2", "g2"))
```

![](static_dbn_files/figure-html/convergence-1.png)

## 4 Model summary and parameter inspection

[`summary()`](https://rdrr.io/r/base/summary.html) shows posterior means
and credible intervals for the variance parameters.
[`param_summary()`](https://netify-dev.github.io/dbn/reference/param_summary.md)
returns a tidy data frame for further analysis or plotting.

``` r
summary(fit_static)
param_summary(fit_static)
#>   parameter         mean           sd           q5          q50          q95
#> 1        s2   1.00000000 0.000000e+00   1.00000000   1.00000000   1.00000000
#> 2        t2   0.03411929 9.043186e-03   0.02222235   0.03283632   0.05078364
#> 3        g2 501.75461287 2.744300e+02 147.71915496 451.92929611 980.56082492
```

## 5 Latent mean structure

The baseline mean M captures the average relational tendency for each
dyad. Extract and summarize M across posterior draws:

``` r
M_summary = latent_summary(fit_static, fun = mean)
head(M_summary)
#>   i j rel      value
#> 1 1 1   1   0.000000
#> 2 2 1   1   3.113400
#> 3 3 1   1 -40.374945
#> 4 4 1   1 -33.221025
#> 5 5 1   1 -18.398423
#> 6 6 1   1  -6.575227
```

## 6 Gaussian and binary families

The static model supports all three outcome families. Switch the
`family` argument depending on your data: use `"gaussian"` for
continuous outcomes, `"binary"` for 0/1 data with a probit link, or
`"ordinal"` when you only trust the rank ordering of ties.

**Gaussian data** (continuous outcomes):

``` r
# simulate continuous data
set.seed(10)
n = 10; p = 1; time = 5
Z_gauss = array(rnorm(n * n * p * time), dim = c(n, n, p, time))
for (t in 1:time) diag(Z_gauss[,,1,t]) = NA

fit_gauss = dbn(Z_gauss, model = "static", family = "gaussian",
                 nscan = 300, burn = 150, verbose = FALSE)
summary(fit_gauss)
```

**Binary data** (0/1 outcomes via probit link):

``` r
set.seed(42)
n_bin = 8
Y_bin = array(rbinom(n_bin * n_bin * p * time, 1, 0.5),
               dim = c(n_bin, n_bin, p, time))
for (t in 1:time) diag(Y_bin[,,1,t]) = NA

fit_bin = dbn(Y_bin, model = "static", family = "binary",
               nscan = 300, burn = 150, verbose = FALSE)
summary(fit_bin)
```

## 7 Wrap-up

- We simulated ordinal network data,
- fit a static bilinear model with MCMC,
- checked convergence,
- inspected parameter posteriors and latent means,
- showed how to switch between ordinal, gaussian, and binary families.

For models with time-varying parameters, see
[`vignette("dynamic_dbn")`](https://netify-dev.github.io/dbn/articles/dynamic_dbn.md).
