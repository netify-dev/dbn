# Dynamic bilinear network models

## 1 The model

The dynamic bilinear model estimates how past network interactions
predict future interactions, and how that predictive mapping itself
changes over time. The core equation is:

$$\Theta_{t} = A_{t}\,\Theta_{t - 1}\, B_{t}^{\top} + M + \varepsilon_{t}$$

where $\Theta_{t}$ is the latent interaction state at time $t$ (an
$n \times n$ matrix), $A_{t}$ and $B_{t}$ are the time-varying sender
and receiver influence matrices, and $M$ is a baseline mean capturing
persistent dyad-specific tendencies.

**What each parameter tells you:**

- **$A_{t}$ (sender influence):** An $n \times n$ matrix where entry
  $a_{i,i\prime,t}$ measures how strongly actor $i\prime$’s past
  behavior as a sender predicts actor $i$’s current sending behavior.
  Large positive values indicate that $i$ follows or echoes $i\prime$;
  negative values indicate opposing patterns.
- **$B_{t}$ (receiver influence):** Entry $b_{j,j\prime,t}$ captures how
  past targeting of actor $j\prime$ predicts current targeting of $j$.
  This reveals shifting alliances and targeting patterns.
- **$M$ (baseline mean):** Captures stable, time-invariant dyad-specific
  tendencies. A large $M_{ij}$ indicates a persistent tendency for actor
  $i$ to interact with actor $j$, regardless of dynamics.
- **$\sigma^{2}$ (process variance):** How much unexplained variation
  remains in the latent network after accounting for the bilinear
  dynamics.
- **$\tau_{A}^{2}$, $\tau_{B}^{2}$ (innovation variances):** How rapidly
  the influence matrices change from one period to the next. Large
  values mean the influence structure is reorganizing quickly; small
  values mean it is relatively stable. If $\tau_{A}^{2}$ is near zero,
  the data favor a static influence structure.
- **$\rho_{A}$, $\rho_{B}$ (AR(1) coefficients, optional):** How
  persistent the influence matrices are. Values near 1 indicate strong
  persistence; values near 0 indicate rapid mean reversion.

See the [methodology
page](https://netify-dev.github.io/dbn/articles/methodology.md) and the
[working
paper](https://netify-dev.github.io/dbn/salau_minhas_when_shocks_cascade.pdf)
for the full mathematical treatment.

## 2 Simulate a dynamic data set

We generate a small network with 12 actors observed over 20 time
periods. The simulation creates time-varying sender and receiver
influence matrices.

``` r
set.seed(6886)
sim = simulate_dynamic_dbn(n = 12,   # network size: 12 actors
                           p   =  1,
                           time  = 20)   # total time periods: 20
Y = sim$Y
dim(Y)
#> [1] 12 12  1 20
```

## 3 Fit the dynamic model

We fit the model with `ar1 = TRUE` to allow autoregressive persistence
in the influence matrices, and `update_rho = TRUE` so the AR(1)
coefficients are estimated from the data. The `time_thin` argument
controls how many time points are stored in the output; set it to 1 to
keep all of them.

``` r
fit_dyn = dbn(Y,
               family = "ordinal",  # ordinal data via rank likelihood
               model = 'dynamic',   # time-varying A_t, B_t
               nscan = 800,
               burn = 400,
               odens = 4,           # save every 4th iteration
               time_thin = 1,       # save every time period
               ar1 = TRUE,          # enable autoregressive dynamics
               update_rho = TRUE,
               verbose = TRUE)
```

## 4 Convergence diagnostics

Always check convergence before interpreting results. The trace plots
show the MCMC chains for each scalar parameter:

- **sigma2**: process variance — should stabilize, indicating the
  sampler has found the right noise level
- **tau_A2 / tau_B2**: innovation variances — these tell you how much
  the influence matrices change between periods
- **rho_A / rho_B**: AR(1) persistence — values near 1 indicate strong
  temporal persistence in the influence structure

``` r
check_convergence(fit_dyn)
#>    sigma2    tau_A2    tau_B2        g2 
#>  2.066464 61.515167 36.456775 51.913459
#> 
#> Fraction in 1st window = 0.1
#> Fraction in 2nd window = 0.5 
#> 
#>  sigma2  tau_A2  tau_B2      g2 
#> -8.2616 -0.2645 -2.5402 -3.6347
```

![](dynamic_dbn_files/figure-html/conv-dyn-1.png)

``` r
plot_trace(fit_dyn, pars = c("sigma2", "tau_A2", "tau_B2", "rho_A", "rho_B"))
```

![](dynamic_dbn_files/figure-html/conv-dyn-2.png)

## 5 Forecasting

Generate forecasts 6 time steps ahead. The model propagates the
estimated dynamics forward:

``` r
Theta_forecast = predict(fit_dyn, H = 6, S = 200, summary = "mean")
dim(Theta_forecast)   # n_row x n_col x p x H
#> [1] 12 12  1  6
```

## 6 Dyad trajectory

Visualize how a specific bilateral relationship evolves over time. The
ribbon shows the posterior credible interval — wider ribbons indicate
more uncertainty about that dyad’s trajectory:

``` r
dyad_path(fit_dyn, i = 2, j = 7)
```

![](dynamic_dbn_files/figure-html/dyad-path-1.png)

## 7 Group-level influence

Track how a group of senders’ aggregate influence changes through time.
This is useful for identifying when a coalition of actors becomes more
or less structurally important in the network:

``` r
plot_group_influence(fit_dyn,
                     group = c(1, 3, 5),
                     type  = "sender",
                     measure = "rowsum",
                     cred = 0.9)
```

![](dynamic_dbn_files/figure-html/group-influence-1.png)

## 8 Parameter summaries

[`param_summary()`](https://netify-dev.github.io/dbn/reference/param_summary.md)
returns posterior means and credible intervals for all scalar parameters
in a tidy data frame.

``` r
param_summary(fit_dyn)
#>   parameter         mean           sd           q5          q50          q95
#> 1    sigma2 3.805361e+03 1.305657e+03 1.870896e+03 3.817688e+03 6.043277e+03
#> 2    tau_A2 7.416453e-02 5.577535e-03 6.387336e-02 7.464214e-02 8.341511e-02
#> 3    tau_B2 9.957576e-02 1.053630e-02 8.458207e-02 9.884585e-02 1.197036e-01
#> 4        g2 1.252089e-01 3.077106e-02 8.612315e-02 1.185952e-01 1.854716e-01
#> 5      rhoA 7.907962e-01 2.107282e-02 7.570238e-01 7.909403e-01 8.268467e-01
#> 6      rhoB 7.608620e-01 2.727771e-02 7.125525e-01 7.623608e-01 8.050317e-01
```

## 9 Credible intervals and network statistics

Compute dyad-level credible intervals over time. This gives you, for
each dyad at each time point, the posterior mean and credible interval
for the latent interaction intensity:

``` r
tc = theta_credible(fit_dyn, i = 1:6, j = 1:6, time = 1:20)
head(tc)
#>   i j rel time        mean     lower      median     upper
#> 1 1 1   1    1   2.0504579 -22.18684   3.3534039  21.58073
#> 2 2 1   1    1  29.1609381  14.82135  28.6266136  42.72343
#> 3 3 1   1    1 -34.1252481 -51.40412 -33.0155181 -24.89329
#> 4 4 1   1    1  36.4255985  14.29147  37.6021427  53.06615
#> 5 5 1   1    1  -0.6984207 -17.40425  -0.6098354  12.21810
#> 6 6 1   1    1 -60.2127657 -75.67646 -59.7203886 -44.05003
```

You can use this output to flag dyads whose credible intervals exclude
zero at a given time point, indicating a reliably non-zero relationship.

Track network-level statistics over time. Network density summarizes how
active the overall network is at each time point, with posterior
uncertainty:

``` r
ns = network_summary(fit_dyn, stat = "density")
head(ns)
#>   time      mean     lower     upper
#> 1    1 0.4722727 0.4469697 0.5000000
#> 2    2 0.5009470 0.4696970 0.5303030
#> 3    3 0.5054545 0.4768939 0.5378788
#> 4    4 0.4954924 0.4621212 0.5303030
#> 5    5 0.5085227 0.4545455 0.5609848
#> 6    6 0.4986364 0.4621212 0.5303030
```

Look for time points where the credible interval is narrow versus wide –
narrow intervals mean the model is confident about the overall level of
activity at that point.

Compute the posterior probability that each edge is positive at a given
time. Values near 1 indicate edges that are consistently positive across
posterior draws; values near 0.5 indicate high uncertainty:

``` r
ep = edge_prob(fit_dyn, rel = 1, time = 20)
ep[1:5, 1:5]
#>      [,1] [,2]  [,3]  [,4]  [,5]
#> [1,]   NA 0.00 1.000 1.000 0.000
#> [2,]    1   NA 1.000 0.025 0.000
#> [3,]    0 0.02    NA 1.000 0.995
#> [4,]    1 0.00 1.000    NA 0.025
#> [5,]    1 1.00 0.015 1.000    NA
```

You can threshold these probabilities (e.g., \> 0.95) to construct a
binary network of edges that are reliably positive at a given time
point.

## 10 Bipartite (rectangular) networks

The package automatically detects bipartite networks when senders and
receivers differ. Simply pass a rectangular array:

``` r
sim_bp = simulate_dynamic_dbn(n = 10, n_col = 8, p = 1, time = 15, seed = 77)
dim(sim_bp$Y)   # 10 senders x 8 receivers x 1 relation x 15 time points
#> [1] 10  8  1 15

fit_bp = dbn(sim_bp$Y, model = "dynamic", family = "ordinal",
              nscan = 400, burn = 200, odens = 2, verbose = FALSE)
summary(fit_bp)
```

When `n_row != n_col`, $A$ is `n_row x n_row` (sender dynamics) and $B$
is `n_col x n_col` (receiver dynamics). Self-loops are retained since
senders and receivers are distinct.

## 11 Parallel computation

The dynamic model parallelizes the row-wise A and B FFBS updates using
OpenMP. Set the number of threads before fitting:

``` r
# use physical cores (not logical/hyperthreaded)
set_dbn_threads(parallel::detectCores(logical = FALSE))
get_dbn_threads()
#> [1] 4
```

Parallelization benefits scale with network size: for networks with 15+
actors, expect 2–4x speedup with 4–8 cores. For small networks (n \<
10), the overhead of thread management may exceed the benefit, so
single-threaded execution is fine.

OpenMP is available by default on Linux and Windows. On macOS, install
it via `brew install libomp` (see the README for details). If OpenMP is
not available, the package falls back to single-threaded execution
automatically.

To reset to single-threaded:

``` r
set_dbn_threads(1)
```

## 12 Scaling to larger networks

The dynamic model can handle networks with many actors. The main
constraint is RAM, not computation. Use
[`estimate_memory()`](https://netify-dev.github.io/dbn/reference/estimate_memory.md)
to check requirements before fitting, and adjust `time_thin` and `odens`
to manage memory:

``` r
# a 50-actor network
estimate_memory(n_row = 50, n_col = 50, p = 1, Tt = 30,
                nscan = 5000, burn = 1000, odens = 5,
                time_thin = 2, family = "ordinal")
#> Dynamic DBN memory estimate:
#>   Network:   50 x 50, 1 relation(s), 30 time points
#>   MCMC:      1000 draws (nscan=5000, odens=5)
#>   Time thin: 2 (keeping 15 of 30 time points)
#>   --------------------------------
#>   Theta:    0.28 GB
#>   Z:        0.28 GB
#>   A:        0.28 GB
#>   B:        0.28 GB
#>   M:        0.02 GB
#>   --------------------------------
#>   TOTAL:    1.14 GB
```

For networks with 100+ actors, you may need 16–32 GB of RAM depending on
the number of time points and MCMC settings. Increasing `time_thin`
(save every $k$-th time point) and `odens` (save every $k$-th MCMC draw)
are the most effective ways to reduce memory usage.
