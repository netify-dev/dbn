# Low-rank bilinear network models

------------------------------------------------------------------------

**NOTE: This vignette is under construction.** The low-rank model
variant is implemented and functional, but the documentation and
examples below are preliminary and subject to change.

------------------------------------------------------------------------

## 1 Overview

The low-rank model parameterizes the sender dynamics matrix as
$A_{t} = U\,\text{diag}\left( \alpha_{t} \right)\, U\prime$, where $U$
is an orthogonal matrix on the Stiefel manifold and $\alpha_{t}$ is a
vector of time-varying factor strengths. This reduces the number of free
parameters from $m^{2}$ to $m \times r + r$, making it **scalable to
larger networks** (50+ nodes).

Use `model = "lowrank"` when you have many actors and want an
interpretable factor decomposition of sender dynamics.

## 2 Simulate a low-rank network

We start by simulating a 20-node network over 25 time periods with `r`
set to 3. This gives us ground truth to check recovery of the factor
structure.

``` r
sim = simulate_lowrank_dbn(
  n = 20, p = 1, time = 25,
  r = 3,
  sigma2 = 0.3,
  tau_alpha2 = 0.05,
  tauB2 = 0.03,
  seed = 42
)

dim(sim$Y)
#> [1] 20 20  1 25
```

## 3 Memory estimation

Low-rank models store the Stiefel matrix $U$ and factor paths at each
MCMC iteration, so memory usage can grow quickly with network size.
Before fitting, you can estimate RAM requirements:

``` r
estimate_memory(
  n_row = 20, n_col = 20, p = 1, Tt = 25,
  nscan = 1000, burn = 500, odens = 2,
  family = "ordinal"
)
#> Dynamic DBN memory estimate:
#>   Network:   20 x 20, 1 relation(s), 25 time points
#>   MCMC:      500 draws (nscan=1000, odens=2)
#>   Time thin: 1 (keeping 25 of 25 time points)
#>   --------------------------------
#>   Theta:    0.04 GB
#>   Z:        0.04 GB
#>   A:        0.04 GB
#>   B:        0.04 GB
#>   M:        0.00 GB
#>   --------------------------------
#>   TOTAL:    0.15 GB
```

## 4 Fit the low-rank model

The key argument is `r`, which sets the rank of the factor
decomposition. Setting `r` much smaller than `n` is what makes this
model tractable for larger networks.

``` r
fit_lr = dbn(
  sim$Y,
  model   = "lowrank",
  family  = "ordinal",
  r       = 3,
  nscan   = 1000,
  burn    = 500,
  odens   = 2,
  verbose = FALSE
)
```

## 5 Model summary

The summary shows posterior means for the variance parameters and the
estimated rank structure.

``` r
summary(fit_lr)
#> Low-rank Dynamic Bilinear Network model
#>   nodes     : 20 
#>   relations : 1 
#>   time pts  : 25 
#>   rank      : 3 
#> 
#>            mean          2.5%          97.5%        
#> sigma2     54210.442     36011.429     78915.593    
#> tau_alpha2 0.225         0.168         0.294        
#> tau_B2     100000000.000 100000000.000 100000000.000
#> g2         0.146         0.122         0.172
```

## 6 Diagnostics

The low-rank plot shows MCMC traces, factor trajectories ($\alpha_{k}$),
and the posterior mean of the node loading matrix $U$:

``` r
plot(fit_lr)
```

![](lowrank_dbn_files/figure-html/plot-lr-1.png)

Scalar parameter traces:

``` r
plot_trace(fit_lr)
```

![](lowrank_dbn_files/figure-html/trace-lr-1.png)

## 7 Factor trajectories

The factor trajectories $\alpha_{t}$ capture how each latent dimension’s
influence on sender dynamics evolves over time. A flat trajectory means
that factor’s contribution is stable; a trending or volatile trajectory
suggests structural change along that dimension. Extract factor paths in
tidy format for custom plotting:

``` r
alpha_df = tidy_dbn_lowrank(fit_lr)
head(alpha_df)
#>   time          mean           lo          hi factor
#> 1    1  0.0053000041 -0.001305681 0.053818449      1
#> 2    2 -0.0088905063 -0.088534301 0.001363029      1
#> 3    3  0.0002627765 -0.001980068 0.002014024      1
#> 4    4  0.0070274876 -0.001512614 0.070152050      1
#> 5    5 -0.0004314414 -0.005152276 0.002132315      1
#> 6    6  0.0186692588 -0.002072614 0.189898130      1
```

``` r
ggplot(alpha_df, aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line() +
  facet_wrap(~factor, scales = "free_y", labeller = label_bquote(alpha[.(factor)])) +
  labs(x = "Time", y = expression(alpha)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text.x = element_text(color = "white", hjust = 0)
  )
```

![](lowrank_dbn_files/figure-html/factor-plot-1.png)

## 8 Network summaries

You can compute summary statistics across the posterior distribution of
fitted networks.
[`network_summary()`](https://netify-dev.github.io/dbn/reference/network_summary.md)
returns time-varying means or densities, and
[`edge_prob()`](https://netify-dev.github.io/dbn/reference/edge_prob.md)
gives the posterior probability that each edge exceeds a threshold.

``` r
network_summary(fit_lr, stat = "mean")
#>    time         mean        lower      upper
#> 1     1 -0.020421509 -0.142362807 0.09360024
#> 2     2 -0.071045382 -0.186990100 0.04609783
#> 3     3  0.035263908 -0.091393392 0.16153756
#> 4     4 -0.036851925 -0.161480314 0.07648696
#> 5     5  0.070041315 -0.044095176 0.17557397
#> 6     6  0.033265357 -0.089773516 0.14524146
#> 7     7  0.056062345 -0.058339393 0.17897512
#> 8     8 -0.048022658 -0.159074743 0.06882154
#> 9     9 -0.012539343 -0.133540748 0.10885173
#> 10   10 -0.008325071 -0.129606574 0.10150096
#> 11   11  0.073922513 -0.060109515 0.19726443
#> 12   12  0.041775746 -0.081816831 0.17374382
#> 13   13  0.069037324 -0.052511000 0.18676551
#> 14   14  0.109312886 -0.007400952 0.22620013
#> 15   15 -0.071518234 -0.190651452 0.04793111
#> 16   16 -0.095027310 -0.198798167 0.01367304
#> 17   17 -0.093121103 -0.213892633 0.02694589
#> 18   18  0.009555333 -0.105015798 0.11314555
#> 19   19 -0.059402213 -0.184659441 0.05792613
#> 20   20 -0.059519496 -0.186489963 0.06539211
#> 21   21  0.016639845 -0.105926703 0.13637307
#> 22   22 -0.043626115 -0.171100829 0.08539066
#> 23   23  0.024727449 -0.090846983 0.15473508
#> 24   24  0.011549380 -0.108555771 0.13511772
#> 25   25  0.066646083 -0.049914737 0.18069156
network_summary(fit_lr, stat = "density")
#>    time      mean     lower     upper
#> 1     1 0.5006158 0.4657895 0.5315789
#> 2     2 0.4866789 0.4552632 0.5184211
#> 3     3 0.5068000 0.4736842 0.5421053
#> 4     4 0.4950421 0.4605263 0.5264474
#> 5     5 0.5168316 0.4815789 0.5500000
#> 6     6 0.5072947 0.4736842 0.5396053
#> 7     7 0.5131579 0.4789474 0.5500000
#> 8     8 0.4958789 0.4631579 0.5289474
#> 9     9 0.5002684 0.4657895 0.5368421
#> 10   10 0.4965263 0.4631579 0.5290789
#> 11   11 0.5182474 0.4842105 0.5526316
#> 12   12 0.5093684 0.4736842 0.5421053
#> 13   13 0.5075579 0.4763158 0.5421053
#> 14   14 0.5232000 0.4894737 0.5552632
#> 15   15 0.4867105 0.4526316 0.5236842
#> 16   16 0.4825368 0.4526316 0.5131579
#> 17   17 0.4841737 0.4473684 0.5210526
#> 18   18 0.5040526 0.4684211 0.5368421
#> 19   19 0.4870684 0.4526316 0.5238158
#> 20   20 0.4932684 0.4605263 0.5289474
#> 21   21 0.5014737 0.4684211 0.5368421
#> 22   22 0.4928421 0.4526316 0.5290789
#> 23   23 0.5057579 0.4684211 0.5394737
#> 24   24 0.4916263 0.4578947 0.5263158
#> 25   25 0.5110211 0.4763158 0.5447368
edge_prob(fit_lr, threshold = 0)
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,]    NA 0.586 0.574 0.024 0.000 0.386 0.282 0.382 0.126 0.516 0.396 0.652
#>  [2,] 0.202    NA 0.674 0.088 0.392 0.406 0.120 0.776 0.040 0.522 0.694 0.812
#>  [3,] 0.150 0.188    NA 0.500 0.780 0.728 0.774 0.858 0.088 0.628 0.526 0.828
#>  [4,] 0.080 0.468 0.910    NA 0.910 0.110 0.064 0.686 0.492 0.106 0.536 0.240
#>  [5,] 0.096 0.154 0.514 0.014    NA 0.388 0.842 0.478 0.154 0.198 0.858 0.290
#>  [6,] 0.094 0.876 0.484 0.860 1.000    NA 0.174 0.538 0.986 0.060 0.376 0.672
#>  [7,] 0.588 0.132 0.270 0.592 0.898 0.880    NA 0.910 0.166 0.248 0.788 0.658
#>  [8,] 0.210 0.960 0.988 0.770 0.882 0.884 0.726    NA 0.532 0.562 0.424 0.348
#>  [9,] 0.692 0.202 0.718 0.638 0.594 0.408 0.938 0.314    NA 0.054 0.998 0.168
#> [10,] 0.520 0.710 0.738 0.240 0.996 0.880 0.776 0.828 0.972    NA 0.784 0.562
#> [11,] 0.040 0.728 0.036 0.620 0.678 0.668 0.724 0.702 0.936 0.930    NA 0.810
#> [12,] 0.668 0.614 0.212 0.096 0.260 0.432 0.210 0.744 0.822 0.532 0.492    NA
#> [13,] 0.156 0.648 0.780 0.316 0.902 0.122 0.328 0.206 0.108 0.988 0.886 0.882
#> [14,] 0.182 0.740 0.224 0.812 0.076 0.772 0.906 0.684 0.752 0.652 0.078 0.840
#> [15,] 0.836 0.500 0.066 0.710 0.456 0.188 0.280 0.344 0.142 0.262 0.878 0.232
#> [16,] 0.552 0.204 0.560 0.968 0.404 0.890 0.924 0.064 0.450 0.242 0.866 0.096
#> [17,] 0.322 0.050 0.070 0.168 0.142 0.188 0.376 0.062 0.218 0.558 0.078 0.208
#> [18,] 0.834 0.616 0.644 0.492 0.350 0.780 0.900 0.232 0.036 0.654 0.192 0.780
#> [19,] 0.990 0.906 0.096 0.974 0.042 0.668 0.134 0.980 0.448 0.308 0.338 0.620
#> [20,] 0.302 0.716 0.188 0.776 0.798 0.472 0.152 0.414 0.932 0.648 0.552 0.262
#>       [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
#>  [1,] 0.492 0.754 0.338 0.216 0.954 0.794 0.816 0.570
#>  [2,] 0.754 0.202 0.232 0.090 0.420 0.756 0.820 0.192
#>  [3,] 0.016 0.466 0.212 0.744 0.174 0.826 0.900 0.052
#>  [4,] 0.916 0.434 0.880 0.932 0.432 0.774 0.234 0.362
#>  [5,] 0.582 0.544 0.254 0.492 0.988 0.290 0.916 0.142
#>  [6,] 0.386 0.682 0.654 0.696 0.950 0.324 0.908 0.236
#>  [7,] 0.520 0.102 0.856 0.926 0.262 0.298 0.290 0.440
#>  [8,] 0.822 0.432 0.654 0.200 0.926 0.944 0.918 0.166
#>  [9,] 0.514 0.684 0.878 0.954 0.094 0.726 0.074 0.146
#> [10,] 0.238 0.470 0.416 0.006 0.752 0.870 0.704 0.254
#> [11,] 0.604 0.028 0.956 0.672 0.354 0.798 0.248 0.684
#> [12,] 0.328 0.882 0.836 0.250 0.268 0.846 0.696 0.804
#> [13,]    NA 0.708 0.412 0.818 0.874 0.634 0.364 0.050
#> [14,] 0.186    NA 0.548 0.712 0.218 0.658 0.650 0.250
#> [15,] 0.242 0.240    NA 0.338 0.122 0.110 0.370 0.958
#> [16,] 0.880 0.954 0.222    NA 0.100 0.872 0.806 0.292
#> [17,] 0.148 0.392 0.810 0.084    NA 0.654 0.938 0.634
#> [18,] 0.922 0.294 0.428 0.596 0.434    NA 0.926 0.150
#> [19,] 0.418 0.894 0.344 0.078 0.024 0.704    NA 0.682
#> [20,] 0.680 0.026 0.298 0.726 0.882 0.556 0.964    NA
```

## 9 Dyad trajectory

Track a specific sender-receiver pair over time with posterior credible
intervals.

``` r
dyad_path(fit_lr, i = 2, j = 7)
```

![](lowrank_dbn_files/figure-html/dyad-lr-1.png)

## 10 Forecasting

[`predict()`](https://rdrr.io/r/stats/predict.html) propagates the
low-rank dynamics forward by `H` steps, sampling from the posterior to
produce forecast distributions.

``` r
pred = predict(fit_lr, H = 3, draws = 50, summary = "mean")
dim(pred)
#> [1] 20 20  1  3
```

## 11 When to use low-rank

The low-rank model is appropriate when:

- The network has **many actors** (50+) making the full dynamic model
  computationally expensive
- You want an **interpretable factor structure** where $U$ reveals actor
  groupings and $\alpha_{t}$ captures time-varying factor strengths
- The rank $r$ is much smaller than the number of actors

For smaller networks, `model = "dynamic"` provides more flexible
dynamics. For regime-switching behavior, consider `model = "hmm"`.
