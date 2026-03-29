# Applied Example: Foreign Policy Alignment Among Major Powers

## 1 Overview

This vignette walks through a complete `dbn` workflow on real data,
fitting a dynamic model to foreign policy alignment among 12 major
powers (G7 + BRICS) and using impulse response functions (IRFs) to study
how shocks propagate through the network.

We use UNGA voting similarity scores from the `peacesciencer` package.
The DBN model captures both the time-varying latent structure of the
alignment network and the transition dynamics that govern how it
evolves. IRFs then allow us to ask counterfactual questions that
standard dyadic approaches cannot address, such as what happens to the
broader network if USA-Russia alignment suddenly improves.

**Requirements:** This vignette requires the `peacesciencer` package and
its external data. Install with:

``` r
install.packages("peacesciencer")
peacesciencer::download_extdata()
```

The code below is shown but not evaluated during the vignette build. To
reproduce the results, copy the code blocks into an R session with
`peacesciencer` installed.

## 2 Data preparation

We use the `kappavv` measure from
[`peacesciencer::add_fpsim()`](https://rdrr.io/pkg/peacesciencer/man/add_fpsim.html),
which captures pairwise agreement in UNGA voting. It is continuous,
symmetric ($\kappa(i,j) = \kappa(j,i)$), and ranges from $- 1$
(systematic disagreement) to $+ 1$ (perfect agreement).

``` r
library(dbn)
library(ggplot2)
library(peacesciencer)

# G7 + BRICS: 12 actors
actor_codes = c(
  USA = 2, Canada = 20, UK = 200, France = 220, Germany = 255,
  Italy = 325, Japan = 740, Brazil = 140, Russia = 365,
  India = 750, China = 710, SouthAfrica = 560
)
actor_labels = names(actor_codes)
n = length(actor_codes)

usa_idx    = which(actor_labels == "USA")
russia_idx = which(actor_labels == "Russia")
china_idx  = which(actor_labels == "China")

# build dyadic panel 2005-2015
dyads = create_dyadyears(subset_years = 2005:2015) |>
  add_fpsim()

dyads_sub = dyads[dyads$ccode1 %in% actor_codes &
                    dyads$ccode2 %in% actor_codes, ]

years = sort(unique(dyads_sub$year))
Tt = length(years)

# build [n, n, 1, T] array
code_to_idx = setNames(seq_along(actor_codes), actor_codes)
Y = array(NA_real_, dim = c(n, n, 1, Tt),
          dimnames = list(actor_labels, actor_labels, "unga",
                          as.character(years)))

for (i in seq_len(nrow(dyads_sub))) {
  row_i = code_to_idx[as.character(dyads_sub$ccode1[i])]
  col_j = code_to_idx[as.character(dyads_sub$ccode2[i])]
  t_idx = which(years == dyads_sub$year[i])
  Y[row_i, col_j, 1, t_idx] = dyads_sub$kappavv[i]
}
for (t in 1:Tt) diag(Y[, , 1, t]) = NA
```

## 3 Fit the dynamic model

Since `kappavv` is symmetric, we set `symmetric = TRUE` to constrain
$B_{t} = A_{t}$, halving the parameter count. This is appropriate
because the sender-receiver distinction is not meaningful for undirected
voting similarity.

``` r
fit = dbn(
  data      = Y,
  model     = "dynamic",
  family    = "gaussian",
  symmetric = TRUE,
  nscan     = 10000,
  burn      = 5000,
  odens     = 5,
  verbose   = FALSE
)
```

### Convergence

Always check convergence before interpreting results. Well-behaved
traces show stationary fluctuation around a stable mean. Under the
symmetric constraint, `tau_A2` and `tau_B2` should be identical.

``` r
check_convergence(fit)
plot_trace(fit)
```

## 4 Posterior alignment scores

The model estimates a time-varying $\Theta_{t}$ and a static baseline
$M$. Summing them gives smoothed alignment trajectories for each dyad,
filtering out observation noise to reveal the underlying relational
dynamics.

``` r
Theta_mean = apply(fit$Theta, 1:4, mean)
M_mean = apply(fit$M, 1:3, mean)

alignment = array(NA, dim = c(n, n, 1, length(fit$time_kept)))
for (t_idx in seq_along(fit$time_kept)) {
  alignment[, , 1, t_idx] = Theta_mean[, , 1, t_idx] + M_mean[, , 1]
}

key_dyads = list(
  "USA - UK"       = c("USA", "UK"),
  "USA - China"    = c("USA", "China"),
  "USA - Russia"   = c("USA", "Russia"),
  "China - Russia" = c("China", "Russia"),
  "India - China"  = c("India", "China"),
  "Brazil - USA"   = c("Brazil", "USA")
)

stored_years = years[fit$time_kept]
df_align = do.call(rbind, lapply(names(key_dyads), function(nm) {
  pair = key_dyads[[nm]]
  i = which(actor_labels == pair[1])
  j = which(actor_labels == pair[2])
  data.frame(dyad = nm, year = stored_years,
             alignment = alignment[i, j, 1, ])
}))

ggplot(df_align, aes(x = year, y = alignment, colour = dyad)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Posterior Mean Foreign Policy Alignment",
    subtitle = "G7 + BRICS, UNGA Voting Similarity (2005-2015)",
    x = "Year", y = "Estimated Alignment", colour = "Dyad"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    legend.position = "bottom"
  )
```

The model should recover patterns consistent with well-known
geopolitical structures. China-Russia alignment is expected to be
positive and stable, while USA-China and USA-Russia should remain
negative, reflecting systematic disagreement on UNGA votes. Western
allied dyads (USA-UK) should cluster at positive values.

## 5 Impulse response: edge shock

How does a bilateral shock propagate through the broader network? We
simulate a positive shock to USA-Russia alignment and track network
density. Because the data are symmetric, the shock matrix is symmetrized
to maintain consistency.

``` r
S_edge = build_shock(m = n, type = "unit_edge",
                     i = usa_idx, j = russia_idx, magnitude = 0.5)
S_edge = S_edge + t(S_edge)

irf_density = compute_irf(
  fit,
  shock    = S_edge,
  H        = 5,
  t0       = 3,
  stat_fun = dbn::stat_density
)

plot(irf_density, ci_level = 0.90,
     title = "USA-Russia Alignment Shock (+0.5) -> Network Density")
print(irf_density)
```

At horizon 0, the effect is exact: $2 \times 0.5/132 \approx 0.0076$. At
subsequent horizons, credible intervals widen reflecting posterior
uncertainty in $A_{t}$.

## 6 Impulse response: node shock

What if China shifts alignment toward all partners at once? We apply a
node-out shock and track both network density and China’s out-degree to
see whether the effect remains localized or propagates.

``` r
S_node = build_shock(m = n, type = "node_out",
                     i = china_idx, magnitude = 0.3)
diag(S_node) = 0
S_node = S_node + t(S_node)

irf_china_dens = compute_irf(
  fit, shock = S_node, H = 5, t0 = 3,
  stat_fun = dbn::stat_density
)

irf_china_outdeg = compute_irf(
  fit, shock = S_node, H = 5, t0 = 3,
  stat_fun = function(X) dbn::stat_out_degree(X)[china_idx]
)

df_irfs = rbind(
  data.frame(irf_china_dens, stat = "Network Density"),
  data.frame(irf_china_outdeg, stat = "China Out-Degree")
)

ggplot(df_irfs, aes(x = horizon, y = mean)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ stat, scales = "free_y") +
  labs(
    title = "Node Shock: China Alignment +0.3 (All Partners)",
    x = "Horizon (Years)", y = "IRF (Change From Baseline)"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text.x = element_text(color = "white", hjust = 0)
  )
```

## 7 Comparative IRF: G7 node shocks

Which country’s engagement shift produces the largest ripple? We apply
identical outward shocks from each G7 member and compare the resulting
density trajectories. All countries produce the same effect at horizon 0
by construction; differences at subsequent horizons reflect each
country’s structural position in the estimated dynamics.

``` r
g7 = c("USA", "Canada", "UK", "France", "Germany", "Italy", "Japan")

irf_list = lapply(setNames(g7, g7), function(country) {
  idx = which(actor_labels == country)
  S = build_shock(m = n, type = "node_out", i = idx, magnitude = 0.3)
  diag(S) = 0
  S = S + t(S)
  compute_irf(fit, shock = S, H = 5, t0 = 3,
              stat_fun = dbn::stat_density)
})

df_compare = do.call(rbind, lapply(names(irf_list), function(nm) {
  d = irf_list[[nm]]
  data.frame(horizon = d$horizon, mean = d$mean,
             q10 = d$q10, q90 = d$q90, country = nm)
}))

ggplot(df_compare, aes(x = horizon, y = mean, colour = country)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Comparative IRF: G7 Node Shocks",
    subtitle = "Same Magnitude; Differences Reflect Structural Position",
    x = "Horizon (Years)", y = "Change in Network Density",
    colour = "Shock Origin"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    legend.position = "bottom"
  )
```

Posterior uncertainty (omitted for visual clarity) is typically large
relative to between-country differences, so these comparisons should be
interpreted as suggestive rather than definitive without longer chains.

## 8 Reciprocity analysis

Nonlinear statistics like reciprocity interact with the baseline
network, producing posterior uncertainty even at horizon 0.

``` r
S_usachn = build_shock(m = n, type = "unit_edge",
                       i = usa_idx, j = china_idx, magnitude = 0.5)
S_usachn = S_usachn + t(S_usachn)

irf_recip = compute_irf(
  fit, shock = S_usachn, H = 5, t0 = 3,
  stat_fun = dbn::stat_reciprocity
)

plot(irf_recip, ci_level = 0.90,
     title = "USA-China Alignment Thaw (+0.5) -> Network Reciprocity")
```

## 9 Practical guidance

### Model settings

For publication-quality results, use longer chains to ensure stable
credible intervals:

``` r
fit = dbn(Y, model = "dynamic", family = "gaussian",
          symmetric = TRUE, nscan = 10000, burn = 5000, odens = 5)
```

### Runtime expectations

| Actors | Time Points | Approx. Time (8 Cores)  |
|--------|-------------|-------------------------|
| 10-15  | 10-20       | 5-15 minutes            |
| 15-25  | 10-20       | 15-60 minutes           |
| 25-50  | 10-20       | 1-4 hours               |
| 50+    | any         | Use `model = "lowrank"` |

### Interpreting IRFs

Horizon 0 captures the mechanical consequence of the shock. For linear
statistics (density, degree), the horizon-0 effect is exact; for
nonlinear statistics (reciprocity, transitivity), it is uncertain even
at horizon 0. The decay rate is governed by the spectral radius of $A$,
and the credible intervals reflect posterior uncertainty in model
parameters. IRFs are counterfactual scenarios that assume transition
dynamics are invariant to the shock.

### Symmetric networks

For undirected data, two choices must be made consistently. First, set
`symmetric = TRUE` in
[`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md) to
constrain $B_{t} = A_{t}$. Second, symmetrize any shock matrix after
construction: `S = S + t(S)`.

### Time-thinned models

When `time_thin > 1`, stored arrays correspond to a subset of the time
series. Map back to original years with `years[fit$time_kept]`.
