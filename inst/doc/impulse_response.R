## ----setup, include = FALSE---------------------------------------------------
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.align = "center",
  fig.width = 6, fig.height = 4,
  message   = FALSE, warning = FALSE,
  eval      = NOT_CRAN
)
library(dbn)
library(ggplot2)


## ----fit, cache = TRUE--------------------------------------------------------
sim = simulate_dynamic_dbn(
  n = 10, p = 1, time = 15,
  sigma2 = 0.3, tauA2 = 0.03, tauB2 = 0.03,
  seed = 6886
)

fit = dbn(
  sim$Z,
  model   = "dynamic",
  family  = "gaussian",
  nscan   = 1000,
  burn    = 500,
  odens   = 2,
  verbose = FALSE
)


## ----check-conv---------------------------------------------------------------
check_convergence(fit)


## ----shocks-------------------------------------------------------------------
m = 10

# shock a single directed edge (actor 1 -> actor 2)
S_edge = build_shock(m, type = "unit_edge", i = 1, j = 2, magnitude = 2)

# shock all outgoing edges from actor 1
S_node = build_shock(m, type = "node_out", i = 1, magnitude = 1)

# uniform shock to all edges
S_dens = build_shock(m, type = "density", magnitude = 0.5)

# custom block shock (e.g., a group-level perturbation)
S_custom = matrix(0, m, m)
S_custom[1:3, 4:6] = 1.0


## ----irf-density, cache = TRUE, fig.height = 3.5------------------------------
irf_dens = compute_irf(
  fit,
  shock    = S_edge,
  H        = 5,
  t0       = 3,
  stat_fun = dbn::stat_density,
  n_draws  = 100
)

plot(irf_dens, title = "Edge Shock (1->2, +2) -> Network Density")
print(irf_dens)


## ----irf-degree, cache = TRUE, fig.height = 3.5-------------------------------
irf_deg = compute_irf(
  fit,
  shock    = S_node,
  H        = 5,
  t0       = 3,
  stat_fun = function(X) stat_out_degree(X)[1],
  n_draws  = 100
)

plot(irf_deg, title = "Node-Out Shock (Actor 1, +1) -> Actor 1 Out-Degree")


## ----irf-recip, cache = TRUE, fig.height = 3.5--------------------------------
irf_recip = compute_irf(
  fit,
  shock    = S_edge,
  H        = 5,
  t0       = 3,
  stat_fun = stat_reciprocity,
  n_draws  = 100
)

plot(irf_recip, title = "Edge Shock (1->2, +2) -> Network Reciprocity")


## ----compare-irfs, fig.width = 7, fig.height = 4------------------------------
irfs = list(
  "Edge -> Density"     = irf_dens,
  "Node -> Out-Degree"  = irf_deg,
  "Edge -> Recipr."     = irf_recip
)

df_all = do.call(rbind, lapply(names(irfs), function(nm) {
  d = irfs[[nm]]
  data.frame(
    horizon = d$horizon,
    mean    = d$mean,
    q025    = d$q025,
    q975    = d$q975,
    shock   = nm
  )
}))

ggplot(df_all, aes(x = horizon, y = mean)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ shock, scales = "free_y") +
  labs(x = "Horizon", y = "IRF", title = "Impulse Response Comparison") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text.x = element_text(color = "white", hjust = 0)
  )


## ----pub-settings, eval = FALSE-----------------------------------------------
# fit = dbn(Y, model = "dynamic", family = "gaussian",
#           nscan = 10000, burn = 5000, odens = 5)
# 
# irf = compute_irf(fit, shock = S, H = 10, t0 = 5,
#                   stat_fun = dbn::stat_density, n_draws = 500)

