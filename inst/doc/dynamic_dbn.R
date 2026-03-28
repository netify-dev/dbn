## ----setup, include = FALSE---------------------------------------------------
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse  = TRUE, comment = "#>",
  fig.align = "center",
  fig.width = 6, fig.height = 4,
  message   = FALSE, warning = FALSE,
  eval      = NOT_CRAN
)
library(dbn)
library(ggplot2)


## ----simulate-----------------------------------------------------------------
sim = simulate_dynamic_dbn(
  n      = 8,
  p      = 1,
  time   = 15,
  sigma2 = 0.3,
  tauA2  = 0.04,
  tauB2  = 0.04,
  ar1    = TRUE,
  rhoA   = 0.9,
  rhoB   = 0.9,
  seed   = 6886
)

Y = sim$Z
dim(Y)


## ----fit, cache = TRUE--------------------------------------------------------
fit = dbn(
  Y,
  family     = "gaussian",
  model      = "dynamic",
  nscan      = 1000,
  burn       = 500,
  odens      = 2,
  ar1        = TRUE,
  update_rho = TRUE,
  verbose    = FALSE
)


## ----convergence--------------------------------------------------------------
check_convergence(fit)


## ----trace, fig.height = 5----------------------------------------------------
plot_trace(fit, pars = c("sigma2", "tau_A2", "tau_B2", "rho_A", "rho_B"))


## ----params-------------------------------------------------------------------
param_summary(fit)


## ----theta-recovery, fig.height = 4.5-----------------------------------------
Theta_mean = apply(fit$Theta, 1:4, mean)
M_mean = apply(fit$M, 1:3, mean)

n_act = dim(sim$Theta)[1]
Tt = dim(sim$Theta)[4]

true_vals = numeric(0)
est_vals  = numeric(0)

for (t in seq_len(Tt)) {
  true_mat = sim$Theta[, , 1, t]
  est_mat  = Theta_mean[, , 1, t] + M_mean[, , 1]
  mask = !is.na(true_mat)
  true_vals = c(true_vals, true_mat[mask])
  est_vals  = c(est_vals, est_mat[mask])
}

df_theta = data.frame(truth = true_vals, estimate = est_vals)
r_sq = round(cor(df_theta$truth, df_theta$estimate)^2, 3)

ggplot(df_theta, aes(x = truth, y = estimate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.1, size = 0.8) +
  labs(
    title    = "Recovery of Latent Network State",
    subtitle = paste0("R-Squared = ", r_sq, " (All Dyads x Time Points)"),
    x = expression("True " * Theta[t]),
    y = expression("Posterior Mean " * Theta[t] + M)
  ) +
  theme_bw() +
  theme(panel.border = element_blank())


## ----dyad-path, fig.height = 3.5----------------------------------------------
dyad_path(fit, i = 1, j = 2)


## ----forecast-----------------------------------------------------------------
fc = predict(fit, H = 5, S = 200, summary = "mean")
dim(fc)


## ----network-stats------------------------------------------------------------
ns = network_summary(fit, stat = "density")
head(ns)


## ----edge-prob----------------------------------------------------------------
ep = edge_prob(fit, rel = 1, time = dim(Y)[4])
round(ep[1:5, 1:5], 2)


## ----bipartite, cache = TRUE--------------------------------------------------
sim_bp = simulate_dynamic_dbn(
  n = 8, n_col = 6, p = 1, time = 12, seed = 6886
)

fit_bp = dbn(
  sim_bp$Z,
  model   = "dynamic",
  family  = "gaussian",
  nscan   = 1000,
  burn    = 500,
  odens   = 2,
  verbose = FALSE
)
summary(fit_bp)


## ----threads, eval = FALSE----------------------------------------------------
# set_dbn_threads(parallel::detectCores(logical = FALSE))


## ----memory-est---------------------------------------------------------------
estimate_memory(
  n_row = 50, n_col = 50, p = 1, Tt = 30,
  nscan = 1000, burn = 500, odens = 2,
  time_thin = 2, family = "gaussian"
)

