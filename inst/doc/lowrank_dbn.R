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


## ----memory-est---------------------------------------------------------------
estimate_memory(
  n_row = 50, n_col = 50, p = 1, Tt = 30,
  nscan = 5000, burn = 2000, odens = 5,
  family = "ordinal"
)


## ----simulate, cache = TRUE---------------------------------------------------
sim = simulate_lowrank_dbn(
  n          = 10,
  p          = 1,
  time       = 15,
  r          = 2,
  sigma2     = 0.3,
  tau_alpha2 = 0.1,
  tauB2      = 0.04,
  seed       = 6886
)

dim(sim$Z)


## ----fit, cache = TRUE--------------------------------------------------------
fit = dbn(
  sim$Z,
  model   = "lowrank",
  family  = "gaussian",
  r       = 2,
  nscan   = 1000,
  burn    = 500,
  odens   = 2,
  verbose = FALSE
)


## ----plot-model, fig.height = 8-----------------------------------------------
plot(fit)


## ----factor-traj, fig.height = 3.5--------------------------------------------
alpha_df = tidy_dbn_lowrank(fit)

ggplot(alpha_df, aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  facet_wrap(
    ~ factor, scales = "free_y",
    labeller = label_bquote(alpha[.(factor)])
  ) +
  labs(
    title = "Estimated Factor Trajectories",
    x = "Time", y = expression(alpha)
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text.x = element_text(color = "white", hjust = 0)
  )


## ----summary------------------------------------------------------------------
summary(fit)


## ----forecast-----------------------------------------------------------------
pred = predict(fit, H = 3, draws = 100, summary = "mean")
dim(pred)

