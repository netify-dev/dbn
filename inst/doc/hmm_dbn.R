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


## ----simulate, cache = TRUE---------------------------------------------------
sim = simulate_hmm_dbn(
  n     = 12,
  p     = 1,
  time  = 30,
  R     = 2,
  sigma2        = 0.3,
  tau_A2        = 0.8,
  tau_B2        = 0.8,
  transition_prob = 0.9,
  seed  = 6886
)

dim(sim$Y)

# true regime sequence and frequencies
sim$S
table(sim$S)


## ----fit, cache = TRUE--------------------------------------------------------
fit = dbn(
  sim$Y,
  model   = "hmm",
  family  = "ordinal",
  R       = 2,
  nscan   = 1000,
  burn    = 500,
  odens   = 2,
  verbose = FALSE
)


## ----summary------------------------------------------------------------------
summary(fit)


## ----regime-probs, fig.height = 3.5-------------------------------------------
plot_regime_probs(fit)


## ----regime-table-------------------------------------------------------------
probs = regime_probs(fit)
head(probs, 10)


## ----plot-hmm, fig.height = 8-------------------------------------------------
plot(fit)


## ----dyad, cache = TRUE, fig.height = 3.5-------------------------------------
dyad_path(fit, i = 2, j = 5)


## ----forecast-----------------------------------------------------------------
pred = predict(fit, H = 5, draws = 100, summary = "mean")
dim(pred)

