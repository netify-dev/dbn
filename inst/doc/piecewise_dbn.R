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
sim = simulate_piecewise_dbn(
  n      = 8,
  time   = 30,
  blocks = c(15, 30),
  p      = 1,
  sigma2 = 0.5,
  tau2   = 0.3,
  seed   = 6886
)

dim(sim$Y)

# block structure
sim$block_info$K
sim$block_info$boundaries
sim$block_info$lengths

# average element-wise difference between regime A matrices
round(mean(abs(sim$true_A[[1]] - sim$true_A[[2]])), 3)


## ----fit, cache = TRUE--------------------------------------------------------
fit = dbn(
  sim$Y,
  model   = "piecewise",
  family  = "ordinal",
  blocks  = c(15, 30),
  nscan   = 1000,
  burn    = 500,
  odens   = 2,
  verbose = FALSE
)

summary(fit)


## ----convergence, fig.height = 5----------------------------------------------
check_convergence(fit)
plot_trace(fit, pars = c("s2", "t2", "g2"))


## ----compare-blocks-----------------------------------------------------------
regime_diffs = compare_blocks(fit)


## ----regime-matrices----------------------------------------------------------
A1_est = fit$A_blocks[[1]]
A2_est = fit$A_blocks[[2]]


## ----influence-change, fig.height = 5-----------------------------------------
diff_A = A2_est - A1_est
n = nrow(diff_A)

df_heatmap = expand.grid(
  receiver = seq_len(n),
  sender   = seq_len(n)
)
df_heatmap$value = as.vector(diff_A)

ggplot(df_heatmap, aes(x = sender, y = receiver, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    name = expression(Delta * A)
  ) +
  labs(
    title = "Change in Sender Influence: Regime 2 vs. Regime 1",
    x = "Sender (Influenced By)", y = "Sender (Influences)"
  ) +
  coord_equal() +
  theme_bw() +
  theme(panel.border = element_blank())


## ----timing, cache = TRUE-----------------------------------------------------
sim_t = simulate_piecewise_dbn(n = 10, time = 30, blocks = 3, seed = 6886)

t_pw = system.time({
  fit_pw = dbn(sim_t$Y, model = "piecewise", blocks = 3,
               nscan = 1000, burn = 500, verbose = FALSE)
})

t_dyn = system.time({
  fit_dyn = dbn(sim_t$Y, model = "dynamic",
                nscan = 1000, burn = 500, verbose = FALSE)
})

data.frame(
  model   = c("Piecewise", "Dynamic"),
  seconds = round(c(t_pw["elapsed"], t_dyn["elapsed"]), 1)
)


## ----memory-compare-----------------------------------------------------------
fit_full = dbn(sim$Y, model = "piecewise", blocks = c(15, 30),
               nscan = 200, burn = 100, verbose = FALSE, store_theta = TRUE)
fit_lean = dbn(sim$Y, model = "piecewise", blocks = c(15, 30),
               nscan = 200, burn = 100, verbose = FALSE, store_theta = FALSE)

data.frame(
  storage = c("With Theta", "Without Theta"),
  mb      = round(c(object.size(fit_full), object.size(fit_lean)) / 1e6, 2)
)


## ----block-specs, eval = FALSE------------------------------------------------
# # integer: equal-sized regimes
# fit = dbn(Y, model = "piecewise", blocks = 4)
# 
# # vector of endpoints
# fit = dbn(Y, model = "piecewise", blocks = c(10, 25, 40))
# 
# # named vector for interpretability
# fit = dbn(Y, model = "piecewise",
#           blocks = c(pre_war = 15, war = 30, post_war = 45))


## ----gaussian, cache = TRUE---------------------------------------------------
fit_gauss = dbn(
  sim$Y_continuous,
  model   = "piecewise",
  family  = "gaussian",
  blocks  = c(15, 30),
  nscan   = 1000,
  burn    = 500,
  odens   = 2,
  verbose = FALSE
)

summary(fit_gauss)


## ----large-settings, eval = FALSE---------------------------------------------
# fit_large = dbn(
#   Y_large,
#   model       = "piecewise",
#   family      = "gaussian",
#   blocks      = c(25, 50),
#   store_theta = FALSE,
#   nscan       = 5000,
#   burn        = 2000,
#   odens       = 5,
#   verbose     = TRUE
# )

