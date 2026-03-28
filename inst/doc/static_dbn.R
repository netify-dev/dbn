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


## ----simulate-----------------------------------------------------------------
sim = simulate_static_dbn(
  n      = 12,   # 12 actors
  p      = 1,    # 1 relation type
  time   = 15,   # 15 time periods
  sigma2 = 0.5,  # process noise variance
  tau2   = 0.1,  # prior variance for A and B
  seed   = 6886
)

Y = sim$Z
dim(Y)


## ----fit, cache = TRUE--------------------------------------------------------
fit = dbn(
  Y,
  model   = "static",
  family  = "gaussian",
  nscan   = 1000,
  burn    = 500,
  odens   = 2,
  verbose = FALSE
)


## ----convergence, fig.height = 5----------------------------------------------
check_convergence(fit)
plot_trace(fit)


## ----param-recovery-----------------------------------------------------------
param_summary(fit)


## ----m-recovery, fig.height = 4.5---------------------------------------------
M_est = latent_summary(fit, fun = mean)

# exclude diagonal entries from both truth and estimate
n_act = 12
off_diag = which(M_est$i != M_est$j)
est_M_vec = M_est$value[off_diag]

true_M_mat = sim$M[, , 1]
diag(true_M_mat) = NA
true_M_vec = as.vector(true_M_mat)
true_M_vec = true_M_vec[!is.na(true_M_vec)]

df_M = data.frame(truth = true_M_vec, estimate = est_M_vec)
r_sq = round(cor(df_M$truth, df_M$estimate)^2, 3)

ggplot(df_M, aes(x = truth, y = estimate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  labs(
    title = "Recovery of Baseline Mean M",
    subtitle = paste0("R-Squared = ", r_sq),
    x = "True M", y = "Posterior Mean M"
  ) +
  theme_bw() +
  theme(panel.border = element_blank())


## ----summary------------------------------------------------------------------
summary(fit)


## ----ppc, fig.height = 4------------------------------------------------------
ppd = posterior_predict_dbn(fit, ndraws = 20, seed = 6886)
plot_ppc_ecdf(fit, ppd, Y_obs = Y, rel = 1)


## ----other-families, cache = TRUE---------------------------------------------
# ordinal: for ranked data (e.g., event counts, Likert items)
fit_ord = dbn(sim$Y, model = "static", family = "ordinal",
              nscan = 1000, burn = 500, odens = 2, verbose = FALSE)
summary(fit_ord)

# binary: for 0/1 tie presence/absence
Y_binary = array(as.integer(sim$Y[, , , ] > 3), dim = dim(sim$Y))
for (t in seq_len(dim(Y_binary)[4])) diag(Y_binary[, , 1, t]) = NA
fit_bin = dbn(Y_binary, model = "static", family = "binary",
              nscan = 1000, burn = 500, odens = 2, verbose = FALSE)
summary(fit_bin)

