# synthetic-data generator that exactly matches the model used by the
# covariate sampler:
#   R_t = A_t R_{t-1} B_t' + state_noise       (state, A_t = B_t = I)
#   Y_t = R_t + L_t + obs_noise                (observation)
.sim_covariate_dynamic <- function(n_row, n_col, Tt, X_list, beta_true,
                                    sigma2_state = 0.05, sigma2_obs = 0.5,
                                    seed = 1L) {
	set.seed(seed)
	nc <- n_row * n_col
	K <- length(beta_true)
	L_fixed <- matrix(0, n_row, n_col)
	for (k in seq_len(K)) L_fixed <- L_fixed + beta_true[k] * X_list[[k]]
	R_arr <- array(0, dim = c(n_row, n_col, Tt))
	R_arr[, , 1] <- matrix(stats::rnorm(nc, 0, sqrt(sigma2_state)), n_row, n_col)
	for (t in 2:Tt) {
		R_arr[, , t] <- R_arr[, , t - 1] +
			matrix(stats::rnorm(nc, 0, sqrt(sigma2_state)), n_row, n_col)
	}
	Y_arr <- array(0, dim = c(n_row, n_col, 1, Tt))
	for (t in seq_len(Tt)) {
		Y_arr[, , 1, t] <- R_arr[, , t] + L_fixed +
			matrix(stats::rnorm(nc, 0, sqrt(sigma2_obs)), n_row, n_col)
	}
	if (n_row == n_col) for (t in seq_len(Tt)) diag(Y_arr[, , 1, t]) <- NA
	list(Y = Y_arr, R = R_arr, L = L_fixed, beta_true = beta_true,
		sigma2_state = sigma2_state, sigma2_obs = sigma2_obs)
}

# joint gibbs recovery: shift-trick FFBS + observation-residual beta update
test_that("covariate engine recovers beta on data matching the model", {
	skip_on_cran()
	set.seed(42L)
	n_row <- 8L
	n_col <- 8L
	Tt    <- 30L
	X1 <- matrix(stats::rbinom(n_row * n_col, 1, 0.5), n_row, n_col)
	X2 <- matrix(stats::rnorm(n_row * n_col), n_row, n_col)
	beta_true <- c(0.7, -0.4)
	sim <- .sim_covariate_dynamic(n_row, n_col, Tt,
		X_list = list(X1, X2), beta_true = beta_true,
		sigma2_state = 0.05, sigma2_obs = 0.5, seed = 13L)
	cv <- dbn_covariates(
		dyad = list(x1 = X1, x2 = X2),
		standardize = FALSE,
		coef_by = "shared"
	)
	fit <- suppressMessages(suppressWarnings(
		dbn:::.dbn_with_covariates(
			sim$Y, cv,
			family = "gaussian",
			model = "dynamic",
			nscan = 2000L,
			burn = 1000L,
			odens = 2L,
			seed = 99L,
			verbose = FALSE
		)
	))
	beta_mean <- colMeans(fit$beta)
	beta_sd   <- apply(fit$beta, 2, stats::sd)
	beta_q025 <- apply(fit$beta, 2, stats::quantile, 0.025)
	beta_q975 <- apply(fit$beta, 2, stats::quantile, 0.975)
	# truth must lie in the 95% credible interval and within ~2 sd of mean
	expect_gt(beta_true[1], beta_q025[1])
	expect_lt(beta_true[1], beta_q975[1])
	expect_gt(beta_true[2], beta_q025[2])
	expect_lt(beta_true[2], beta_q975[2])
	expect_lt(abs(beta_mean[1] - beta_true[1]) / beta_sd[1], 2)
	expect_lt(abs(beta_mean[2] - beta_true[2]) / beta_sd[2], 2)
	# fit object has expected structure
	expect_s3_class(fit, "dbn")
	expect_s3_class(fit, "dbn_covariates_fit")
	expect_equal(dim(fit$Theta), c(n_row, n_col, 1L, Tt, 500L))
	expect_equal(length(fit$A), 500L)
})

# two-chain consistency: posterior means should be within a few SD of each other
test_that("two independent chains agree on beta posterior", {
	skip_on_cran()
	set.seed(2L)
	n_row <- 6L
	n_col <- 6L
	Tt    <- 25L
	X1 <- matrix(stats::rnorm(n_row * n_col), n_row, n_col)
	beta_true <- 0.6
	sim <- .sim_covariate_dynamic(n_row, n_col, Tt,
		X_list = list(X1), beta_true = beta_true,
		sigma2_state = 0.05, sigma2_obs = 0.5, seed = 71L)
	cv <- dbn_covariates(
		dyad = list(x1 = X1),
		standardize = FALSE,
		coef_by = "shared"
	)
	make_fit <- function(seed) {
		suppressMessages(suppressWarnings(
			dbn:::.dbn_with_covariates(
				sim$Y, cv,
				family = "gaussian",
				model = "dynamic",
				nscan = 1500L,
				burn = 750L,
				odens = 2L,
				seed = seed,
				verbose = FALSE
			)
		))
	}
	fit_a <- make_fit(13L)
	fit_b <- make_fit(17L)
	mu_a <- mean(fit_a$beta[, 1])
	mu_b <- mean(fit_b$beta[, 1])
	sd_a <- stats::sd(fit_a$beta[, 1])
	expect_lt(abs(mu_a - mu_b) / sd_a, 2)
})

# binary path: probit observation with augmented Z sampled via truncated normal
test_that("binary path runs and recovers sign of beta", {
	skip_on_cran()
	set.seed(99L)
	n_row <- 8L
	n_col <- 8L
	Tt    <- 25L
	nc    <- n_row * n_col
	# simulate latent gaussian theta + L; threshold at 0 for binary
	X1 <- matrix(stats::rnorm(n_row * n_col), n_row, n_col)
	beta_true <- 1.2
	L_fixed <- X1 * beta_true
	sigma2_state <- 0.05
	R_arr <- array(0, dim = c(n_row, n_col, Tt))
	R_arr[, , 1] <- matrix(stats::rnorm(nc, 0, sqrt(sigma2_state)), n_row, n_col)
	for (t in 2:Tt) {
		R_arr[, , t] <- R_arr[, , t - 1] +
			matrix(stats::rnorm(nc, 0, sqrt(sigma2_state)), n_row, n_col)
	}
	# observation: probit Y = 1 if (R + L + N(0,1)) > 0
	Y_arr <- array(0L, dim = c(n_row, n_col, 1, Tt))
	for (t in seq_len(Tt)) {
		Z_t <- R_arr[, , t] + L_fixed +
			matrix(stats::rnorm(nc, 0, 1), n_row, n_col)
		Y_arr[, , 1, t] <- as.integer(Z_t > 0)
	}
	for (t in seq_len(Tt)) diag(Y_arr[, , 1, t]) <- NA
	cv <- dbn_covariates(dyad = list(x1 = X1), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn:::.dbn_with_covariates(
			Y_arr, cv,
			family = "binary",
			nscan = 1200L,
			burn = 600L,
			odens = 2L,
			seed = 7L,
			verbose = FALSE
		)
	))
	beta_mean <- mean(fit$beta[, 1])
	# binary identifiability requires fixed observation variance (sigma2_obs = 1)
	# so the beta scale is identified; sign should match
	expect_gt(beta_mean, 0)
	expect_s3_class(fit, "dbn")
	expect_equal(dim(fit$beta), c(300L, 1L))
})

# public dbn() entry routes through .dbn_with_covariates when covariates set
test_that("dbn() public entry dispatches to covariate sampler", {
	skip_on_cran()
	set.seed(7L)
	n_row <- 6L
	Tt    <- 20L
	X <- matrix(stats::rnorm(n_row * n_row), n_row, n_row)
	sim <- .sim_covariate_dynamic(n_row, n_row, Tt,
		X_list = list(X), beta_true = 0.6,
		sigma2_state = 0.05, sigma2_obs = 0.5, seed = 11L)
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn(sim$Y, model = "dynamic", family = "gaussian",
			covariates = cv,
			nscan = 600L, burn = 300L, odens = 2L,
			seed = 17L, verbose = FALSE)
	))
	expect_s3_class(fit, "dbn_covariates_fit")
	expect_true(!is.null(fit$beta))
	expect_equal(dim(fit$beta)[2], 1L)
})

# covariates now accepted on static / piecewise / hmm / lowrank via the
# block-Gibbs shift-trick driver; this test confirms the static path runs
test_that("dbn() accepts covariates on static models via shift-trick driver", {
	skip_on_cran()
	n_row <- 5L; Tt <- 6L
	X <- matrix(stats::rnorm(n_row * n_row), n_row, n_row)
	Y <- array(0, c(n_row, n_row, 1, Tt))
	for (t in seq_len(Tt))
		Y[, , 1, t] <- 0.5 * X + matrix(stats::rnorm(n_row * n_row, 0, 0.5),
			n_row, n_row)
	for (t in seq_len(Tt)) diag(Y[, , 1, t]) <- NA
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn(Y, model = "static", family = "gaussian", covariates = cv,
			nscan = 60L, burn = 20L, verbose = FALSE)
	))
	expect_s3_class(fit, "dbn_covariates_fit")
	expect_true(!is.null(fit$beta))
	expect_equal(ncol(fit$beta), 1L)
})

# accessor methods: coef / tidy / summary / vcov / confint
test_that("accessors expose beta posterior on dbn_covariates_fit", {
	skip_on_cran()
	set.seed(5L)
	n_row <- 5L; Tt <- 18L
	X1 <- matrix(stats::rnorm(n_row * n_row), n_row, n_row)
	X2 <- matrix(stats::rbinom(n_row * n_row, 1, 0.5), n_row, n_row)
	sim <- .sim_covariate_dynamic(n_row, n_row, Tt,
		X_list = list(X1, X2), beta_true = c(0.6, -0.3),
		sigma2_state = 0.05, sigma2_obs = 0.4, seed = 41L)
	cv <- dbn_covariates(dyad = list(x1 = X1, x2 = X2), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn(sim$Y, model = "dynamic", family = "gaussian",
			covariates = cv,
			nscan = 600L, burn = 300L, odens = 2L,
			seed = 19L, verbose = FALSE)
	))
	# coef(fit, "beta") returns named vector
	cb <- coef(fit, what = "beta")
	expect_equal(length(cb), 2L)
	expect_equal(names(cb), c("x1", "x2"))
	# coef(fit, "all") returns named list with beta + A + B + M
	ca <- coef(fit, what = "all")
	expect_true("beta" %in% names(ca))
	expect_true("A" %in% names(ca))
	# tidy returns long data frame with parameter == "beta" rows
	td <- tidy(fit)
	beta_rows <- td[td$parameter == "beta", ]
	expect_equal(nrow(beta_rows), 2L)
	expect_true(all(c("estimate", "std.error", "conf.low", "conf.high") %in% names(beta_rows)))
	expect_equal(sort(beta_rows$term), c("x1", "x2"))
	# vcov returns K x K cov of beta draws
	V <- vcov(fit)
	expect_equal(dim(V), c(2L, 2L))
	expect_equal(rownames(V), c("x1", "x2"))
	# confint(what = "beta") returns K x 2 quantile-based intervals
	ci_beta <- confint(fit, level = 0.95, what = "beta")
	expect_equal(dim(ci_beta), c(2L, 2L))
	expect_equal(rownames(ci_beta), c("x1", "x2"))
	# confint(what = "all") (default) returns a list with $beta appended
	ci_all <- confint(fit, level = 0.95)
	expect_true(is.list(ci_all))
	expect_true("beta" %in% names(ci_all))
	# summary appends a beta posterior table
	smry <- summary(fit)
	expect_true("beta" %in% names(smry))
	expect_equal(nrow(smry$beta), 2L)
})

# H-step covariate-aware forecast
test_that("forecast_with_covariates produces correct shape", {
	skip_on_cran()
	set.seed(99L)
	n_row <- 5L; Tt <- 18L
	X <- matrix(stats::rnorm(n_row * n_row), n_row, n_row)
	sim <- .sim_covariate_dynamic(n_row, n_row, Tt,
		X_list = list(X), beta_true = 0.6,
		sigma2_state = 0.05, sigma2_obs = 0.5, seed = 23L)
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn(sim$Y, model = "dynamic", family = "gaussian",
			covariates = cv,
			nscan = 400L, burn = 200L, odens = 2L,
			seed = 31L, verbose = FALSE)
	))
	# constant design over the H future steps
	cv_future <- dbn_covariates(
		dyad = list(x = array(X, dim = c(n_row, n_row, 5L))),
		standardize = FALSE
	)
	# none-summary returns full draws
	pred_draws <- forecast_with_covariates(fit, H = 5L,
		new_covariates = cv_future, draws = 40L, seed = 1L)
	expect_equal(dim(pred_draws), c(n_row, n_row, 1L, 5L, 40L))
	# mean-summary collapses draws dimension
	pred_mean <- forecast_with_covariates(fit, H = 5L,
		new_covariates = cv_future, draws = 40L,
		summary = "mean", seed = 1L)
	expect_equal(dim(pred_mean), c(n_row, n_row, 1L, 5L))
	# ci-summary returns named list
	pred_ci <- forecast_with_covariates(fit, H = 5L,
		new_covariates = cv_future, draws = 40L,
		summary = "ci", level = 0.9, seed = 1L)
	expect_true(all(c("mean", "lower", "upper", "level") %in% names(pred_ci)))
	expect_equal(pred_ci$level, 0.9)
	expect_equal(dim(pred_ci$mean), c(n_row, n_row, 1L, 5L))
	# rejects mismatched terms
	cv_wrong <- dbn_covariates(
		dyad = list(other = X),
		standardize = FALSE
	)
	expect_error(
		forecast_with_covariates(fit, H = 3L, new_covariates = cv_wrong, draws = 10L),
		"terms must match"
	)
})

# synthetic generator with random actor (sender + receiver) effects
.sim_actor_effects <- function(n, Tt, beta_true, u_row_true, u_col_true,
                                 X, sigma2_state = 0.05, sigma2_obs = 0.5,
                                 seed = 1L) {
	set.seed(seed)
	L_fixed <- X * beta_true
	u_row_mat <- matrix(u_row_true, n, n, byrow = FALSE)
	u_col_mat <- matrix(u_col_true, n, n, byrow = TRUE)
	R <- array(0, c(n, n, Tt))
	R[, , 1] <- matrix(stats::rnorm(n * n, 0, sqrt(sigma2_state)), n, n)
	for (t in 2:Tt) {
		R[, , t] <- R[, , t - 1] +
			matrix(stats::rnorm(n * n, 0, sqrt(sigma2_state)), n, n)
	}
	Y <- array(0, c(n, n, 1, Tt))
	for (t in seq_len(Tt)) {
		Y[, , 1, t] <- R[, , t] + L_fixed + u_row_mat + u_col_mat +
			matrix(stats::rnorm(n * n, 0, sqrt(sigma2_obs)), n, n)
	}
	for (t in seq_len(Tt)) diag(Y[, , 1, t]) <- NA
	list(Y = Y, u_row_true = u_row_true, u_col_true = u_col_true)
}

test_that("actor_effects = 'both' recovers centered u_row and u_col", {
	skip_on_cran()
	set.seed(123L)
	n  <- 8L
	Tt <- 30L
	u_row_true <- stats::rnorm(n, 0, 0.5)
	u_row_true <- u_row_true - mean(u_row_true)  # sum to zero
	u_col_true <- stats::rnorm(n, 0, 0.5)
	u_col_true <- u_col_true - mean(u_col_true)
	X <- matrix(stats::rnorm(n * n), n, n)
	beta_true <- 0.5
	sim <- .sim_actor_effects(n, Tt, beta_true, u_row_true, u_col_true, X,
		sigma2_state = 0.05, sigma2_obs = 0.4, seed = 41L)
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn:::.dbn_with_covariates(
			sim$Y, cv,
			family = "gaussian",
			actor_effects = "both",
			nscan = 1500L, burn = 750L, odens = 3L,
			seed = 17L, verbose = FALSE
		)
	))
	expect_true(!is.null(fit$u_row))
	expect_equal(dim(fit$u_row)[2], n)
	expect_true(!is.null(fit$u_col))
	expect_true(!is.null(fit$sigma_u_row2))
	# posterior means within ~1.5 sd of truth for each actor
	u_row_mean <- colMeans(fit$u_row)
	u_col_mean <- colMeans(fit$u_col)
	u_row_sd <- apply(fit$u_row, 2, stats::sd)
	u_col_sd <- apply(fit$u_col, 2, stats::sd)
	# at least 60% of actors should have truth within 2 sd of posterior mean
	hits_row <- mean(abs(u_row_mean - u_row_true) < 2 * u_row_sd)
	hits_col <- mean(abs(u_col_mean - u_col_true) < 2 * u_col_sd)
	expect_gt(hits_row, 0.6)
	expect_gt(hits_col, 0.6)
	# sum-to-zero invariant should be maintained per draw (after centering)
	row_means_per_draw <- rowMeans(fit$u_row)
	col_means_per_draw <- rowMeans(fit$u_col)
	expect_true(all(abs(row_means_per_draw) < 1e-8))
	expect_true(all(abs(col_means_per_draw) < 1e-8))
})

# synthetic generator with time-varying beta_t (random walk on a scalar
# coefficient applied to a time-invariant dyadic covariate X)
.sim_tv_beta <- function(n, Tt, X, beta_path,
                          sigma2_state = 0.05, sigma2_obs = 0.4, seed = 1L) {
	set.seed(seed)
	R <- array(0, c(n, n, Tt))
	R[, , 1] <- matrix(stats::rnorm(n * n, 0, sqrt(sigma2_state)), n, n)
	for (t in 2:Tt) {
		R[, , t] <- R[, , t - 1] +
			matrix(stats::rnorm(n * n, 0, sqrt(sigma2_state)), n, n)
	}
	Y <- array(0, c(n, n, 1, Tt))
	for (t in seq_len(Tt)) {
		Y[, , 1, t] <- R[, , t] + X * beta_path[t] +
			matrix(stats::rnorm(n * n, 0, sqrt(sigma2_obs)), n, n)
	}
	for (t in seq_len(Tt)) diag(Y[, , 1, t]) <- NA
	list(Y = Y, beta_path = beta_path)
}

test_that("accessors handle 3D beta from time_varying_beta = TRUE", {
	skip_on_cran()
	set.seed(13L)
	n <- 5L; Tt <- 14L
	X <- matrix(stats::rnorm(n * n), n, n)
	beta_path <- seq(0.3, 0.8, length.out = Tt)
	sim <- .sim_tv_beta(n, Tt, X, beta_path,
		sigma2_state = 0.05, sigma2_obs = 0.3, seed = 23L)
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn:::.dbn_with_covariates(
			sim$Y, cv,
			family = "gaussian",
			time_varying_beta = TRUE,
			nscan = 300L, burn = 150L, odens = 2L,
			seed = 9L, verbose = FALSE
		)
	))
	# coef returns K x Tt matrix for time-varying fits
	cb <- coef(fit, what = "beta")
	expect_equal(dim(cb), c(1L, Tt))
	# tidy emits K * Tt beta rows
	td <- tidy(fit)
	beta_rows <- td[td$parameter == "beta", ]
	expect_equal(nrow(beta_rows), Tt)
	# vcov returns K x K averaged covariance
	V <- vcov(fit)
	expect_equal(dim(V), c(1L, 1L))
	# confint(what="beta") returns K x Tt x 2
	ci <- confint(fit, what = "beta")
	expect_equal(length(dim(ci)), 3L)
	expect_equal(dim(ci)[2], Tt)
	expect_equal(dim(ci)[3], 2L)
})

test_that("time_varying_beta recovers a smooth beta_t path", {
	skip_on_cran()
	set.seed(99L)
	n <- 6L; Tt <- 30L
	X <- matrix(stats::rnorm(n * n), n, n)
	# beta drifts smoothly from 0.2 to 1.0
	beta_path <- seq(0.2, 1.0, length.out = Tt)
	sim <- .sim_tv_beta(n, Tt, X, beta_path,
		sigma2_state = 0.05, sigma2_obs = 0.3, seed = 31L)
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn:::.dbn_with_covariates(
			sim$Y, cv,
			family = "gaussian",
			time_varying_beta = TRUE,
			tau_beta2_init = 0.05,
			nscan = 1000L, burn = 500L, odens = 2L,
			seed = 7L, verbose = FALSE
		)
	))
	# beta_store has shape n_keep x K x Tt
	expect_equal(length(dim(fit$beta)), 3L)
	expect_equal(dim(fit$beta)[3], Tt)
	# posterior mean per t
	beta_mean <- apply(fit$beta, c(2, 3), mean)  # K x Tt
	# the recovered path should correlate strongly with the truth
	corr <- stats::cor(beta_mean[1, ], beta_path)
	expect_gt(corr, 0.7)
	# innovation variance positive
	expect_true(all(fit$tau_beta2 > 0))
})

test_that("actor_effects = 'none' leaves u_row / u_col NULL", {
	skip_on_cran()
	set.seed(2L)
	n <- 5L; Tt <- 12L
	X <- matrix(stats::rnorm(n * n), n, n)
	sim <- .sim_covariate_dynamic(n, n, Tt, list(X), 0.5,
		sigma2_state = 0.05, sigma2_obs = 0.4, seed = 7L)
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn:::.dbn_with_covariates(
			sim$Y, cv,
			family = "gaussian",
			actor_effects = "none",
			nscan = 200L, burn = 100L, odens = 2L,
			seed = 5L, verbose = FALSE
		)
	))
	expect_null(fit$u_row)
	expect_null(fit$u_col)
})

# variance components remain finite and positive
test_that("variance components stay positive and finite", {
	skip_on_cran()
	set.seed(31L)
	n_row <- 5L
	Tt    <- 15L
	X <- matrix(stats::rnorm(n_row * n_row), n_row, n_row)
	sim <- .sim_covariate_dynamic(n_row, n_row, Tt,
		X_list = list(X), beta_true = 0.5,
		sigma2_state = 0.05, sigma2_obs = 0.5, seed = 23L)
	cv <- dbn_covariates(dyad = list(x = X), standardize = FALSE)
	fit <- suppressMessages(suppressWarnings(
		dbn:::.dbn_with_covariates(
			sim$Y, cv,
			family = "gaussian",
			nscan = 800L,
			burn = 400L,
			odens = 2L,
			seed = 41L,
			verbose = FALSE
		)
	))
	expect_true(all(fit$tau_A2 > 0))
	expect_true(all(fit$tau_B2 > 0))
	expect_true(all(is.finite(fit$sigma2)))
	expect_true(all(is.finite(fit$sigma2_obs)))
	expect_true(all(fit$sigma2 > 0))
	expect_true(all(fit$sigma2_obs > 0))
})
