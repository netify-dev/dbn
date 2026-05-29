####
# dbn samplers with exogenous covariates
#
# implements the unified covariate engine documented in covariates.R for
# the dynamic asymmetric path. the math:
#   theta_t = L_t + A_t (theta_{t-1} - L_{t-1}) B_t' + eps_t
#   L_t     = D_t %*% beta
# the shift trick: pass Z' = Z - L to the existing samplers with M = 0
# and define R_t = theta_t - L_t. then:
#   R_t  = A_t R_{t-1} B_t' + eps_t           (standard state eq)
#   Z'_t = R_t + obs_noise                    (standard observation eq)
# so the existing batch_ffbs_all_relations and update_AB_batch_extended
# code paths apply unchanged; beta is sampled in a separate conjugate
# gaussian step inside the gibbs loop.
####

####
#' Fit a dbn with exogenous covariates (dynamic, asymmetric)
#'
#' @description Internal entry point. Called by \code{\link{dbn}} when
#'   \code{covariates} is non-NULL. Users should not call this directly.
#' @keywords internal
#' @noRd
.dbn_with_covariates <- function(Y, covariates,
                                  family = "gaussian",
                                  model = "dynamic",
                                  nscan = 2000L,
                                  burn = 500L,
                                  odens = 1L,
                                  symmetric = FALSE,
                                  prior_beta_scale = 2.5,
                                  prior_beta_mean = NULL,
                                  actor_effects = c("none", "row", "col", "both"),
                                  time_varying_beta = FALSE,
                                  tau_beta2_init = 0.05,
                                  a_beta = 0.5,
                                  b_beta = 0.5,
                                  ar1 = FALSE,
                                  rhoA = 0,
                                  rhoB = 0,
                                  prior_kind = c("rw", "iid"),
                                  init = NULL,
                                  seed = NULL,
                                  verbose = TRUE,
                                  shrink_rho = 0.9,
                                  ...) {
	actor_effects <- match.arg(actor_effects)
	# covariates can be NULL when the user only wants actor effects
	if (!is.null(covariates) && !inherits(covariates, "dbn_covariates"))
		cli::cli_abort("{.arg covariates} must be a {.cls dbn_covariates} object.")
	if (is.null(covariates) && identical(actor_effects, "none"))
		cli::cli_abort("Either {.arg covariates} or {.arg actor_effects} must be set.")
	if (!identical(model, "dynamic"))
		cli::cli_abort("Covariate / actor-effects support currently requires {.code model = \"dynamic\"}.")
	if (isTRUE(symmetric))
		cli::cli_abort("Covariate / actor-effects support for {.code symmetric = TRUE} is on the roadmap.")
	family <- match.arg(family, c("gaussian", "ordinal", "binary"))
	prior_kind <- match.arg(prior_kind)

	dY <- dim(Y)
	if (length(dY) != 4L)
		cli::cli_abort("{.arg Y} must be a 4D array {.code [n_row, n_col, p, T]}.")
	n_row <- dY[1L]
	n_col <- dY[2L]
	p     <- dY[3L]
	Tt    <- dY[4L]
	nc    <- n_row * n_col
	# contractivity ridge precision on A_t, B_t, as in dbn_dynamic. per-entry
	# prior variance shrink_rho^2 / n; null shrink_rho = off.
	if (!is.null(shrink_rho) && (length(shrink_rho) != 1L || !is.numeric(shrink_rho) ||
		!is.finite(shrink_rho) || shrink_rho <= 0 || shrink_rho > 1))
		cli::cli_abort("{.arg shrink_rho} must be a scalar in (0, 1] (or NULL).")
	kappaA_inv <- if (is.null(shrink_rho)) 0 else n_row / (shrink_rho^2)
	kappaB_inv <- if (is.null(shrink_rho)) 0 else n_col / (shrink_rho^2)
	if (p != 1L)
		cli::cli_abort(c(
			"Covariate / actor-effects support currently requires {.code p = 1} (single relation).",
			"i" = "Multi-relation covariates are on the roadmap."
		))

	if (is.null(covariates)) {
		# actor-effects-only path: build a placeholder covars with K = 0
		covars <- list(
			version = 1L, entries = list(),
			terms = data.frame(name = character(0), component = character(0),
				time_varying = logical(0), mean = numeric(0), sd = numeric(0),
				stringsAsFactors = FALSE),
			standardize = FALSE, coef_by = "shared",
			validated = TRUE,
			n_row = n_row, n_col = n_col, p = p, Tt = Tt, K = 0L
		)
		class(covars) <- c("dbn_covariates", "list")
	} else {
		covars <- .dbn_validate_covariates(covariates, n_row, n_col, p, Tt,
			symmetric = symmetric)
	}
	K <- covars$K

	if (K > 0L) {
		if (is.null(prior_beta_mean)) prior_beta_mean <- rep(0, K)
		if (length(prior_beta_mean) != K)
			cli::cli_abort("{.arg prior_beta_mean} must have length K = {K}.")
		if (length(prior_beta_scale) == 1L) prior_beta_scale <- rep(prior_beta_scale, K)
		if (length(prior_beta_scale) != K)
			cli::cli_abort("{.arg prior_beta_scale} must be scalar or length K = {K}.")
		prior_prec_vec <- 1 / prior_beta_scale ^ 2
	} else {
		prior_beta_mean <- numeric(0)
		prior_beta_scale <- numeric(0)
		prior_prec_vec <- numeric(0)
	}

	if (!is.null(seed)) set.seed(seed)

	# delegate Y → Z conversion and NA handling to the standard preprocessor.
	# for gaussian: Z = Y with non-finite cells set to 0 (matches the
	# observed-zero convention used by batch_ffbs_all_relations and
	# update_variances_dynamic in the main dbn() path)
	pre <- shared_preprocess(Y, family = family)
	Z   <- pre$Z
	Omega <- pre$Omega

	if (K > 0L) {
		# initial OLS estimate of beta on the preprocessed Z (post-NA-fill).
		# this matches the data passed to FFBS so the residual frame is
		# consistent from iteration 0 onward
		HtH_init <- matrix(0, K, K)
		Htc_init <- numeric(K)
		off_diag_mask <- if (n_row == n_col) as.numeric(diag(n_row) == 0) else
			rep(1, nc)
		for (t in seq_len(Tt)) {
			D_t <- .dbn_build_design_t(covars, t)
			z_t <- as.numeric(Z[, , 1, t])
			keep <- (off_diag_mask > 0)
			HtH_init <- HtH_init + crossprod(D_t[keep, , drop = FALSE])
			Htc_init <- Htc_init + crossprod(D_t[keep, , drop = FALSE], z_t[keep])
		}
		beta_init <- as.numeric(solve(diag(prior_prec_vec, K) + HtH_init,
			prior_prec_vec * prior_beta_mean + Htc_init))
	} else {
		beta_init <- numeric(0)
	}

	# beta storage: K-vector when static, K x Tt matrix when time-varying
	# (one column per t; warm-start each column at the static OLS init)
	if (time_varying_beta) {
		beta_mat <- matrix(beta_init, nrow = K, ncol = Tt)
		tau_beta2 <- rep(tau_beta2_init, K)
	} else {
		beta_mat <- matrix(beta_init, nrow = K, ncol = 1L)
	}

	# build L given beta vector / matrix
	# static: L_t = D_t %*% beta (column 1 of beta_mat repeated)
	# tv:     L_t = D_t %*% beta_mat[, t]
	# K == 0: actor-effects-only path; L from beta is identically zero
	build_L_local <- function(beta_mat) {
		L <- array(0, dim = c(n_row, n_col, Tt))
		if (K == 0L) return(L)
		for (t in seq_len(Tt)) {
			b_t <- if (ncol(beta_mat) == 1L) beta_mat[, 1L] else beta_mat[, t]
			L[, , t] <- .dbn_build_linear_predictor_t(covars, b_t, t)
		}
		L
	}
	build_L <- build_L_local
	L <- build_L(beta_mat)

	# initial residual state R = Theta - L; warm-start at Z - L
	R <- array(0, dim = c(n_row, n_col, Tt))
	for (t in seq_len(Tt)) R[, , t] <- Z[, , 1, t] - L[, , t]

	# operator + variance state. honors `init` like the no-covariate path:
	# init = "smart" (default if NULL) runs a brief ALS on the residual R
	# to seed A_t / B_t / sigma2 near the posterior mode; init = "default"
	# falls back to the identity. ALS init blocked on symmetric=TRUE since
	# the covariate path errors on symmetric anyway.
	init_mode <- if (is.character(init) && length(init) == 1L) init else "smart"
	smart_pkg <- NULL
	if (identical(init_mode, "smart") || is.null(init)) {
		R_for_init <- array(R, dim = c(n_row, n_col, 1L, Tt))
		smart_pkg <- tryCatch(
			.dbn_smart_init(R_for_init, family = "gaussian", model = "dynamic",
				symmetric = FALSE, Tt = Tt,
				n_row = n_row, n_col = n_col, p = 1L,
				verbose = FALSE),
			error = function(e) NULL
		)
	}
	if (!is.null(smart_pkg)) {
		Aarray <- smart_pkg$A_init
		Barray <- smart_pkg$B_init
		sigma2     <- smart_pkg$sigma2_init     %||% 1
		sigma2_obs <- smart_pkg$sigma2_obs_init %||% 1
		tauA2      <- smart_pkg$tauA2_init      %||% 1
		tauB2      <- smart_pkg$tauB2_init      %||% 1
	} else {
		Aarray <- array(0, dim = c(n_row, n_row, Tt))
		Barray <- array(0, dim = c(n_col, n_col, Tt))
		for (t in seq_len(Tt)) {
			Aarray[, , t] <- diag(n_row)
			Barray[, , t] <- diag(n_col)
		}
		sigma2     <- 1
		sigma2_obs <- 1
		tauA2 <- 1
		tauB2 <- 1
	}

	# weakly-informative IG(0.5, 0.5) priors for sigma2 / tauA2 / tauB2
	a_sig <- 0.5
	b_sig <- 0.5
	a_tau <- 0.5
	b_tau <- 0.5

	# random actor effects: latent sender (row) and receiver (col) effects
	# with sum-to-zero centering enforced post-sample. priors:
	#   u_row[i] ~ N(0, sigma_u_row2), sigma_u_row2 ~ IG(0.5, 0.5)
	use_u_row <- actor_effects %in% c("row", "both")
	use_u_col <- actor_effects %in% c("col", "both")
	u_row  <- numeric(n_row)
	u_col  <- numeric(n_col)
	sigma_u_row2 <- 1
	sigma_u_col2 <- 1
	a_u <- 0.5
	b_u <- 0.5

	# total linear shift L_total = L + u_row 1' + 1 u_col' (time-invariant
	# actor effects are broadcast across t)
	build_L_total <- function(L, u_row, u_col) {
		L_total <- L
		if (use_u_row) {
			for (t in seq_len(dim(L)[3])) L_total[, , t] <- L_total[, , t] +
				matrix(u_row, n_row, n_col, byrow = FALSE)
		}
		if (use_u_col) {
			for (t in seq_len(dim(L)[3])) L_total[, , t] <- L_total[, , t] +
				matrix(u_col, n_row, n_col, byrow = TRUE)
		}
		L_total
	}
	L_total <- build_L_total(L, u_row, u_col)

	# storage
	n_keep <- floor((nscan - burn) / odens)
	if (n_keep <= 0L)
		cli::cli_abort("nscan - burn must yield at least one stored draw.")
	# beta_store: n_keep x K (static) or n_keep x K x Tt (time-varying);
	# NULL when K == 0 (actor-effects-only path)
	if (K == 0L) {
		beta_store <- NULL
		tau_beta2_store <- NULL
	} else if (time_varying_beta) {
		beta_store <- array(NA_real_, dim = c(n_keep, K, Tt))
		dimnames(beta_store) <- list(NULL, covars$terms$name, NULL)
		tau_beta2_store <- matrix(NA_real_, n_keep, K)
		colnames(tau_beta2_store) <- covars$terms$name
	} else {
		beta_store <- matrix(NA_real_, n_keep, K)
		colnames(beta_store) <- covars$terms$name
		tau_beta2_store <- NULL
	}
	sigma2_store     <- numeric(n_keep)
	sigma2_obs_store <- numeric(n_keep)
	tauA2_store      <- numeric(n_keep)
	tauB2_store      <- numeric(n_keep)
	A_store     <- vector("list", n_keep)
	B_store     <- vector("list", n_keep)
	Theta_store <- array(NA_real_, dim = c(n_row, n_col, p, Tt, n_keep))
	L_store     <- array(NA_real_, dim = c(n_row, n_col, Tt, n_keep))
	u_row_store <- if (use_u_row) matrix(NA_real_, n_keep, n_row) else NULL
	u_col_store <- if (use_u_col) matrix(NA_real_, n_keep, n_col) else NULL
	sigma_u_row2_store <- if (use_u_row) numeric(n_keep) else NULL
	sigma_u_col2_store <- if (use_u_col) numeric(n_keep) else NULL
	keep_idx <- 0L

	verbose_int <- if (is.logical(verbose)) as.integer(verbose) else
		as.integer(verbose)
	if (verbose_int > 0L)
		cli::cli_progress_bar("dbn_with_covariates", total = nscan)

	# precompute zero-mean M for the shift trick
	M_zero <- array(0, dim = c(n_row, n_col, p))

	for (iter in seq_len(nscan)) {
		# 1. update augmented Z given current theta = R + L (binary / ordinal)
		if (family != "gaussian") {
			Theta_curr <- array(0, dim = c(n_row, n_col, Tt))
			for (t in seq_len(Tt)) Theta_curr[, , t] <- R[, , t] + L[, , t]
			Y_3d <- array(Y[, , 1, ], dim = c(n_row, n_col, Tt))
			if (family == "binary") {
				Z[, , 1, ] <- .dbn_sample_z_binary(Y_3d, Theta_curr, n_row, n_col, Tt)
			} else {
				Z[, , 1, ] <- .dbn_sample_z_ordinal(Y_3d, Theta_curr, n_row, n_col, Tt)
			}
		}

		# 2. shifted observation Z' = Z - L_total (as a 4D matrix nc x p*Tt)
		# L_total absorbs covariate L plus any random actor effects
		Zp <- Z
		for (t in seq_len(Tt)) Zp[, , 1, t] <- Z[, , 1, t] - L_total[, , t]
		Zp_4d <- matrix(Zp, nrow = nc, ncol = p * Tt)

		# 3. sample R via batch_ffbs_all_relations with M = 0
		# returns n_row*n_col x p*Tt; reshape into n_row x n_col x Tt
		R_cube <- batch_ffbs_all_relations(Zp_4d, M_zero, Aarray, Barray,
			sigma2, n_row, n_col, p, Tt)
		R_4d <- matrix(R_cube, nrow = nc, ncol = p * Tt)
		R[] <- array(R_cube, dim = c(n_row, n_col, Tt))

		# 4. sample A_t, B_t via update_AB_batch_extended on the residual state
		AB <- update_AB_batch_extended(
			R_4d, Aarray, Barray,
			sigma2, tauA2, tauB2,
			ar1, rhoA, rhoB,
			n_row, n_col, p, Tt,
			prior_kind = prior_kind,
			kappaA_inv = kappaA_inv, kappaB_inv = kappaB_inv
		)
		Aarray <- AB$Aarray
		Barray <- AB$Barray
		if (!all(is.finite(Aarray)) || !all(is.finite(Barray)))
			cli::cli_abort(c(
				"Sampler diverged at iteration {iter}: A_t or B_t contain non-finite values.",
				"i" = "Shorten the chain or rescale the design toward unit variance."
			))

		# 5. sample sigma2 (state-eq) and sigma2_obs (obs-eq, gaussian only)
		# via update_variances_dynamic with the shifted observation
		var_res <- update_variances_dynamic(
			R_4d, Zp_4d, M_zero, Aarray, Barray,
			a_sig, b_sig, n_row, n_col, p, Tt,
			is_gaussian = (family == "gaussian"),
			exclude_diagonal = (n_row == n_col)
		)
		sigma2 <- var_res$sigma2
		if (family == "gaussian") sigma2_obs <- var_res$sigma2_obs

		# 6. sample tauA2 / tauB2 from conjugate IG over A_t - A_{t-1} innovations
		# under prior_kind = "rw" the anchor is A_0 = B_0 = I
		if (prior_kind == "rw") {
			A0 <- diag(n_row); B0 <- diag(n_col)
			ssA <- sum((Aarray[, , 1] - A0)^2)
			ssB <- sum((Barray[, , 1] - B0)^2)
			if (Tt >= 2) {
				for (t in 2:Tt) {
					ssA <- ssA + sum((Aarray[, , t] - Aarray[, , t - 1])^2)
					ssB <- ssB + sum((Barray[, , t] - Barray[, , t - 1])^2)
				}
			}
			shape_A <- a_tau + (n_row * n_row * Tt) / 2
			shape_B <- a_tau + (n_col * n_col * Tt) / 2
			tauA2 <- 1 / stats::rgamma(1, shape = shape_A, rate = b_tau + ssA / 2)
			tauB2 <- 1 / stats::rgamma(1, shape = shape_B, rate = b_tau + ssB / 2)
		} else {
			# iid prior: innovations relative to identity at every t
			A0 <- diag(n_row); B0 <- diag(n_col)
			ssA <- 0; ssB <- 0
			for (t in seq_len(Tt)) {
				ssA <- ssA + sum((Aarray[, , t] - A0)^2)
				ssB <- ssB + sum((Barray[, , t] - B0)^2)
			}
			shape_A <- a_tau + (n_row * n_row * Tt) / 2
			shape_B <- a_tau + (n_col * n_col * Tt) / 2
			tauA2 <- 1 / stats::rgamma(1, shape = shape_A, rate = b_tau + ssA / 2)
			tauB2 <- 1 / stats::rgamma(1, shape = shape_B, rate = b_tau + ssB / 2)
		}

		# 7. update beta from the observation residual Z - R - u_row - u_col
		Z_3d <- array(Z[, , 1, ], dim = c(n_row, n_col, Tt))
		# subtract actor effects (broadcast across time) from Z before regressing
		Z_minus_u <- Z_3d
		if (use_u_row || use_u_col) {
			for (t in seq_len(Tt)) {
				Z_minus_u[, , t] <- Z_3d[, , t]
				if (use_u_row) Z_minus_u[, , t] <- Z_minus_u[, , t] -
					matrix(u_row, n_row, n_col, byrow = FALSE)
				if (use_u_col) Z_minus_u[, , t] <- Z_minus_u[, , t] -
					matrix(u_col, n_row, n_col, byrow = TRUE)
			}
		}
		if (K == 0L) {
			# actor-effects-only path: no beta to update
		} else if (time_varying_beta) {
			# FFBS over beta_{1:T} with state eq beta_t = beta_{t-1} + N(0, diag(tau_beta2))
			tv <- .dbn_update_beta_tv(covars, Z_minus_u, R, sigma2_obs,
				tau_beta2 = tau_beta2,
				prior_mean = prior_beta_mean,
				prior_sd   = prior_beta_scale[1])
			beta_mat <- tv$beta_mat
			# update tau_beta2_k from RW innovations (conjugate IG per term)
			for (k in seq_len(K)) {
				if (Tt >= 2) {
					innov <- diff(beta_mat[k, ])
					ss <- sum(innov^2)
				} else {
					ss <- 0
				}
				tau_beta2[k] <- 1 / stats::rgamma(1,
					shape = a_beta + (Tt - 1L) / 2,
					rate  = b_beta + ss / 2)
			}
		} else {
			beta_update <- .dbn_update_beta_obs(covars, Z_minus_u, R, sigma2_obs,
				prior_mean = prior_beta_mean, prior_sd = prior_beta_scale[1])
			beta_mat[, 1L] <- beta_update$beta
		}
		L <- build_L(beta_mat)

		# 8. update random actor effects (sender u_row and receiver u_col)
		# full conditional given Z, R, L, sigma2_obs, sigma_u_row2 / sigma_u_col2
		if (use_u_row) {
			# residual after removing R and L (and current u_col)
			# Y_resid[i, j, t] = Z[i, j, t] - R[i, j, t] - L[i, j, t] - u_col[j]
			# u_row[i] | rest ~ N(mu_i, V_i) with:
			#   V_i = (n_eff / sigma2_obs + 1 / sigma_u_row2)^{-1}
			#   mu_i = V_i * sum_{j != i, t} resid[i, j, t] / sigma2_obs
			# where n_eff = number of off-diagonal cells per actor per time, summed over t
			cell_sums <- numeric(n_row)
			off_diag <- if (n_row == n_col) seq_len(n_col) else seq_len(n_col)
			n_eff <- if (n_row == n_col) (n_col - 1L) * Tt else n_col * Tt
			for (i in seq_len(n_row)) {
				js <- if (n_row == n_col) setdiff(seq_len(n_col), i) else seq_len(n_col)
				ssum <- 0
				for (t in seq_len(Tt)) {
					rij <- Z_3d[i, js, t] - R[i, js, t] - L[i, js, t]
					if (use_u_col) rij <- rij - u_col[js]
					ssum <- ssum + sum(rij)
				}
				cell_sums[i] <- ssum
			}
			V <- 1 / (n_eff / sigma2_obs + 1 / sigma_u_row2)
			mu <- V * cell_sums / sigma2_obs
			u_row <- mu + sqrt(V) * stats::rnorm(n_row)
			# sum-to-zero centering
			u_row <- u_row - mean(u_row)
			# update sigma_u_row2 via conjugate IG
			ss_u_row <- sum(u_row^2)
			sigma_u_row2 <- 1 / stats::rgamma(1,
				shape = a_u + n_row / 2,
				rate = b_u + ss_u_row / 2)
		}
		if (use_u_col) {
			n_eff <- if (n_row == n_col) (n_row - 1L) * Tt else n_row * Tt
			cell_sums <- numeric(n_col)
			for (j in seq_len(n_col)) {
				is <- if (n_row == n_col) setdiff(seq_len(n_row), j) else seq_len(n_row)
				ssum <- 0
				for (t in seq_len(Tt)) {
					rij <- Z_3d[is, j, t] - R[is, j, t] - L[is, j, t]
					if (use_u_row) rij <- rij - u_row[is]
					ssum <- ssum + sum(rij)
				}
				cell_sums[j] <- ssum
			}
			V <- 1 / (n_eff / sigma2_obs + 1 / sigma_u_col2)
			mu <- V * cell_sums / sigma2_obs
			u_col <- mu + sqrt(V) * stats::rnorm(n_col)
			u_col <- u_col - mean(u_col)
			ss_u_col <- sum(u_col^2)
			sigma_u_col2 <- 1 / stats::rgamma(1,
				shape = a_u + n_col / 2,
				rate = b_u + ss_u_col / 2)
		}

		# rebuild L_total with new beta and actor effects
		L_total <- build_L_total(L, u_row, u_col)
		# also stash theta in original frame for storage / downstream
		Theta_curr <- array(0, dim = c(n_row, n_col, Tt))
		for (t in seq_len(Tt)) Theta_curr[, , t] <- R[, , t] + L_total[, , t]

		# store
		if (iter > burn && (iter - burn) %% odens == 0L) {
			keep_idx <- keep_idx + 1L
			if (K == 0L) {
				# actor-effects-only: nothing to store for beta
			} else if (time_varying_beta) {
				beta_store[keep_idx, , ] <- beta_mat
				tau_beta2_store[keep_idx, ] <- tau_beta2
			} else {
				beta_store[keep_idx, ]   <- beta_mat[, 1L]
			}
			sigma2_store[keep_idx]    <- sigma2
			sigma2_obs_store[keep_idx] <- sigma2_obs
			tauA2_store[keep_idx]     <- tauA2
			tauB2_store[keep_idx]     <- tauB2
			A_store[[keep_idx]]       <- Aarray
			B_store[[keep_idx]]       <- Barray
			Theta_store[, , 1, , keep_idx] <- Theta_curr
			L_store[, , , keep_idx]   <- L_total
			if (use_u_row) {
				u_row_store[keep_idx, ] <- u_row
				sigma_u_row2_store[keep_idx] <- sigma_u_row2
			}
			if (use_u_col) {
				u_col_store[keep_idx, ] <- u_col
				sigma_u_col2_store[keep_idx] <- sigma_u_col2
			}
		}

		if (verbose_int > 0L && (iter %% max(1L, floor(nscan / 50)) == 0L))
			cli::cli_progress_update()
	}
	if (verbose_int > 0L) cli::cli_progress_done()

	# assemble fit object compatible with downstream consumers
	dims <- list(n_row = n_row, n_col = n_col, p = p, Tt = Tt,
		is_bipartite = (n_row != n_col), is_symmetric = FALSE)
	params_df <- data.frame(
		sigma2     = sigma2_store,
		sigma2_obs = sigma2_obs_store,
		tau_A2     = tauA2_store,
		tau_B2     = tauB2_store
	)
	out <- list(
		model = "dynamic",
		family = family,
		Y = Y,
		dims = dims,
		settings = list(nscan = nscan, burn = burn, odens = odens, draws = n_keep),
		params = params_df,
		A = A_store,
		B = B_store,
		Theta = Theta_store,
		sigma2     = sigma2_store,
		sigma2_obs = sigma2_obs_store,
		tau_A2     = tauA2_store,
		tau_B2     = tauB2_store,
		covariates = covars,
		beta  = beta_store,
		tau_beta2 = tau_beta2_store,
		time_varying_beta = time_varying_beta,
		L     = L_store,
		u_row = u_row_store,
		u_col = u_col_store,
		sigma_u_row2 = sigma_u_row2_store,
		sigma_u_col2 = sigma_u_col2_store,
		actor_effects = actor_effects,
		draws = list(
			beta = beta_store,
			pars = params_df
		),
		symmetric = FALSE,
		meta = list(
			sampler_used = "ffbs_with_covariates",
			uncertainty_available = TRUE,
			prior_kind = prior_kind
		)
	)
	class(out) <- c("dbn_covariates_fit", "dbn", "list")
	out
}

####
# initialize augmented Z for ordinal / binary based on Y
# binary: sign-coded probit init; ordinal: rank-to-normal-quantile
.dbn_init_z_aug <- function(Y, family, n_row, n_col, Tt) {
	Z <- array(0, dim = c(n_row, n_col, Tt))
	if (family == "binary") {
		for (t in seq_len(Tt)) {
			Y_t <- Y[, , t]
			Z[, , t] <- ifelse(is.na(Y_t), 0,
				ifelse(Y_t == 1, 0.5, -0.5))
		}
	} else {
		# ordinal: rough rank-to-normal-quantile init per time
		for (t in seq_len(Tt)) {
			Y_t <- Y[, , t]
			n_obs <- sum(!is.na(Y_t))
			if (n_obs > 0L) {
				Z[, , t] <- stats::qnorm(
					rank(Y_t, na.last = "keep") / (n_obs + 1L))
			}
		}
	}
	Z
}

####
# truncated normal helpers (univariate, N(mu, 1) constrained)
.rtruncnorm_above <- function(mu, lower) {
	p_lower <- stats::pnorm(lower - mu)
	u <- stats::runif(1, p_lower, 1)
	mu + stats::qnorm(u)
}
.rtruncnorm_below <- function(mu, upper) {
	p_upper <- stats::pnorm(upper - mu)
	u <- stats::runif(1, 0, p_upper)
	mu + stats::qnorm(u)
}

####
# sample binary augmented Z element-wise given current theta and Y
.dbn_sample_z_binary <- function(Y_all, Theta_all, n_row, n_col, Tt) {
	Z <- array(0, dim = c(n_row, n_col, Tt))
	bipartite <- (n_row != n_col)
	for (t in seq_len(Tt)) {
		Y_t  <- Y_all[, , t]
		mu_t <- Theta_all[, , t]
		for (i in seq_len(n_row)) {
			for (j in seq_len(n_col)) {
				if (!bipartite && i == j) next
				y_ij <- Y_t[i, j]
				m_ij <- mu_t[i, j]
				if (is.na(y_ij)) {
					Z[i, j, t] <- stats::rnorm(1, m_ij, 1)
				} else if (y_ij == 1) {
					Z[i, j, t] <- .rtruncnorm_above(m_ij, 0)
				} else {
					Z[i, j, t] <- .rtruncnorm_below(m_ij, 0)
				}
			}
		}
	}
	Z
}

####
# sample ordinal augmented Z via rank likelihood (hoff 2007) per time
# preserves order of Z within each time slice consistent with Y
.dbn_sample_z_ordinal <- function(Y_all, Theta_all, n_row, n_col, Tt) {
	Z <- array(0, dim = c(n_row, n_col, Tt))
	for (t in seq_len(Tt)) {
		Z[, , t] <- .rz_rank_likelihood(Y_all[, , t], Theta_all[, , t],
			n_row, n_col)
	}
	Z
}

####
# rank-likelihood Z sampler given current theta and ordered Y
# initializes Z at theta then resamples each observed cell from a normal
# centered at theta truncated to the rank-induced interval
.rz_rank_likelihood <- function(Y_t, Theta_t, n_row, n_col) {
	Z <- Theta_t
	ranks_Y <- rank(Y_t, ties.method = "average", na.last = "keep")
	idx <- which(!is.na(Y_t))
	if (length(idx) == 0L) return(Z)
	ord <- order(ranks_Y[idx])
	idx_sorted <- idx[ord]
	prev_z <- -Inf
	for (k in seq_along(idx_sorted)) {
		i <- ((idx_sorted[k] - 1L) %% n_row) + 1L
		j <- ((idx_sorted[k] - 1L) %/% n_row) + 1L
		mu <- Theta_t[i, j]
		next_z <- if (k < length(idx_sorted)) {
			i_n <- ((idx_sorted[k + 1L] - 1L) %% n_row) + 1L
			j_n <- ((idx_sorted[k + 1L] - 1L) %/% n_row) + 1L
			Z[i_n, j_n]
		} else Inf
		if (prev_z >= next_z) {
			Z[i, j] <- mu
		} else {
			p_lo <- stats::pnorm(prev_z - mu)
			p_hi <- stats::pnorm(next_z - mu)
			u <- stats::runif(1, p_lo, p_hi)
			Z[i, j] <- mu + stats::qnorm(u)
		}
		prev_z <- Z[i, j]
	}
	Z
}
