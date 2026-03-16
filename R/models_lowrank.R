####
#' Low-rank helper functions
#' @name lowrank-helpers
#' @keywords internal
NULL

#' Initialize low-rank model
#' @param Y Data array
#' @param r Rank
#' @return List with U and alpha
#' @keywords internal
init_lowrank <- function(Y, r) {
	m <- dim(Y)[1]
	Tt <- dim(Y)[4]

	# SVD of time-averaged adjacency gives starting directions
	Ybar <- apply(Y, c(1, 2), mean, na.rm = TRUE)
	Ybar[is.na(Ybar)] <- 0
	sv <- svd(Ybar, nu = r, nv = r)
	U <- sv$u[, 1:r, drop = FALSE]
	alpha <- matrix(0, r, Tt)

	list(U = U, alpha = alpha)
}

#' Log-likelihood for U
#' @keywords internal
loglik_U <- function(U, alpha, Theta, Barray, sigma2) {
	m <- nrow(U)
	Tt <- ncol(alpha)
	ll <- 0

	for (t in 2:Tt) {
		alpha_t <- alpha[, t, drop = TRUE]
		if (length(alpha_t) == 1) {
			A_t <- alpha_t * U %*% t(U)
		} else {
			A_t <- U %*% diag(alpha_t) %*% t(U)
		}
		pred <- A_t %*% Theta[, , t - 1] %*% t(Barray[, , t])
		resid <- Theta[, , t] - pred
		ll <- ll - 0.5 * sum(resid^2) / sigma2
	}
	ll
}
####

####
#' DBN Low-rank Model (Accurate)
#'
#' @description Fits DBN with low-rank sender effects using batch FFBS
#' @param Y Data array (nodes x nodes x relations x time)
#' @param r Rank for low-rank factorization. A good starting point is
#'   \code{ceiling(log2(n))} where \code{n} is the number of nodes. Increase
#'   if posterior predictive checks show poor fit. Default: 2.
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param ar1_alpha Use AR(1) for alpha dynamics
#' @param update_rho_alpha Update AR coefficient for alpha
#' @param ar1_B Use AR(1) for B dynamics
#' @param update_rho_B Update AR coefficient for B
#' @param seed Random seed
#' @param verbose Logical or numeric. TRUE prints every 100 iterations, numeric prints every n iterations, FALSE suppresses output.
#' @param time_thin Save every nth time point
#' @param family Data family
#' @param previous Previous results to continue from
#' @param init Initial values
#' @param symmetric Logical. Not supported for low-rank models (will error). Default: FALSE.
#' @param ... Additional arguments (currently unused)
#' @return List containing MCMC results
#' @seealso \code{\link{dbn}} for the main dispatcher, \code{\link{param_summary}} for posterior summaries
#' @examples
#' \donttest{
#' sim <- simulate_lowrank_dbn(n = 8, time = 5, r = 2, seed = 1)
#' fit <- dbn_lowrank(sim$Y, r = 2, nscan = 200, burn = 100, verbose = FALSE)
#' }
#' @export
dbn_lowrank_accurate <- function(Y,
						family = c("ordinal", "gaussian", "binary"),
						r = 2,
						nscan = 10000,
						burn = 1000,
						odens = 1,
						ar1_alpha = TRUE,
						update_rho_alpha = FALSE,
						ar1_B = FALSE,
						update_rho_B = FALSE,
						seed = 6886,
						verbose = TRUE,
						time_thin = 1,
						previous = NULL,
						init = NULL,
						symmetric = FALSE,
						...) {

	if (symmetric) {
		cli::cli_abort(c(
			"Symmetric networks are not yet supported for low-rank models.",
			"i" = "Use {.code model = \"dynamic\"} with {.code symmetric = TRUE} instead."
		))
	}

	# family setup
	family <- match.arg(family)
	if (isTRUE(verbose)) verbose <- 100L
	if (isFALSE(verbose)) verbose <- 0L
	FAM <- switch(family,
		ordinal = family_ordinal(),
		gaussian = family_gaussian(),
		binary = family_binary()
	)

	set.seed(seed)

	if (dim(Y)[4] < 2) {
		cli::cli_abort("Low-rank model requires at least 2 time points. Use {.code model = \"static\"} for cross-sectional data.")
	}
	####

	# preprocessing
	pre <- shared_preprocess(Y, family = family)
	dims <- pre$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	is_bipartite <- dims$is_bipartite
	m <- n_row
	p <- dims$p
	Tt <- dims$Tt
	nc <- n_row * n_col

	# bipartite networks not yet supported for low-rank
	if (is_bipartite) {
		cli::cli_abort(c(
			"Low-rank model does not yet support bipartite (rectangular) networks.",
			"i" = "Your data has {n_row} senders and {n_col} receivers.",
			"i" = "Use {.code model = \"dynamic\"} for bipartite networks."
		))
	}

	# check rank bounds
	if (r < 1) cli::cli_abort("{.arg r} must be at least 1.")
	if (r > min(n_row, n_col)) cli::cli_abort("{.arg r} ({r}) cannot exceed min(n_row, n_col) = {min(n_row, n_col)}.")
	####

	# initialization
	if (!is.null(previous)) {
		last_idx <- length(previous$sigma2_proc)
		sigma2_proc <- previous$sigma2_proc[last_idx]
		tau_alpha2 <- previous$tau_alpha2[last_idx]
		tau_B2 <- previous$tau_B2[last_idx]
		g2 <- previous$g2[last_idx]
		U <- previous$U[[last_idx]]
		alpha <- previous$alpha[[last_idx]]

		# expand alpha to full time grid if thinned
		if (ncol(alpha) < Tt) {
			alpha_full <- matrix(0, r, Tt)
			prev_times <- seq(1, Tt, by = previous$settings$time_thin)
			for (t in 1:Tt) {
				nearest <- which.min(abs(prev_times - t))
				alpha_full[, t] <- alpha[, nearest]
			}
			alpha <- alpha_full
		}

		Barray <- previous$B[[last_idx]]
		if (dim(Barray)[3] < Tt) {
			Barray_full <- array(0, dim = c(n_col, n_col, Tt))
			prev_times <- seq(1, Tt, by = previous$settings$time_thin)
			for (t in 1:Tt) {
				nearest <- which.min(abs(prev_times - t))
				Barray_full[, , t] <- Barray[, , nearest]
			}
			Barray <- Barray_full
		}

		rho_alpha <- if (ar1_alpha && !is.null(previous$rho_alpha)) previous$rho_alpha[last_idx] else if (ar1_alpha) 0.9 else 0
		rho_B <- if (ar1_B && !is.null(previous$rho_B)) previous$rho_B[last_idx] else if (ar1_B) 0.9 else 0
		sigma2_obs <- if (FAM$name == "gaussian" && !is.null(previous$sigma2_obs)) previous$sigma2_obs[last_idx] else FAM$init_pars$sigma2_obs %||% 1

		if (verbose) cli::cli_inform("Continuing from previous run")
	} else if (!is.null(init)) {
		U <- init$U %||% init_lowrank(Y, r)$U
		alpha <- init$alpha %||% init_lowrank(Y, r)$alpha
		Barray <- init$B %||% {
			B <- array(0, dim = c(n_col, n_col, Tt))
			for (t in 1:Tt) B[, , t] <- diag(n_col)
			B
		}
		sigma2_proc <- init$sigma2_proc %||% 0.5
		sigma2_obs <- init$sigma2_obs %||% (FAM$init_pars$sigma2_obs %||% 1)
		tau_alpha2 <- init$tau_alpha2 %||% 1
		tau_B2 <- init$tau_B2 %||% 1
		g2 <- init$g2 %||% 1
		rho_alpha <- if (ar1_alpha) (init$rho_alpha %||% 0.9) else 0
		rho_B <- if (ar1_B) (init$rho_B %||% 0.9) else 0

		if (verbose) cli::cli_inform("Using provided initial values")
	} else {
		lr <- init_lowrank(Y, r)
		U <- lr$U
		alpha <- lr$alpha
		tau_alpha2 <- 1
		rho_alpha <- if (ar1_alpha) 0.9 else 0
		sigma2_proc <- 0.5
		sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
		g2 <- 1

		Barray <- array(0, dim = c(n_col, n_col, Tt))
		for (t in 1:Tt) Barray[, , t] <- diag(n_col)
		tau_B2 <- 1
		rho_B <- if (ar1_B) 0.9 else 0
	}
	####

	# MCMC storage
	n_iter <- burn + nscan
	keep_idx <- ((burn + 1):n_iter)[((burn + 1):n_iter) %% odens == 0]
	S <- length(keep_idx)

	Usave <- vector("list", S)
	alphasave <- vector("list", S)
	Bsave <- vector("list", S)
	Msave <- vector("list", S)
	Thetasave <- vector("list", S)
	sigmasave <- numeric(S)
	tau_alpha_save <- numeric(S)
	tau_B_save <- numeric(S)
	g2save <- numeric(S)
	if (FAM$name == "gaussian") sigma2_obs_save <- numeric(S)
	if (ar1_alpha && update_rho_alpha) rho_alpha_save <- numeric(S)
	if (ar1_B && update_rho_B) rho_B_save <- numeric(S)

	time_keep <- seq(1, Tt, by = time_thin)

	if (verbose) {
		cli::cli_alert_info("Running low-rank DBN MCMC ({n_iter} iterations)")
		cli::cli_progress_bar("MCMC iterations", total = n_iter)
	}

	s <- 0
	epsilon <- 0.005
	accept_U <- 0

	a_tau <- 10
	b_tau <- 10

	Aarray <- compute_all_A_lowrank(U, alpha, Tt)

	# flatten arrays for batch operations
	total_slices <- p * Tt
	Z_flat <- array(0, dim = c(n_row, n_col, total_slices))
	Theta_flat <- array(0, dim = c(n_row, n_col, total_slices))
	M_flat <- array(0, dim = c(n_row, n_col, p))

	for (j in 1:p) {
		M_flat[, , j] <- pre$M[, , j]
		slice_idx <- ((j-1)*Tt + 1):(j*Tt)
		for (i in seq_along(slice_idx)) {
			Z_flat[, , slice_idx[i]] <- pre$Z[, , j, i]
			Theta_flat[, , slice_idx[i]] <- pre$Theta[, , j, i]
		}
	}
	####

	# main MCMC loop
	for (g in seq_len(n_iter)) {
		# update Z for ordinal/binary
		if (FAM$name == "ordinal") {
			pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "ordinal")
			tmp <- aperm(pre$Z, c(1L, 2L, 4L, 3L))
			dim(tmp) <- c(n_row, n_col, p * Tt)
			Z_flat[] <- tmp
		} else if (FAM$name == "binary") {
			pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "binary")
			tmp <- aperm(pre$Z, c(1L, 2L, 4L, 3L))
			dim(tmp) <- c(n_row, n_col, p * Tt)
			Z_flat[] <- tmp
		}

		# M update (vectorized on flat arrays)
		for (j in 1:p) {
			offset <- (j-1) * Tt
			resid_sum <- Z_flat[, , offset + seq_len(Tt), drop = FALSE] -
						 Theta_flat[, , offset + seq_len(Tt), drop = FALSE]
			mu_sum <- rowSums(resid_sum, dims = 2L)
			mu_var <- 1 / (Tt + 1 / g2)
			M_flat[, , j] <- mu_sum * mu_var + matrix(sqrt(mu_var) * rnorm(nc), n_row, n_col)
		}
		M_ss <- sum(M_flat^2)
		g2 <- safe_rinv_gamma(10 + length(M_flat)/2, 10 + M_ss/2)

		# Theta update via FFBS
		Theta_flat <- ffbs_theta_all_relations(Z_flat, M_flat, Aarray, Barray,
											  sigma2_proc, sigma2_obs, m, p, Tt)

		# alpha update (blocked)
		alpha <- update_alpha_batch(Theta_flat, U, Barray,
								   sigma2_proc,
								   tau_alpha2,
								   ar1_alpha, rho_alpha, m, p, Tt, r)

		Aarray <- compute_all_A_lowrank(U, alpha, Tt)

		# B update (parallelized)
		Barray <- update_B_parallel(Theta_flat, Aarray,
										 sigma2_proc,
										 tau_B2,
										 ar1_B, rho_B, m, p, Tt)

		# update tau_alpha2
		if (ar1_alpha) {
			innov <- alpha[, -1, drop = FALSE] - rho_alpha * alpha[, -Tt, drop = FALSE]
		} else {
			innov <- alpha[, -1, drop = FALSE] - alpha[, -Tt, drop = FALSE]
		}
		innov_ss <- sum(innov^2)
		tau_alpha2 <- safe_rinv_gamma(a_tau + 0.5 * r * (Tt - 1), b_tau + 0.5 * innov_ss)

		# rho_alpha MH step (optional)
		if (ar1_alpha && update_rho_alpha) {
			rho_prop <- rho_alpha + rnorm(1, 0, 0.05)
			if (abs(rho_prop) < 1) {
				ll_prop <- -sum((alpha[, -1] - rho_prop * alpha[, -Tt])^2) / (2 * tau_alpha2)
				ll_curr <- -sum((alpha[, -1] - rho_alpha * alpha[, -Tt])^2) / (2 * tau_alpha2)
				if (log(runif(1)) < ll_prop - ll_curr) {
					rho_alpha <- rho_prop
				}
			}
		}

		# U via Metropolis on the Stiefel manifold
		u_update_freq <- max(5, floor(m / 10))
		if (g %% u_update_freq == 0) {
			W_vec <- rnorm(m * (m - 1) / 2)
			W <- matrix(0, m, m)
			W[lower.tri(W)] <- W_vec
			W <- W - t(W)

			norm_cap <- sqrt(m) * epsilon
			W_norm <- sqrt(sum(W^2))
			if (W_norm > norm_cap) {
				W <- W * (norm_cap / W_norm)
			}

			I_m <- diag(m)
			half_W <- epsilon * W / 2
			U_prop <- U + 2 * half_W %*% solve(I_m - half_W %*% half_W) %*% U

			# average Theta across relations for the likelihood
			if (p == 1L) {
				Theta_avg <- Theta_flat[, , 1:Tt, drop = FALSE]
			} else {
				idx <- as.vector(outer(seq_len(Tt), (seq_len(p) - 1L) * Tt, "+"))
				dim_tf <- dim(Theta_flat)
				arr <- array(Theta_flat[, , idx], c(dim_tf[1] * dim_tf[2], Tt, p))
				Theta_avg <- array(rowMeans(arr, dims = 2L), c(n_row, n_col, Tt))
			}

			log_acc <- loglik_U(U_prop, alpha, Theta_avg, Barray, sigma2_proc) -
					   loglik_U(U, alpha, Theta_avg, Barray, sigma2_proc)

			if (!is.na(log_acc) && is.finite(log_acc) && log(runif(1)) < log_acc) {
				U <- U_prop
				accept_U <- accept_U + 1
				Aarray <- compute_all_A_lowrank(U, alpha, Tt)
			}
		}

		# update tau_B2
		if (ar1_B) {
			innovB <- Barray[, , 2:Tt] - rho_B * Barray[, , 1:(Tt - 1)]
		} else {
			innovB <- Barray[, , 2:Tt] - Barray[, , 1:(Tt - 1)]
		}
		innovB_ss <- sum(innovB^2)
		tau_B2 <- safe_rinv_gamma(a_tau + n_col * n_col * (Tt - 1) / 2, b_tau + innovB_ss / 2)

		# rho_B MH step (optional)
		if (ar1_B && update_rho_B) {
			rho_prop <- rho_B + rnorm(1, 0, 0.05)
			if (abs(rho_prop) < 1) {
				innovB_prop <- Barray[, , 2:Tt] - rho_prop * Barray[, , 1:(Tt - 1)]
				ll_prop <- -sum(innovB_prop^2) / (2 * tau_B2)
				innovB_curr <- Barray[, , 2:Tt] - rho_B * Barray[, , 1:(Tt - 1)]
				ll_curr <- -sum(innovB_curr^2) / (2 * tau_B2)
				if (log(runif(1)) < ll_prop - ll_curr) {
					rho_B <- rho_prop
				}
			}
		}

		# sigma2_proc via batch C++ residuals
		resid_ss <- compute_sigma2_lowrank_batch(Theta_flat, Aarray, Barray, m, p, Tt)
		sigma2_proc <- safe_rinv_gamma(a_tau + nc * p * (Tt - 1) / 2, b_tau + resid_ss / 2)

		# observation variance for gaussian
		if (FAM$name == "gaussian") {
			resid_obs <- compute_gaussian_obs_residuals_flat(Z_flat, Theta_flat, M_flat, m, p, Tt)
			sigma2_obs <- safe_rinv_gamma(a_tau + nc * p * Tt / 2, b_tau + resid_obs / 2)
		}

		# store draws
		if (g %in% keep_idx) {
			s <- s + 1
			Usave[[s]] <- U
			alphasave[[s]] <- alpha[, time_keep, drop = FALSE]
			Bsave[[s]] <- Barray[, , time_keep, drop = FALSE]
			sigmasave[s] <- sigma2_proc
			tau_alpha_save[s] <- tau_alpha2
			tau_B_save[s] <- tau_B2
			g2save[s] <- g2

			if (FAM$name == "gaussian") sigma2_obs_save[s] <- sigma2_obs
			if (ar1_alpha && update_rho_alpha) rho_alpha_save[s] <- rho_alpha
			if (ar1_B && update_rho_B) rho_B_save[s] <- rho_B

			Msave[[s]] <- M_flat
			Theta_4d <- array(NA, c(n_row, n_col, p, length(time_keep)))
			for (j_save in 1:p) {
				for (i_save in seq_along(time_keep)) {
					idx_save <- (j_save - 1) * Tt + time_keep[i_save]
					Theta_4d[, , j_save, i_save] <- Theta_flat[, , idx_save]
				}
			}
			Thetasave[[s]] <- Theta_4d
		}

		# adapt step size for U proposal
		if (g %% (u_update_freq * 20) == 0) {
			accept_rate <- accept_U / 20
			if (accept_rate < 0.2) {
				epsilon <- max(epsilon * 0.9, 1e-4)
			} else if (accept_rate > 0.3) {
				epsilon <- min(epsilon * 1.1, 0.1)
			}
			accept_U <- 0
		}

		if (verbose) cli::cli_progress_update()

		if (is.numeric(verbose) && verbose > 0 && g %% verbose == 0) {
			msg <- "iter {g}: sigma^2={round(sigma2_proc, 3)} tauA^2={round(tau_alpha2, 3)} tauB^2={round(tau_B2, 3)}"
			if (ar1_alpha && update_rho_alpha) msg <- paste0(msg, " rhoA={round(rho_alpha, 3)}")
			if (ar1_B && update_rho_B) msg <- paste0(msg, " rhoB={round(rho_B, 3)}")
			cli::cli_inform(msg)
		}
	}
	####

	# assemble output
	if (verbose) cli::cli_progress_done()

	# unflatten arrays back to 4D
	for (j in 1:p) {
		pre$M[, , j] <- M_flat[, , j]
		slice_idx <- ((j-1)*Tt + 1):(j*Tt)
		for (i in seq_along(slice_idx)) {
			pre$Z[, , j, i] <- Z_flat[, , slice_idx[i]]
			pre$Theta[, , j, i] <- Theta_flat[, , slice_idx[i]]
		}
	}

	# reconstruct A arrays from U and alpha draws
	Asave <- vector("list", S)
	for (s in 1:S) {
		A_s <- array(0, dim = c(n_row, n_row, length(time_keep)))
		for (i in seq_along(time_keep)) {
			alpha_i <- alphasave[[s]][, i, drop = TRUE]
			if (length(alpha_i) == 1) {
				A_s[, , i] <- alpha_i * Usave[[s]] %*% t(Usave[[s]])
			} else {
				A_s[, , i] <- Usave[[s]] %*% diag(alpha_i) %*% t(Usave[[s]])
			}
		}
		Asave[[s]] <- A_s
	}

	pars_df <- data.frame(
		sigma2_proc = sigmasave,
		tau_alpha2 = tau_alpha_save,
		tau_B2 = tau_B_save,
		g2 = g2save
	)

	if (FAM$name == "gaussian") pars_df$sigma2_obs <- sigma2_obs_save
	if (ar1_alpha && update_rho_alpha) pars_df$rho_alpha <- rho_alpha_save
	if (ar1_B && update_rho_B) pars_df$rho_B <- rho_B_save

	out <- structure(
		list(
			draws = list(
				theta = Thetasave,
				z = NULL,
				pars = pars_df,
				misc = list(
					U = Usave,
					alpha = alphasave,
					A = Asave,
					B = Bsave,
					M = Msave
				)
			),
			meta = list(
				dims = dims,
				draws = S,
				thin = odens,
				time_thin = time_thin,
				time_kept = time_keep,
				model = "lowrank"
			),
			family = FAM,
			model = "lowrank",
			U = Usave,
			alpha = alphasave,
			A = Asave,
			B = Bsave,
			M = Msave,
			sigma2_proc = sigmasave,
			tau_alpha2 = tau_alpha_save,
			tau_B2 = tau_B_save,
			g2 = g2save,
			Y = Y,
			dims = dims,
			settings = list(
				r = r,
				nscan = nscan,
				burn = burn,
				odens = odens,
				ar1_alpha = ar1_alpha,
				ar1_B = ar1_B,
				time_thin = time_thin,
				family = family
			)
		),
		class = "dbn"
	)

	if (FAM$name == "gaussian") {
		out$sigma2_obs <- sigma2_obs_save
	}

	if (ar1_alpha && update_rho_alpha) {
		out$rho_alpha <- rho_alpha_save
	}
	if (ar1_B && update_rho_B) {
		out$rho_B <- rho_B_save
	}

	out
}
####

####
#' DBN Low-rank Model
#'
#' @description Fits DBN with low-rank sender effects
#' @inheritParams dbn_lowrank_accurate
#' @param ... Additional arguments (currently unused)
#' @seealso \code{\link{dbn}} for the main dispatcher, \code{\link{param_summary}} for posterior summaries
#' @examples
#' \donttest{
#' sim <- simulate_lowrank_dbn(n = 8, time = 5, r = 2, seed = 1)
#' fit <- dbn_lowrank(sim$Y, r = 2, nscan = 200, burn = 100, verbose = FALSE)
#' }
#' @export
dbn_lowrank <- function(Y,
						family = c("ordinal", "gaussian", "binary"),
						r = 2,
						nscan = 10000,
						burn = 1000,
						odens = 1,
						ar1_alpha = TRUE,
						update_rho_alpha = FALSE,
						ar1_B = FALSE,
						update_rho_B = FALSE,
						seed = 6886,
						verbose = TRUE,
						time_thin = 1,
						previous = NULL,
						init = NULL,
						symmetric = FALSE,
						...) {

	if (symmetric) {
		cli::cli_abort(c(
			"Symmetric networks are not yet supported for low-rank models.",
			"i" = "Use {.code model = \"dynamic\"} with {.code symmetric = TRUE} instead."
		))
	}

	# family setup
	family <- match.arg(family)
	if (isTRUE(verbose)) verbose <- 100L
	if (isFALSE(verbose)) verbose <- 0L
	FAM <- switch(family,
		ordinal = family_ordinal(),
		gaussian = family_gaussian(),
		binary = family_binary()
	)

	set.seed(seed)

	if (dim(Y)[4] < 2) {
		cli::cli_abort("Low-rank model requires at least 2 time points. Use {.code model = \"static\"} for cross-sectional data.")
	}
	####

	# preprocessing
	pre <- shared_preprocess(Y, family = family)
	dims <- pre$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	is_bipartite <- dims$is_bipartite
	m <- n_row
	p <- dims$p
	Tt <- dims$Tt
	nc <- n_row * n_col

	# bipartite networks not yet supported for low-rank
	if (is_bipartite) {
		cli::cli_abort(c(
			"Low-rank model does not yet support bipartite (rectangular) networks.",
			"i" = "Your data has {n_row} senders and {n_col} receivers.",
			"i" = "Use {.code model = \"dynamic\"} for bipartite networks."
		))
	}

	# check rank bounds
	if (r < 1) cli::cli_abort("{.arg r} must be at least 1.")
	if (r > min(n_row, n_col)) cli::cli_abort("{.arg r} ({r}) cannot exceed min(n_row, n_col) = {min(n_row, n_col)}.")
	####

	# initialization
	if (!is.null(previous)) {
		last_idx <- length(previous$sigma2_proc)
		sigma2_proc <- previous$sigma2_proc[last_idx]
		tau_alpha2 <- previous$tau_alpha2[last_idx]
		tau_B2 <- previous$tau_B2[last_idx]
		g2 <- previous$g2[last_idx]
		U <- previous$U[[last_idx]]
		alpha <- previous$alpha[[last_idx]]

		# expand alpha to full time grid if thinned
		if (ncol(alpha) < Tt) {
			alpha_full <- matrix(0, r, Tt)
			prev_times <- seq(1, Tt, by = previous$settings$time_thin)
			for (t in 1:Tt) {
				nearest <- which.min(abs(prev_times - t))
				alpha_full[, t] <- alpha[, nearest]
			}
			alpha <- alpha_full
		}

		Barray <- previous$B[[last_idx]]
		if (dim(Barray)[3] < Tt) {
			Barray_full <- array(0, dim = c(n_col, n_col, Tt))
			prev_times <- seq(1, Tt, by = previous$settings$time_thin)
			for (t in 1:Tt) {
				nearest <- which.min(abs(prev_times - t))
				Barray_full[, , t] <- Barray[, , nearest]
			}
			Barray <- Barray_full
		}

		rho_alpha <- if (ar1_alpha && !is.null(previous$rho_alpha)) previous$rho_alpha[last_idx] else if (ar1_alpha) 0.9 else 0
		rho_B <- if (ar1_B && !is.null(previous$rho_B)) previous$rho_B[last_idx] else if (ar1_B) 0.9 else 0
		sigma2_obs <- if (FAM$name == "gaussian" && !is.null(previous$sigma2_obs)) previous$sigma2_obs[last_idx] else FAM$init_pars$sigma2_obs %||% 1

		if (verbose) cli::cli_inform("Continuing from previous run")
	} else if (!is.null(init)) {
		U <- init$U %||% init_lowrank(Y, r)$U
		alpha <- init$alpha %||% init_lowrank(Y, r)$alpha
		Barray <- init$B %||% {
			B <- array(0, dim = c(n_col, n_col, Tt))
			for (t in 1:Tt) B[, , t] <- diag(n_col)
			B
		}
		sigma2_proc <- init$sigma2_proc %||% 0.5
		sigma2_obs <- init$sigma2_obs %||% (FAM$init_pars$sigma2_obs %||% 1)
		tau_alpha2 <- init$tau_alpha2 %||% 1
		tau_B2 <- init$tau_B2 %||% 1
		g2 <- init$g2 %||% 1
		rho_alpha <- if (ar1_alpha) (init$rho_alpha %||% 0.9) else 0
		rho_B <- if (ar1_B) (init$rho_B %||% 0.9) else 0

		if (verbose) cli::cli_inform("Using provided initial values")
	} else {
		lr <- init_lowrank(Y, r)
		U <- lr$U
		alpha <- lr$alpha
		tau_alpha2 <- 1
		rho_alpha <- if (ar1_alpha) 0.9 else 0
		sigma2_proc <- 0.5
		sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
		g2 <- 1

		Barray <- array(0, dim = c(n_col, n_col, Tt))
		for (t in 1:Tt) Barray[, , t] <- diag(n_col)
		tau_B2 <- 1
		rho_B <- if (ar1_B) 0.9 else 0
	}
	####

	# MCMC storage
	n_iter <- burn + nscan
	keep_idx <- ((burn + 1):n_iter)[((burn + 1):n_iter) %% odens == 0]
	S <- length(keep_idx)

	Usave <- vector("list", S)
	alphasave <- vector("list", S)
	Bsave <- vector("list", S)
	Msave <- vector("list", S)
	Thetasave <- vector("list", S)
	sigmasave <- numeric(S)
	tau_alpha_save <- numeric(S)
	tau_B_save <- numeric(S)
	g2save <- numeric(S)
	if (FAM$name == "gaussian") sigma2_obs_save <- numeric(S)
	if (ar1_alpha && update_rho_alpha) rho_alpha_save <- numeric(S)
	if (ar1_B && update_rho_B) rho_B_save <- numeric(S)

	time_keep <- seq(1, Tt, by = time_thin)

	if (verbose) {
		cli::cli_alert_info("Running low-rank DBN MCMC ({n_iter} iterations)")
		cli::cli_progress_bar("MCMC iterations", total = n_iter)
	}

	s <- 0
	epsilon <- 0.005
	accept_U <- 0

	a_tau <- 10
	b_tau <- 10

	Aarray <- compute_all_A_lowrank(U, alpha, Tt)

	# flatten arrays for batch operations
	total_slices <- p * Tt
	Z_flat <- array(0, dim = c(n_row, n_col, total_slices))
	Theta_flat <- array(0, dim = c(n_row, n_col, total_slices))
	M_flat <- array(0, dim = c(n_row, n_col, p))

	for (j in 1:p) {
		M_flat[, , j] <- pre$M[, , j]
		for (t in 1:Tt) {
			idx <- (j-1)*Tt + t
			Z_flat[, , idx] <- pre$Z[, , j, t]
			Theta_flat[, , idx] <- pre$Theta[, , j, t]
		}
	}
	####

	# main MCMC loop
	for (g in seq_len(n_iter)) {
		# update Z for ordinal/binary
		if (FAM$name == "ordinal") {
			pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "ordinal")
			tmp <- aperm(pre$Z, c(1L, 2L, 4L, 3L))
			dim(tmp) <- c(n_row, n_col, p * Tt)
			Z_flat[] <- tmp
		} else if (FAM$name == "binary") {
			pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "binary")
			tmp <- aperm(pre$Z, c(1L, 2L, 4L, 3L))
			dim(tmp) <- c(n_row, n_col, p * Tt)
			Z_flat[] <- tmp
		}

		# M update (vectorized on flat arrays)
		for (j in 1:p) {
			offset <- (j-1) * Tt
			resid_sum <- Z_flat[, , offset + seq_len(Tt), drop = FALSE] -
						 Theta_flat[, , offset + seq_len(Tt), drop = FALSE]
			mu_sum <- rowSums(resid_sum, dims = 2L)
			mu_var <- 1 / (Tt + 1 / g2)
			M_flat[, , j] <- mu_sum * mu_var + matrix(sqrt(mu_var) * rnorm(nc), n_row, n_col)
		}
		M_ss <- sum(M_flat^2)
		g2 <- safe_rinv_gamma(10 + length(M_flat)/2, 10 + M_ss/2)

		# Theta update via blocked FFBS
		Theta_flat <- ffbs_theta_blocked(Z_flat, M_flat, Aarray, Barray,
										sigma2_proc, sigma2_obs, m, p, Tt)

		# alpha update (optimized)
		alpha <- update_alpha_optimized(Theta_flat, U, Barray,
									   sigma2_proc,
									   tau_alpha2,
									   ar1_alpha, rho_alpha, m, p, Tt, r)

		Aarray <- compute_all_A_lowrank(U, alpha, Tt)

		# B update (parallelized)
		Barray <- update_B_parallel(Theta_flat, Aarray,
									sigma2_proc,
									tau_B2,
									ar1_B, rho_B, m, p, Tt)

		# update tau_alpha2
		if (ar1_alpha) {
			innov <- alpha[, -1, drop = FALSE] - rho_alpha * alpha[, -Tt, drop = FALSE]
		} else {
			innov <- alpha[, -1, drop = FALSE] - alpha[, -Tt, drop = FALSE]
		}
		innov_ss <- sum(innov^2)
		tau_alpha2 <- safe_rinv_gamma(a_tau + 0.5 * r * (Tt - 1), b_tau + 0.5 * innov_ss)

		# rho_alpha MH step (optional)
		if (ar1_alpha && update_rho_alpha) {
			rho_prop <- rho_alpha + rnorm(1, 0, 0.05)
			if (abs(rho_prop) < 1) {
				ll_prop <- -sum((alpha[, -1] - rho_prop * alpha[, -Tt])^2) / (2 * tau_alpha2)
				ll_curr <- -sum((alpha[, -1] - rho_alpha * alpha[, -Tt])^2) / (2 * tau_alpha2)
				if (log(runif(1)) < ll_prop - ll_curr) {
					rho_alpha <- rho_prop
				}
			}
		}

		# U via Stiefel manifold Metropolis step
		u_update_freq <- max(5, floor(m / 10))
		if (g %% u_update_freq == 0) {
			W <- generate_skew_proposal(m, sqrt(m) * epsilon)

			# Cayley transform proposal
			U_prop <- cayley_transform(U, W, epsilon)

			# re-orthogonalize if numerical drift exceeds tolerance
			orth_err <- norm(crossprod(U_prop) - diag(ncol(U_prop)), "F")
			if (orth_err > 1e-10) {
				U_prop <- qr.Q(qr(U_prop))
			}

			# average Theta across relations for the likelihood
			Theta_avg <- array(0, dim = c(n_row, n_col, Tt))
			if (p == 1) {
				Theta_avg <- Theta_flat[, , 1:Tt]
			} else {
				for (t in 1:Tt) {
					for (j in 1:p) {
						Theta_avg[, , t] <- Theta_avg[, , t] + Theta_flat[, , (j-1)*Tt + t]
					}
					Theta_avg[, , t] <- Theta_avg[, , t] / p
				}
			}

			log_acc <- loglik_U(U_prop, alpha, Theta_avg, Barray, sigma2_proc) -
					   loglik_U(U, alpha, Theta_avg, Barray, sigma2_proc)

			if (!is.na(log_acc) && is.finite(log_acc) && log(runif(1)) < log_acc) {
				U <- U_prop
				accept_U <- accept_U + 1
				Aarray <- compute_all_A_lowrank(U, alpha, Tt)
			}
		}

		# update tau_B2
		if (ar1_B) {
			innovB <- Barray[, , 2:Tt] - rho_B * Barray[, , 1:(Tt - 1)]
		} else {
			innovB <- Barray[, , 2:Tt] - Barray[, , 1:(Tt - 1)]
		}
		innovB_ss <- sum(innovB^2)
		tau_B2 <- safe_rinv_gamma(a_tau + n_col * n_col * (Tt - 1) / 2, b_tau + innovB_ss / 2)

		# rho_B MH step (optional)
		if (ar1_B && update_rho_B) {
			rho_prop <- rho_B + rnorm(1, 0, 0.05)
			if (abs(rho_prop) < 1) {
				innovB_prop <- Barray[, , 2:Tt] - rho_prop * Barray[, , 1:(Tt - 1)]
				ll_prop <- -sum(innovB_prop^2) / (2 * tau_B2)
				innovB_curr <- Barray[, , 2:Tt] - rho_B * Barray[, , 1:(Tt - 1)]
				ll_curr <- -sum(innovB_curr^2) / (2 * tau_B2)
				if (log(runif(1)) < ll_prop - ll_curr) {
					rho_B <- rho_prop
				}
			}
		}

		# sigma2_proc via SIMD residuals
		resid_ss <- compute_sigma2_simd(Theta_flat, U, alpha, Barray, m, p, Tt, r)
		sigma2_proc <- safe_rinv_gamma(a_tau + nc * p * (Tt - 1) / 2, b_tau + resid_ss / 2)

		# observation variance for gaussian
		if (FAM$name == "gaussian") {
			resid_obs <- compute_gaussian_obs_residuals_flat(Z_flat, Theta_flat, M_flat, m, p, Tt)
			sigma2_obs <- safe_rinv_gamma(a_tau + nc * p * Tt / 2, b_tau + resid_obs / 2)
		}

		# store draws
		if (g %in% keep_idx) {
			s <- s + 1
			Usave[[s]] <- U
			alphasave[[s]] <- alpha[, time_keep, drop = FALSE]
			Bsave[[s]] <- Barray[, , time_keep, drop = FALSE]
			sigmasave[s] <- sigma2_proc
			tau_alpha_save[s] <- tau_alpha2
			tau_B_save[s] <- tau_B2
			g2save[s] <- g2

			if (FAM$name == "gaussian") sigma2_obs_save[s] <- sigma2_obs
			if (ar1_alpha && update_rho_alpha) rho_alpha_save[s] <- rho_alpha
			if (ar1_B && update_rho_B) rho_B_save[s] <- rho_B

			Msave[[s]] <- M_flat
			Theta_4d <- array(NA, c(n_row, n_col, p, length(time_keep)))
			for (j_save in 1:p) {
				for (i_save in seq_along(time_keep)) {
					idx_save <- (j_save - 1) * Tt + time_keep[i_save]
					Theta_4d[, , j_save, i_save] <- Theta_flat[, , idx_save]
				}
			}
			Thetasave[[s]] <- Theta_4d
		}

		# adapt step size for U proposal
		if (g %% (u_update_freq * 20) == 0) {
			accept_rate <- accept_U / 20
			if (accept_rate < 0.2) {
				epsilon <- max(epsilon * 0.9, 1e-4)
			} else if (accept_rate > 0.3) {
				epsilon <- min(epsilon * 1.1, 0.1)
			}
			accept_U <- 0
		}

		if (verbose) cli::cli_progress_update()

		if (is.numeric(verbose) && verbose > 0 && g %% verbose == 0) {
			msg <- "iter {g}: sigma^2={round(sigma2_proc, 3)} tauA^2={round(tau_alpha2, 3)} tauB^2={round(tau_B2, 3)}"
			if (ar1_alpha && update_rho_alpha) msg <- paste0(msg, " rhoA={round(rho_alpha, 3)}")
			if (ar1_B && update_rho_B) msg <- paste0(msg, " rhoB={round(rho_B, 3)}")
			cli::cli_inform(msg)
		}
	}
	####

	# assemble output
	if (verbose) cli::cli_progress_done()

	# unflatten arrays back to 4D
	for (j in 1:p) {
		pre$M[, , j] <- M_flat[, , j]
		for (t in 1:Tt) {
			idx <- (j-1)*Tt + t
			pre$Z[, , j, t] <- Z_flat[, , idx]
			pre$Theta[, , j, t] <- Theta_flat[, , idx]
		}
	}

	# reconstruct A arrays from U and alpha draws
	Asave <- vector("list", S)
	for (s in 1:S) {
		A_s <- array(0, dim = c(n_row, n_row, length(time_keep)))
		for (i in seq_along(time_keep)) {
			alpha_i <- alphasave[[s]][, i, drop = TRUE]
			A_s[, , i] <- Usave[[s]] %*% (alpha_i * t(Usave[[s]]))
		}
		Asave[[s]] <- A_s
	}

	pars_df <- data.frame(
		sigma2_proc = sigmasave,
		tau_alpha2 = tau_alpha_save,
		tau_B2 = tau_B_save,
		g2 = g2save
	)

	if (FAM$name == "gaussian") pars_df$sigma2_obs <- sigma2_obs_save
	if (ar1_alpha && update_rho_alpha) pars_df$rho_alpha <- rho_alpha_save
	if (ar1_B && update_rho_B) pars_df$rho_B <- rho_B_save

	out <- structure(
		list(
			draws = list(
				theta = Thetasave,
				z = NULL,
				pars = pars_df,
				misc = list(
					U = Usave,
					alpha = alphasave,
					A = Asave,
					B = Bsave,
					M = Msave
				)
			),
			meta = list(
				dims = dims,
				draws = S,
				thin = odens,
				time_thin = time_thin,
				time_kept = time_keep,
				model = "lowrank"
			),
			family = FAM,
			model = "lowrank",
			U = Usave,
			alpha = alphasave,
			A = Asave,
			B = Bsave,
			M = Msave,
			sigma2_proc = sigmasave,
			tau_alpha2 = tau_alpha_save,
			tau_B2 = tau_B_save,
			g2 = g2save,
			Y = Y,
			dims = dims,
			settings = list(
				r = r,
				nscan = nscan,
				burn = burn,
				odens = odens,
				ar1_alpha = ar1_alpha,
				ar1_B = ar1_B,
				time_thin = time_thin,
				family = family
			)
		),
		class = "dbn"
	)

	if (FAM$name == "gaussian") {
		out$sigma2_obs <- sigma2_obs_save
	}

	if (ar1_alpha && update_rho_alpha) {
		out$rho_alpha <- rho_alpha_save
	}
	if (ar1_B && update_rho_B) {
		out$rho_B <- rho_B_save
	}

	out
}
####
