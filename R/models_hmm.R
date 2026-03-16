####
#' DBN HMM Model
#'
#' @description Regime-switching DBN for large-scale networks (100+ actors, 50+ time steps)
#' @param Y Data array (nodes x nodes x relations x time)
#' @param R Number of regimes
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param delta Dirichlet prior for transition matrix
#' @param seed Random seed
#' @param verbose Print progress
#' @param time_thin Save every nth time point, when 1, save all time points
#' @param family Data family (ordinal, gaussian, or binary)
#' @param previous Previous dbn_hmm results to continue from (optional)
#' @param init List of initial values: S, A_list, B_list, Pi, sigma2_proc, tau_A2, tau_B2, g2, pi0 (optional)
#' @param symmetric Logical. If TRUE, enforce B = A for each regime. Default: FALSE.
#' @param ... Additional arguments
#' @return List containing MCMC results
#' @seealso \code{\link{dbn}} for the main dispatcher, \code{\link{param_summary}} for posterior summaries
#' @examples
#' \donttest{
#' sim <- simulate_hmm_dbn(n = 8, time = 10, R = 2, seed = 1)
#' fit <- dbn_hmm(sim$Y, R = 2, nscan = 200, burn = 100, verbose = FALSE)
#' }
#' @export
####
dbn_hmm <- function(Y,
					family = c("ordinal", "gaussian", "binary"),
					R = 3,
					nscan = 10000,
					burn = 1000,
					odens = 1,
					delta = rep(1, R),
					seed = 6886,
					verbose = TRUE,
					time_thin = 1,
					previous = NULL,
					init = NULL,
					symmetric = FALSE,
					...) {

	# family setup

	family <- match.arg(family)
	if (isTRUE(verbose)) verbose <- 100L
	if (isFALSE(verbose)) verbose <- 0L
	FAM <- switch(family,
		ordinal = family_ordinal(),
		gaussian = family_gaussian(),
		binary = family_binary()
	)

	# need at least 2 time points
	if (dim(Y)[4] < 2) {
		cli::cli_abort("HMM model requires at least 2 time points. Use {.code model = \"static\"} for cross-sectional data.")
	}

	# check for zero variance and infinite values
	check_edge_cases_hmm <- function(Y) {
		has_zero_var <- TRUE
		has_inf <- FALSE
		n_inf <- 0

		for (rel in 1:dim(Y)[3]) {
			Y_rel <- Y[, , rel, , drop = FALSE]
			if (length(dim(Y_rel)) == 4) {
				Y_rel <- array(Y_rel, c(dim(Y)[1], dim(Y)[2], dim(Y)[4]))
			}
			edge_cases <- check_edge_cases(Y_rel)
			if (edge_cases$variance > 0) has_zero_var <- FALSE
			if (edge_cases$has_infinite) {
				has_inf <- TRUE
				n_inf <- n_inf + edge_cases$n_infinite
			}
		}

		list(
			variance = if (has_zero_var) 0 else 1,
			has_infinite = has_inf,
			n_infinite = n_inf,
			n_obs = prod(dim(Y))
		)
	}

	edge_cases <- check_edge_cases_hmm(Y)
	if (edge_cases$variance == 0) {
		cli::cli_abort(c(
			"Input data has zero variance",
			"i" = "All values in the data array are constant",
			"x" = "DBN models require variation in the data to estimate parameters"
		))
	}
	if (edge_cases$has_infinite) {
		cli::cli_abort(c(
			"Input data contains infinite values",
			"i" = "Found {edge_cases$n_infinite} infinite values out of {edge_cases$n_obs} total observations",
			"x" = "Please check your data for invalid entries"
		))
	}

	# preprocessing

	set.seed(seed)

	pre <- shared_preprocess(Y, family = family)
	dims <- pre$dims
	dims$is_symmetric <- symmetric
	n_row <- dims$n_row
	n_col <- dims$n_col
	is_bipartite <- dims$is_bipartite
	p <- dims$p
	Tt <- dims$Tt

	# bipartite networks not yet supported for HMM
	if (is_bipartite) {
		cli::cli_abort(c(
			"HMM model does not yet support bipartite (rectangular) networks.",
			"i" = "Your data has {n_row} senders and {n_col} receivers.",
			"i" = "Use {.code model = \"dynamic\"} for bipartite networks."
		))
	}

	# tuning based on problem size
	m_max <- max(n_row, n_col)
	large_scale <- (m_max > 50 || Tt > 100)
	use_beam_search <- (Tt > 200 && R > 5)
	beam_width <- if(use_beam_search) min(R, max(3, R/2)) else 0

	# total cells per time slice
	nc <- n_row * n_col

	# initialization

	if (!is.null(previous)) {
		# continue from a previous run
		if (!inherits(previous, "dbn") || previous$model != "hmm") {
			cli::cli_abort("previous must be results from dbn_hmm()")
		}

		last_idx <- length(previous$sigma2_proc)
		sigma2_proc <- previous$sigma2_proc[last_idx]
		tau_A2 <- previous$tau_A2[last_idx]
		tau_B2 <- previous$tau_B2[last_idx]
		g2 <- previous$g2[last_idx]
		S <- previous$S[[last_idx]]

		# expand states to full time grid if thinned
		if (length(S) < Tt) {
			S_full <- integer(Tt)
			prev_times <- seq(1, Tt, by = previous$settings$time_thin)
			for (t in 1:Tt) {
				nearest <- which.min(abs(prev_times - t))
				S_full[t] <- S[nearest]
			}
			S <- S_full
		}

		A_list <- previous$A[[last_idx]]
		B_list <- previous$B[[last_idx]]
		Pi <- previous$Pi[[last_idx]]
		pi0 <- if (!is.null(previous$pi0)) previous$pi0 else rep(1 / R, R)
		sigma2_obs <- if (FAM$name == "gaussian" && !is.null(previous$sigma2_obs)) {
			previous$sigma2_obs[last_idx]
		} else {
			FAM$init_pars$sigma2_obs %||% 1
		}

		if (verbose) cli::cli_inform("Continuing from previous run with sigma2_proc={round(sigma2_proc,3)}, tau_A2={round(tau_A2,3)}, tau_B2={round(tau_B2,3)}")
	} else if (!is.null(init)) {
		# user-supplied initial values
		S <- init$S %||% sample(R, Tt, replace = TRUE)
		Pi <- init$Pi %||% matrix(1/R, R, R)
		diag(Pi) <- diag(Pi) + 0.5
		Pi <- Pi / rowSums(Pi)
		pi0 <- init$pi0 %||% rep(1 / R, R)
		A_list <- init$A_list %||% replicate(R, diag(n_row) + matrix(rnorm(n_row * n_row, 0, 0.05), n_row, n_row), simplify = FALSE)
		B_list <- init$B_list %||% replicate(R, diag(n_col) + matrix(rnorm(n_col * n_col, 0, 0.05), n_col, n_col), simplify = FALSE)
		sigma2_proc <- init$sigma2_proc %||% 1
		sigma2_obs <- init$sigma2_obs %||% (FAM$init_pars$sigma2_obs %||% 1)
		tau_A2 <- init$tau_A2 %||% 1
		tau_B2 <- init$tau_B2 %||% 1
		g2 <- init$g2 %||% 1

		if (verbose) cli::cli_inform("Using initial values: sigma2_proc={round(sigma2_proc,3)}, tau_A2={round(tau_A2,3)}, tau_B2={round(tau_B2,3)}")
	} else {
		# default initialization, scaled to the data
		pi0 <- rep(1 / R, R)

		# k-means initialization when there are enough time points
		if (Tt >= R * 2) {
			n_init <- min(20, Tt)
			Y_early <- apply(Y[, , , 1:n_init, drop = FALSE], 4, function(x) c(x))
			Y_early_clean <- Y_early

			# impute NAs with row means
			for (i in 1:nrow(Y_early_clean)) {
				na_idx <- is.na(Y_early_clean[i, ])
				if (any(na_idx)) {
					row_mean <- mean(Y_early_clean[i, !na_idx])
					Y_early_clean[i, na_idx] <- if (is.finite(row_mean)) row_mean else 0
				}
			}

			if (any(is.finite(Y_early_clean))) {
				km <- tryCatch({
					kmeans(t(Y_early_clean), centers = R, nstart = 10, iter.max = 50)
				}, error = function(e) {
					list(cluster = sample(R, ncol(Y_early_clean), replace = TRUE))
				})

				S <- integer(Tt)
				S[1:n_init] <- km$cluster

				# fill remaining time points with persistence
				if (Tt > n_init) {
					for (t in (n_init + 1):Tt) {
						if (runif(1) < 0.85) {
							S[t] <- S[t-1]
						} else {
							S[t] <- sample(R, 1)
						}
					}
				}
			} else {
				S <- sample(R, Tt, replace = TRUE)
			}
		} else {
			S <- sample(R, Tt, replace = TRUE)
		}

		# transition matrix biased toward self-transitions
		Pi <- matrix(delta[1] / sum(delta), R, R)
		diag(Pi) <- diag(Pi) + 2
		Pi <- Pi / rowSums(Pi)

		# A and B matrices near identity
		A_list <- vector("list", R)
		B_list <- vector("list", R)

		Y_vals <- as.vector(Y[!is.na(Y)])
		Y_scale <- if(length(Y_vals) > 0) sd(Y_vals) else 1
		noise_scale <- min(0.1, Y_scale / (5 * sqrt(m_max)))

		for (r in 1:R) {
			A_list[[r]] <- diag(n_row) + matrix(rnorm(n_row * n_row, 0, noise_scale), n_row, n_row)
			B_list[[r]] <- diag(n_col) + matrix(rnorm(n_col * n_col, 0, noise_scale), n_col, n_col)

			A_list[[r]] <- stabilize_spectral_radius(A_list[[r]], 0.95)
			B_list[[r]] <- stabilize_spectral_radius(B_list[[r]], 0.95)
		}

		# initial variance parameters
		tau_A2 <- tau_B2 <- 1
		sigma2_proc <- max(0.01, Y_scale^2 / 100)
		sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
		g2 <- 1
	}

	# MCMC storage

	n_iter <- burn + nscan
	keep_idx <- ((burn + 1):n_iter)[((burn + 1):n_iter) %% odens == 0]
	n_keep <- length(keep_idx)

	# allocate output arrays
	Ssave <- vector("list", n_keep)
	Asave <- vector("list", n_keep)
	Bsave <- vector("list", n_keep)
	Pisave <- vector("list", n_keep)
	sigmasave <- numeric(n_keep)
	tau_A_save <- numeric(n_keep)
	tau_B_save <- numeric(n_keep)
	g2save <- numeric(n_keep)
	if (FAM$name == "gaussian") sigma2_obs_save <- numeric(n_keep)

	Msave <- vector("list", n_keep)
	Thetasave <- if(time_thin > 1 || m_max > 100) NULL else vector("list", n_keep)
	if (FAM$name %in% c("ordinal", "binary")) {
		Zsave <- if(time_thin > 1 || m_max > 100) NULL else vector("list", n_keep)
	}

	# which time points to keep
	time_keep <- seq(1, Tt, by = time_thin)

	# progress bar
	if (verbose) {
		cli::cli_alert_info("Running HMM DBN MCMC ({n_iter} iterations)")
		cli::cli_progress_bar("MCMC iterations", total = n_iter)
	}

	s <- 0

	# preallocated constants and work arrays
	I_nrow <- diag(n_row)
	I_ncol <- diag(n_col)
	Theta_avg <- array(0, dim = c(n_row, n_col, Tt))
	regime_counts <- matrix(0, n_iter, R)
	log_2pi_sigma2 <- -0.5 * nc * log(2 * pi)

	# main MCMC loop

	for (g in seq_len(n_iter)) {

		# update Z and M

		regime_arrays <- build_regime_arrays(S, A_list, B_list, n_row, Tt)
		Aarray <- regime_arrays$Aarray
		Barray <- regime_arrays$Barray

		# M update
		mu_var <- 1 / (pre$dims$Tt + 1 / g2)
		resid <- pre$Z - pre$Theta
		mu_hat <- mu_var * rowSums(resid, dims = 3L)
		pre$M <- mu_hat + sqrt(mu_var) * rsan(dim(pre$M))

		# Z update
		if (FAM$name == "ordinal") {
			pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "ordinal")
		} else if (FAM$name == "binary") {
			pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "binary")
		}

		# theta via ffbs

		for (rel in 1:p) {
			pre$Theta[, , rel, ] <- FAM$ffbs_wrapper(
				pre$Z[, , rel, ], pre$M[, , rel],
				Aarray, Barray,
				sigma2_proc = sigma2_proc,
				sigma2_obs = sigma2_obs
			)
		}

		# average Theta across relations for state updates
		if (p == 1L) {
			Theta_avg[] <- pre$Theta[, , 1, ]
		} else {
			Theta_avg[] <- rowMeans(aperm(pre$Theta, c(1L, 2L, 4L, 3L)), dims = 3L)
		}

		# forward-backward sample of regime sequence

		if (use_beam_search) {
			log_alpha <- forward_hmm_fast(Theta_avg, A_list, B_list, Pi, sigma2_proc, pi0, beam_width)
			S <- backward_sample_fast(log_alpha, Pi)
		} else {
			log_alpha <- forward_hmm(Theta_avg, A_list, B_list, Pi, sigma2_proc, pi0)
			S <- backward_sample(log_alpha, Pi)
		}

		# record regime occupancy
		regime_counts[g, ] <- tabulate(S, nbins = R)

		# update initial-state probabilities after burn-in
		if (g > burn) {
			pi0 <- (pi0 * (g - burn - 1) + as.numeric(S[1] == 1:R)) / (g - burn)
		}

		# transition matrix pi

		n_ij <- count_transitions(S, R)
		# diagonal pseudo-counts to prevent empty rows
		n_ij <- n_ij + diag(R)
		for (i in 1:R) {
			Pi[i, ] <- rdirichlet(delta + n_ij[i, ])
		}

		# regime-specific A, B

		for (r in 1:R) {
			regime_data <- collect_regime_thetas(Theta_avg, S, r, n_row)
			if (regime_data$n_obs > 0) {
				upd <- update_AB_static_cpp(
					regime_data$Th_prev, regime_data$Th_curr,
					B_list[[r]], tau_A2, tau_B2, sigma2_proc
				)
				A_list[[r]] <- upd$A
				B_list[[r]] <- upd$B
			} else {
				# empty regime: draw from prior
				A_list[[r]] <- diag(n_row) + matrix(rnorm(n_row * n_row, 0, sqrt(tau_A2 / 10)), n_row, n_row)
				B_list[[r]] <- diag(n_col) + matrix(rnorm(n_col * n_col, 0, sqrt(tau_B2 / 10)), n_col, n_col)
			}
		}
		if (symmetric) for (r in 1:R) B_list[[r]] <- A_list[[r]]

		# hyperparameters pooled across regimes

		# tau_A2
		A_resid <- compute_regime_residuals(A_list, I_nrow, R, n_row)
		n_active <- sum(regime_counts[g, ] > 0)
		shape0_A <- 10 + n_row * n_row
		rate0 <- 10
		tau_A2 <- safe_rinv_gamma(shape0_A + n_row * n_row * n_active / 2, rate0 + A_resid / 2)

		# tau_B2
		B_resid <- compute_regime_residuals(B_list, I_ncol, R, n_col)
		shape0_B <- 10 + n_col * n_col
		tau_B2 <- safe_rinv_gamma(shape0_B + n_col * n_col * n_active / 2, rate0 + B_resid / 2)
		if (symmetric) tau_B2 <- tau_A2

		# sigma2_proc

		shape0_proc <- 10 + nc
		if (large_scale) {
			resid <- compute_bilinear_residuals_fast(pre$Theta, Aarray, Barray, n_row, p, Tt)
		} else {
			Theta_flat <- array(pre$Theta, dim = c(n_row, n_col, p * Tt))
			resid <- compute_bilinear_residuals(Theta_flat, Aarray, Barray, n_row, p, Tt)
		}
		sigma2_proc <- safe_rinv_gamma(shape0_proc + nc * p * (Tt - 1) / 2, rate0 + resid / 2)

		# sigma2_obs for gaussian family
		if (FAM$name == "gaussian") {
			resid_obs <- compute_gaussian_obs_residuals(pre$Z, pre$Theta, pre$M)
			sigma2_obs <- safe_rinv_gamma(10 + nc + length(pre$Z) / 2, rate0 + resid_obs / 2)
		}

		# g2

		g2 <- safe_rinv_gamma(10 + nc + prod(dim(pre$M)) / 2, (rate0 + sum(pre$M^2)) / 2)

		# store this iteration

		if (g %in% keep_idx) {
			s <- s + 1
			Ssave[[s]] <- S[time_keep]
			Asave[[s]] <- A_list
			Bsave[[s]] <- B_list
			Pisave[[s]] <- Pi
			sigmasave[s] <- sigma2_proc
			tau_A_save[s] <- tau_A2
			tau_B_save[s] <- tau_B2
			g2save[s] <- g2

			if (FAM$name == "gaussian") {
				sigma2_obs_save[s] <- sigma2_obs
			}

			Msave[[s]] <- pre$M
			if (!is.null(Thetasave)) {
				if (time_thin > 1) {
					Thetasave[[s]] <- pre$Theta[, , , time_keep, drop = FALSE]
				} else {
					Thetasave[[s]] <- pre$Theta
				}
			}

			if (FAM$name %in% c("ordinal", "binary") && !is.null(Zsave)) {
				if (time_thin > 1) {
					Zsave[[s]] <- pre$Z[, , , time_keep, drop = FALSE]
				} else {
					Zsave[[s]] <- pre$Z
				}
			}
		}

		if (verbose > 0) cli::cli_progress_update()

		if (is.numeric(verbose) && verbose > 0 && g %% verbose == 0) {
			regime_table <- table(factor(S, levels = 1:R))
			msg <- "iter {g}: sigma^2={round(sigma2_proc, 3)} tauA^2={round(tau_A2, 3)} tauB^2={round(tau_B2, 3)}"
			msg <- paste0(msg, " regimes=", paste(regime_table, collapse = "/"))
			cli::cli_inform(msg)
		}
	}

	# assemble output

	if (verbose) cli::cli_progress_done()

	# regime occupancy summary
	if (verbose) {
		cli::cli_h3("Regime occupancy summary")
		regime_props <- colMeans(regime_counts[(burn + 1):n_iter, ]) / Tt
		for (r in 1:R) {
			cli::cli_inform("  Regime {r}: {round(regime_props[r] * 100, 1)}% ({round(regime_props[r] * Tt)} time points on average)")
		}

		if (any(regime_props < 0.05)) {
			cli::cli_warn("Some regimes are rarely occupied (<5%). Consider reducing R.")
		}
	}

	# reshape A and B lists into arrays
	Asave_array <- vector("list", n_keep)
	Bsave_array <- vector("list", n_keep)

	for (s in 1:n_keep) {
		A_s <- array(0, dim = c(n_row, n_row, R))
		B_s <- array(0, dim = c(n_col, n_col, R))
		for (r in 1:R) {
			A_s[, , r] <- Asave[[s]][[r]]
			B_s[, , r] <- Bsave[[s]][[r]]
		}
		Asave_array[[s]] <- A_s
		Bsave_array[[s]] <- B_s
	}

	# scalar parameters as data frame
	pars_df <- data.frame(
		sigma2_proc = sigmasave,
		tau_A2 = tau_A_save,
		tau_B2 = tau_B_save,
		g2 = g2save
	)

	if (FAM$name == "gaussian") {
		pars_df$sigma2_obs <- sigma2_obs_save
	}

	# collected draws
	draws <- list(
		theta = Thetasave,
		z = if (FAM$name %in% c("ordinal", "binary")) Zsave else NULL,
		pars = pars_df,
		misc = list(
			M = Msave,
			S = Ssave,
			A = Asave_array,
			B = Bsave_array,
			Pi = Pisave
		)
	)

	# model metadata
	meta <- list(
		dims = dims,
		draws = n_keep,
		thin = odens,
		time_thin = time_thin,
		time_kept = time_keep,
		model = "hmm",
		R = R,
		pi0 = pi0
	)

	# final output list
	out <- list(
		draws = draws,
		meta = meta,
		family = FAM,
		model = "hmm",
		S = Ssave,
		A = Asave_array,
		B = Bsave_array,
		Pi = Pisave,
		pi0 = pi0,
		sigma2_proc = sigmasave,
		tau_A2 = tau_A_save,
		tau_B2 = tau_B_save,
		g2 = g2save,
		Y = Y,
		n_iter = n_iter,
		burn = burn,
		thin = odens,
		dims = dims,
		settings = list(
			R = R,
			nscan = nscan,
			burn = burn,
			odens = odens,
			delta = delta,
			time_thin = time_thin,
			family = family
		),
		regime_counts = regime_counts
	)

	if (FAM$name == "gaussian") {
		out$sigma2_obs <- sigma2_obs_save
	}

	# accumulate iteration counts for continued runs
	if (!is.null(previous)) {
		prev_iter <- if (!is.null(previous$total_iter)) {
			previous$total_iter
		} else if (!is.null(previous$settings$nscan)) {
			previous$settings$burn + previous$settings$nscan
		} else {
			previous$n_iter
		}
		out$total_iter <- prev_iter + nscan
		out$continued_from <- prev_iter
	}

	class(out) <- "dbn"
	out
}

####
#' Step Log-likelihood
#'
#' @description Computes log-likelihood for one time step given regime
#' @param A Sender effects matrix
#' @param B Receiver effects matrix
#' @param Theta_prev Previous Theta
#' @param Theta_curr Current Theta
#' @param sigma2_proc Innovation variance
#' @return Log-likelihood value
#' @keywords internal
####
step_ll <- function(A, B, Theta_prev, Theta_curr, sigma2_proc) {
	# clamp sigma2_proc away from zero
	sigma2_proc <- max(sigma2_proc, 1e-6)

	# residual from bilinear prediction A * Theta_{t-1} * B'
	resid <- Theta_curr - bilinear_step(A, Theta_prev, B)

	# Gaussian log-likelihood
	-0.5 * sum(resid^2) / sigma2_proc - 0.5 * length(resid) * log(2 * pi * sigma2_proc)
}

####
#' Sample Dirichlet
#'
#' @description Samples from Dirichlet distribution
#' @param alpha Parameter vector
#' @return Sample vector
#' @keywords internal
####
rdirichlet <- function(alpha) {
	x <- rgamma(length(alpha), alpha, 1)
	x / sum(x)
}
