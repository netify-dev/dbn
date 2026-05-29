####
#' Block Utilities for Piecewise DBN Models
#'
#' @description parsing and validation utilities for piecewise-static models
#' @name utils_blocks
#' @keywords internal
NULL
####

####
#' Parse Blocks Parameter
#'
#' @description converts user-specified blocks into standardized format
#' @param blocks block specification (integer, vector, or "auto")
#' @param Tt total number of time points
#' @param n number of actors (for minimum block length)
#' @return list with K, boundaries, lengths, block_indices, names
#' @keywords internal
parse_blocks <- function(blocks, Tt, n) {

	# minimum block length depends on network size
	min_block_length <- get_min_block_length(n)

	if (is.numeric(blocks) && length(blocks) == 1) {
		# integer: number of equal-sized blocks
		K <- as.integer(blocks)
		if (K < 1) {
			cli::cli_abort("Number of blocks must be at least 1.")
		}
		if (K > Tt) {
			cli::cli_abort("Number of blocks ({K}) cannot exceed time points ({Tt}).")
		}
		boundaries <- round(seq(0, Tt, length.out = K + 1))
		block_names <- paste0("block_", 1:K)

	} else if (is.numeric(blocks) && length(blocks) > 1) {
		# vector: explicit boundaries
		boundaries <- sort(unique(c(0, blocks)))
		if (max(boundaries) < Tt) {
			boundaries <- c(boundaries, Tt)
		}
		# remove 0 if it was duplicated
		boundaries <- unique(boundaries)
		K <- length(boundaries) - 1

		# use names if provided
		if (!is.null(names(blocks))) {
			block_names <- names(blocks)
			if (length(block_names) < K) {
				block_names <- c(block_names, paste0("block_", (length(block_names) + 1):K))
			}
		} else {
			block_names <- paste0("block_", 1:K)
		}

	} else if (identical(blocks, "auto")) {
		# auto: handled by caller before this function
		cli::cli_abort("blocks = 'auto' should be resolved before calling parse_blocks")

	} else {
		cli::cli_abort(c(
			"Invalid blocks specification.",
			"i" = "Use an integer (number of blocks), a numeric vector (boundaries), or 'auto'."
		))
	}

	# compute block lengths
	lengths <- diff(boundaries)

	# validate minimum block lengths. blocks shorter than 2 cannot identify
	# the per-block operator at all (no transition pair); abort directively
	# instead of letting the C++ block update return a singular matrix.
	infeasible <- which(lengths < 2L)
	if (length(infeasible) > 0L) {
		cli::cli_abort(c(
			"Block specification yields blocks of length < 2.",
			"x" = "Blocks {paste(infeasible, collapse = ', ')} have length {paste(lengths[infeasible], collapse = ', ')} (need at least 2 for a transition).",
			"i" = "For {.code Tt = {Tt}}, the maximum number of equal-length blocks is {.val {Tt %/% 2}}.",
			"i" = "Pass a smaller integer K or a breakpoint vector with at least 2 time points per block."
		))
	}
	short_blocks <- which(lengths < min_block_length)
	if (length(short_blocks) > 0) {
		cli::cli_warn(c(
			"Some blocks shorter than recommended minimum ({min_block_length}).",
			"i" = "Short blocks: {paste(short_blocks, collapse = ', ')} with lengths {paste(lengths[short_blocks], collapse = ', ')}.",
			"i" = "This affects precision of actor-specific influence estimates (A[i,j]).",
			"i" = "Detecting overall structural change (compare_blocks) is typically robust to short blocks.",
			"i" = "Consider longer time spans or fewer actors if precise individual estimates are needed."
		))
	}

	# compute time indices for each block
	block_indices <- lapply(1:K, function(k) {
		(boundaries[k] + 1):boundaries[k + 1]
	})
	names(block_indices) <- block_names

	list(
		K = K,
		boundaries = boundaries,
		lengths = lengths,
		block_indices = block_indices,
		names = block_names
	)
}
####

####
#' Get Minimum Block Length
#'
#' @description determines minimum block length based on network size
#' @param n number of actors
#' @return minimum recommended block length
#' @keywords internal
get_min_block_length <- function(n) {
	if (n <= 20) {
		5L
	} else if (n <= 50) {
		8L
	} else if (n <= 100) {
		10L
	} else {
		15L
	}
}
####

####
#' Validate Block Structure
#'
#' @description checks that block structure is valid for data
#' @param block_info parsed block information
#' @param Tt total time points
#' @param verbose show validation messages
#' @return TRUE if valid, aborts otherwise
#' @keywords internal
validate_blocks <- function(block_info, Tt, verbose = TRUE) {
	# check boundaries cover full time range
	if (min(block_info$boundaries) != 0) {
		cli::cli_abort("Block boundaries must start at 0.")
	}
	if (max(block_info$boundaries) != Tt) {
		cli::cli_abort("Block boundaries must end at {Tt} (total time points).")
	}

	# check no empty blocks
	if (any(block_info$lengths <= 0)) {
		cli::cli_abort("All blocks must have positive length.")
	}

	# check total coverage
	total_covered <- sum(block_info$lengths)
	if (total_covered != Tt) {
		cli::cli_abort("Block lengths ({total_covered}) must sum to total time ({Tt}).")
	}

	if (verbose) {
		cli::cli_alert_success("Block structure validated: {block_info$K} blocks")
	}

	TRUE
}
####

####
#' Structural Stability Scan
#'
#' @description rolling window analysis to detect structural changes
#' @param Y data array (n_row x n_col x p x Tt)
#' @param family distribution family
#' @param window_size rolling window size
#' @param step_size step between windows
#' @param nscan MCMC iterations per window
#' @param verbose show progress
#' @return list with change magnitudes and suggested boundaries
#' @keywords internal
structural_stability <- function(Y, family = "ordinal",
								 window_size = NULL, step_size = NULL,
								 nscan = 500, verbose = TRUE) {
	dims <- dim(Y)
	n_row <- dims[1]
	n_col <- dims[2]
	p <- dims[3]
	Tt <- dims[4]

	# default window and step sizes
	if (is.null(window_size)) {
		min_block <- get_min_block_length(n_row)
		window_size <- max(min_block, ceiling(Tt / 10))
	}
	if (is.null(step_size)) {
		step_size <- max(1L, window_size %/% 2)
	}

	# window start positions
	starts <- seq(1, Tt - window_size + 1, by = step_size)
	n_windows <- length(starts)

	if (n_windows < 2) {
		cli::cli_warn("Too few windows for stability analysis. Using single block.")
		return(list(
			change_magnitude = numeric(0),
			window_centers = numeric(0),
			suggested_boundaries = integer(0)
		))
	}

	if (verbose) {
		cli::cli_alert_info("Scanning {n_windows} windows (size {window_size}, step {step_size})")
	}

	# fit static model in each window
	B_windows <- list()
	for (i in seq_along(starts)) {
		t_start <- starts[i]
		t_end <- min(t_start + window_size - 1, Tt)
		Y_window <- Y[, , , t_start:t_end, drop = FALSE]

		if (verbose) {
			cli::cli_progress_step("Window {i}/{n_windows}: t = {t_start}:{t_end}")
		}

		# fit a quick static model. A single window failing is tolerable
		# (we fall back to identity for that window's B), but a systematic
		# failure makes the breakpoint detector unreliable -- so we surface
		# the error class via cli_warn rather than swallowing silently.
		fit_window <- tryCatch({
			dbn_static(Y_window, family = family,
					   nscan = nscan, burn = nscan %/% 2, odens = 1,
					   verbose = FALSE, seed = 12345 + i)
		}, error = function(e) {
			cli::cli_warn(c(
				"Window {i}/{n_windows} (t = {t_start}:{t_end}) static fit failed.",
				"x" = "{conditionMessage(e)}",
				"i" = "Using identity B for this window; breakpoint detection may be unreliable if this happens for many windows."
			))
			NULL
		})

		if (!is.null(fit_window) && !is.null(fit_window$B)) {
			# extract posterior mean of B[[1]]
			B_windows[[i]] <- apply(fit_window$B[[1]], c(1, 2), mean)
		} else {
			B_windows[[i]] <- diag(n_row)
		}
	}

	# compute change magnitude between adjacent windows
	change_magnitude <- numeric(n_windows - 1)
	window_centers <- numeric(n_windows - 1)

	for (i in 1:(n_windows - 1)) {
		delta_B <- B_windows[[i + 1]] - B_windows[[i]]
		change_magnitude[i] <- sqrt(sum(delta_B^2))
		window_centers[i] <- (starts[i] + starts[i + 1]) / 2 + window_size / 2
	}

	# identify peaks (potential change points)
	if (length(change_magnitude) >= 3) {
		# find local maxima
		peaks <- which(diff(sign(diff(change_magnitude))) == -2) + 1

		# filter by threshold (above median)
		threshold <- median(change_magnitude)
		significant_peaks <- peaks[change_magnitude[peaks] > threshold]

		# convert to time boundaries
		suggested_boundaries <- round(window_centers[significant_peaks])
		suggested_boundaries <- suggested_boundaries[suggested_boundaries > 0 & suggested_boundaries < Tt]
	} else {
		suggested_boundaries <- integer(0)
	}

	if (verbose) {
		cli::cli_progress_done()
		if (length(suggested_boundaries) > 0) {
			cli::cli_alert_info("Detected potential change points at t = {paste(suggested_boundaries, collapse = ', ')}")
		} else {
			cli::cli_alert_info("No significant structural changes detected")
		}
	}

	list(
		change_magnitude = change_magnitude,
		window_centers = window_centers,
		suggested_boundaries = suggested_boundaries,
		window_size = window_size,
		B_windows = B_windows
	)
}
####

####
#' Automatic Block Selection
#'
#' @description two-stage algorithm for selecting number of blocks
#' @param Y data array (n_row x n_col x p x Tt)
#' @param family distribution family
#' @param K_min minimum blocks to consider
#' @param K_max maximum blocks to consider
#' @param nscan_stage1 MCMC iterations for stage 1
#' @param nscan_stage2 MCMC iterations for stage 2
#' @param verbose show progress
#' @param ... additional arguments passed to dbn_piecewise
#' @return list with selected_K, selected_boundaries, comparison results
#' @keywords internal
select_K_auto <- function(Y, family = "ordinal",
						  K_min = 1L, K_max = NULL,
						  nscan_stage1 = 500, nscan_stage2 = 2000,
						  verbose = TRUE, ...) {
	dims <- dim(Y)
	n <- dims[1]
	Tt <- dims[4]

	# determine K_max based on data
	min_block <- get_min_block_length(n)
	if (is.null(K_max)) {
		K_max <- max(2L, min(8L, floor(Tt / min_block)))
	}

	if (verbose) {
		cli::cli_h3("Automatic Block Selection")
		cli::cli_alert_info("K range: {K_min} to {K_max}")
	}

	# stage 1: quick structural stability scan
	if (verbose) cli::cli_h3("Stage 1: Structural Stability Scan")

	stab_result <- structural_stability(Y, family = family,
										nscan = nscan_stage1,
										verbose = verbose)

	# suggested K based on detected change points
	K_suggested <- length(stab_result$suggested_boundaries) + 1
	K_suggested <- max(K_min, min(K_max, K_suggested))

	if (verbose) {
		cli::cli_alert_info("Stage 1 suggests K = {K_suggested}")
	}

	# stage 2: formal model comparison around suggested K
	if (verbose) cli::cli_h3("Stage 2: Model Comparison via WAIC")

	K_range <- max(K_min, K_suggested - 1):min(K_max, K_suggested + 1)

	fits <- list()
	waics <- numeric(length(K_range))

	for (i in seq_along(K_range)) {
		K <- K_range[i]
		if (verbose) cli::cli_alert("Fitting K = {K}...")

		if (K == 1) {
			# use static model for K=1
			fits[[i]] <- dbn_static(Y, family = family,
									nscan = nscan_stage2,
									burn = nscan_stage2 %/% 2,
									verbose = FALSE, ...)
		} else {
			# use piecewise model
			fits[[i]] <- dbn_piecewise(Y, family = family,
									   blocks = K,
									   nscan = nscan_stage2,
									   burn = nscan_stage2 %/% 2,
									   verbose = FALSE, ...)
		}

		# compute WAIC
		waics[i] <- compute_waic_dbn(fits[[i]])$waic
	}

	# select best K
	best_idx <- which.min(waics)
	selected_K <- K_range[best_idx]

	if (verbose) {
		cli::cli_alert_success("Selected K = {selected_K} (WAIC = {round(waics[best_idx], 2)})")
		for (i in seq_along(K_range)) {
			marker <- if (i == best_idx) "*" else " "
			cli::cli_alert("{marker} K = {K_range[i]}: WAIC = {round(waics[i], 2)}")
		}
	}

	# determine boundaries for selected K
	if (selected_K == 1) {
		selected_boundaries <- c(Tt)
	} else if (selected_K == K_suggested && length(stab_result$suggested_boundaries) > 0) {
		# use stability-detected boundaries
		selected_boundaries <- c(sort(stab_result$suggested_boundaries)[1:(selected_K - 1)], Tt)
	} else {
		# equal-sized blocks
		selected_boundaries <- round(seq(Tt / selected_K, Tt, length.out = selected_K))
	}

	list(
		selected_K = selected_K,
		selected_boundaries = selected_boundaries,
		K_range = K_range,
		waics = waics,
		stage1_suggestion = K_suggested,
		stage1_boundaries = stab_result$suggested_boundaries,
		stability_result = stab_result
	)
}
####

####
#' Compute WAIC for a dbn fit
#'
#' Computes the Watanabe-Akaike Information Criterion (WAIC) for model
#' comparison. Returns a list compatible with `loo::loo_compare()`: the
#' `estimates` slot has rows `elpd_waic`, `p_waic`, `waic` with `Estimate`
#' and `SE` columns, and the object is `class = c("waic", "loo")`.
#'
#' @param fit A `dbn` fit object.
#' @param model_name Optional character label attached to the returned WAIC
#'   object so [loo::loo_compare()] can distinguish competing fits in its
#'   output table.
#' @return A list with `estimates`, `pointwise`, `lppd`, `p_waic`, `waic`,
#'   and `class = c("waic", "loo")` so `loo::loo_compare(w1, w2)` dispatches.
#' @author Tosin Salau and Shahryar Minhas
#' @export
compute_waic_dbn <- function(fit, model_name = NULL) {
	# accept model_name = so loo_compare rows are distinguishable; two
	# competing specs would otherwise both label as their model type.
	if (!is.null(model_name)) {
		if (length(model_name) != 1L || !is.character(model_name)) {
			cli::cli_abort("{.arg model_name} must be a single string.")
		}
	}
	# ALS fits don't have a posterior over parameters:
	# - without bootstrap, n_draws = 1 -> p_waic = 0 and WAIC reduces to
	#   -2 * lppd, which trivially "wins" against any MCMC competitor.
	# - with bootstrap, the variance is empirical-Bayes flavoured and
	#   under-counts model complexity.
	# refuse outright and point users at the right tool.
	su <- fit$meta$sampler_used %||% ""
	is_als <- su %in% c("als", "als_tv", "als_piecewise")
	if (is_als) {
		cli::cli_abort(c(
			"WAIC is not a valid criterion for ALS fits.",
			"i" = "ALS produces a point estimate (no posterior over parameters), so {.code p_waic = 0} and the WAIC reduces to {.code -2*lppd}, which trivially favors ALS over any MCMC competitor.",
			"i" = "For ALS-vs-MCMC comparison, use {.fn compute_irf} stability across resamples or refit the candidate with MCMC ({.fn dbn}, {.code model = \"dynamic\"})."
		))
	}
	# extract relevant info based on model type
	if (fit$model == "static") {
		waic <- compute_waic_static(fit)
	} else if (fit$model == "piecewise") {
		waic <- compute_waic_piecewise(fit)
	} else if (fit$model == "dynamic") {
		waic <- compute_waic_dynamic(fit)
	} else {
		# fallback: use DIC-based approximation
		if (!is.null(fit$diagnostics$dic) && !is.na(fit$diagnostics$dic)) {
			waic <- list(waic = fit$diagnostics$dic, lppd = NA_real_, p_waic = NA_real_)
		} else {
			cli::cli_warn("WAIC not implemented for model type {fit$model}. Using NA.")
			waic <- list(waic = NA_real_, lppd = NA_real_, p_waic = NA_real_)
		}
	}
	# repackage as a loo-shaped object (estimates matrix + class) so that
	# loo::loo_compare() dispatches.
	lppd_v <- as.numeric(waic$lppd %||% NA_real_)
	pwaic_v <- as.numeric(waic$p_waic %||% NA_real_)
	waic_v <- as.numeric(waic$waic %||% NA_real_)
	elpd_waic <- lppd_v - pwaic_v
	# if pointwise contributions are available, compute SEs honestly
	pointwise <- waic$pointwise
	if (is.null(pointwise) || !is.matrix(pointwise)) {
		se_elpd <- NA_real_
		se_pwaic <- NA_real_
		se_waic <- NA_real_
	} else {
		n_obs <- nrow(pointwise)
		se_elpd  <- sqrt(n_obs * var(pointwise[, "elpd_waic"]))
		se_pwaic <- sqrt(n_obs * var(pointwise[, "p_waic"]))
		se_waic  <- 2 * se_elpd
	}
	# warn when p_waic dwarfs the number of observations. p_waic > N is a
	# red flag for under-mixed chains: the per-observation log-lik
	# variance is dominating the criterion, so WAIC differences are not
	# interpretable until the sampler is run longer.
	n_obs_total <- tryCatch({
		Y <- fit$Y
		if (is.null(Y)) NA_real_
		else sum(!is.na(Y))
	}, error = function(e) NA_real_)
	if (is.finite(n_obs_total) && is.finite(pwaic_v) && n_obs_total > 0 &&
	    pwaic_v > n_obs_total) {
		cli::cli_warn(c(
			"{.field p_waic} ({sprintf('%.3g', pwaic_v)}) is {sprintf('%.1fx', pwaic_v / n_obs_total)} the number of observations ({n_obs_total}).",
			"i" = "This typically signals the chain has not mixed long enough for stable WAIC, not an over-parameterized model.",
			"i" = "Increase {.arg nscan} (10x is a good starting heuristic) or check {.fn check_convergence}; do not interpret WAIC differences in this regime."
		))
	}
	estimates <- matrix(
		c(elpd_waic, se_elpd, pwaic_v, se_pwaic, waic_v, se_waic),
		nrow = 3, ncol = 2, byrow = TRUE,
		dimnames = list(c("elpd_waic", "p_waic", "waic"), c("Estimate", "SE"))
	)
	out <- list(
		estimates = estimates,
		pointwise = pointwise,
		lppd      = lppd_v,
		p_waic    = pwaic_v,
		waic      = waic_v,
		elpd_waic = elpd_waic,
		n_obs     = n_obs_total,
		model_name = model_name %||% fit$model %||% NA_character_
	)
	class(out) <- c("waic", "loo")
	out
}
####

####
#' Compute WAIC for Static Model
#'
#' @description WAIC computation for static DBN
#' @param fit static dbn fit
#' @return list with waic, lppd, p_waic
#' @keywords internal
#' @noRd
compute_waic_static <- function(fit) {
	Y <- fit$Y
	family <- fit$family
	n_draws <- fit$settings$draws

	# get M samples
	M_samples <- fit$draws$misc$M
	if (is.null(M_samples)) {
		return(list(waic = NA, lppd = NA, p_waic = NA))
	}

	# compute log-likelihood for each draw
	log_lik_samples <- matrix(NA, nrow = n_draws, ncol = prod(dim(Y)))

	for (s in seq_len(n_draws)) {
		M_s <- M_samples[[s]]
		s2_s <- fit$draws$pars$s2[s]

		# expand M to match Y dimensions
		M_expanded <- array(M_s, dim = dim(Y))

		if (family == "gaussian") {
			log_lik_samples[s, ] <- dnorm(c(Y), c(M_expanded), sqrt(s2_s), log = TRUE)
		} else if (family == "ordinal" || family == "binary") {
			# approximate with probit likelihood
			log_lik_samples[s, ] <- dnorm(c(Y), c(M_expanded), 1, log = TRUE)
		}
	}

	# handle NAs
	log_lik_samples[!is.finite(log_lik_samples)] <- NA

	# per-observation lppd and p_waic; loo::loo_compare needs the pointwise
	# matrix to compute SE on differences. drop cells with <2 finite draws.
	col_ok <- apply(log_lik_samples, 2, function(x) sum(is.finite(x)) >= 2L)
	pointwise <- matrix(NA_real_, nrow = sum(col_ok), ncol = 2,
		dimnames = list(NULL, c("elpd_waic", "p_waic")))
	if (sum(col_ok) > 0L) {
		ll_ok <- log_lik_samples[, col_ok, drop = FALSE]
		lse <- apply(ll_ok, 2, function(x) {
			x <- x[is.finite(x)]
			if (length(x) == 0L) return(NA_real_)
			log(mean(exp(x - max(x)))) + max(x)
		})
		var_per_col <- apply(ll_ok, 2, function(x) {
			x <- x[is.finite(x)]
			if (length(x) < 2L) return(0)
			var(x)
		})
		pointwise[, "elpd_waic"] <- lse - var_per_col
		pointwise[, "p_waic"]    <- var_per_col
	}
	lppd <- sum(apply(log_lik_samples, 2, function(x) {
		x <- x[is.finite(x)]
		if (length(x) == 0) return(NA)
		log(mean(exp(x - max(x)))) + max(x)
	}), na.rm = TRUE)
	p_waic <- sum(apply(log_lik_samples, 2, function(x) {
		x <- x[is.finite(x)]
		if (length(x) < 2) return(0)
		var(x)
	}), na.rm = TRUE)
	waic <- -2 * (lppd - p_waic)

	list(waic = waic, lppd = lppd, p_waic = p_waic, pointwise = pointwise)
}
####

####
#' Compute WAIC for Piecewise Model
#'
#' @description WAIC computation for piecewise DBN. Evaluates the
#'   per-observation Gaussian (or probit) log-likelihood under the
#'   block-specific predictive distribution
#'   `Y_t | Theta_t, M, sigma^2 ~ N(Theta_t + M, sigma^2)`, where `Theta_t`
#'   is governed by the block-specific operators `A_k(t), B_k(t)`. Iterates
#'   over posterior draws to compute lppd (mean of likelihood, log) and
#'   p_waic (variance of log-likelihood).
#' @param fit piecewise dbn fit
#' @return list with waic, lppd, p_waic
#' @keywords internal
#' @noRd
compute_waic_piecewise <- function(fit) {
	Y <- fit$Y
	family <- fit$family
	n_draws <- fit$settings$draws %||% 0L
	M_samples <- fit$draws$misc$M
	A_samples <- fit$draws$misc$A
	B_samples <- fit$draws$misc$B
	Theta_samples <- fit$draws$misc$Theta
	s2_vec <- fit$draws$pars$s2

	if (is.null(M_samples) || is.null(A_samples) || is.null(B_samples) ||
	    n_draws < 2L || is.null(s2_vec)) {
		# fall back to the static approximation if the piecewise draws are
		# missing; this typically only happens when keep dropped them
		return(compute_waic_static(fit))
	}

	dims <- fit$dims
	n_row <- dims$n_row; n_col <- dims$n_col; p <- dims$p %||% 1L; Tt <- dims$Tt

	# block_indices[[k]] is the vector of time indices assigned to block k
	block_indices <- fit$blocks$block_indices
	K_blocks <- length(block_indices)
	t2k <- integer(Tt)
	for (k in seq_len(K_blocks)) t2k[block_indices[[k]]] <- k

	# subsample draws to keep the cost manageable
	S_max <- 100L
	s_idx <- if (n_draws > S_max) {
		as.integer(round(seq(1, n_draws, length.out = S_max)))
	} else {
		seq_len(n_draws)
	}
	S <- length(s_idx)

	n_obs <- prod(dim(Y))
	Y_vec <- as.numeric(Y)
	have_theta <- !is.null(Theta_samples)

	log_lik_samples <- matrix(NA_real_, nrow = S, ncol = n_obs)
	for (ii in seq_len(S)) {
		s <- s_idx[ii]
		s2_s <- s2_vec[s]
		if (!is.finite(s2_s) || s2_s <= 0) next
		M_s <- M_samples[[s]]
		if (is.null(M_s)) next

		if (have_theta) {
			Theta_s <- Theta_samples[[s]]
		} else {
			# reconstruct Theta from block operators starting at Theta_1 = M
			# (deterministic point trajectory; ignores process innovations)
			A_s <- A_samples[[s]]
			B_s <- B_samples[[s]]
			Theta_s <- array(0, dim = c(n_row, n_col, p, Tt))
			Theta_s[,,, 1] <- 0
			for (t in 2:Tt) {
				k <- t2k[t]
				Ak <- A_s[[k]]; Bk <- B_s[[k]]
				for (r in seq_len(p)) {
					Theta_s[,, r, t] <- Ak %*% Theta_s[,, r, t - 1L] %*% t(Bk)
				}
			}
		}

		# predicted mean for Y: Theta_t + M (broadcast over r and t)
		M_arr <- array(M_s, dim = c(n_row, n_col, p))
		mu <- array(0, dim = c(n_row, n_col, p, Tt))
		for (t in seq_len(Tt)) {
			for (r in seq_len(p)) {
				mu[,, r, t] <- Theta_s[,, r, t] + M_arr[,, r]
			}
		}

		if (family == "gaussian") {
			log_lik_samples[ii, ] <- dnorm(Y_vec, as.numeric(mu),
			                                sqrt(s2_s), log = TRUE)
		} else {
			# ordinal / binary: probit-scale approximation matching the
			# static path's convention. Sigma fixed at 1 on the latent scale.
			log_lik_samples[ii, ] <- dnorm(Y_vec, as.numeric(mu), 1,
			                                log = TRUE)
		}
	}
	log_lik_samples[!is.finite(log_lik_samples)] <- NA_real_

	# per-observation lppd and p_waic for loo::loo_compare SE
	col_ok <- apply(log_lik_samples, 2, function(x) sum(is.finite(x)) >= 2L)
	pointwise <- matrix(NA_real_, nrow = sum(col_ok), ncol = 2,
		dimnames = list(NULL, c("elpd_waic", "p_waic")))
	if (sum(col_ok) > 0L) {
		ll_ok <- log_lik_samples[, col_ok, drop = FALSE]
		lse <- apply(ll_ok, 2, function(x) {
			x <- x[is.finite(x)]
			if (length(x) == 0L) return(NA_real_)
			mx <- max(x); log(mean(exp(x - mx))) + mx
		})
		var_per_col <- apply(ll_ok, 2, function(x) {
			x <- x[is.finite(x)]
			if (length(x) < 2L) return(0)
			var(x)
		})
		pointwise[, "elpd_waic"] <- lse - var_per_col
		pointwise[, "p_waic"]    <- var_per_col
	}
	lppd <- sum(apply(log_lik_samples, 2, function(x) {
		x <- x[is.finite(x)]
		if (length(x) == 0) return(NA_real_)
		mx <- max(x)
		log(mean(exp(x - mx))) + mx
	}), na.rm = TRUE)
	p_waic <- sum(apply(log_lik_samples, 2, function(x) {
		x <- x[is.finite(x)]
		if (length(x) < 2) return(0)
		var(x)
	}), na.rm = TRUE)
	waic <- -2 * (lppd - p_waic)
	list(waic = waic, lppd = lppd, p_waic = p_waic, pointwise = pointwise)
}
####

####
#' Compute WAIC for Dynamic Model
#'
#' @description WAIC computation for dynamic DBN
#' @param fit dynamic dbn fit
#' @return list with waic, lppd, p_waic
#' @keywords internal
#' @noRd
compute_waic_dynamic <- function(fit) {
	# simplified WAIC using Theta samples
	Y <- fit$Y
	family <- fit$family
	Theta <- fit$Theta

	if (is.null(Theta)) {
		return(list(waic = NA, lppd = NA, p_waic = NA))
	}

	dims <- dim(Theta)
	n_draws <- dims[5]

	# subsample if too many draws
	if (n_draws > 100) {
		sample_idx <- seq(1, n_draws, length.out = 100)
	} else {
		sample_idx <- 1:n_draws
	}

	# compute approximate WAIC
	log_lik_samples <- matrix(NA, nrow = length(sample_idx), ncol = prod(dim(Y)))

	# use the per-draw latent state Theta_s as the location parameter so
	# WAIC reflects the dynamic model fit (rather than a constant-mean
	# proxy), and populate pointwise contributions for loo::loo_compare.
	for (i in seq_along(sample_idx)) {
		s <- sample_idx[i]
		Theta_s <- Theta[, , , , s]  # [n_row, n_col, p, Tt]

		if (family == "gaussian") {
			sigma2_s <- fit$sigma2_obs[s] %||% fit$sigma2[s]
			log_lik_samples[i, ] <- dnorm(c(Y), c(Theta_s), sqrt(sigma2_s), log = TRUE)
		} else {
			# binary / ordinal: Theta_s is the latent location on the probit
			# / rank scale; use standard-normal kernel as a coarse proxy.
			log_lik_samples[i, ] <- dnorm(c(Y), c(Theta_s), 1, log = TRUE)
		}
	}

	# NA cells in Y produce NA log-lik; drop them from pointwise summaries.
	log_lik_samples[!is.finite(log_lik_samples)] <- NA

	# pointwise contributions: per-observation lppd and p_waic, drop NA cols.
	col_ok <- apply(log_lik_samples, 2, function(x) sum(is.finite(x)) >= 2L)
	pointwise <- matrix(NA_real_, nrow = sum(col_ok), ncol = 2,
		dimnames = list(NULL, c("elpd_waic", "p_waic")))
	if (sum(col_ok) > 0L) {
		ll_ok <- log_lik_samples[, col_ok, drop = FALSE]
		lse <- apply(ll_ok, 2, function(x) {
			x <- x[is.finite(x)]
			if (length(x) == 0L) return(NA_real_)
			log(mean(exp(x - max(x)))) + max(x)
		})
		var_per_col <- apply(ll_ok, 2, function(x) {
			x <- x[is.finite(x)]
			if (length(x) < 2L) return(0)
			var(x)
		})
		pointwise[, "elpd_waic"] <- lse - var_per_col
		pointwise[, "p_waic"]    <- var_per_col
	}

	lppd <- sum(apply(log_lik_samples, 2, function(x) {
		x <- x[is.finite(x)]
		if (length(x) == 0) return(NA)
		log(mean(exp(x - max(x)))) + max(x)
	}), na.rm = TRUE)

	p_waic <- sum(apply(log_lik_samples, 2, function(x) {
		x <- x[is.finite(x)]
		if (length(x) < 2) return(0)
		var(x)
	}), na.rm = TRUE)

	waic <- -2 * (lppd - p_waic)

	list(waic = waic, lppd = lppd, p_waic = p_waic, pointwise = pointwise)
}
####
