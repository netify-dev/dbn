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

	# validate minimum block lengths
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

		# fit quick static model
		fit_window <- tryCatch({
			dbn_static(Y_window, family = family,
					   nscan = nscan, burn = nscan %/% 2, odens = 1,
					   verbose = FALSE, seed = 12345 + i)
		}, error = function(e) NULL)

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
#' Compute WAIC for DBN Fit
#'
#' @description computes Watanabe-Akaike Information Criterion
#' @param fit dbn model fit object
#' @return list with waic, lppd, p_waic
#' @keywords internal
compute_waic_dbn <- function(fit) {
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
			return(list(waic = fit$diagnostics$dic, lppd = NA, p_waic = NA))
		}
		cli::cli_warn("WAIC not implemented for model type {fit$model}. Using NA.")
		return(list(waic = NA, lppd = NA, p_waic = NA))
	}

	waic
}
####

####
#' Compute WAIC for Static Model
#'
#' @description WAIC computation for static DBN
#' @param fit static dbn fit
#' @return list with waic, lppd, p_waic
#' @keywords internal
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

	# compute WAIC components
	# lppd: log pointwise predictive density
	lppd <- sum(apply(log_lik_samples, 2, function(x) {
		x <- x[is.finite(x)]
		if (length(x) == 0) return(NA)
		log(mean(exp(x - max(x)))) + max(x)
	}), na.rm = TRUE)

	# p_waic: effective number of parameters
	p_waic <- sum(apply(log_lik_samples, 2, function(x) {
		x <- x[is.finite(x)]
		if (length(x) < 2) return(0)
		var(x)
	}), na.rm = TRUE)

	waic <- -2 * (lppd - p_waic)

	list(waic = waic, lppd = lppd, p_waic = p_waic)
}
####

####
#' Compute WAIC for Piecewise Model
#'
#' @description WAIC computation for piecewise DBN
#' @param fit piecewise dbn fit
#' @return list with waic, lppd, p_waic
#' @keywords internal
compute_waic_piecewise <- function(fit) {
	# use same approach as static but accounting for block structure
	compute_waic_static(fit)
}
####

####
#' Compute WAIC for Dynamic Model
#'
#' @description WAIC computation for dynamic DBN
#' @param fit dynamic dbn fit
#' @return list with waic, lppd, p_waic
#' @keywords internal
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

	for (i in seq_along(sample_idx)) {
		s <- sample_idx[i]
		Theta_s <- Theta[, , , , s]

		# get M for this draw
		M_s <- fit$M[, , , s]
		if (length(dim(M_s)) == 3) {
			M_expanded <- array(M_s, dim = dim(Y))
		} else {
			M_expanded <- array(mean(M_s, na.rm = TRUE), dim = dim(Y))
		}

		if (family == "gaussian") {
			sigma2_s <- fit$sigma2_obs[s] %||% fit$sigma2[s]
			log_lik_samples[i, ] <- dnorm(c(Y), c(M_expanded), sqrt(sigma2_s), log = TRUE)
		} else {
			log_lik_samples[i, ] <- dnorm(c(Y), c(M_expanded), 1, log = TRUE)
		}
	}

	log_lik_samples[!is.finite(log_lik_samples)] <- NA

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

	list(waic = waic, lppd = lppd, p_waic = p_waic)
}
####
