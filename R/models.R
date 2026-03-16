####
#' Dynamic Bilinear Network Analysis
#'
#' @description
#' Main wrapper function for Dynamic Bilinear Network (DBN) analysis. This is the primary
#' interface for fitting DBN models to network data. DBN models capture complex dependencies
#' in network data through bilinear interactions between latent sender and receiver effects.
#'
#' @param data Data array or path to .RData file containing Y array.
#'   Array should be 3-dimensional (actors x actors x time) for single relation
#'   or 4-dimensional (actors x actors x relations x time) for multiple relations
#' @param family Character string specifying the data family/distribution:
#'   - "ordinal": For ordinal/ranked data (e.g., ratings 1-5). Data should be positive integers.
#'   - "gaussian": For continuous data. Data can be any real numbers.
#'   - "binary": For binary data. Data should be 0/1 or logical values.
#' @param model Character string specifying model type:
#'   - "static": Fixed sender/receiver effects across time
#'   - "dynamic": Time-varying sender/receiver effects
#'   - "lowrank": Low-rank factorization of sender effects
#'   - "hmm": Regime-switching model with hidden Markov states
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain (save every odens-th iteration)
#' @param verbose Logical or numeric. If TRUE, show progress. If numeric, print detailed info every n iterations (default: TRUE)
#' @param symmetric Logical. If TRUE, enforce B = A (symmetric/undirected network). Requires square network (n_row == n_col). Not supported for lowrank models. Default: FALSE.
#' @param ... Additional model-specific parameters:
#'   \describe{
#'     \item{\code{r}}{Rank for lowrank model (default: 2)}
#'     \item{\code{R}}{Number of regimes for HMM model (default: 3)}
#'     \item{\code{ar1}}{Use AR(1) dynamics for dynamic model (default: FALSE)}
#'     \item{\code{update_rho}}{Update AR coefficient in dynamic model (default: FALSE)}
#'     \item{\code{seed}}{Random seed for reproducibility (default: 6886)}
#'     \item{\code{previous}}{Previous fit object to continue MCMC from}
#'     \item{\code{init}}{List of initial values for parameters}
#'     \item{\code{time_thin}}{Time thinning factor for dynamic/lowrank/HMM (default: auto for dynamic, 1 for others)}
#'     \item{\code{store_z}}{Store Z draws for dynamic model (default: auto based on memory)}
#'   }
#' @return A list of class "dbn" containing:
#'   \item{B}{List of posterior samples for B matrices (static model)}
#'   \item{A}{List of posterior samples for time-varying A matrices (dynamic model)}
#'   \item{params}{Matrix of parameter traces (static model)}
#'   \item{sigma2}{Vector of sigma^2 samples (dynamic model)}
#'   \item{model}{Character string indicating which model was run}
#'   \item{dims}{List containing data dimensions}
#'   \item{settings}{List of model settings used}
#' @export
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_data)
#'
#' # Run static model with default settings
#' results <- dbn(example_data, model = "static")
#'
#' # Run dynamic model with custom MCMC settings
#' results <- dbn(example_data,
#'     model = "dynamic",
#'     nscan = 5000, burn = 1000, odens = 10
#' )
#'
#' # Run HMM model with 3 regimes
#' results <- dbn(example_data, model = "hmm", R = 3)
#'
#' # Run low-rank model with rank 2
#' results <- dbn(example_data, model = "lowrank", r = 2)
#'
#' # Run quietly without progress output
#' results <- dbn(example_data, model = "static", verbose = FALSE)
#'
#' # Run with detailed output every 100 iterations
#' results <- dbn(example_data, model = "dynamic", verbose = 100)
#' }
####
dbn <- function(data,
				family = c("ordinal", "gaussian", "binary"),
				model = c("static", "dynamic", "lowrank", "hmm"),
				nscan = 10000,
				burn = 1000,
				odens = 1,
				verbose = TRUE,
				symmetric = FALSE,
				...) {
	family <- match.arg(family)
	model <- match.arg(model)

	# load data from file path or use the array directly
	if (is.character(data)) {
		cli::cli_inform("Loading data from: {.path {data}}")
		env <- new.env()
		load(data, envir = env)
		Y <- env$Y
		if (is.null(Y)) cli::cli_abort("Data file must contain object {.var Y}")
	} else {
		Y <- data
	}
	####

	# validate dimensions and reshape 3D to 4D if needed
	if (length(dim(Y)) == 3) {
		dim_orig <- dim(Y)
		Y <- array(Y, dim = c(dim_orig[1], dim_orig[2], 1, dim_orig[3]))
		cli::cli_inform("Converting 3D array to 4D array with single relation")
	} else if (length(dim(Y)) != 4) {
		cli::cli_abort("Data must be a 3D array [actors x actors x time] or 4D array [actors x actors x relations x time]")
	}

	if (nscan <= 0) cli::cli_abort("{.arg nscan} must be positive.")
	if (burn < 0) cli::cli_abort("{.arg burn} must be non-negative.")
	if (odens < 1) cli::cli_abort("{.arg odens} must be at least 1.")
	if (burn >= nscan) cli::cli_abort("{.arg burn} ({burn}) must be less than {.arg nscan} ({nscan}).")
	if (odens > nscan) cli::cli_abort("{.arg odens} ({odens}) too large: no iterations would be saved ({.arg nscan} = {nscan}).")
	####

	# enforce static model for single time point
	Tt <- dim(Y)[4]
	if (Tt < 2) {
		if (model != "static") {
			cli::cli_abort(c(
				"Model {.val {model}} requires at least 2 time points.",
				"i" = "Your data has {Tt} time point{?s}.",
				"i" = "Use {.code model = \"static\"} for cross-sectional data."
			))
		}
		cli::cli_inform("Single time point detected -- using static model.")
	}
	####

	# print dimension summary
	if (verbose) {
		n_row <- dim(Y)[1]
		n_col <- dim(Y)[2]
		cli::cli_h3("Data dimensions")
		if (n_row != n_col) {
			cli::cli_bullets(c(
				" " = "Senders: {n_row}",
				" " = "Receivers: {n_col}",
				" " = "Relations: {dim(Y)[3]}",
				" " = "Time points: {Tt}",
				"i" = "Bipartite network detected"
			))
		} else {
			cli::cli_bullets(c(
				" " = "Nodes: {n_row}",
				" " = "Relations: {dim(Y)[3]}",
				" " = "Time points: {Tt}"
			))
		}
	}
	####

	# validate symmetric constraint
	if (symmetric) {
		n_r <- dim(Y)[1]
		n_c <- dim(Y)[2]
		if (n_r != n_c) {
			cli::cli_abort(c(
				"Symmetric networks require equal sender and receiver dimensions.",
				"i" = "Your data has {n_r} senders and {n_c} receivers.",
				"i" = "Symmetric networks are not compatible with bipartite data."
			))
		}
		if (model == "lowrank") {
			cli::cli_abort(c(
				"Symmetric networks are not yet supported for low-rank models.",
				"i" = "Use {.code model = \"dynamic\"} or {.code model = \"static\"} instead."
			))
		}
		if (verbose) {
			cli::cli_inform(c("i" = "Symmetric constraint active: B will be set equal to A."))
		}
	}
	####

	# suggest lowrank for large networks
	if (model == "dynamic" && max(dim(Y)[1], dim(Y)[2]) > 80 && verbose) {
		n_suggest <- max(dim(Y)[1], dim(Y)[2])
		r_suggest <- min(max(2L, as.integer(ceiling(log2(n_suggest)))), n_suggest - 1L)
		cli::cli_inform(c(
			"i" = "Large network ({n_suggest} nodes). Consider {.code model = \"lowrank\"} for better scalability.",
			"i" = "Suggested rank: {.code r = {r_suggest}} (log2(n), increase if fit is poor)."
		))
	}
	####

	# dispatch to model-specific sampler
	results <- switch(model,
		static = dbn_static(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, symmetric = symmetric, ...),
		dynamic = dbn_dynamic(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, symmetric = symmetric, ...),
		lowrank = dbn_lowrank(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, symmetric = symmetric, ...),
		hmm = dbn_hmm(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, symmetric = symmetric, ...)
	)

	results$model <- model
	class(results) <- "dbn"

	return(results)
	####
}

####
#' Static DBN MCMC
#'
#' @description Fits DBN model with fixed sender/receiver effects
#' @param Y Data array (nodes x nodes x relations x time)
#' @param family Character string specifying the data family/distribution:
#'   - "ordinal": Ordinal data (ordered categories). Data should be positive integers.
#'   - "gaussian": Continuous data with Gaussian errors. Data can be any real numbers.
#'   - "binary": Binary (0/1) data. Data should be 0/1 or logical values.
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param seed Random seed for reproducibility
#' @param verbose Logical or numeric. TRUE prints every 100 iterations, numeric prints every n iterations, FALSE suppresses output.
#' @param previous Previous dbn_static results to continue from (optional)
#' @param init List of initial values for parameters: B, s2, t2, g2, M, Z (optional)
#' @param symmetric Logical. If TRUE, store symmetric flag in output dims. Default: FALSE.
#' @return List containing MCMC results
#' @seealso \code{\link{dbn}} for the main dispatcher, \code{\link{param_summary}} for posterior summaries
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 5, seed = 1)
#' fit <- dbn_static(sim$Y, nscan = 200, burn = 100, verbose = FALSE)
#' }
#' @export
dbn_static <- function(Y, family = c("ordinal", "gaussian", "binary"),
					   nscan = 10000, burn = 1000, odens = 1,
					   seed = 6886, verbose = TRUE,
					   previous = NULL, init = NULL,
					   symmetric = FALSE) {

	# set up family and preprocess data
	family <- match.arg(family)
	if (isTRUE(verbose)) verbose <- 100L
	if (isFALSE(verbose)) verbose <- 0L
	FAM <- switch(family,
		ordinal = family_ordinal(),
		gaussian = family_gaussian(),
		binary = family_binary()
	)

	set.seed(seed)

	pre <- shared_preprocess(Y, family = family)
	Z <- pre$Z
	R <- pre$R
	M <- pre$M
	dims <- pre$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	m <- n_row
	p <- dims$p
	n <- Tt <- dims$Tt
	is_bipartite <- dims$is_bipartite
	nc <- n_row * n_col

	if (n_row == 1 || n_col == 1) {
		cli::cli_abort("Network must have at least 2 nodes on each side.")
	}

	is_large_network <- (max(n_row, n_col) > 15) || (p > 1) || (Tt > 20) || (nc * p * Tt > 10000)

	if (family == "ordinal") {
		R_flat <- array(R, c(n_row, n_col, p * Tt))
		IR <- precompute_rank_structure(R_flat, n_row, n_col, p, Tt)
	} else {
		IR <- pre$IR
	}
	K <- 3
	d <- c(n_row, n_col, p)
	####

	# initialize parameters (warm-start from previous fit if available)
	if (!is.null(previous)) {
		s2 <- tail(previous$draws$pars$s2, 1)
		t2 <- tail(previous$draws$pars$t2, 1)
		g2 <- tail(previous$draws$pars$g2, 1)
		M <- previous$draws$misc$M[[length(previous$draws$misc$M)]]
		B <- lapply(previous$draws$misc$B, function(b) {
			mat <- b[, , dim(b)[3], drop = FALSE]
			dim(mat) <- dim(b)[1:2]
			mat
		})
	} else {
		s2 <- 1
		t2 <- 1
		g2 <- max(0.1, mean(M^2, na.rm = TRUE))
		B <- list(diag(n_row), diag(n_col), diag(p))
	}
	####

	# allocate MCMC storage
	n_iter <- burn + nscan
	keep_idx <- seq(burn + 1, n_iter, by = odens)
	n_keep <- length(keep_idx)

	B_samples <- list()
	for (k in 1:K) B_samples[[k]] <- array(NA, c(d[k], d[k], n_keep))
	param_samples <- matrix(NA, n_keep, 3)
	colnames(param_samples) <- c("s2", "t2", "g2")
	Msave <- vector("list", n_keep)
	if (FAM$name %in% c("ordinal", "binary")) {
		Zsave <- vector("list", n_keep)
	}

	Z_flat <- array(Z, c(n_row, n_col, p * Tt))

	if (verbose) {
		cli::cli_alert_info("Running static DBN MCMC ({n_iter} iterations)")
		cli::cli_progress_bar("MCMC iterations", total = n_iter)
	}

	if (is_large_network) {
		Z_cube <- reshape_Z_to_cube_parallel(Z, n_row, n_col, p, Tt)
	} else {
		Z_cube <- reshape_Z_to_cube(Z, n_row, n_col, p, Tt)
	}

	if (!is.array(M) || length(dim(M)) != 3) {
		M <- array(as.numeric(M), dim = c(n_row, n_col, p))
	}
	storage.mode(Z_cube) <- "double"
	storage.mode(M) <- "double"

	if (is_large_network && family == "ordinal") {
		Z_flat_mat <- matrix(0, nrow = nc, ncol = p * Tt)
		EZ_flat_mat <- matrix(0, nrow = nc, ncol = p * Tt)
	}

	# precompute ordinal flags/arrays once (R doesn't change)
	if (family == "ordinal") {
		static_use_approx <- should_use_gaussian_approximation(R) ||
							(nc * p * Tt > 5000)
		static_R_flat <- array(R, c(n_row, n_col, p * Tt))
		static_has_rz_cpp <- exists("rz_gaussian_approx_cpp", mode = "function")
	}
	####

	# timing diagnostics (activated by verbose >= 200, e.g., verbose = 200L)
	static_do_timing <- is.numeric(verbose) && verbose >= 200L
	if (static_do_timing) {
		stiming <- list(b_update = 0, z_update = 0, m_update = 0, s2_update = 0, storage = 0)
		stiming_iters <- 0L
	}

	# main MCMC loop
	for (iter in 1:n_iter) {

		# update B matrices
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (is_large_network) {
			B[[1]] <- update_B_static_tiled(Z_cube, M, s2, t2, n_row, n_col, p, Tt)
		} else {
			B[[1]] <- update_B_static(Z_cube, M, s2, t2, n_row, n_col, p, Tt)
		}

		# update t2 (B precision)
		sse <- compute_diagonal_sse(B, K)
		t2 <- safe_rinv_gamma((sum(d) + 1) / 2, (sse + 1) / 2)
		if (static_do_timing) stiming$b_update <- stiming$b_update + (proc.time()[[3]] - .t0)
		####

		# update latent Z (ordinal/binary families)
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (FAM$name == "binary") {
			for (j in 1:p) {
				for (t in 1:Tt) {
					eta <- M[, , j]
					Y_jt <- R[, , j, t]
					pos <- which(Y_jt == 1)
					neg <- which(Y_jt == 0)
					if (length(pos) > 0) {
						Z[, , j, t][pos] <- truncnorm::rtruncnorm(
							length(pos), a = 0, b = Inf, mean = eta[pos], sd = 1)
					}
					if (length(neg) > 0) {
						Z[, , j, t][neg] <- truncnorm::rtruncnorm(
							length(neg), a = -Inf, b = 0, mean = eta[neg], sd = 1)
					}
				}
			}
			if (is_large_network) {
				Z_cube <- reshape_Z_to_cube_parallel(Z, n_row, n_col, p, Tt)
			} else {
				Z_cube <- reshape_Z_to_cube(Z, n_row, n_col, p, Tt)
			}
		} else if (FAM$name == "ordinal") {
			Z_flat <- array(Z, c(n_row, n_col, p * Tt))
			EZ_cube <- broadcast_M_and_compute_EZ(M, s2, n_row, n_col, p, Tt)

			if (static_use_approx) {
				if (static_has_rz_cpp) {
					Z_flat <- rz_gaussian_approx_cpp(static_R_flat, Z_flat, EZ_cube)
				} else {
					Z_flat <- rz_gaussian_approx(static_R_flat, Z_flat, EZ_cube)
				}
			} else {
				Z_flat <- rz_fc_batch(static_R_flat, Z_flat, EZ_cube, IR, n_row, n_col, p, Tt)
			}

			Z <- array(Z_flat, c(n_row, n_col, p, Tt))

			if (is_large_network) {
				Z_cube <- reshape_Z_to_cube_parallel(Z, n_row, n_col, p, Tt)
			} else {
				Z_cube <- reshape_Z_to_cube(Z, n_row, n_col, p, Tt)
			}
		}
		if (static_do_timing) stiming$z_update <- stiming$z_update + (proc.time()[[3]] - .t0)
		####

		# update baseline mean M
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (FAM$name == "ordinal" || iter %% 5 == 0) {
			if (FAM$name != "ordinal") {
				Z_flat <- array(Z, c(n_row, n_col, p * Tt))
			}
			Z_flat_mat <- matrix(Z_flat, nrow = nc, ncol = p * Tt)

			if (is_large_network) {
				M <- compute_M_static_blocked(Z_flat_mat, n_row, n_col, p, Tt)
			} else {
				M <- compute_M_static(Z_flat_mat, n_row, n_col, p, Tt)
			}

			M_sum_sq <- sum(M^2)
			g2 <- (1 + M_sum_sq) / (2 * rgamma(1, shape = (1 + nc * p) / 2, rate = 1))
		}
		if (static_do_timing) stiming$m_update <- stiming$m_update + (proc.time()[[3]] - .t0)
		####

		# update observation variance s2 (gaussian only)
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (FAM$name == "gaussian") {
			if (is_large_network) {
				rss <- compute_rss_static_parallel(Z_cube, M, n_row, n_col, p, Tt)
			} else {
				rss <- compute_rss_static(Z, M, n_row, n_col, p, Tt)
			}
			s2 <- safe_rinv_gamma(1 + nc * p * Tt / 2, 1 + rss / 2)
		} else {
			s2 <- 1
		}
		if (static_do_timing) stiming$s2_update <- stiming$s2_update + (proc.time()[[3]] - .t0)
		####

		# store thinned samples
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (iter %in% keep_idx) {
			idx <- which(keep_idx == iter)
			B_samples[[1]][, , idx] <- B[[1]]
			if (K > 1) {
				if (length(B) >= 2 && !is.null(B[[2]])) {
					B_samples[[2]][, , idx] <- B[[2]]
				} else {
					B_samples[[2]][, , idx] <- diag(d[2])
				}
			}
			if (K > 2) {
				if (length(B) >= 3 && !is.null(B[[3]])) {
					B_samples[[3]][, , idx] <- B[[3]]
				} else {
					id_mat <- diag(d[3])
					if (d[3] == 1) {
						id_mat <- matrix(1, 1, 1)
					}
					B_samples[[3]][, , idx] <- id_mat
				}
			}
			param_samples[idx, ] <- c(s2, t2, g2)
			Msave[[idx]] <- M
			if (FAM$name %in% c("ordinal", "binary")) {
				Zsave[[idx]] <- Z
			}
		}
		if (static_do_timing) stiming$storage <- stiming$storage + (proc.time()[[3]] - .t0)
		####

		if (static_do_timing) {
			stiming_iters <- stiming_iters + 1L
			if (iter == 10L) {
				total <- Reduce(`+`, stiming)
				if (total > 0) {
					cli::cli_alert_info("Static MCMC timing (first 10 iterations, avg per iter):")
					for (nm in names(stiming)) {
						pct <- round(100 * stiming[[nm]] / total, 1)
						ms <- round(1000 * stiming[[nm]] / stiming_iters, 1)
						cli::cli_alert("{nm}: {ms}ms ({pct}%)")
					}
					cli::cli_alert("Total: {round(1000 * total / stiming_iters, 1)}ms/iter")
				}
			}
		}

		if (verbose && iter %% verbose == 0) {
			cli::cli_progress_update()
		}
	}
	####

	if (verbose) cli::cli_progress_done()

	# assemble output list
	draws <- list(
		theta = NULL,
		z = if (FAM$name %in% c("ordinal", "binary")) Zsave else NULL,
		pars = data.frame(
			s2 = param_samples[, "s2"],
			t2 = param_samples[, "t2"],
			g2 = param_samples[, "g2"]
		),
		misc = list(
			M = Msave,
			B = B_samples
		)
	)

	dic <- NA
	pd <- NA
	deviance_mean <- NA
	if (FAM$name == "gaussian") {
		dev_samples <- numeric(n_keep)
		for (idx in seq_len(n_keep)) {
			M_idx <- Msave[[idx]]
			s2_idx <- param_samples[idx, "s2"]
			dev <- 0
			for (j in 1:p) {
				for (t in 1:Tt) {
					resid <- Y[, , j, t] - M_idx[, , j]
					dev <- dev + sum(resid^2, na.rm = TRUE)
				}
			}
			n_obs <- sum(!is.na(Y))
			dev_samples[idx] <- n_obs * log(2 * pi * s2_idx) + dev / s2_idx
		}
		deviance_mean <- mean(dev_samples)
		pd <- var(dev_samples) / 2
		dic <- deviance_mean + pd
	}

	out <- list(
		model = "static",
		family = family,
		Y = Y,
		R = R,
		dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, n = n,
					 is_bipartite = is_bipartite, is_symmetric = symmetric),
		settings = list(
			nscan = nscan,
			burn = burn,
			odens = odens,
			draws = n_keep
		),
		meta = list(
			family = family,
			dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
						is_bipartite = is_bipartite, is_symmetric = symmetric),
			draws = n_keep,
			settings = list(nscan = nscan, burn = burn, odens = odens)
		),
		params = param_samples,
		M = M,
		B = B_samples,
		draws = draws,
		diagnostics = list(
			deviance = deviance_mean,
			pd = pd,
			dic = dic
		),
		sigma2_obs = if (FAM$name == "gaussian") param_samples[, "s2"] else NULL,
		Z_final = if (FAM$name %in% c("ordinal", "binary")) Z else NULL,
		M_final = M
	)

	if (!is.null(previous)) {
		prev_total <- previous$total_iter %||% (previous$settings$burn + previous$settings$nscan)
		out$total_iter <- prev_total + nscan
		out$continued_from <- prev_total
	}

	class(out) <- "dbn"
	return(out)
	####
}
####

####
#' Dynamic DBN MCMC
#'
#' @description Fits DBN model with time-varying sender/receiver effects
#' @param Y Data array (nodes x nodes x relations x time)
#' @param family Character string specifying the data family/distribution:
#'   - "ordinal": Ordinal data (ordered categories). Data should be positive integers.
#'   - "gaussian": Continuous data with Gaussian errors. Data can be any real numbers.
#'   - "binary": Binary (0/1) data. Data should be 0/1 or logical values.
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param ar1 Use AR(1) dynamics instead of random walk (default: FALSE)
#' @param update_rho Update AR coefficients (default: FALSE)
#' @param seed Random seed
#' @param verbose Logical or numeric. TRUE prints every 100 iterations, numeric prints every n iterations, FALSE suppresses output.
#' @param time_thin Save every nth time point to reduce memory (default: NULL = auto)
#' @param store_z Whether to store Z draws (default: NULL = auto based on memory)
#' @param previous Previous dbn_dynamic results to continue from (optional)
#' @param init List of initial values: A, B, sigma2, tau_A2, tau_B2, g2, rho_A, rho_B, Theta, M, Z (optional)
#' @param symmetric Logical. If TRUE, enforce B = A after each MCMC update. Default: FALSE.
#' @return List containing MCMC results
#' @seealso \code{\link{dbn}} for the main dispatcher, \code{\link{param_summary}} for posterior summaries
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn_dynamic(sim$Y, nscan = 200, burn = 100, verbose = FALSE)
#' }
#' @export
dbn_dynamic <- function(Y,
						family = c("ordinal", "gaussian", "binary"),
						nscan = 10000,
						burn = 1000,
						odens = 1,
						ar1 = FALSE,
						update_rho = FALSE,
						seed = 6886,
						verbose = TRUE,
						time_thin = NULL,
						store_z = NULL,
						previous = NULL,
						init = NULL,
						symmetric = FALSE) {

	# set up family and preprocess data
	set.seed(seed)
	family <- match.arg(family)
	if (isTRUE(verbose)) verbose <- 100L
	if (isFALSE(verbose)) verbose <- 0L
	FAM <- switch(family,
		ordinal = family_ordinal(),
		gaussian = family_gaussian(),
		binary = family_binary()
	)

	if (dim(Y)[4] < 2) {
		cli::cli_abort("Dynamic model requires at least 2 time points. Use {.code model = \"static\"} for cross-sectional data.")
	}

	# small binary networks can be numerically unstable
	if (family == "binary" && min(dim(Y)[1], dim(Y)[2]) < 15) {
		cli::cli_warn(c(
			"Dynamic binary models with small networks (n < 15) may encounter numerical singularities.",
			"i" = "Consider using {.code model = \"static\"} or a larger network.",
			"i" = "The model will attempt to run, but may produce unreliable results."
		))
	}

	pre <- shared_preprocess(Y, family = family)
	Z <- pre$Z
	R <- pre$R
	IR <- pre$IR
	M <- pre$M
	dims <- pre$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	m <- n_row
	p <- dims$p
	Tt <- dims$Tt
	is_bipartite <- dims$is_bipartite
	nc <- n_row * n_col
	d <- nc
	####

	# configure time-thinning and estimate memory usage
	if (is.null(time_thin)) {
		time_thin <- max(1L, Tt %/% 20L)
		if (time_thin > 1 && isTRUE(verbose > 0)) {
			cli::cli_inform(c(
				"i" = "Auto time-thinning: storing every {time_thin}th time point ({length(seq(1, Tt, by = time_thin))} of {Tt}).",
				"i" = "Override with {.code time_thin = 1} to store all."
			))
		}
	}

	if (is.null(store_z)) {
		n_keep_est <- floor(nscan / odens)
		n_time_est <- length(seq(1, Tt, by = time_thin))
		z_gb <- nc * p * n_time_est * n_keep_est * 8 / 1024^3
		store_z <- (z_gb < 2)
		if (!store_z && family %in% c("ordinal", "binary") && isTRUE(verbose > 0)) {
			cli::cli_inform(c(
				"i" = "Z array would be {round(z_gb, 1)} GB. Skipping Z storage to save memory.",
				"i" = "Override with {.code store_z = TRUE} to force storage."
			))
		}
	}

	if (isTRUE(verbose > 0)) {
		n_keep_est <- floor(nscan / odens)
		n_time_est <- length(seq(1, Tt, by = time_thin))
		est_gb <- estimate_memory(n_row, n_col, p, Tt, nscan, burn, odens,
								  time_thin, family, quiet = TRUE)
		if (est_gb > 4) {
			cli::cli_warn(c(
				"!" = "Estimated output size is {round(est_gb, 1)} GB.",
				"i" = "Consider increasing {.code time_thin} or {.code odens}, or using {.code model = \"lowrank\"}."
			))
		}
	}

	is_large_network <- (max(n_row, n_col) > 15) || (p > 1) || (Tt > 20) || (nc * p * Tt > 10000)

	if (family == "ordinal" && !is.null(IR)) {
		IR_time_indices <- precompute_time_indices(IR, n_row, n_col, p, Tt)
	}
	####

	# flatten 4D arrays to matrices for vectorized operations
	Z_4d <- matrix(Z, nrow = nc, ncol = p * Tt)
	R_4d <- matrix(R, nrow = nc, ncol = p * Tt)
	####

	# initialize parameters (warm-start from previous fit if available)
	if (!is.null(previous)) {
		if (is.null(previous$A) || is.null(previous$B) || is.null(previous$sigma2)) {
			cli::cli_abort("{.arg previous} must be results from {.fun dbn_dynamic}.")
		}
		last_idx <- length(previous$sigma2)

		if (!is.null(previous$Theta)) {
			n_prev_iter <- dim(previous$Theta)[5]
			last_theta_idx <- n_prev_iter
			prev_time_thin <- previous$time_thin %||% 1
			n_time_stored <- dim(previous$Theta)[4]

			if (n_time_stored < Tt) {
				time_indices <- seq(1, Tt, by = prev_time_thin)
				Theta_all <- array(0, dim = c(n_row, n_col, p, Tt))
				for (i in 1:min(length(time_indices), n_time_stored)) {
					if (time_indices[i] <= Tt) {
						Theta_all[,,,time_indices[i]] <- previous$Theta[,,,i,last_theta_idx]
					}
				}
				for (t in 1:Tt) {
					if (all(Theta_all[,,,t] == 0)) {
						available_times <- time_indices[time_indices <= Tt]
						nearest_idx <- which.min(abs(available_times - t))
						nearest_stored <- nearest_idx
						if (nearest_stored <= n_time_stored) {
							Theta_all[,,,t] <- previous$Theta[,,,nearest_stored,last_theta_idx]
						}
					}
				}
			} else {
				Theta_all <- previous$Theta[,,,1:Tt,last_theta_idx]
			}
		} else {
			Theta_all <- pre$Theta
		}

		sigma2 <- previous$sigma2[last_idx]
		sigma2_obs <- if (!is.null(previous$sigma2_obs)) previous$sigma2_obs[last_idx] else 1
		tauA2 <- previous$tau_A2[last_idx]
		tauB2 <- previous$tau_B2[last_idx]
		g2 <- previous$g2[last_idx]
		rhoA <- if (!is.null(previous$rhoA)) previous$rhoA[last_idx] else 0
		rhoB <- if (!is.null(previous$rhoB)) previous$rhoB[last_idx] else 0

		last_A <- previous$A[[last_idx]]
		last_B <- previous$B[[last_idx]]

		if (dim(last_A)[3] < Tt) {
			prev_time_thin <- previous$time_thin %||% 1
			time_indices <- seq(1, Tt, by = prev_time_thin)
			Aarray <- array(0, dim = c(n_row, n_row, Tt))
			Barray <- array(0, dim = c(n_col, n_col, Tt))
			for (i in 1:length(time_indices)) {
				if (time_indices[i] <= Tt && i <= dim(last_A)[3]) {
					Aarray[,,time_indices[i]] <- last_A[,,i]
					Barray[,,time_indices[i]] <- last_B[,,i]
				}
			}
			for (t in 1:Tt) {
				if (all(Aarray[,,t] == 0)) {
					available_times <- time_indices[time_indices <= Tt]
					nearest_idx <- which.min(abs(available_times - t))
					nearest_time <- available_times[nearest_idx]
					nearest_stored <- which(time_indices == nearest_time)
					Aarray[,,t] <- last_A[,,nearest_stored]
					Barray[,,t] <- last_B[,,nearest_stored]
				}
			}
		} else {
			Aarray <- last_A
			Barray <- last_B
		}
	} else {
		Theta_all <- pre$Theta
		Aarray <- array(0, dim = c(n_row, n_row, Tt))
		Barray <- array(0, dim = c(n_col, n_col, Tt))
		for (t in 1:Tt) {
			Aarray[, , t] <- diag(n_row)
			Barray[, , t] <- diag(n_col)
		}
		sigma2 <- 1
		sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
		tauA2 <- tauB2 <- 1
		g2 <- 1
		rhoA <- rhoB <- 0
	}

	Theta_4d <- matrix(Theta_all, nrow = nc, ncol = p * Tt)
	####

	# allocate MCMC storage
	n_iter <- burn + nscan
	keep <- seq(burn + 1, n_iter, by = odens)
	n_keep <- length(keep)
	time_keep <- seq(1, Tt, by = time_thin)
	n_time_keep <- length(time_keep)

	Theta_store <- array(NA, dim = c(n_row, n_col, p, n_time_keep, n_keep))
	if (FAM$name %in% c("ordinal", "binary") && store_z) {
		Z_store <- array(NA, dim = c(n_row, n_col, p, n_time_keep, n_keep))
	}

	A_store <- B_store <- vector("list", n_keep)
	M_store <- array(NA, dim = c(n_row, n_col, p, n_keep))

	sigma2_store <- numeric(n_keep)
	sigma2_obs_store <- if (FAM$name == "gaussian") numeric(n_keep) else NULL
	tau_A2_store <- tau_B2_store <- g2_store <- numeric(n_keep)
	rhoA_store <- rhoB_store <- if (ar1) numeric(n_keep) else NULL

	keep_id <- 0
	a_sig <- b_sig <- 2
	a_g <- b_g <- 2
	eye_nr <- diag(n_row)
	eye_nc <- diag(n_col)

	if (is_large_network) {
		shape_tauA <- (1 + n_row * n_row * (Tt - 1)) / 2
		shape_tauB <- (1 + n_col * n_col * (Tt - 1)) / 2
		shape_sigma_proc <- (a_sig + nc * (Tt - 1) * p) / 2.0
		shape_sigma_obs <- (1.0 + nc * Tt * p) / 2.0
	}

	use_approx <- FALSE
	if (FAM$name == "ordinal") {
		use_approx <- should_use_gaussian_approximation(R) ||
					 (nc * p * Tt > 5000)
		if (use_approx) {
			EZ_cube <- array(0, c(n_row, n_col, p * Tt))
			Z_cube  <- array(0, c(n_row, n_col, p * Tt))
			R_cube  <- array(R_4d, c(n_row, n_col, p * Tt))
			M_expanded <- rep(as.vector(M), times = Tt)
		}
	}

	# precompute function availability flags (avoid per-iteration exists() calls)
	has_rz_gaussian_approx_cpp <- exists("rz_gaussian_approx_cpp", mode = "function")
	has_IR_time_indices <- exists("IR_time_indices")
	has_shape_sigma_proc <- exists("shape_sigma_proc")

	if (verbose) {
		cli::cli_alert_info("Running dynamic DBN MCMC ({n_iter} iterations)")
		cli::cli_progress_bar("MCMC iterations", total = n_iter)
	}
	####

	# timing diagnostics (activated by verbose >= 200, e.g., verbose = 200L)
	do_timing <- is.numeric(verbose) && verbose >= 200L
	if (do_timing) {
		timing_accum <- list(z_update = 0, mu_update = 0, ffbs = 0, ab_update = 0,
							 tau_update = 0, var_update = 0, rho_update = 0, storage = 0)
		timing_iters <- 0L
	}

	# main MCMC loop
	for (g in 1:n_iter) {

		# update latent Z (ordinal/binary families)
		if (do_timing) .t0 <- proc.time()[[3]]
		if (FAM$name == "ordinal") {
			if (use_approx) {
				EZ_cube[] <- Theta_4d + M_expanded
				Z_cube[]  <- Z_4d
				R_cube[]  <- R_4d

				if (has_rz_gaussian_approx_cpp) {
					Z_cube <- rz_gaussian_approx_cpp(R_cube, Z_cube, EZ_cube)
				} else {
					Z_cube <- rz_gaussian_approx(R_cube, Z_cube, EZ_cube)
				}
				Z_4d <- matrix(Z_cube, nrow = nc, ncol = p * Tt)
			} else {
				if (has_IR_time_indices) {
					Z_4d <- batch_update_Z_ordinal_fast(R_4d, Z_4d, Theta_4d, M, IR, IR_time_indices, n_row, n_col, p, Tt)
				} else {
					Z_4d <- batch_update_Z_ordinal(R_4d, Z_4d, Theta_4d, M, IR, n_row, n_col, p, Tt)
				}
			}
			if (g %in% keep) {
				Z <- array(Z_4d, dim = c(n_row, n_col, p, Tt))
			}
		} else if (FAM$name == "binary") {
			Z <- update_Z_optimized(R, Z, Theta_all, M, IR = NULL, family = "binary")
			Z_4d <- matrix(Z, nrow = nc, ncol = p * Tt)
		}
		if (do_timing) timing_accum$z_update <- timing_accum$z_update + (proc.time()[[3]] - .t0)
		####

		# update baseline mean M and g2
		if (do_timing) .t0 <- proc.time()[[3]]
		mu_result <- update_mu_dynamic(Z_4d, Theta_4d, g2, a_g, b_g, n_row, n_col, p, Tt)
		M <- mu_result$M
		g2 <- mu_result$g2
		if (use_approx) M_expanded <- rep(as.vector(M), times = Tt)
		if (do_timing) timing_accum$mu_update <- timing_accum$mu_update + (proc.time()[[3]] - .t0)
		####

		# FFBS for Theta
		if (do_timing) .t0 <- proc.time()[[3]]
		if (is_large_network && max(n_row, n_col) > 100) {
			Theta_cube <- batch_ffbs_all_relations_blocked(Z_4d, M, Aarray, Barray, sigma2, n_row, n_col, p, Tt)
		} else {
			Theta_cube <- batch_ffbs_all_relations(Z_4d, M, Aarray, Barray, sigma2, n_row, n_col, p, Tt)
		}
		if (g %in% keep) {
			Theta_all <- array(Theta_cube, dim = c(n_row, n_col, p, Tt))
		}
		Theta_4d <- matrix(Theta_cube, nrow = nc, ncol = p * Tt)
		if (do_timing) timing_accum$ffbs <- timing_accum$ffbs + (proc.time()[[3]] - .t0)
		####

		# update A and B matrices
		if (do_timing) .t0 <- proc.time()[[3]]
		if (is_large_network && max(n_row, n_col) > 100) {
			AB_result <- update_AB_batch_large(
				Theta_4d, Aarray, Barray,
				sigma2, tauA2, tauB2,
				ar1, rhoA, rhoB,
				n_row, n_col, p, Tt
			)
		} else {
			AB_result <- update_AB_batch_extended(
				Theta_4d, Aarray, Barray,
				sigma2, tauA2, tauB2,
				ar1, rhoA, rhoB,
				n_row, n_col, p, Tt
			)
		}
		Aarray <- AB_result$Aarray
		Barray <- AB_result$Barray
		if (symmetric) Barray <- Aarray
		if (do_timing) timing_accum$ab_update <- timing_accum$ab_update + (proc.time()[[3]] - .t0)
		####

		# update innovation variances tauA2, tauB2
		if (do_timing) .t0 <- proc.time()[[3]]
		if (ar1) {
			innovA_ss <- compute_ar1_innovation_ss_cpp(Aarray, rhoA, n_row, Tt)
			innovB_ss <- compute_ar1_innovation_ss_cpp(Barray, rhoB, n_col, Tt)
			if (is_large_network) {
				tauA2 <- safe_rinv_gamma(shape_tauA, (1 + innovA_ss)/2)
				tauB2 <- safe_rinv_gamma(shape_tauB, (1 + innovB_ss)/2)
			} else {
				tauA2 <- safe_rinv_gamma((1 + n_row * n_row * (Tt-1))/2, (1 + innovA_ss)/2)
				tauB2 <- safe_rinv_gamma((1 + n_col * n_col * (Tt-1))/2, (1 + innovB_ss)/2)
			}
		} else {
			A_sum <- compute_deviation_sum(Aarray, n_row, Tt)
			B_sum <- compute_deviation_sum(Barray, n_col, Tt)
			if (is_large_network) {
				tauA2 <- safe_rinv_gamma(shape_tauA, (1 + A_sum)/2)
				tauB2 <- safe_rinv_gamma(shape_tauB, (1 + B_sum)/2)
			} else {
				tauA2 <- safe_rinv_gamma((1 + n_row * n_row * (Tt-1))/2, (1 + A_sum)/2)
				tauB2 <- safe_rinv_gamma((1 + n_col * n_col * (Tt-1))/2, (1 + B_sum)/2)
			}
		}
		if (symmetric) tauB2 <- tauA2
		if (do_timing) timing_accum$tau_update <- timing_accum$tau_update + (proc.time()[[3]] - .t0)
		####

		# update process and observation variances
		if (do_timing) .t0 <- proc.time()[[3]]
		if (is_large_network && max(n_row, n_col) > 100) {
			proc_rss <- compute_process_variance_blocked(Theta_4d, Aarray, Barray, n_row, n_col, p, Tt)
			if (has_shape_sigma_proc) {
				sigma2 <- (b_sig + proc_rss / 2.0) / rgamma(1, shape = shape_sigma_proc, rate = 1)
			} else {
				sigma2 <- (b_sig + proc_rss / 2.0) / rgamma(1, shape = (a_sig + nc * (Tt - 1) * p) / 2.0, rate = 1)
			}

			if (FAM$name == "gaussian") {
				obs_rss <- compute_gaussian_obs_residuals_dynamic_cpp(Z_4d, Theta_4d, M, n_row, n_col, p, Tt)
				sigma2_obs <- (1 + obs_rss / 2) / rgamma(1, shape = (1 + nc * Tt * p) / 2, rate = 1)
			}
		} else {
			var_result <- update_variances_dynamic(
				Theta_4d, Z_4d, M, Aarray, Barray,
				a_sig, b_sig, n_row, n_col, p, Tt,
				is_gaussian = (FAM$name == "gaussian")
			)
			sigma2 <- var_result$sigma2
			if (FAM$name == "gaussian") {
				sigma2_obs <- var_result$sigma2_obs
			}
		}
		if (do_timing) timing_accum$var_update <- timing_accum$var_update + (proc.time()[[3]] - .t0)
		####

		# update AR(1) coefficients rhoA, rhoB
		if (do_timing) .t0 <- proc.time()[[3]]
		if (ar1 && update_rho) {
			rhoA_result <- compute_rho_update_cpp(Aarray, n_row, Tt)
			rho_mean <- rhoA_result$num / (rhoA_result$denom + 1e-10)
			rho_var  <- tauA2 / (rhoA_result$denom + 1e-10)
			rhoA <- truncnorm::rtruncnorm(1, a = -0.99, b = 0.99, mean = rho_mean, sd = sqrt(rho_var))

			rhoB_result <- compute_rho_update_cpp(Barray, n_col, Tt)
			rho_mean <- rhoB_result$num / (rhoB_result$denom + 1e-10)
			rho_var  <- tauB2 / (rhoB_result$denom + 1e-10)
			rhoB <- truncnorm::rtruncnorm(1, a = -0.99, b = 0.99, mean = rho_mean, sd = sqrt(rho_var))
			if (symmetric) rhoB <- rhoA
		}
		if (do_timing) timing_accum$rho_update <- timing_accum$rho_update + (proc.time()[[3]] - .t0)
		####

		# store thinned samples
		if (do_timing) .t0 <- proc.time()[[3]]
		if (g %in% keep) {
			keep_id <- keep_id + 1

			if (!exists("Theta_all") || !is.array(Theta_all) || length(dim(Theta_all)) != 4) {
				Theta_all <- array(Theta_cube, dim = c(n_row, n_col, p, Tt))
			}
			Theta_store[,,,,keep_id] <- Theta_all[,,,time_keep, drop = FALSE]

			if (FAM$name %in% c("ordinal", "binary") && store_z) {
				if (!is.array(Z) || length(dim(Z)) != 4) {
					Z <- array(Z_4d, dim = c(n_row, n_col, p, Tt))
				}
				Z_store[,,,,keep_id] <- Z[,,,time_keep, drop = FALSE]
			}

			A_store[[keep_id]] <- Aarray[,,time_keep, drop = FALSE]
			B_store[[keep_id]] <- Barray[,,time_keep, drop = FALSE]

			M_store[,,,keep_id] <- M
			sigma2_store[keep_id] <- sigma2
			if (FAM$name == "gaussian") {
				sigma2_obs_store[keep_id] <- sigma2_obs
			}

			tau_A2_store[keep_id] <- tauA2
			tau_B2_store[keep_id] <- tauB2
			g2_store[keep_id] <- g2

			if (ar1) {
				rhoA_store[keep_id] <- rhoA
				rhoB_store[keep_id] <- rhoB
			}
		}
		if (do_timing) timing_accum$storage <- timing_accum$storage + (proc.time()[[3]] - .t0)
		####

		if (do_timing) {
			timing_iters <- timing_iters + 1L
			if (g == 10L) {
				total <- Reduce(`+`, timing_accum)
				if (total > 0) {
					cli::cli_alert_info("MCMC timing (first 10 iterations, avg per iter):")
					for (nm in names(timing_accum)) {
						pct <- round(100 * timing_accum[[nm]] / total, 1)
						ms <- round(1000 * timing_accum[[nm]] / timing_iters, 1)
						cli::cli_alert("{nm}: {ms}ms ({pct}%)")
					}
					cli::cli_alert("Total: {round(1000 * total / timing_iters, 1)}ms/iter")
				}
			}
		}

		if (verbose && (g %% verbose == 0)) cli::cli_progress_update()
	}
	####

	if (verbose) cli::cli_progress_done()

	# assemble output list
	out <- list(
		model = "dynamic",
		Y = Y,
		Theta = Theta_store,
		A = A_store,
		B = B_store,
		M = M_store,
		sigma2 = sigma2_store,
		tau_A2 = tau_A2_store,
		tau_B2 = tau_B2_store,
		g2 = g2_store,
		n_iter = n_iter,
		burn = burn,
		thin = odens,
		time_thin = time_thin,
		dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
					is_bipartite = is_bipartite, is_symmetric = symmetric),
		family = family,
		time_kept = time_keep,
		settings = list(time_thin = time_thin, store_z = store_z)
	)

	if (FAM$name == "gaussian") {
		out$sigma2_obs <- sigma2_obs_store
	}

	if (FAM$name %in% c("ordinal", "binary") && store_z) {
		out$Z <- Z_store
	}

	if (ar1) {
		out$rhoA <- rhoA_store
		out$rhoB <- rhoB_store
		out$ar1 <- TRUE
	}

	if (!is.null(previous)) {
		prev_total <- previous$n_iter %||% (previous$burn + length(previous$sigma2))
		out$total_iter <- prev_total + nscan
		out$continued_from <- prev_total
	}

	class(out) <- "dbn"
	out
	####
}
####
