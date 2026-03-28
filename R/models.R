####
#' Dynamic Bilinear Network Analysis
#'
#' @description
#' Main entry point for fitting Dynamic Bilinear Network (DBN) models. These
#' models estimate how past network interactions predict future interactions,
#' recovering time-varying influence structures from temporal relational data.
#'
#' The core model is: \eqn{\Theta_t = A_t \Theta_{t-1} B_t' + M + \varepsilon_t},
#' where \eqn{A_t} captures sender influence, \eqn{B_t} captures receiver
#' influence, and \eqn{M} captures stable dyad-specific tendencies.
#'
#' @details
#' **Choosing a family:**
#' - `"gaussian"`: Use when your relational data are continuous measurements
#'   (e.g., trade volumes, similarity scores). This is the simplest family and
#'   converges fastest.
#' - `"ordinal"`: Use when your data are ordered categories or counts (e.g.,
#'   conflict severity 1-5, event counts) and you trust the ordering but not
#'   the exact values. Uses a rank likelihood.
#' - `"binary"`: Use for presence/absence data (0/1), e.g., whether a tie
#'   exists. Uses a probit link with data augmentation.
#'
#' **Choosing a model:**
#' - `"static"`: Simplest. Influence structure is fixed over time. Good
#'   starting point and for short time series.
#' - `"dynamic"`: Influence structure changes over time. Use when you expect
#'   shifting alliances, evolving trade patterns, etc.
#' - `"piecewise"`: Influence is constant within known regimes but differs
#'   across them. Use when you know when structural breaks occurred (e.g.,
#'   before/after a crisis).
#' - `"hmm"`: Like piecewise but discovers regimes from data. Use when breaks
#'   are unknown.
#' - `"lowrank"`: Like dynamic but with dimensionality reduction for large
#'   networks (50+ actors).
#'
#' **MCMC settings:**
#' The sampler draws `nscan` posterior samples after discarding the first
#' `burn` as warm-up. Setting `odens > 1` thins the output by saving every
#' k-th sample. For initial exploration, `nscan = 5000, burn = 2000, odens = 5`
#' is a reasonable starting point. For publication, use longer chains
#' (`nscan = 10000+`) and verify convergence with [check_convergence()].
#'
#' @param data Numeric array of network data, or a file path to an `.RData`
#'   file that contains an object named `Y`.  The array should be
#'   3-dimensional `[actors, actors, time]` for a single relation type, or
#'   4-dimensional `[actors, actors, relations, time]` for multiple relation
#'   types.  Diagonal entries (self-ties) should be `NA` for unipartite
#'   networks.  For bipartite networks, pass a rectangular array where the
#'   first dimension (senders) differs from the second (receivers).
#' @param family Character string specifying the outcome distribution.
#'   See **Details** for guidance on choosing:
#'   \itemize{
#'     \item `"ordinal"`: For ordinal/ranked data (positive integers)
#'     \item `"gaussian"`: For continuous data (any real numbers)
#'     \item `"binary"`: For binary data (0/1 or logical)
#'   }
#' @param model Character string specifying the model type.
#'   See **Details** for guidance on choosing:
#'   \itemize{
#'     \item `"static"`: Fixed sender/receiver effects across time
#'     \item `"dynamic"`: Time-varying sender/receiver effects
#'     \item `"lowrank"`: Low-rank factorization of sender effects (large networks)
#'     \item `"hmm"`: Regime-switching with data-driven regime discovery
#'     \item `"piecewise"`: Block-constant influence with known break points
#'   }
#' @param nscan Number of posterior samples to draw after burn-in
#' @param burn Number of initial MCMC samples to discard (warm-up period)
#' @param odens Thinning interval: save every odens-th sample (reduces autocorrelation and memory)
#' @param verbose Logical or numeric. If TRUE, show progress. If numeric, print detailed info every n iterations (default: TRUE)
#' @param symmetric Logical. If TRUE, enforce B = A (symmetric/undirected network). Requires square network (n_row == n_col). Not supported for lowrank models. Default: FALSE.
#' @param ... Additional model-specific parameters:
#'   \describe{
#'     \item{\code{r}}{Rank for lowrank model (default: 2)}
#'     \item{\code{R}}{Number of regimes for HMM model (default: 3)}
#'     \item{\code{blocks}}{Block specification for piecewise model: integer (number of equal blocks),
#'       numeric vector (block boundaries), named vector (labeled boundaries), or "auto" for automatic selection}
#'     \item{\code{ar1}}{Use AR(1) dynamics for dynamic model (default: FALSE)}
#'     \item{\code{update_rho}}{Update AR coefficient in dynamic model (default: FALSE)}
#'     \item{\code{seed}}{Random seed for reproducibility (default: 6886)}
#'     \item{\code{previous}}{Previous fit object to continue MCMC from}
#'     \item{\code{init}}{List of initial values for parameters}
#'     \item{\code{time_thin}}{Time thinning factor for dynamic/lowrank/HMM (default: auto for dynamic, 1 for others)}
#'     \item{\code{store_z}}{Store Z draws for dynamic model (default: auto based on memory)}
#'     \item{\code{store_theta}}{Store full Theta trajectory draws for piecewise model (default: TRUE).
#'       \strong{Critical for large networks:} Set to FALSE for networks with 100+ actors to avoid
#'       memory issues. Theta storage scales as O(n^2 * T * draws) -- a 200-actor network with 50 time
#'       points and 500 draws requires ~40 GB. With \code{store_theta = FALSE}, you retain posterior
#'       draws for A, B, M and variance parameters, \code{compare_blocks()} functionality, and
#'       convergence diagnostics. You lose full posterior uncertainty on individual Theta entries
#'       and \code{posterior_predict_dbn()} with uncertainty propagation.}
#'   }
#' @return A list of class `"dbn"` with model-specific contents. Common elements:
#'   \item{model}{Character string indicating which model was fit}
#'   \item{family}{Character string indicating the outcome family}
#'   \item{dims}{List of data dimensions (n_row, n_col, p, Tt)}
#'   \item{settings}{List of MCMC settings used}
#'   \item{Y}{Original data array}
#'   \item{M}{Posterior draws for baseline mean M}
#'   \item{Theta}{Posterior draws for latent network state}
#'
#'   Model-specific elements include:
#'   \item{A}{Posterior draws for sender influence matrices}
#'   \item{B}{Posterior draws for receiver influence matrices}
#'   \item{sigma2, tau_A2, tau_B2, g2}{Posterior draws for variance parameters}
#'   \item{rhoA, rhoB}{AR(1) persistence parameters (dynamic model with ar1=TRUE)}
#'   \item{A_blocks}{List of regime-specific posterior mean A matrices (piecewise)}
#'   \item{time_kept}{Which time indices are stored (dynamic/lowrank/HMM)}
#'
#'   Use [summary()], [plot()], [param_summary()], and [check_convergence()]
#'   to inspect results. See model-specific vignettes for full workflows.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
#'
#' # static model with gaussian family
#' fit <- dbn(sim$Z, model = "static", family = "gaussian",
#'     nscan = 200, burn = 100, verbose = FALSE)
#'
#' # dynamic model
#' fit_dyn <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'     nscan = 200, burn = 100, verbose = FALSE)
#'
#' # piecewise model with 2 blocks
#' fit_pw <- dbn(sim$Y, model = "piecewise", blocks = 2,
#'     nscan = 200, burn = 100, verbose = FALSE)
#' }
####
dbn <- function(data,
				family = c("ordinal", "gaussian", "binary"),
				model = c("static", "dynamic", "lowrank", "hmm", "piecewise"),
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

	# piecewise model requires sufficient time points
	if (model == "piecewise" && Tt < 4) {
		cli::cli_abort(c(
			"Piecewise model requires at least 4 time points.",
			"i" = "Your data has {Tt} time point{?s}.",
			"i" = "Use {.code model = \"static\"} for short time series."
		))
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
	dots <- list(...)

	results <- switch(model,
		static = dbn_static(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, symmetric = symmetric, ...),
		dynamic = dbn_dynamic(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, symmetric = symmetric, ...),
		lowrank = dbn_lowrank(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, symmetric = symmetric, ...),
		hmm = dbn_hmm(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, symmetric = symmetric, ...),
		piecewise = {
			blocks <- dots$blocks
			if (is.null(blocks)) {
				cli::cli_abort(c(
					"Piecewise model requires {.arg blocks} parameter.",
					"i" = "Use {.code blocks = 4} for 4 equal blocks,",
					"i" = "or {.code blocks = c(25, 50, 75)} for custom boundaries,",
					"i" = "or {.code blocks = \"auto\"} for automatic selection."
				))
			}

			# handle automatic block selection
			auto_K_result <- NULL
			if (is.character(blocks) && blocks == "auto") {
				if (verbose) cli::cli_h2("Automatic Block Selection")
				auto_K_result <- select_K_auto(Y, family = family,
											   K_min = dots$K_min %||% 1L,
											   K_max = dots$K_max %||% NULL,
											   verbose = (verbose > 0))
				blocks <- auto_K_result$selected_boundaries
			}

			# filter out blocks parameter from dots
			piecewise_dots <- dots[!names(dots) %in% c("blocks", "K_min", "K_max")]

			fit <- do.call(dbn_piecewise, c(
				list(Y = Y, family = family, blocks = blocks,
					 nscan = nscan, burn = burn, odens = odens,
					 verbose = verbose, symmetric = symmetric),
				piecewise_dots
			))

			if (!is.null(auto_K_result)) {
				fit$auto_K <- auto_K_result
			}
			fit
		}
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
							(nc * p * Tt > 500)
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

			M_sum_sq <- sum(M^2, na.rm = TRUE)
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

	# convert arrays to lists for unified draws structure
	Theta_list <- lapply(seq_len(n_keep), function(i) Theta_store[, , , , i, drop = FALSE])
	for (i in seq_along(Theta_list)) {
		dim(Theta_list[[i]]) <- dim(Theta_store)[1:4]
	}
	M_list <- lapply(seq_len(n_keep), function(i) M_store[, , , i, drop = FALSE])
	for (i in seq_along(M_list)) {
		dim(M_list[[i]]) <- dim(M_store)[1:3]
	}

	# build unified draws structure (matching static/piecewise)
	draws <- list(
		theta = Theta_list,
		z = NULL,
		pars = data.frame(
			sigma2 = sigma2_store,
			tau_A2 = tau_A2_store,
			tau_B2 = tau_B2_store,
			g2 = g2_store
		),
		misc = list(
			A = A_store,
			B = B_store,
			M = M_list
		)
	)

	# add Z to draws if stored
	if (FAM$name %in% c("ordinal", "binary") && store_z) {
		Z_list <- lapply(seq_len(n_keep), function(i) Z_store[, , , , i, drop = FALSE])
		for (i in seq_along(Z_list)) {
			dim(Z_list[[i]]) <- dim(Z_store)[1:4]
		}
		draws$z <- Z_list
	}

	# add AR1 parameters if used
	if (ar1) {
		draws$pars$rhoA <- rhoA_store
		draws$pars$rhoB <- rhoB_store
	}

	# add sigma2_obs for gaussian
	if (FAM$name == "gaussian") {
		draws$pars$sigma2_obs <- sigma2_obs_store
	}

	# params data.frame for API consistency
	params <- draws$pars

	# assemble output list with unified structure
	out <- list(
		model = "dynamic",
		family = family,
		Y = Y,
		dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
					is_bipartite = is_bipartite, is_symmetric = symmetric),
		settings = list(
			nscan = nscan,
			burn = burn,
			odens = odens,
			time_thin = time_thin,
			store_z = store_z,
			draws = n_keep
		),
		meta = list(
			family = family,
			dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
						is_bipartite = is_bipartite, is_symmetric = symmetric),
			draws = n_keep,
			settings = list(nscan = nscan, burn = burn, odens = odens, time_thin = time_thin)
		),
		params = params,
		draws = draws,
		time_kept = time_keep,

		# top-level accessors for convenience
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
		time_thin = time_thin
	)

	# gaussian-specific outputs
	if (FAM$name == "gaussian") {
		out$sigma2_obs <- sigma2_obs_store
	}
	if (FAM$name %in% c("ordinal", "binary") && store_z) {
		out$Z <- Z_store
		out$Z_final <- Z_4d  # match static/piecewise naming
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
