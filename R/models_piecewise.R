####
#' Piecewise-Static DBN Model
#'
#' @description block-constant influence matrices for structural change detection
#' @name models_piecewise
#' @keywords internal
NULL
####

####
#' Piecewise DBN MCMC
#'
#' @description fits DBN model with block-constant influence matrices
#' @param Y data array (n_row x n_col x p x Tt)
#' @param family distribution family: "ordinal", "gaussian", or "binary"
#' @param blocks block specification: integer, vector, or parsed block_info
#' @param nscan number of MCMC iterations after burn-in
#' @param burn burn-in iterations
#' @param odens output density (save every odens-th iteration)
#' @param seed random seed
#' @param verbose show progress
#' @param symmetric enforce B = A (symmetric networks)
#' @param previous previous fit to continue from
#' @param store_theta Logical. Store full Theta trajectory draws (default TRUE).
#'   \strong{For large networks (100+ actors), set to FALSE.} Theta storage scales as
#'   O(n^2 * T * draws) and becomes prohibitive for large networks (e.g., 200 actors
#'   with 50 time points and 500 draws requires ~40 GB). With FALSE, you retain A, B,
#'   M draws and \code{compare_blocks()} functionality but lose \code{theta_slice()},
#'   \code{theta_summary()}, and \code{posterior_predict_dbn()} with uncertainty.
#' @return list containing MCMC results
#' @seealso \code{\link{dbn}} for main dispatcher
#' @keywords internal
dbn_piecewise <- function(Y,
						  family = c("ordinal", "gaussian", "binary"),
						  blocks = 4,
						  nscan = 10000,
						  burn = 1000,
						  odens = 1,
						  seed = 6886,
						  verbose = TRUE,
						  symmetric = FALSE,
						  previous = NULL,
						  store_theta = TRUE) {

	# set up family and preprocess
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
	p <- dims$p
	Tt <- dims$Tt
	is_bipartite <- dims$is_bipartite
	nc <- n_row * n_col
	####

	# parse block structure
	if (is.list(blocks) && !is.null(blocks$K)) {
		# already parsed
		block_info <- blocks
	} else {
		block_info <- parse_blocks(blocks, Tt, n_row)
	}
	K_blocks <- block_info$K
	validate_blocks(block_info, Tt, verbose = (verbose > 0))

	if (verbose) {
		cli::cli_h3("Piecewise Model Configuration")
		cli::cli_bullets(c(
			" " = "Blocks: {K_blocks}",
			" " = "Boundaries: {paste(block_info$boundaries, collapse = ', ')}"
		))
	}
	####

	# set up ordinal sampling
	if (family == "ordinal") {
		R_flat <- array(R, c(n_row, n_col, p * Tt))
		IR <- precompute_rank_structure(R_flat, n_row, n_col, p, Tt)
		# use Gaussian approximation for rank likelihood - much faster with minimal accuracy loss
		use_approx <- should_use_gaussian_approximation(R) || (nc * p * Tt > 500)
		has_rz_cpp <- exists("rz_gaussian_approx_cpp", mode = "function")
	} else {
		IR <- pre$IR
		use_approx <- FALSE
		has_rz_cpp <- FALSE
	}
	####

	# initialize parameters
	if (!is.null(previous)) {
		# warm start from previous fit
		s2 <- tail(previous$draws$pars$s2, 1)
		t2 <- tail(previous$draws$pars$t2, 1)
		g2 <- tail(previous$draws$pars$g2, 1)
		M <- previous$draws$misc$M[[length(previous$draws$misc$M)]]

		# initialize block-specific A, B
		A_blocks <- lapply(1:K_blocks, function(k) diag(n_row))
		B_blocks <- lapply(1:K_blocks, function(k) diag(n_col))
	} else {
		s2 <- 1
		t2 <- 1
		g2 <- max(0.1, mean(M^2, na.rm = TRUE))
		A_blocks <- lapply(1:K_blocks, function(k) diag(n_row))
		B_blocks <- lapply(1:K_blocks, function(k) diag(n_col))
	}

	# initialize Theta
	Theta <- sweep(Z, c(1, 2, 3), M, "-") + rsan(dim(Z))
	####

	# allocate MCMC storage
	n_iter <- burn + nscan
	keep_idx <- seq(burn + 1, n_iter, by = odens)
	n_keep <- length(keep_idx)

	# parameter storage
	param_samples <- matrix(NA, n_keep, 3)
	colnames(param_samples) <- c("s2", "t2", "g2")
	Msave <- vector("list", n_keep)

	# block-specific A, B storage
	A_store <- vector("list", n_keep)
	B_store <- vector("list", n_keep)

	# Theta storage for actor position inference (optional for memory savings)
	if (store_theta) {
		Theta_store <- vector("list", n_keep)
	}

	if (isTRUE(verbose > 0)) {
		pw_notes <- character()
		if (.dbn_data_looks_undirected(Y, n_row, n_col) && !isTRUE(symmetric)) {
			pw_notes <- c(pw_notes,
				"Data looks undirected (symmetric adjacency) but {.code symmetric = FALSE}; {.code symmetric = TRUE} is better identified for undirected networks.")
		}
		.dbn_preflight("piecewise", n_row, n_col, p, Tt, burn, nscan, odens,
					   n_keep, notes = pw_notes)
	}

	t_start <- proc.time()[[3]]
	if (verbose) {
		cli::cli_alert_info("Running piecewise DBN MCMC ({n_iter} iterations)")
		cli::cli_progress_bar("MCMC iterations", total = n_iter)
	}
	####

	# check for C++ block update function
	has_block_update_cpp <- p == 1 && exists("piecewise_block_update_cpp", mode = "function")

	# precompute block boundaries for C++ (0-indexed)
	if (has_block_update_cpp) {
		block_starts <- as.integer(sapply(block_info$block_indices, function(x) x[1] - 1))
		block_ends <- as.integer(sapply(block_info$block_indices, function(x) x[length(x)] - 1))
	}

	# main MCMC loop
	for (iter in 1:n_iter) {

		# update latent Z (ordinal/binary families)
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
		} else if (FAM$name == "ordinal") {
			Z_flat <- array(Z, c(n_row, n_col, p * Tt))
			EZ_cube <- array(rep(as.vector(M), Tt), c(n_row, n_col, p * Tt))

			if (use_approx) {
				if (has_rz_cpp) {
					Z_flat <- rz_gaussian_approx_cpp(R_flat, Z_flat, EZ_cube)
				} else {
					Z_flat <- rz_gaussian_approx(R_flat, Z_flat, EZ_cube)
				}
			} else {
				Z_flat <- rz_fc_batch(R_flat, Z_flat, EZ_cube, IR, n_row, n_col, p, Tt)
			}
			Z <- array(Z_flat, c(n_row, n_col, p, Tt))
		}
		####

		# update baseline M
		M_sum <- apply(Z - Theta, c(1, 2, 3), sum)
		M <- M_sum / (Tt + 1 / g2) + rsan(c(n_row, n_col, p)) / sqrt(Tt + 1 / g2)

		# update g2 (M prior variance)
		g2 <- 1 / rgamma(1, (1 + nc * p) / 2, (1 + sum(M^2, na.rm = TRUE)) / 2)
		####

		# update each block: Theta, A_k, B_k
		# use C++ for single relation case (p=1), R fallback for p>1
		if (has_block_update_cpp) {
			# reshape for C++ (drop p dimension since p=1)
			Z_3d <- array(Z[, , 1, ], c(n_row, n_col, Tt))
			Theta_3d <- array(Theta[, , 1, ], c(n_row, n_col, Tt))
			M_2d <- M[, , 1]

			# call C++ block update
			cpp_result <- piecewise_block_update_cpp(
				Z_3d, Theta_3d, M_2d, A_blocks, B_blocks,
				block_starts, block_ends, s2, t2, symmetric
			)

			# extract results
			Theta[, , 1, ] <- cpp_result$Theta
			A_blocks <- cpp_result$A_blocks
			B_blocks <- cpp_result$B_blocks
			sse_total <- cpp_result$sse
		} else {
			# R fallback for p > 1 or when C++ unavailable
			sse_total <- 0
			for (k in 1:K_blocks) {
				t_idx <- block_info$block_indices[[k]]
				block_len <- length(t_idx)

				# extract block data
				Z_block <- Z[, , , t_idx, drop = FALSE]
				Theta_block <- Theta[, , , t_idx, drop = FALSE]

				# compute Y (centered Z) and X (lagged Theta) for this block
				Y_block <- sweep(Z_block, c(1, 2, 3), M, "-")

				if (block_len > 1) {
					# X = lagged Theta (with first time point from previous block or noise)
					if (k == 1) {
						X_first <- rsan(c(n_row, n_col, p, 1)) * 0.1
					} else {
						prev_block_end <- block_info$block_indices[[k - 1]]
						last_t <- prev_block_end[length(prev_block_end)]
						Theta_last <- Theta[, , , last_t, drop = FALSE]
						X_first <- sweep(Theta_last, c(1, 2, 3), M, "-")
					}
					X_block <- abind::abind(X_first, Y_block[, , , -block_len, drop = FALSE], along = 4)
				} else {
					# single time point: no dynamics
					X_block <- rsan(c(n_row, n_col, p, 1)) * 0.1
				}

				# update Theta within block
				A_k <- A_blocks[[k]]
				B_k <- B_blocks[[k]]

				# compute expected Theta from AR dynamics
				for (t_local in 1:block_len) {
					if (t_local == 1 && k == 1) {
						# first time point: prior centered at M
						Theta[, , , t_idx[t_local]] <- M + rsan(c(n_row, n_col, p))
					} else {
						# AR update: Theta_t = A * Theta_{t-1} * B' + noise
						if (t_local == 1) {
							prev_t <- block_info$block_indices[[k - 1]][block_info$lengths[k - 1]]
						} else {
							prev_t <- t_idx[t_local - 1]
						}

						# apply A and B via matrix operations for each relation
						for (j in 1:p) {
							Theta_prev_j <- Theta[, , j, prev_t]
							M_j <- M[, , j]
							E_Theta <- A_k %*% (Theta_prev_j - M_j) %*% t(B_k) + M_j

							# sample from posterior
							resid <- Z[, , j, t_idx[t_local]] - E_Theta
							post_var <- s2 / (1 + s2)
							post_mean <- E_Theta + post_var * resid
							Theta[, , j, t_idx[t_local]] <- post_mean + sqrt(post_var) * rsan(c(n_row, n_col))
						}
					}
				}

				# update A_k (sender influence)
				if (block_len > 1) {
					# collect sufficient statistics across block
					sum_YX <- matrix(0, n_row, n_row)
					sum_XX <- matrix(0, n_row, n_row)

					for (t_local in 2:block_len) {
						for (j in 1:p) {
							Y_tj <- Theta[, , j, t_idx[t_local]] - M[, , j]
							X_tj <- Theta[, , j, t_idx[t_local - 1]] - M[, , j]
							X_tj_B <- X_tj %*% t(B_k)

							sum_YX <- sum_YX + Y_tj %*% t(X_tj_B)
							sum_XX <- sum_XX + X_tj_B %*% t(X_tj_B)
						}
					}

					# posterior for A_k
					V_A <- tryCatch(
						solve(sum_XX / s2 + diag(n_row) / t2),
						error = function(e) solve(sum_XX / s2 + diag(n_row) / t2 + diag(n_row) * 1e-6)
					)
					E_A <- diag(n_row) * rnorm(1, mean(diag(A_k)), sqrt(1 / (n_row / t2 + 1)))
					A_new <- (sum_YX / s2 + E_A / t2) %*% V_A + rsan(c(n_row, n_row)) %*% chol(V_A)
					A_blocks[[k]] <- A_new
					sse_total <- sse_total + sum((A_new - E_A)^2)

					# update B_k (receiver influence)
					if (!symmetric) {
						sum_YX_B <- matrix(0, n_col, n_col)
						sum_XX_B <- matrix(0, n_col, n_col)

						for (t_local in 2:block_len) {
							for (j in 1:p) {
								Y_tj <- t(Theta[, , j, t_idx[t_local]] - M[, , j])
								X_tj <- t(Theta[, , j, t_idx[t_local - 1]] - M[, , j])
								X_tj_A <- X_tj %*% t(A_new)

								sum_YX_B <- sum_YX_B + Y_tj %*% t(X_tj_A)
								sum_XX_B <- sum_XX_B + X_tj_A %*% t(X_tj_A)
							}
						}

						V_B <- tryCatch(
							solve(sum_XX_B / s2 + diag(n_col) / t2),
							error = function(e) solve(sum_XX_B / s2 + diag(n_col) / t2 + diag(n_col) * 1e-6)
						)
						E_B <- diag(n_col) * rnorm(1, mean(diag(B_k)), sqrt(1 / (n_col / t2 + 1)))
						B_new <- (sum_YX_B / s2 + E_B / t2) %*% V_B + rsan(c(n_col, n_col)) %*% chol(V_B)
						B_blocks[[k]] <- B_new
						sse_total <- sse_total + sum((B_new - E_B)^2)
					} else {
						B_blocks[[k]] <- A_new
						sse_total <- sse_total + sum((A_new - E_A)^2)
					}
				}
			}
		}

		# update t2 (A, B precision)
		df_AB <- K_blocks * (n_row^2 + n_col^2)
		t2 <- 1 / rgamma(1, (1 + df_AB) / 2, (1 + sse_total) / 2)
		####

		# update s2 (observation variance for gaussian)
		if (FAM$name == "gaussian") {
			rss <- sum((Z - Theta - rep(M, each = Tt))^2, na.rm = TRUE)
			s2 <- 1 / rgamma(1, (1 + nc * p * Tt) / 2, (1 + rss) / 2)
		} else {
			s2 <- 1
		}
		####

		# store thinned samples
		if (iter %in% keep_idx) {
			idx <- which(keep_idx == iter)
			param_samples[idx, ] <- c(s2, t2, g2)
			Msave[[idx]] <- M

			# store block A, B matrices (deep copy to avoid reference issues)
			A_store[[idx]] <- lapply(A_blocks, function(m) m + 0)
			B_store[[idx]] <- lapply(B_blocks, function(m) m + 0)

			# store Theta trajectory for actor position inference
			if (store_theta) {
				Theta_store[[idx]] <- Theta
			}
		}
		####

		if (verbose) {
			if (iter %% verbose == 0 || iter == n_iter) {
				cli::cli_progress_update(set = iter)
			}
			.dbn_heartbeat(iter, n_iter, t_start, verbose)
		}
	}

	if (verbose) cli::cli_progress_done()
	if (verbose) {
		.dbn_closing("piecewise", t_start,
					 conv = .dbn_mixing_note(param_samples[, "s2"], "s2"))
	}
	####

	# assemble output
	misc_list <- list(
		M = Msave,
		A = A_store,
		B = B_store
	)
	if (store_theta) {
		misc_list$Theta <- Theta_store
	}

	draws <- list(
		pars = data.frame(
			s2 = param_samples[, "s2"],
			t2 = param_samples[, "t2"],
			g2 = param_samples[, "g2"]
		),
		misc = misc_list
	)

	out <- list(
		model = "piecewise",
		family = family,
		Y = Y,
		R = R,
		dims = list(
			m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
			is_bipartite = is_bipartite, is_symmetric = symmetric
		),
		settings = list(
			nscan = nscan,
			burn = burn,
			odens = odens,
			draws = n_keep
		),
		blocks = block_info,
		params = param_samples,
		M = M,
		Theta = Theta,
		A_blocks = A_blocks,
		B_blocks = B_blocks,
		draws = draws,
		Z_final = if (FAM$name %in% c("ordinal", "binary")) Z else NULL,
		M_final = M
	)

	if (!is.null(previous)) {
		prev_total <- previous$total_iter %||% (previous$settings$burn + previous$settings$nscan)
		out$total_iter <- prev_total + nscan
		out$continued_from <- prev_total
	}

	class(out) <- "dbn"
	out
	####
}
####

####
#' Simulate Data from Piecewise DBN
#'
#' @description generates synthetic data from piecewise-static model
#' @param n number of actors
#' @param time number of time points
#' @param blocks block specification
#' @param p number of relations
#' @param sigma2 observation variance
#' @param tau2 A/B prior variance
#' @param seed random seed
#' @return list with Y, true_A, true_B, true_M, block_info
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_piecewise_dbn(n = 10, time = 40, blocks = c(15, 30, 40))
#' dim(sim$Y)
#' }
simulate_piecewise_dbn <- function(n = 20, time = 50, blocks = 4,
								   p = 1, sigma2 = 1, tau2 = 0.5,
								   seed = NULL) {
	if (!is.null(seed)) set.seed(seed)

	# parse blocks
	block_info <- parse_blocks(blocks, time, n)
	K_blocks <- block_info$K

	# generate block-specific A, B matrices
	true_A <- lapply(1:K_blocks, function(k) {
		A <- diag(n) + matrix(rnorm(n * n, 0, sqrt(tau2 / n)), n, n)
		# ensure stability
		eig <- eigen(A)$values
		if (max(Mod(eig)) > 0.95) {
			A <- A * 0.9 / max(Mod(eig))
		}
		A
	})

	true_B <- lapply(1:K_blocks, function(k) {
		B <- diag(n) + matrix(rnorm(n * n, 0, sqrt(tau2 / n)), n, n)
		eig <- eigen(B)$values
		if (max(Mod(eig)) > 0.95) {
			B <- B * 0.9 / max(Mod(eig))
		}
		B
	})

	# generate baseline M
	true_M <- array(rnorm(n * n * p), dim = c(n, n, p))

	# generate Theta trajectory
	Theta <- array(0, dim = c(n, n, p, time))
	Theta[, , , 1] <- true_M + rsan(c(n, n, p))

	for (k in 1:K_blocks) {
		t_idx <- block_info$block_indices[[k]]
		A_k <- true_A[[k]]
		B_k <- true_B[[k]]

		for (t in t_idx) {
			if (t == 1) next
			for (j in 1:p) {
				Theta_prev <- Theta[, , j, t - 1] - true_M[, , j]
				Theta[, , j, t] <- A_k %*% Theta_prev %*% t(B_k) + true_M[, , j] +
					sqrt(sigma2) * rsan(c(n, n))
			}
		}
	}

	# generate observed Y
	Y <- Theta + sqrt(sigma2) * rsan(dim(Theta))

	# convert to ordinal (integer ranks)
	Y_ord <- array(0, dim = dim(Y))
	for (j in 1:p) {
		vals <- c(Y[, , j, ])
		Y_ord[, , j, ] <- array(
			as.integer(cut(vals, breaks = quantile(vals, seq(0, 1, 0.2)), include.lowest = TRUE)),
			dim = c(n, n, time)
		)
	}

	list(
		Y = Y_ord,
		Y_continuous = Y,
		Theta = Theta,
		true_A = true_A,
		true_B = true_B,
		true_M = true_M,
		block_info = block_info,
		sigma2 = sigma2,
		tau2 = tau2
	)
}
####
