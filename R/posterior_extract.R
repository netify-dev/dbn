####
#' Posterior Extraction and Summary Functions
#'
#' @description Core functions for extracting and summarizing posterior samples
#' @name posterior_extract
#' @keywords internal
NULL
####

####
#' Extract Theta slices from posterior draws
#'
#' @description Extract specific slices of Theta arrays from posterior draws
#' @param fit A dbn model fit object
#' @param draws Integer vector of draw indices to extract
#' @param i Row indices (sender nodes)
#' @param j Column indices (receiver nodes)
#' @param rel Relation indices
#' @param time Time indices
#' @return List of Theta slices
#' @seealso \code{\link{theta_summary}}, \code{\link{theta_credible}}, \code{\link{param_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' ts <- theta_slice(fit, time = 1)
#' }
#' @export
theta_slice <- function(fit, draws = NULL, i = NULL, j = NULL, rel = NULL, time = NULL) {

	# find theta draws - different models store them differently
	theta_draws <- NULL
	if (!is.null(fit$draws$theta)) {
		theta_draws <- fit$draws$theta
	} else if (!is.null(fit$draws$misc$Theta)) {
		# piecewise model stores at draws$misc$Theta
		theta_draws <- fit$draws$misc$Theta
	}

	if (is.null(theta_draws)) {
		warning("Model fit does not contain theta draws")
		return(NULL)
	}

	dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims

	if (is.null(i)) i <- seq_len(dims$n_row)
	if (is.null(j)) j <- seq_len(dims$n_col)
	if (is.null(rel)) rel <- seq_len(dims$p)
	if (is.null(time)) {
		time_dim <- dims$Tt
		time <- seq_len(time_dim)
	}

	if (is.null(draws)) {
		n_draws <- length(theta_draws)
		draws <- seq_len(n_draws)
	}

	out <- lapply(theta_draws[draws], function(th) {
		th[i, j, rel, time, drop = FALSE]
	})

	# single value: return as a vector
	if (length(i) == 1 && length(j) == 1 && length(rel) == 1 && length(time) == 1) {
		return(unlist(out))
	}

	# single dyad/rel with multiple times: return as a matrix
	if (length(i) == 1 && length(j) == 1 && length(rel) == 1 && length(time) > 1) {
		return(do.call(rbind, lapply(out, as.vector)))
	}

	out
}
####

####
#' Summarize Theta over posterior draws
#'
#' @description Compute summary statistics for Theta parameters
#' @param fit A dbn model fit object
#' @param fun Summary function (default: mean)
#' @param draws Draw indices (default: all)
#' @param i Row indices
#' @param j Column indices
#' @param rel Relation indices
#' @param time Time indices
#' @param chunk Chunk size for memory-efficient processing
#' @return Data frame with summarized values
#' @seealso \code{\link{theta_slice}}, \code{\link{theta_credible}}, \code{\link{plot_theta}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' ts <- theta_summary(fit)
#' }
#' @export
theta_summary <- function(fit, fun = mean,
						  draws = NULL, i = NULL, j = NULL,
						  rel = NULL, time = NULL, chunk = 20) {
	if (is.null(draws)) {
		# find theta draws - different models store them differently
		if (!is.null(fit$draws$theta)) {
			n_draws <- length(fit$draws$theta)
		} else if (!is.null(fit$draws$misc$Theta)) {
			n_draws <- length(fit$draws$misc$Theta)
		} else {
			n_draws <- fit$meta$draws %||% 0
		}
		draws <- seq_len(n_draws)
	}

	# non-linear functions need all draws at once
	if (!identical(fun, mean) && chunk < length(draws)) {
		chunk <- length(draws)
	}

	dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims

	out <- NULL
	for (block in split(draws, ceiling(seq_along(draws) / chunk))) {
		slices <- theta_slice(fit, block, i, j, rel, time)

		if (is.null(slices)) {
			return(NULL)
		}

		# single value
		if (is.vector(slices) && !is.list(slices)) {
			res <- fun(slices)
			if (length(res) == 1) {
				df <- data.frame(
					i = i, j = j, rel = rel, time = time,
					value = res, .w_chunk = length(block)
				)
			} else {
				df <- data.frame(
					i = rep(i, length(res)),
					j = rep(j, length(res)),
					rel = rep(rel, length(res)),
					time = rep(time, length(res)),
					value = res, .w_chunk = length(block)
				)
			}
			out <- rbind(out, df)
			next
		}
		####

		# single dyad/rel with multiple times
		if (is.matrix(slices)) {
			vals <- apply(slices, 2, fun)
			time_seq <- if (is.null(time)) seq_len(length(vals)) else time
			df <- data.frame(
				i = rep(i, length(vals)),
				j = rep(j, length(vals)),
				rel = rep(rel, length(vals)),
				time = time_seq,
				value = vals, .w_chunk = length(block)
			)
			out <- rbind(out, df)
			next
		}
		####

		# general array
		if (requireNamespace("abind", quietly = TRUE)) {
			arr <- do.call(abind::abind, c(slices, list(along = 5)))
		} else {
			d <- dim(slices[[1]])
			arr <- array(unlist(slices), c(d, length(slices)))
		}

		res <- apply(arr, 1:4, fun)

		if (is.matrix(res) || (is.array(res) && length(dim(res)) > 4)) {
			df_list <- list()
			n_vals <- if (is.matrix(res)) ncol(res) else dim(res)[5]

			for (v in 1:n_vals) {
				if (is.matrix(res)) {
					res_v <- res[, v]
					dim(res_v) <- dim(arr)[1:4]
				} else {
					res_v <- res[, , , , v]
				}
				df_v <- suppressWarnings(as.data.frame.table(res_v, responseName = "value"))
				dim_names <- c("i", "j", "rel", "time")[1:length(dim(res_v))]
				names(df_v) <- c(dim_names, "value")
				for (col in dim_names) {
					if (is.factor(df_v[[col]])) {
						df_v[[col]] <- as.numeric(df_v[[col]])
					}
				}
				if (!"i" %in% names(df_v)) df_v$i <- i
				if (!"j" %in% names(df_v)) df_v$j <- j
				if (!"rel" %in% names(df_v)) df_v$rel <- rel
				if (!"time" %in% names(df_v)) df_v$time <- time
				df_list[[v]] <- df_v
			}
			df <- do.call(rbind, df_list)
		} else {
			df <- suppressWarnings(as.data.frame.table(res, responseName = "value"))
			names(df) <- c("i", "j", "rel", "time", "value")
			for (col in c("i", "j", "rel", "time")) {
				if (col %in% names(df) && is.factor(df[[col]])) {
					df[[col]] <- as.numeric(df[[col]])
				}
			}
		}

		df$.w_chunk <- length(block)
		out <- rbind(out, df)
		####
	}

	# combine chunked results
	if (!is.null(out) && nrow(out) > 0) {
		key_cols <- intersect(c("i", "j", "rel", "time"), names(out))

		if (length(key_cols) > 0 && "value" %in% names(out) && chunk < length(draws)) {
			if (length(key_cols) == 1) {
				keys <- out[[key_cols[1]]]
			} else {
				keys <- do.call(paste, c(out[key_cols], sep = "_"))
			}

			if (anyDuplicated(keys)) {
				if (identical(fun, mean)) {
					out$num <- out$value * out$.w_chunk
					out$den <- out$.w_chunk
					agg <- aggregate(cbind(num, den) ~ .,
						data = out[c(key_cols, "num", "den")],
						FUN = sum, na.action = na.pass
					)
					agg$value <- agg$num / agg$den
					out <- agg[c(key_cols, "value")]
				} else {
					cli::cli_abort("Chunked processing with non-linear `fun` not supported; set `chunk = length(draws)`.")
				}
			}
		}
	}
	####

	if (!is.null(out) && ".w_chunk" %in% names(out)) {
		out$.w_chunk <- NULL
	}

	out
}
####

####
#' Summarize Scalar Parameters
#'
#' @description Returns posterior mean, standard deviation, and quantiles
#'   for the scalar variance parameters estimated by the model. These
#'   typically include:
#'   - `sigma2` or `s2`: process noise variance
#'   - `tau_A2` / `tau_B2` or `t2`: innovation variance for A/B
#'   - `g2`: latent variance
#'   - `rhoA` / `rhoB`: AR(1) persistence (dynamic model with `ar1 = TRUE`)
#'   - `sigma2_obs`: observation variance (gaussian family)
#' @param fit A dbn model fit object (output from [dbn()])
#' @param probs Quantile probabilities (default: 5th, 50th, 95th percentiles)
#' @return Data frame with columns: `parameter`, `mean`, `sd`, and one
#'   column per requested quantile
#' @seealso \code{\link{theta_summary}}, \code{\link{plot_trace}}, \code{\link{derive_draws}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' ps <- param_summary(fit)
#' }
#' @export
param_summary <- function(fit, probs = c(0.025, 0.5, 0.975)) {
	if (!is.null(fit$draws$pars)) {
		pars <- fit$draws$pars
	} else {
		pars <- NULL

		scalar_pars <- c(
			"sigma2", "sigma2_proc", "sigma2_obs",
			"tau_A2", "tau_B2", "g2", "rho_A", "rho_B"
		)

		for (par in scalar_pars) {
			if (!is.null(fit[[par]])) {
				if (is.null(pars)) {
					pars <- data.frame(fit[[par]])
					names(pars) <- par
				} else {
					pars[[par]] <- fit[[par]]
				}
			}
		}
	}

	if (is.null(pars)) {
		warning("No scalar parameters found in model fit")
		return(NULL)
	}

	quants <- as.data.frame(t(apply(pars, 2, quantile, probs = probs, na.rm = TRUE)))
	quant_names <- paste0("q", round(100 * probs, 1))
	names(quants) <- quant_names

	quants$parameter <- rownames(quants)
	rownames(quants) <- NULL
	quants$mean <- apply(pars, 2, mean, na.rm = TRUE)
	quants$sd <- apply(pars, 2, sd, na.rm = TRUE)
	quants <- quants[, c("parameter", "mean", "sd", quant_names)]

	quants
}
####

####
#' Summarize Baseline Mean M
#'
#' @description Computes posterior summaries of the baseline mean matrix M,
#'   which captures persistent dyad-specific tendencies (e.g., stable
#'   alliances or rivalries) that do not change over time. Returns a
#'   data frame with one row per dyad, suitable for plotting or comparison
#'   against known ground truth in simulation studies.
#' @param fit A fitted `dbn` object returned by [dbn()].
#' @param fun Summary function applied across posterior draws (default:
#'   [mean]). Use [median] for a robust alternative, or a custom function.
#' @param draws Integer vector of posterior draw indices to use. If `NULL`
#'   (default), all available draws are used.
#' @param rel Integer vector of relation indices to summarize. If `NULL`
#'   (default), all relations are included.
#' @param chunk Integer controlling memory-efficient processing. Draws are
#'   processed in blocks of this size. Only relevant for very large fits.
#' @return Data frame with columns `i` (sender), `j` (receiver), `rel`,
#'   and `value` (the summary statistic).
#' @seealso \code{\link{param_summary}}, \code{\link{theta_summary}}, \code{\link{derive_draws}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' ls <- latent_summary(fit)
#' }
#' @export
latent_summary <- function(fit, fun = mean, draws = NULL, rel = NULL, chunk = 20) {
	if (is.null(fit$draws$misc$M) && is.null(fit$M)) {
		warning("No M arrays found in model fit")
		return(NULL)
	}

	if (is.null(draws)) {
		n_draws <- fit$meta$draws %||% length(fit$draws$misc$M %||% 1)
		draws <- seq_len(n_draws)
	}

	# non-linear functions need all draws at once
	if (!identical(fun, mean) && chunk < length(draws)) {
		chunk <- length(draws)
	}

	if (!is.null(fit$draws$misc$M)) {
		M_list <- fit$draws$misc$M
	} else {
		M_list <- replicate(length(draws), fit$M, simplify = FALSE)
	}

	dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims

	out <- NULL
	blocks <- split(draws, ceiling(seq_along(draws) / chunk))

	for (idx in seq_along(blocks)) {
		block <- blocks[[idx]]
		m_block <- M_list[block]

		if (length(m_block) == 0) {
			next
		}

		if (requireNamespace("abind", quietly = TRUE)) {
			arr <- do.call(abind::abind, c(m_block, list(along = 4)))
		} else {
			d <- dim(m_block[[1]])
			arr <- array(unlist(m_block), c(d, length(m_block)))
		}

		if (!is.null(rel)) {
			arr <- arr[, , rel, , drop = FALSE]
			res <- apply(arr, 1:2, fun)
		} else {
			res <- apply(arr, 1:3, fun)
		}

		if (is.null(res) || length(res) == 0) {
			next
		}

		df <- NULL
		tryCatch(
			{
				df <- as.data.frame.table(res, responseName = "value")
			},
			error = function(e) {
				warning("Failed to convert array to data frame: ", e$message)
				df <- NULL
			}
		)

		if (is.null(df) || ncol(df) == 0) {
			next
		}

		n_vars <- ncol(df) - 1

		if (!is.null(rel)) {
			if (n_vars == 2) {
				names(df) <- c("i", "j", "value")
				df$rel <- rel
			} else {
				var_names <- paste0("Var", 1:n_vars)
				names(df) <- c(var_names, "value")
			}
		} else {
			if (n_vars == 3) {
				names(df) <- c("i", "j", "rel", "value")
			} else if (n_vars == 2) {
				names(df) <- c("i", "j", "value")
				df$rel <- 1
			} else {
				var_names <- paste0("Var", 1:n_vars)
				names(df) <- c(var_names, "value")
			}
		}

		for (col in c("i", "j", "rel")) {
			if (col %in% names(df) && is.factor(df[[col]])) {
				df[[col]] <- as.numeric(df[[col]])
			}
		}

		df$.w_chunk <- length(block)
		out <- rbind(out, df)
	}

	# combine chunked results
	if (!is.null(out) && nrow(out) > 0) {
		key_cols <- intersect(c("i", "j", "rel"), names(out))

		if (length(key_cols) > 0 && "value" %in% names(out) && chunk < length(draws)) {
			if (length(key_cols) == 1) {
				keys <- out[[key_cols[1]]]
			} else {
				keys <- do.call(paste, c(out[key_cols], sep = "_"))
			}

			if (anyDuplicated(keys)) {
				if (identical(fun, mean)) {
					out$num <- out$value * out$.w_chunk
					out$den <- out$.w_chunk
					agg <- aggregate(cbind(num, den) ~ .,
						data = out[c(key_cols, "num", "den")],
						FUN = sum, na.action = na.pass
					)
					agg$value <- agg$num / agg$den
					out <- agg[c(key_cols, "value")]
				} else {
					cli::cli_abort("Chunked processing with non-linear `fun` not supported; set `chunk = length(draws)`.")
				}
			}
		}
	}
	####

	if (!is.null(out) && ".w_chunk" %in% names(out)) {
		out$.w_chunk <- NULL
	}

	out
}
####

####
#' Extract regime probabilities for HMM models
#'
#' @description Compute posterior probabilities of regime assignments
#' @param fit A dbn_hmm model fit object
#' @return Matrix of regime probabilities (T x R) or NULL
#' @seealso \code{\link{plot_regime_probs}}, \code{\link{param_summary}}, \code{\link{dbn_hmm}}
#' @examples
#' \donttest{
#' sim <- simulate_hmm_dbn(n = 6, time = 10, R = 2, seed = 1)
#' fit <- dbn(sim$Y, model = "hmm", R = 2, nscan = 200, burn = 100, verbose = FALSE)
#' rp <- regime_probs(fit)
#' }
#' @export
regime_probs <- function(fit) {
	if (fit$model != "hmm") {
		return(NULL)
	}

	if (!is.null(fit$draws$misc$S)) {
		S_list <- fit$draws$misc$S
	} else if (!is.null(fit$S)) {
		S_list <- fit$S
	} else {
		return(NULL)
	}

	if (length(S_list) == 0) {
		return(NULL)
	}

	S_mat <- do.call(cbind, S_list)

	Tt <- nrow(S_mat)
	R <- fit$meta$R %||% fit$settings$R

	if (is.null(R)) {
		R <- max(unlist(S_list))
	}

	probs <- matrix(0, nrow = Tt, ncol = R)
	colnames(probs) <- paste0("Regime", 1:R)
	rownames(probs) <- paste0("Time", 1:Tt)

	for (r in 1:R) {
		probs[, r] <- rowMeans(S_mat == r)
	}

	probs
}
####

####
#' Derive new quantities from posterior draws
#'
#' @description Apply a transformation function `fun` to each posterior draw.
#'   The function receives a list with one element per stored draw component
#'   (`theta`, `z`, `pars`, and any `misc` slots), and should return the
#'   derived quantity for that draw.
#' @param fit A `dbn` model fit object
#' @param fun Function applied to each posterior draw; receives a per-draw
#'   list of components and returns the derived quantity.
#' @param draws Integer vector of draw indices to process. `NULL` (default)
#'   processes all stored draws.
#' @param chunk Chunk size for memory efficiency
#' @param name Name for the derived quantity
#' @return List of derived quantities with class `"dbn_derived"`.
#' @seealso \code{\link{param_summary}}, \code{\link{theta_slice}}, \code{\link{theta_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' dd <- derive_draws(fit, function(d) sum(d$theta^2))
#' length(dd)
#' }
#' @export
derive_draws <- function(fit, fun, draws = NULL, chunk = 20, name = "derived") {
	if (is.null(fit$draws)) {
		warning("Model fit does not contain posterior draws in expected format")
		return(NULL)
	}

	if (is.null(draws)) {
		# find n_draws from various possible storage locations
		n_draws <- 0
		if (!is.null(fit$meta$draws) && fit$meta$draws > 0) {
			n_draws <- fit$meta$draws
		} else if (length(fit$draws$theta) > 0) {
			n_draws <- length(fit$draws$theta)
		} else if (length(fit$draws$misc$Theta) > 0) {
			n_draws <- length(fit$draws$misc$Theta)
		} else if (length(fit$draws$misc$A) > 0) {
			n_draws <- length(fit$draws$misc$A)
		} else if (length(fit$draws$misc$B) > 0) {
			n_draws <- length(fit$draws$misc$B)
		}
		if (n_draws == 0) {
			warning("Could not determine number of draws")
			return(NULL)
		}
		draws <- seq_len(n_draws)
	}

	out <- vector("list", length(draws))

	for (block in split(draws, ceiling(seq_along(draws) / chunk))) {
		for (k_idx in seq_along(block)) {
			draw_idx <- block[k_idx]
			k <- which(draws == draw_idx)

			draw <- list()
			if (!is.null(fit$draws$theta)) draw$theta <- fit$draws$theta[[draw_idx]]
			if (!is.null(fit$draws$z)) draw$z <- fit$draws$z[[draw_idx]]
			if (!is.null(fit$draws$pars)) draw$pars <- fit$draws$pars[draw_idx, ]
			if (!is.null(fit$draws$misc)) {
				for (nm in names(fit$draws$misc)) {
					if (length(fit$draws$misc[[nm]]) >= draw_idx) {
						draw[[nm]] <- fit$draws$misc[[nm]][[draw_idx]]
					}
				}
			}

			out[[k]] <- fun(draw)
		}
	}

	# simplify if all elements are numeric vectors of the same length
	if (length(out) > 0 && is.numeric(out[[1]])) {
		lens <- sapply(out, length)
		if (all(lens == lens[1])) {
			out <- do.call(rbind, out)
		}
	}

	structure(out, class = "dbn_derived", name = name)
}
####

####
#' Print method for derived quantities
#' @param x An object of class "dbn_derived"
#' @param ... Additional arguments (currently unused)
#' @export
print.dbn_derived <- function(x, ...) {
	name <- attr(x, "name")
	cat("Derived posterior quantity:", name, "\n")
	cat("Number of draws:", length(x), "\n")

	if (length(x) > 0) {
		cat("First draw dimensions:", paste(dim(x[[1]]), collapse = " x "), "\n")
	}

	invisible(x)
}
####

####
#' Posterior credible intervals for Theta
#'
#' @description Compute pointwise credible intervals for each dyad over time
#' @param fit A dbn model fit object
#' @param probs Probability levels (default: 90 percent interval + median)
#' @param i Row indices (sender nodes, default: all)
#' @param j Column indices (receiver nodes, default: all)
#' @param rel Relation index (default: 1)
#' @param time Time indices (default: all)
#' @return Data frame with columns: i, j, rel, time, mean, lower, median, upper
#' @seealso \code{\link{theta_summary}}, \code{\link{theta_slice}}, \code{\link{edge_prob}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' tc <- theta_credible(fit)
#' }
#' @export
theta_credible <- function(fit, probs = c(0.05, 0.5, 0.95),
						   i = NULL, j = NULL, rel = 1, time = NULL) {
	if (is.null(fit$draws$theta)) {
		warning("Model fit does not contain theta draws")
		return(NULL)
	}

	dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims
	n_row <- dims$n_row
	n_col <- dims$n_col

	if (is.null(i)) i <- seq_len(n_row)
	if (is.null(j)) j <- seq_len(n_col)
	if (is.null(time)) {
		time_dim <- dims$Tt
		time <- seq_len(time_dim)
	}

	n_draws <- length(fit$draws$theta)
	grid <- expand.grid(i = i, j = j, time = time)
	n_cells <- nrow(grid)

	res_mean  <- numeric(n_cells)
	res_lower <- numeric(n_cells)
	res_med   <- numeric(n_cells)
	res_upper <- numeric(n_cells)

	for (k in seq_len(n_cells)) {
		ii <- grid$i[k]
		jj <- grid$j[k]
		tt <- grid$time[k]
		vals <- vapply(fit$draws$theta, function(th) th[ii, jj, rel, tt],
					   numeric(1))
		q <- quantile(vals, probs = probs, na.rm = TRUE)
		res_mean[k]  <- mean(vals, na.rm = TRUE)
		res_lower[k] <- q[1]
		res_med[k]   <- q[2]
		res_upper[k] <- q[3]
	}

	data.frame(
		i      = grid$i,
		j      = grid$j,
		rel    = rel,
		time   = grid$time,
		mean   = res_mean,
		lower  = res_lower,
		median = res_med,
		upper  = res_upper
	)
}
####

####
#' Network-level posterior summary
#'
#' @description Compute network-level statistics at each time point across
#'   posterior draws (e.g., mean edge weight, density of positive edges)
#' @param fit A dbn model fit object
#' @param stat Statistic to compute: "mean" (average theta), "density"
#'   (fraction of positive values), or "strength" (mean absolute theta)
#' @param rel Relation index (default: 1)
#' @param time Time indices (default: all)
#' @return Data frame with columns: time, mean, lower, upper
#' @seealso \code{\link{edge_prob}}, \code{\link{theta_summary}}, \code{\link{param_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' ns <- network_summary(fit)
#' }
#' @export
network_summary <- function(fit, stat = c("mean", "density", "strength"),
							rel = 1, time = NULL) {
	stat <- match.arg(stat)

	if (is.null(fit$draws$theta)) {
		warning("Model fit does not contain theta draws")
		return(NULL)
	}

	dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	is_bipartite <- isTRUE(dims$is_bipartite)

	if (is.null(time)) {
		time_dim <- dims$Tt
		time <- seq_len(time_dim)
	}

	n_draws <- length(fit$draws$theta)

	stat_fn <- switch(stat,
		mean = function(mat) mean(mat, na.rm = TRUE),
		density = function(mat) mean(mat > 0, na.rm = TRUE),
		strength = function(mat) mean(abs(mat), na.rm = TRUE)
	)

	result <- data.frame(time = time, mean = NA_real_,
						 lower = NA_real_, upper = NA_real_)

	for (ti in seq_along(time)) {
		tt <- time[ti]
		vals <- vapply(fit$draws$theta, function(th) {
			mat <- th[, , rel, tt]
			if (!is_bipartite) diag(mat) <- NA
			stat_fn(mat)
		}, numeric(1))
		result$mean[ti]  <- mean(vals)
		result$lower[ti] <- quantile(vals, 0.05)
		result$upper[ti] <- quantile(vals, 0.95)
	}

	result
}
####

####
#' Posterior edge probability
#'
#' @description For each dyad, compute the posterior probability that
#'   the latent theta is positive (or exceeds a threshold)
#' @param fit A dbn model fit object
#' @param threshold Threshold value (default: 0)
#' @param rel Relation index (default: 1)
#' @param time Time index (default: last time point)
#' @return Matrix of posterior probabilities (n_row x n_col)
#' @seealso \code{\link{network_summary}}, \code{\link{theta_credible}}, \code{\link{theta_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' ep <- edge_prob(fit)
#' }
#' @export
edge_prob <- function(fit, threshold = 0, rel = 1, time = NULL) {
	if (is.null(fit$draws$theta)) {
		warning("Model fit does not contain theta draws")
		return(NULL)
	}

	dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims
	n_row <- dims$n_row
	n_col <- dims$n_col

	if (is.null(time)) {
		time_dim <- dims$Tt
		time <- time_dim
	}

	n_draws <- length(fit$draws$theta)
	prob_mat <- matrix(0, n_row, n_col)

	for (s in seq_len(n_draws)) {
		prob_mat <- prob_mat + (fit$draws$theta[[s]][, , rel, time] > threshold)
	}
	prob_mat <- prob_mat / n_draws

	if (!isTRUE(dims$is_bipartite)) diag(prob_mat) <- NA

	prob_mat
}
####
