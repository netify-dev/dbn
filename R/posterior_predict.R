####
#' Posterior Predictive Functions
#'
#' @description Functions for generating posterior predictive samples and checks
#' @name posterior_predict
#' @keywords internal
NULL
####

####
#' Generate posterior predictive samples
#'
#' @description Generate new observations from the posterior predictive distribution
#' @param fit A dbn model fit object
#' @param ndraws Number of posterior draws to use (default: 100)
#' @param seed Random seed for reproducibility
#' @param draws Specific draw indices to use (overrides ndraws)
#' @return List of predicted observations with class "dbn_ppd"
#' @seealso \code{\link{plot_ppc_ecdf}}, \code{\link{plot_ppc_density}}, \code{\link{param_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' ppd <- posterior_predict_dbn(fit, ndraws = 5)
#' }
#' @export
posterior_predict_dbn <- function(fit, ndraws = 100, seed = NULL, draws = NULL) {
	if (!is.null(seed)) set.seed(seed)

	fam <- get_family(fit)
	if (is.null(fam)) {
		cli::cli_abort("Model fit does not contain family information.")
	}

	# determine which draws to use
	if (is.null(draws)) {
		if (fit$model == "dynamic") {
			n_total_draws <- length(fit$sigma2)
		} else if (fit$model == "static") {
			n_total_draws <- nrow(fit$params)
		} else if (fit$model == "hmm") {
			n_total_draws <- length(fit$draws$pars$sigma2_proc)
		} else if (fit$model == "lowrank" || fit$model == "lowrank_accurate") {
			n_total_draws <- length(fit$sigma2)
		} else {
			n_total_draws <- fit$meta$draws %||% length(fit$draws$theta)
		}

		if (ndraws >= n_total_draws) {
			draws <- seq_len(n_total_draws)
		} else {
			draws <- sample(seq_len(n_total_draws), ndraws, replace = TRUE)
		}
	} else {
		ndraws <- length(draws)
	}
	####

	Yrep <- vector("list", ndraws)

	for (d in seq_along(draws)) {
		draw_idx <- draws[d]

		# extract theta based on model type
		if (fit$model == "dynamic") {

			# use stored FFBS theta draws (includes absorbed baseline mean)
			# only if Theta has the full time dimension (not time-thinned)
			Tt_orig <- fit$dims$Tt
			Theta_has_full_time <- !is.null(fit$Theta) && is.array(fit$Theta) &&
				length(dim(fit$Theta)) == 5 && dim(fit$Theta)[4] == Tt_orig
			if (Theta_has_full_time) {
				th <- fit$Theta[, , , , draw_idx, drop = FALSE]
				dim(th) <- dim(fit$Theta)[1:4]
			} else {
				# fallback: reconstruct from A/B/M
				if (is.null(fit$A) || is.null(fit$B) || is.null(fit$M)) {
					cli::cli_abort(c(
						"x" = "Cannot compute theta: Theta, A, B, or M not found in fit object.",
						"i" = "Ensure the model was fitted properly."
					))
				}

				A_draw <- fit$A[[draw_idx]]
				B_draw <- fit$B[[draw_idx]]

				if (is.array(fit$M) && length(dim(fit$M)) == 4) {
					M_draw <- fit$M[, , , draw_idx, drop = FALSE]
					dim(M_draw) <- dim(fit$M)[1:3]
				} else if (is.list(fit$M)) {
					M_draw <- fit$M[[draw_idx]]
				} else {
					cli::cli_abort("Unexpected M structure in dynamic model fit")
				}

				dims <- fit$dims
				n_row <- dims$n_row
				n_col <- dims$n_col
				p <- dims$p
				Tt <- dims$Tt

				th <- array(0, dim = c(n_row, n_col, p, Tt))

				for (j in 1:p) {
					th[, , j, 1] <- M_draw[, , j]
					for (t in 2:Tt) {
						t_idx <- if (dim(A_draw)[3] < Tt) {
							time_thin <- fit$time_thin %||% 1
							ceiling(t / time_thin)
						} else {
							t
						}
						th[, , j, t] <- A_draw[, , t_idx] %*% th[, , j, t-1] %*% t(B_draw[, , t_idx]) + M_draw[, , j]
					}
				}
			}

		} else if (fit$model == "static") {

			if (is.null(fit$draws$misc$B) || is.null(fit$draws$misc$M)) {
				cli::cli_abort(c(
					"x" = "Cannot compute theta: B matrices or M not found in fit object.",
					"i" = "Ensure the model was fitted properly."
				))
			}

			B <- fit$draws$misc$B[[1]][, , draw_idx]
			M <- fit$draws$misc$M[[draw_idx]]

			if (!is.null(fit$draws$z) && length(fit$draws$z) >= draw_idx) {
				Z <- fit$draws$z[[draw_idx]]
			} else if (fam$name == "gaussian") {
				Z <- fit$Y
			} else {
				cli::cli_abort(c(
					"x" = "Cannot determine Z values for theta computation.",
					"i" = "Model family not supported for on-demand theta computation."
				))
			}

			th <- compute_theta_static(B, Z, M)

		} else if (fit$model == "hmm" && !is.null(fit$draws$theta)) {

			th <- fit$draws$theta[[draw_idx]]

		} else if (fit$model == "lowrank" || fit$model == "lowrank_accurate") {

			# use stored FFBS theta draws when available
			if (!is.null(fit$draws$theta) && length(fit$draws$theta) >= draw_idx) {
				th <- fit$draws$theta[[draw_idx]]
			} else {
				# fallback: reconstruct from U/alpha/B
				if (is.null(fit$U) || is.null(fit$alpha) || is.null(fit$B)) {
					cli::cli_abort(c(
						"x" = "Cannot reconstruct theta: U, alpha, or B not found in fit object.",
						"i" = "Ensure the lowrank model was fitted properly."
					))
				}

				dims <- fit$meta$dims %||% fit$dims
				n_row <- dims$n_row
				n_col <- dims$n_col
				p <- dims$p
				T_keep <- dim(fit$alpha[[1]])[2]

				th <- array(0, dim = c(n_row, n_col, p, T_keep))

				U_s <- fit$U[[draw_idx]]
				alpha_s <- fit$alpha[[draw_idx]]
				B_s <- fit$B[[draw_idx]]

				# extract M for initialization if available
				M_draw <- NULL
				if (is.list(fit$M) && length(fit$M) >= draw_idx) {
					M_draw <- fit$M[[draw_idx]]
				} else if (!is.null(fit$draws$misc$M) && length(fit$draws$misc$M) >= draw_idx) {
					M_draw <- fit$draws$misc$M[[draw_idx]]
				}

				for (t in 1:T_keep) {
					if (length(alpha_s[, t]) == 1) {
						A_t <- alpha_s[, t] * U_s %*% t(U_s)
					} else {
						A_t <- U_s %*% diag(alpha_s[, t]) %*% t(U_s)
					}

					B_t <- B_s[, , t]

					if (t == 1) {
						for (rel in 1:p) {
							if (!is.null(M_draw)) {
								th[, , rel, t] <- M_draw[, , rel]
							}
						}
					} else {
						for (rel in 1:p) {
							th[, , rel, t] <- A_t %*% th[, , rel, t-1] %*% t(B_t)
							if (!is.null(M_draw)) {
								th[, , rel, t] <- th[, , rel, t] + M_draw[, , rel]
							}
						}
					}
				}
			}

		} else {

			if (!is.null(fit$draws$theta)) {
				th <- fit$draws$theta[[draw_idx]]
			} else {
				cli::cli_abort(c(
					"x" = "Cannot extract theta for model type '{fit$model}'.",
					"i" = "Model may not support posterior predictive checks yet."
				))
			}
		}
		####

		# extract miscellaneous parameters
		misc <- list()

		if (fit$model == "dynamic" && !is.null(fit$M)) {
			misc$M <- fit$M[, , , draw_idx, drop = FALSE]
			dim(misc$M) <- dim(fit$M)[1:3]
		} else if ((fit$model == "lowrank" || fit$model == "lowrank_accurate") && !is.null(fit$M)) {
			if (is.list(fit$M)) {
				misc$M <- fit$M[[draw_idx]]
			} else {
				misc$M <- fit$M[, , , draw_idx, drop = FALSE]
				dim(misc$M) <- dim(fit$M)[1:3]
			}
		} else if (!is.null(fit$draws$misc$M) && length(fit$draws$misc$M) >= draw_idx) {
			misc$M <- fit$draws$misc$M[[draw_idx]]
		} else if (!is.null(fit$M)) {
			misc$M <- fit$M
		}

		if (fam$name == "gaussian") {
			if ((fit$model == "dynamic" || fit$model == "lowrank" || fit$model == "lowrank_accurate")
				&& !is.null(fit$sigma2_obs)) {
				misc$sigma2_obs <- fit$sigma2_obs[draw_idx]
			} else if (!is.null(fit$draws$pars) && "sigma2_obs" %in% names(fit$draws$pars)) {
				misc$sigma2_obs <- fit$draws$pars$sigma2_obs[draw_idx]
			} else if (!is.null(fit$sigma2_obs) && length(fit$sigma2_obs) >= draw_idx) {
				misc$sigma2_obs <- fit$sigma2_obs[draw_idx]
			}
		}
		####

		Yrep[[d]] <- fam$rgen_obs(th, misc)
	}

	structure(Yrep,
		class = "dbn_ppd",
		family = fam$name,
		dims = fit$meta$dims %||% fit$dims
	)
}
####

####
#' Print method for posterior predictive distribution
#' @param x An object of class "dbn_ppd"
#' @param ... Additional arguments (currently unused)
#' @export
print.dbn_ppd <- function(x, ...) {
	cat("Posterior predictive distribution\n")
	cat("Family:", attr(x, "family"), "\n")
	cat("Number of draws:", length(x), "\n")

	dims <- attr(x, "dims")
	if (!is.null(dims)) {
		nr <- dims$n_row
		nc <- dims$n_col
		if (isTRUE(dims$is_bipartite)) {
			cat(
				"Data dimensions:",
				nr, "senders x", nc, "receivers x",
				dims$p, "relations x",
				dims$Tt, "time points\n"
			)
		} else {
			cat(
				"Data dimensions:",
				nr, "nodes x", dims$p, "relations x",
				dims$Tt, "time points\n"
			)
		}
	}

	invisible(x)
}
####

####
#' Get family object from model fit
#'
#' @description Extract or reconstruct family object from model fit
#' @param fit A dbn model fit object
#' @return A dbn_family object or NULL
#' @keywords internal
get_family <- function(fit) {
	if (!is.null(fit$family) && inherits(fit$family, "dbn_family")) {
		return(fit$family)
	}

	family_name <- fit$settings$family %||% fit$family

	if (is.null(family_name)) {
		if ("IR" %in% names(fit)) {
			family_name <- "ordinal"
		} else if (!is.null(fit$sigma2_obs)) {
			family_name <- "gaussian"
		} else {
			warning("Cannot determine family type from model fit")
			return(NULL)
		}
	}

	switch(family_name,
		ordinal = family_ordinal(),
		gaussian = family_gaussian(),
		binary = family_binary(),
		NULL
	)
}
####

####
#' Posterior predictive ECDF plot
#'
#' @description Plot empirical CDFs of observed vs replicated data
#' @param fit A dbn model fit object
#' @param ppd Posterior predictive samples (from posterior_predict_dbn)
#' @param ndraws_plot Number of replicates to show (default: 20)
#' @param alpha Transparency for replicate lines (default: 0.3)
#' @param rel Relation index to plot (default: 1)
#' @param time Time index to plot (default: all)
#' @param Y_obs Observed data array (required)
#' @return A ggplot object
#' @seealso \code{\link{posterior_predict_dbn}}, \code{\link{plot_ppc_density}}, \code{\link{param_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' ppd <- posterior_predict_dbn(fit, ndraws = 5)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot_ppc_ecdf(fit, ppd, Y_obs = sim$Y)
#' }
#' @export
plot_ppc_ecdf <- function(fit, ppd = NULL, ndraws_plot = 20,
						  alpha = 0.3, rel = 1, time = NULL, Y_obs = NULL) {
	if (is.null(ppd)) {
		ppd <- posterior_predict_dbn(fit, ndraws = ndraws_plot)
	}

	if (is.null(Y_obs)) {
		cli::cli_abort("Observed data {.arg Y_obs} must be provided.")
	}

	dims <- fit$meta$dims %||% fit$dims
	if (is.null(time)) {
		time <- seq_len(dims$Tt)
	}

	obs_vals <- c(Y_obs[, , rel, time])
	obs_vals <- obs_vals[!is.na(obs_vals)]

	if (length(obs_vals) > 5e5) {
		obs_vals <- sample(obs_vals, 5e5)
	}

	vals <- sort(unique(obs_vals))

	df <- data.frame(
		value = vals,
		ecdf = ecdf(obs_vals)(vals),
		type = "Observed",
		set = "Observed"
	)

	draws_show <- sample(seq_along(ppd), min(ndraws_plot, length(ppd)))

	for (d in draws_show) {
		rep_vals <- c(ppd[[d]][, , rel, time])
		rep_vals <- rep_vals[!is.na(rep_vals)]

		if (length(rep_vals) > 5e5) {
			rep_vals <- sample(rep_vals, 5e5)
		}

		df_rep <- data.frame(
			value = vals,
			ecdf = ecdf(rep_vals)(vals),
			type = "Replicated",
			set = paste0("Rep", d)
		)

		df <- rbind(df, df_rep)
	}

	# ggplot2 version
	if (requireNamespace("ggplot2", quietly = TRUE)) {
		p <- ggplot2::ggplot(df, ggplot2::aes(
			x = value, y = ecdf,
			group = set,
			color = type,
			linewidth = type,
			alpha = type
		)) +
			ggplot2::geom_line() +
			ggplot2::scale_color_manual(values = c(
				"Observed" = "black",
				"Replicated" = "grey60"
			)) +
			ggplot2::scale_linewidth_manual(values = c(
				"Observed" = 1.2,
				"Replicated" = 0.8
			)) +
			ggplot2::scale_alpha_manual(values = c(
				"Observed" = 1,
				"Replicated" = alpha
			)) +
			ggplot2::labs(
				title = paste0("Posterior Predictive Check: Relation ", rel),
				x = "Y",
				y = "ECDF",
				color = NULL, linewidth = NULL, alpha = NULL
			) +
			ggplot2::theme_minimal() +
			ggplot2::theme(legend.position = "bottom")

		return(p)
	}
	####

	# base R fallback
	plot(ecdf(obs_vals),
		main = paste0("PPC: Relation ", rel),
		xlab = "Y", ylab = "ECDF", lwd = 2
	)

	for (d in draws_show) {
		rep_vals <- c(ppd[[d]][, , rel, time])
		rep_vals <- rep_vals[!is.na(rep_vals)]
		lines(ecdf(rep_vals), col = adjustcolor("grey60", alpha.f = alpha))
	}

	legend("bottomright", c("Observed", "Replicated"),
		lwd = c(2, 1), col = c("black", "grey60")
	)

	invisible(NULL)
	####
}
####

####
#' Posterior predictive density plot
#'
#' @description Plot density of observed vs replicated data
#' @param fit A dbn model fit object
#' @param ppd Posterior predictive samples
#' @param rel Relation index
#' @param time Time indices
#' @param Y_obs Observed data array (required)
#' @return A ggplot object or NULL
#' @seealso \code{\link{posterior_predict_dbn}}, \code{\link{plot_ppc_ecdf}}, \code{\link{param_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' ppd <- posterior_predict_dbn(fit, ndraws = 5)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot_ppc_density(fit, ppd, Y_obs = sim$Y)
#' }
#' @export
plot_ppc_density <- function(fit, ppd = NULL, rel = 1, time = NULL, Y_obs = NULL) {
	if (is.null(ppd)) {
		ppd <- posterior_predict_dbn(fit, ndraws = 50)
	}

	if (is.null(Y_obs)) {
		cli::cli_abort("Observed data {.arg Y_obs} must be provided.")
	}
	dims <- fit$meta$dims %||% fit$dims

	if (is.null(time)) {
		time <- seq_len(dims$Tt)
	}

	obs_vals <- c(Y_obs[, , rel, time])
	obs_vals <- obs_vals[!is.na(obs_vals)]

	# discrete data dispatches to bar plot
	if (length(unique(obs_vals)) <= 20) {
		return(plot_ppc_bars(fit, ppd, rel, time, Y_obs))
	}

	# ggplot2 version
	if (requireNamespace("ggplot2", quietly = TRUE)) {
		df <- data.frame(value = obs_vals, type = "Observed")

		for (d in seq_along(ppd)) {
			rep_vals <- c(ppd[[d]][, , rel, time])
			rep_vals <- rep_vals[!is.na(rep_vals)]
			df <- rbind(df, data.frame(
				value = rep_vals,
				type = paste0("Rep", d)
			))
		}

		p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = type == "Observed")) +
			ggplot2::geom_density(alpha = 0.5) +
			ggplot2::scale_fill_manual(
				values = c("TRUE" = "black", "FALSE" = "grey60"),
				labels = c("TRUE" = "Observed", "FALSE" = "Replicated")
			) +
			ggplot2::labs(
				title = paste0("Posterior Predictive Density: Relation ", rel),
				x = "Y", y = "Density", fill = NULL
			) +
			ggplot2::theme_minimal()

		return(p)
	}
	####

	# base R fallback
	plot(density(obs_vals),
		main = paste0("PPC Density: Relation ", rel),
		xlab = "Y", ylab = "Density", lwd = 2
	)

	for (d in seq_along(ppd)) {
		rep_vals <- c(ppd[[d]][, , rel, time])
		rep_vals <- rep_vals[!is.na(rep_vals)]
		lines(density(rep_vals), col = adjustcolor("grey60", alpha.f = 0.3))
	}

	invisible(NULL)
	####
}
####

####
#' Posterior predictive bar plot for discrete data
#'
#' @description Bar plot comparison for discrete outcomes
#' @param fit A dbn model fit object
#' @param ppd Posterior predictive samples
#' @param rel Relation index
#' @param time Time indices
#' @param Y_obs Observed data array (required)
#' @return A ggplot object or NULL
#' @keywords internal
plot_ppc_bars <- function(fit, ppd, rel = 1, time = NULL, Y_obs = NULL) {
	if (is.null(Y_obs)) {
		cli::cli_abort("Observed data {.arg Y_obs} must be provided.")
	}
	dims <- fit$meta$dims %||% fit$dims

	if (is.null(time)) {
		time <- seq_len(dims$Tt)
	}

	obs_vals <- c(Y_obs[, , rel, time])
	obs_vals <- obs_vals[!is.na(obs_vals)]
	obs_freq <- table(obs_vals) / length(obs_vals)

	rep_freq_list <- lapply(ppd, function(yrep) {
		rep_vals <- c(yrep[, , rel, time])
		rep_vals <- rep_vals[!is.na(rep_vals)]
		table(rep_vals) / length(rep_vals)
	})

	all_vals <- sort(unique(c(
		names(obs_freq),
		unlist(lapply(rep_freq_list, names))
	)))

	rep_freq_mean <- numeric(length(all_vals))
	rep_freq_lower <- numeric(length(all_vals))
	rep_freq_upper <- numeric(length(all_vals))

	for (i in seq_along(all_vals)) {
		val <- all_vals[i]
		freqs <- sapply(rep_freq_list, function(f) f[val] %||% 0)
		rep_freq_mean[i] <- mean(freqs, na.rm = TRUE)
		rep_freq_lower[i] <- quantile(freqs, 0.05, na.rm = TRUE)
		rep_freq_upper[i] <- quantile(freqs, 0.95, na.rm = TRUE)
	}

	# ggplot2 version
	if (requireNamespace("ggplot2", quietly = TRUE)) {
		df <- data.frame(
			value = rep(all_vals, 2),
			freq = c(as.numeric(obs_freq[all_vals]), rep_freq_mean),
			type = rep(c("Observed", "Replicated"), each = length(all_vals)),
			lower = c(rep(NA, length(all_vals)), rep_freq_lower),
			upper = c(rep(NA, length(all_vals)), rep_freq_upper)
		)
		df$freq[is.na(df$freq)] <- 0

		p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = freq, fill = type)) +
			ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
			ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
				position = ggplot2::position_dodge(0.9),
				width = 0.25, na.rm = TRUE
			) +
			ggplot2::scale_fill_manual(values = c(
				"Observed" = "black",
				"Replicated" = "grey60"
			)) +
			ggplot2::labs(
				title = paste0("Posterior Predictive Check: Relation ", rel),
				x = "Y", y = "Frequency", fill = NULL
			) +
			ggplot2::theme_minimal()

		return(p)
	}
	####

	# base R fallback
	barplot(rbind(as.numeric(obs_freq[all_vals]), rep_freq_mean),
		beside = TRUE, names.arg = all_vals,
		col = c("black", "grey60"),
		main = paste0("PPC: Relation ", rel),
		xlab = "Y", ylab = "Frequency"
	)
	legend("topright", c("Observed", "Replicated"),
		fill = c("black", "grey60")
	)

	invisible(NULL)
	####
}
####
