####
#' Piecewise Model Methods
#'
#' @description plot, summary, and prediction methods for piecewise DBN models
#' @name piecewise_methods
#' @keywords internal
NULL
####

####
#' Plot Piecewise DBN Results
#'
#' @description ggplot-based diagnostic plots for the piecewise-static
#'   model. Selects between four panel types: variance traces, pairwise
#'   block-difference heatmaps, per-block influence-norm dotplot, or a
#'   single-block influence heatmap.
#' @param x piecewise dbn fit object
#' @param type plot type: "trace", "blocks", "stability", "influence"
#' @param block which block to plot (for `type = "influence"`)
#' @param ... additional arguments
#' @return A `ggplot` (or `gridExtra::arrangeGrob`) object.
#' @keywords internal
plot_piecewise <- function(x, type = c("trace", "blocks", "stability", "influence"),
						   block = 1, ...) {
	type <- match.arg(type)
	if (!requireNamespace("ggplot2", quietly = TRUE))
		cli::cli_abort("Package {.pkg ggplot2} is required for piecewise plotting.")
	switch(type,
		trace     = plot_piecewise_trace(x, ...),
		blocks    = plot_piecewise_blocks(x, ...),
		stability = plot_piecewise_stability(x, ...),
		influence = plot_piecewise_influence(x, block = block, ...)
	)
}
####

####
#' Trace Plot for Piecewise Model
#'
#' @description ggplot panel of variance-parameter MCMC traces.
#' @param x piecewise dbn fit
#' @param ... unused
#' @return A `ggplot` object.
#' @keywords internal
plot_piecewise_trace <- function(x, ...) {
	params <- x$draws$pars
	if (is.null(params) || nrow(params) == 0L)
		cli::cli_abort("No scalar parameter draws found on the piecewise fit.")
	n_keep <- nrow(params)
	long <- data.frame(
		iter = rep(seq_len(n_keep), 3L),
		value = c(params$s2, params$t2, params$g2),
		param = factor(rep(c("sigma2", "tau2", "g2"), each = n_keep),
		               levels = c("sigma2", "tau2", "g2"))
	)
	ggplot2::ggplot(long, ggplot2::aes(x = .data$iter, y = .data$value)) +
		ggplot2::geom_line(linewidth = 0.4, colour = "steelblue") +
		ggplot2::facet_wrap(~ .data$param, scales = "free_y", ncol = 1L,
			labeller = ggplot2::labeller(param = c(
				sigma2 = "sigma^2 (observation variance)",
				tau2   = "tau^2 (A/B prior variance)",
				g2     = "g^2 (M prior variance)"
			))) +
		ggplot2::labs(x = "Iteration", y = NULL,
		              title = "Piecewise model: variance parameter traces") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(),
		               strip.background = ggplot2::element_rect(fill = "black", colour = "black"),
		               strip.text = ggplot2::element_text(colour = "white", hjust = 0))
}
####

####
#' Block Comparison Plot
#'
#' @description Per-pair heatmap of posterior-mean differences between
#'   block-specific A matrices.
#' @param x piecewise dbn fit
#' @param ... unused
#' @return A `ggplot` object.
#' @keywords internal
plot_piecewise_blocks <- function(x, ...) {
	K <- x$blocks$K
	if (K < 2) {
		cli::cli_warn("Single block model -- no block comparison to plot.")
		return(invisible(NULL))
	}
	n_row <- x$dims$n_row
	# posterior mean A per block
	A_means <- lapply(seq_len(K), function(k) {
		A_samples <- lapply(x$draws$misc$A, function(a) a[[k]])
		Reduce(`+`, A_samples) / length(A_samples)
	})
	block_names <- x$blocks$names
	pairs <- utils::combn(K, 2L)
	df <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(p) {
		i <- pairs[1L, p]; j <- pairs[2L, p]
		diff_mat <- A_means[[j]] - A_means[[i]]
		data.frame(
			sender = rep(seq_len(n_row), each = n_row),
			receiver = rep(seq_len(n_row), times = n_row),
			diff = as.vector(diff_mat),
			pair = sprintf("A[%s] - A[%s]", block_names[j], block_names[i])
		)
	}))
	lim <- max(abs(df$diff), na.rm = TRUE)
	ggplot2::ggplot(df, ggplot2::aes(x = .data$sender, y = .data$receiver, fill = .data$diff)) +
		ggplot2::geom_tile() +
		ggplot2::scale_y_reverse(breaks = pretty(seq_len(n_row))) +
		ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
		                              midpoint = 0, limits = c(-lim, lim),
		                              name = expression(Delta * A)) +
		ggplot2::facet_wrap(~ .data$pair) +
		ggplot2::coord_equal() +
		ggplot2::labs(x = "Sender", y = "Receiver",
		              title = "Block-to-block differences in posterior mean A") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(),
		               strip.background = ggplot2::element_rect(fill = "black", colour = "black"),
		               strip.text = ggplot2::element_text(colour = "white", hjust = 0))
}
####

####
#' Plot Block Stability
#'
#' @description Per-block Frobenius norm of A with 95% credible
#'   intervals, plotted at the block midpoints.
#' @param x piecewise dbn fit
#' @param ... unused
#' @return A `ggplot` object.
#' @keywords internal
plot_piecewise_stability <- function(x, ...) {
	K <- x$blocks$K
	A_norms <- matrix(NA_real_, nrow = x$settings$draws, ncol = K)
	for (s in seq_len(x$settings$draws)) {
		for (k in seq_len(K)) {
			A_norms[s, k] <- sqrt(sum(x$draws$misc$A[[s]][[k]]^2))
		}
	}
	midpoints <- vapply(x$blocks$block_indices, function(idx) mean(idx), numeric(1))
	df <- data.frame(
		block = factor(x$blocks$names, levels = x$blocks$names),
		midpoint = midpoints,
		mean = colMeans(A_norms),
		lo = apply(A_norms, 2L, stats::quantile, 0.025),
		hi = apply(A_norms, 2L, stats::quantile, 0.975)
	)
	bounds <- x$blocks$boundaries[-c(1L, length(x$blocks$boundaries))]
	p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$midpoint, y = .data$mean)) +
		ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
		                       width = 0.6, colour = "steelblue") +
		ggplot2::geom_point(colour = "steelblue", size = 2.5) +
		ggplot2::geom_line(colour = "steelblue", linewidth = 0.4) +
		ggplot2::geom_text(ggplot2::aes(y = .data$hi, label = .data$block),
		                   vjust = -0.8, size = 3.2) +
		ggplot2::labs(x = "Time (block midpoint)",
		              y = expression("||" * A[k] * "||"[F]),
		              title = "Influence magnitude by block") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank())
	if (length(bounds) > 0L) {
		p <- p + ggplot2::geom_vline(xintercept = bounds,
		                              linetype = "dashed", colour = "grey50")
	}
	p
}
####

####
#' Plot Influence Matrix for Block
#'
#' @description Heatmap of the posterior mean A for a single block.
#' @param x piecewise dbn fit
#' @param block block index or name
#' @param ... unused
#' @return A `ggplot` object.
#' @keywords internal
plot_piecewise_influence <- function(x, block = 1, ...) {
	K <- x$blocks$K
	if (is.character(block)) {
		block_idx <- which(x$blocks$names == block)
		if (length(block_idx) == 0L)
			cli::cli_abort("Block name not found: {block}")
		block <- block_idx
	}
	if (block < 1L || block > K)
		cli::cli_abort("Block index out of range: {block}")
	n_row <- x$dims$n_row
	A_samples <- lapply(x$draws$misc$A, function(a) a[[block]])
	A_mean <- Reduce(`+`, A_samples) / length(A_samples)
	df <- data.frame(
		sender = rep(seq_len(n_row), each = n_row),
		receiver = rep(seq_len(n_row), times = n_row),
		value = as.vector(A_mean)
	)
	lim <- max(abs(df$value), na.rm = TRUE)
	ggplot2::ggplot(df, ggplot2::aes(x = .data$sender, y = .data$receiver, fill = .data$value)) +
		ggplot2::geom_tile() +
		ggplot2::scale_y_reverse(breaks = pretty(seq_len(n_row))) +
		ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
		                              midpoint = 0, limits = c(-lim, lim),
		                              name = "A entry") +
		ggplot2::coord_equal() +
		ggplot2::labs(x = "Sender", y = "Receiver",
		              title = paste0("Influence matrix A -- ", x$blocks$names[block])) +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank())
}
####

####
#' Summary for Piecewise DBN
#'
#' @description generates summary statistics for piecewise model
#' @param object piecewise dbn fit
#' @param ... additional arguments
#' @return summary object
#' @keywords internal
summary_piecewise <- function(object, ...) {
	cat("Piecewise-Static DBN Model Summary\n")
	cat(rep("=", 40), "\n\n", sep = "")
	cat("Data:\n")
	cat("  Nodes: ", object$dims$n_row, "\n", sep = "")
	if (object$dims$is_bipartite) {
		cat("  Receivers: ", object$dims$n_col, "\n", sep = "")
	}
	cat("  Relations: ", object$dims$p, "\n", sep = "")
	cat("  Time points: ", object$dims$Tt, "\n\n", sep = "")
	cat("Block Structure:\n")
	cat("  Number of blocks: ", object$blocks$K, "\n", sep = "")
	cat("  Boundaries: ", paste(object$blocks$boundaries, collapse = " -> "), "\n", sep = "")
	cat("  Block lengths: ", paste(object$blocks$lengths, collapse = ", "), "\n\n", sep = "")
	cat("MCMC:\n")
	cat("  Iterations: ", object$settings$nscan, "\n", sep = "")
	cat("  Burn-in: ", object$settings$burn, "\n", sep = "")
	cat("  Saved draws: ", object$settings$draws, "\n\n", sep = "")
	cat("Parameter Estimates (posterior mean [95% CI]):\n")
	params <- object$draws$pars
	s2_mean <- mean(params$s2); s2_ci <- stats::quantile(params$s2, c(0.025, 0.975))
	cat("  s2: ", round(s2_mean, 4), " [", round(s2_ci[1], 4), ", ", round(s2_ci[2], 4), "]\n", sep = "")
	t2_mean <- mean(params$t2); t2_ci <- stats::quantile(params$t2, c(0.025, 0.975))
	cat("  t2: ", round(t2_mean, 4), " [", round(t2_ci[1], 4), ", ", round(t2_ci[2], 4), "]\n", sep = "")
	g2_mean <- mean(params$g2); g2_ci <- stats::quantile(params$g2, c(0.025, 0.975))
	cat("  g2: ", round(g2_mean, 4), " [", round(g2_ci[1], 4), ", ", round(g2_ci[2], 4), "]\n\n", sep = "")
	cat("Block-Specific Influence (||A_k||_F):\n")
	for (k in seq_len(object$blocks$K)) {
		A_norms <- vapply(object$draws$misc$A, function(a) sqrt(sum(a[[k]]^2)), numeric(1))
		norm_mean <- mean(A_norms); norm_ci <- stats::quantile(A_norms, c(0.025, 0.975))
		cat("  ", object$blocks$names[k], ": ", round(norm_mean, 3),
		    " [", round(norm_ci[1], 3), ", ", round(norm_ci[2], 3), "]\n", sep = "")
	}
	if (!is.null(object$auto_K)) {
		cat("\nAutomatic Block Selection:\n")
		cat("  Stage 1 suggestion: K = ", object$auto_K$stage1_suggestion, "\n", sep = "")
		cat("  Selected K: ", object$auto_K$selected_K, "\n", sep = "")
		cat("  WAIC comparison: ", paste(round(object$auto_K$waics, 2), collapse = ", "), "\n", sep = "")
	}
	invisible(object)
}
####

####

####
#' Simulate from Piecewise Model
#'
#' @description Generates forecasts or simulations from a piecewise DBN model.
#' @param object Piecewise `dbn` fit.
#' @param H Forecast horizon.
#' @param draws Number of posterior draws to use.
#' @param summary `"mean"` (default) or `"none"` (return all draws).
#' @param seed Random seed.
#' @param ... Additional arguments (currently unused).
#' @return Array of simulated values; 4D if `summary = "mean"`, 5D
#'   (with a trailing draws dimension) otherwise.
#' @keywords internal
simulate_piecewise <- function(object, H = 10, draws = 100,
							   summary = c("mean", "none"), seed = NULL, ...) {
	if (!is.null(seed)) set.seed(seed)
	summary <- match.arg(summary)

	n_row <- object$dims$n_row
	n_col <- object$dims$n_col
	p <- object$dims$p
	Tt <- object$dims$Tt
	K <- object$blocks$K

	n_available <- object$settings$draws
	draw_idx <- sample(1:n_available, min(draws, n_available))

	Y_sim <- array(NA, dim = c(n_row, n_col, p, H, length(draw_idx)))

	for (i in seq_along(draw_idx)) {
		s <- draw_idx[i]

		# get last block's A, B
		A_last <- object$draws$misc$A[[s]][[K]]
		B_last <- object$draws$misc$B[[s]][[K]]
		M <- object$draws$misc$M[[s]]
		s2 <- object$draws$pars$s2[s]

		# get last Theta from stored draws
		if (!is.null(object$draws$misc$Theta) && length(object$draws$misc$Theta) >= s) {
			Theta_prev <- object$draws$misc$Theta[[s]][, , , Tt, drop = FALSE]
			dim(Theta_prev) <- c(n_row, n_col, p)
		} else {
			cli::cli_abort(c(
				"Theta draws not available for forecasting.",
				"i" = "Refit with {.code store_theta = TRUE} to enable forecasting."
			))
		}

		for (h in 1:H) {
			# forecast Theta
			Theta_h <- array(0, dim = c(n_row, n_col, p))
			for (j in 1:p) {
				Theta_dev <- Theta_prev[, , j] - M[, , j]
				Theta_h[, , j] <- A_last %*% Theta_dev %*% t(B_last) + M[, , j] +
					sqrt(s2) * rsan(c(n_row, n_col))
			}

			Y_sim[, , , h, i] <- Theta_h
			Theta_prev <- Theta_h
		}
	}

	if (summary == "mean") {
		return(apply(Y_sim, 1:4, mean))
	}
	Y_sim
}
####

####
#' Compare Influence Matrices Across Blocks
#'
#' @description computes posterior probability of difference between blocks
#' @param fit piecewise dbn fit object
#' @param blocks pair of block indices to compare (default: adjacent blocks)
#' @param parameter which parameter to compare: "A" or "B"
#' @param threshold minimum difference to be considered meaningful
#' @return list with comparison results
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_piecewise_dbn(n = 10, time = 40, blocks = c(15, 30, 40))
#' fit <- dbn(sim$Y, model = "piecewise", blocks = 3,
#'            nscan = 200, burn = 100, verbose = FALSE)
#' compare_blocks(fit)
#' }
compare_blocks <- function(fit, blocks = NULL, parameter = "A", threshold = 0.1) {
	if (fit$model != "piecewise") {
		cli::cli_abort("compare_blocks() requires a piecewise model fit")
	}

	K <- fit$blocks$K
	if (K < 2) {
		cli::cli_warn("Single block model - no comparison possible")
		return(NULL)
	}

	# default: all adjacent pairs
	if (is.null(blocks)) {
		pairs <- lapply(1:(K - 1), function(k) c(k, k + 1))
	} else if (is.numeric(blocks) && length(blocks) == 2) {
		pairs <- list(blocks)
	} else {
		cli::cli_abort("blocks must be NULL or a pair of indices")
	}

	# get parameter samples
	if (parameter == "A") {
		param_list <- fit$draws$misc$A
	} else if (parameter == "B") {
		param_list <- fit$draws$misc$B
	} else {
		cli::cli_abort("parameter must be 'A' or 'B'")
	}

	n_draws <- length(param_list)
	results <- list()

	for (pair in pairs) {
		k1 <- pair[1]
		k2 <- pair[2]

		# compute difference norm for each draw
		diff_norms <- numeric(n_draws)
		for (s in 1:n_draws) {
			diff_mat <- param_list[[s]][[k2]] - param_list[[s]][[k1]]
			diff_norms[s] <- sqrt(sum(diff_mat^2))
		}

		# posterior probability of meaningful difference
		prob_diff <- mean(diff_norms > threshold)

		# mean and CI of difference norm
		mean_diff <- mean(diff_norms)
		ci_diff <- quantile(diff_norms, c(0.025, 0.975))

		pair_name <- paste0(fit$blocks$names[k1], "_vs_", fit$blocks$names[k2])
		results[[pair_name]] <- list(
			blocks = c(k1, k2),
			block_names = c(fit$blocks$names[k1], fit$blocks$names[k2]),
			mean_diff = mean_diff,
			ci = ci_diff,
			prob_above_threshold = prob_diff,
			threshold = threshold,
			diff_norms = diff_norms
		)
	}

	# attach class so print method controls output
	class(results) <- "dbn_block_comparison"
	attr(results, "parameter") <- parameter

	# print summary
	print(results)

	invisible(results)
}
####

####
#' @export
print.dbn_block_comparison <- function(x, ...) {
	parameter <- attr(x, "parameter") %||% "A"
	cli::cli_h3("Block Comparison Results ({parameter})")
	for (nm in names(x)) {
		r <- x[[nm]]
		cli::cli_alert_info("{r$block_names[1]} vs {r$block_names[2]}: ||\u0394{parameter}|| = {round(r$mean_diff, 3)} [{round(r$ci[1], 3)}, {round(r$ci[2], 3)}]")
		cli::cli_alert("  P(||\u0394{parameter}|| > {r$threshold}) = {round(r$prob_above_threshold, 3)}")
	}
	invisible(x)
}
####
