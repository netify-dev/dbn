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
#' @description diagnostic plots for piecewise-static model
#' @param x piecewise dbn fit object
#' @param type plot type: "trace", "blocks", "stability", "influence"
#' @param block which block to plot (for block-specific plots)
#' @param ... additional arguments
#' @return plot object (invisibly)
#' @keywords internal
plot_piecewise <- function(x, type = c("trace", "blocks", "stability", "influence"),
						   block = NULL, ...) {
	type <- match.arg(type)

	if (type == "trace") {
		plot_piecewise_trace(x, ...)
	} else if (type == "blocks") {
		plot_piecewise_blocks(x, ...)
	} else if (type == "stability") {
		plot_piecewise_stability(x, ...)
	} else if (type == "influence") {
		plot_piecewise_influence(x, block = block, ...)
	}
}
####

####
#' Trace Plot for Piecewise Model
#'
#' @description plots MCMC traces for variance parameters
#' @param x piecewise dbn fit
#' @param ... additional arguments
#' @keywords internal
plot_piecewise_trace <- function(x, ...) {
	params <- x$draws$pars
	n_keep <- nrow(params)

	old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
	on.exit(par(old_par))

	plot(1:n_keep, params$s2, type = "l", col = "steelblue",
		 xlab = "Iteration", ylab = expression(sigma^2),
		 main = "Observation Variance")

	plot(1:n_keep, params$t2, type = "l", col = "darkgreen",
		 xlab = "Iteration", ylab = expression(tau^2),
		 main = "A/B Prior Variance")

	plot(1:n_keep, params$g2, type = "l", col = "darkred",
		 xlab = "Iteration", ylab = expression(gamma^2),
		 main = "M Prior Variance")

	# autocorrelation of s2
	if (n_keep > 50) {
		acf_vals <- acf(params$s2, lag.max = 50, plot = FALSE)
		plot(acf_vals$lag, acf_vals$acf, type = "h", col = "steelblue",
			 xlab = "Lag", ylab = "ACF", main = "s2 Autocorrelation")
		abline(h = 0, lty = 2)
	} else {
		plot.new()
		text(0.5, 0.5, "Insufficient draws for ACF")
	}

	invisible(x)
}
####

####
#' Block Comparison Plot
#'
#' @description visualizes A matrix differences across blocks
#' @param x piecewise dbn fit
#' @param ... additional arguments
#' @keywords internal
plot_piecewise_blocks <- function(x, ...) {
	K <- x$blocks$K
	n_row <- x$dims$n_row

	if (K < 2) {
		cli::cli_warn("Single block model - no block comparison to plot")
		return(invisible(x))
	}

	# get posterior mean A for each block
	A_means <- lapply(1:K, function(k) {
		A_samples <- lapply(x$draws$misc$A, function(a) a[[k]])
		Reduce(`+`, A_samples) / length(A_samples)
	})

	# compute pairwise differences
	n_pairs <- K * (K - 1) / 2
	n_cols <- min(3, n_pairs)
	n_rows <- ceiling(n_pairs / n_cols)

	old_par <- par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
	on.exit(par(old_par))

	pair_idx <- 1
	for (i in 1:(K - 1)) {
		for (j in (i + 1):K) {
			diff_mat <- A_means[[j]] - A_means[[i]]

			# heatmap of difference
			image(1:n_row, 1:n_row, t(diff_mat)[, n_row:1],
				  col = colorRampPalette(c("blue", "white", "red"))(50),
				  xlab = "Sender", ylab = "Receiver",
				  main = paste0("A[", x$blocks$names[j], "] - A[", x$blocks$names[i], "]"))

			pair_idx <- pair_idx + 1
		}
	}

	invisible(x)
}
####

####
#' Plot Block Stability
#'
#' @description shows Frobenius norm of A across blocks
#' @param x piecewise dbn fit
#' @param ... additional arguments
#' @keywords internal
plot_piecewise_stability <- function(x, ...) {
	K <- x$blocks$K

	# compute posterior mean and CI for ||A_k||_F
	A_norms <- matrix(NA, nrow = x$settings$draws, ncol = K)

	for (s in 1:x$settings$draws) {
		for (k in 1:K) {
			A_norms[s, k] <- sqrt(sum(x$draws$misc$A[[s]][[k]]^2))
		}
	}

	# block midpoints for x-axis
	midpoints <- sapply(x$blocks$block_indices, function(idx) mean(idx))

	# compute means and CIs
	means <- colMeans(A_norms)
	lower <- apply(A_norms, 2, quantile, 0.025)
	upper <- apply(A_norms, 2, quantile, 0.975)

	old_par <- par(mar = c(5, 5, 3, 1))
	on.exit(par(old_par))

	plot(midpoints, means, type = "b", pch = 19, col = "steelblue",
		 ylim = range(c(lower, upper)),
		 xlab = "Time (block midpoint)", ylab = expression(paste("||", A[k], "||"[F])),
		 main = "Influence Magnitude by Block")

	# add CI
	arrows(midpoints, lower, midpoints, upper,
		   length = 0.05, angle = 90, code = 3, col = "steelblue")

	# add block boundaries
	bounds <- x$blocks$boundaries[-c(1, length(x$blocks$boundaries))]
	if (length(bounds) > 0) {
		abline(v = bounds, lty = 2, col = "gray50")
	}

	# add block labels
	for (k in 1:K) {
		text(midpoints[k], upper[k] + 0.1 * diff(range(c(lower, upper))),
			 x$blocks$names[k], cex = 0.8)
	}

	invisible(x)
}
####

####
#' Plot Influence Matrix for Block
#'
#' @description heatmap of A matrix for specified block
#' @param x piecewise dbn fit
#' @param block block index or name
#' @param ... additional arguments
#' @keywords internal
plot_piecewise_influence <- function(x, block = 1, ...) {
	K <- x$blocks$K

	# resolve block index
	if (is.character(block)) {
		block <- which(x$blocks$names == block)
		if (length(block) == 0) {
			cli::cli_abort("Block name not found: {block}")
		}
	}
	if (block < 1 || block > K) {
		cli::cli_abort("Block index out of range: {block}")
	}

	n_row <- x$dims$n_row

	# compute posterior mean A
	A_samples <- lapply(x$draws$misc$A, function(a) a[[block]])
	A_mean <- Reduce(`+`, A_samples) / length(A_samples)

	old_par <- par(mar = c(5, 5, 4, 2))
	on.exit(par(old_par))

	image(1:n_row, 1:n_row, t(A_mean)[, n_row:1],
		  col = colorRampPalette(c("blue", "white", "red"))(50),
		  xlab = "Sender", ylab = "Receiver",
		  main = paste0("Influence Matrix A - ", x$blocks$names[block]))

	invisible(x)
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

	# data dimensions
	cat("Data:\n")
	cat("  Nodes: ", object$dims$n_row, "\n", sep = "")
	if (object$dims$is_bipartite) {
		cat("  Receivers: ", object$dims$n_col, "\n", sep = "")
	}
	cat("  Relations: ", object$dims$p, "\n", sep = "")
	cat("  Time points: ", object$dims$Tt, "\n\n", sep = "")

	# block structure
	cat("Block Structure:\n")
	cat("  Number of blocks: ", object$blocks$K, "\n", sep = "")
	cat("  Boundaries: ", paste(object$blocks$boundaries, collapse = " -> "), "\n", sep = "")
	cat("  Block lengths: ", paste(object$blocks$lengths, collapse = ", "), "\n\n", sep = "")

	# MCMC info
	cat("MCMC:\n")
	cat("  Iterations: ", object$settings$nscan, "\n", sep = "")
	cat("  Burn-in: ", object$settings$burn, "\n", sep = "")
	cat("  Saved draws: ", object$settings$draws, "\n\n", sep = "")

	# parameter summaries
	cat("Parameter Estimates (posterior mean [95% CI]):\n")
	params <- object$draws$pars

	s2_mean <- mean(params$s2)
	s2_ci <- quantile(params$s2, c(0.025, 0.975))
	cat("  s2: ", round(s2_mean, 4), " [", round(s2_ci[1], 4), ", ", round(s2_ci[2], 4), "]\n", sep = "")

	t2_mean <- mean(params$t2)
	t2_ci <- quantile(params$t2, c(0.025, 0.975))
	cat("  t2: ", round(t2_mean, 4), " [", round(t2_ci[1], 4), ", ", round(t2_ci[2], 4), "]\n", sep = "")

	g2_mean <- mean(params$g2)
	g2_ci <- quantile(params$g2, c(0.025, 0.975))
	cat("  g2: ", round(g2_mean, 4), " [", round(g2_ci[1], 4), ", ", round(g2_ci[2], 4), "]\n\n", sep = "")

	# block-specific summaries
	cat("Block-Specific Influence (||A_k||_F):\n")
	for (k in 1:object$blocks$K) {
		A_norms <- sapply(object$draws$misc$A, function(a) sqrt(sum(a[[k]]^2)))
		norm_mean <- mean(A_norms)
		norm_ci <- quantile(A_norms, c(0.025, 0.975))
		cat("  ", object$blocks$names[k], ": ", round(norm_mean, 3),
			" [", round(norm_ci[1], 3), ", ", round(norm_ci[2], 3), "]\n", sep = "")
	}

	# auto-K info if available
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
#' Simulate from Piecewise Model
#'
#' @description generates forecasts or simulations from piecewise model
#' @param object piecewise dbn fit
#' @param H forecast horizon
#' @param ndraws number of posterior draws to use
#' @param seed random seed
#' @param ... additional arguments
#' @return list with simulated values
#' @keywords internal
simulate_piecewise <- function(object, H = 10, ndraws = 100, S = NULL,
							   summary = c("mean", "none"), seed = NULL, ...) {
	if (!is.null(seed)) set.seed(seed)
	if (!is.null(S)) ndraws <- S
	summary <- match.arg(summary)

	n_row <- object$dims$n_row
	n_col <- object$dims$n_col
	p <- object$dims$p
	Tt <- object$dims$Tt
	K <- object$blocks$K

	# use last block's A, B for forecasting
	n_available <- object$settings$draws
	draw_idx <- sample(1:n_available, min(ndraws, n_available))

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
