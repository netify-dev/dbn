####
#' Posterior Visualization Functions
#'
#' @description Basic plotting utilities for posterior analysis
#' @name posterior_plots
#' @keywords internal
NULL
####

####
#' Plot Parameter Trace Plots
#'
#' @description Trace plots with running mean and posterior mean overlay
#' @param fit A dbn model fit object
#' @param pars Character vector of parameter names to plot
#' @param ncol Number of columns for multi-panel plot
#' @return A ggplot object or NULL
#' @seealso \code{\link{param_summary}}, \code{\link{plot_theta}}, \code{\link{plot_regime_probs}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot_trace(fit)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
plot_trace <- function(fit, pars = NULL, ncol = 2) {
	param_df <- param_summary(fit, probs = c(0.025, 0.5, 0.975))
	if (is.null(param_df)) {
		cli::cli_warn("No scalar parameters found to plot.")
		return(invisible(NULL))
	}

	if (is.null(pars)) {
		pars <- param_df$parameter
	} else {
		missing_pars <- setdiff(pars, param_df$parameter)
		pars <- intersect(pars, param_df$parameter)
		if (length(pars) == 0) {
			cli::cli_abort(c(
				"None of the requested parameters were found in this fit.",
				"x" = "Requested: {.val {missing_pars}}.",
				"i" = "Available parameters: {.val {param_df$parameter}}."
			))
		}
		if (length(missing_pars) > 0) {
			cli::cli_warn(c(
				"Skipping parameter{?s} not found in this fit: {.val {missing_pars}}.",
				"i" = "Available parameters: {.val {param_df$parameter}}."
			))
		}
	}

	# get trace data
	if (!is.null(fit$draws$pars)) {
		trace_data <- fit$draws$pars[, pars, drop = FALSE]
	} else {
		trace_data <- data.frame(row.names = seq_len(length(fit[[pars[1]]])))
		for (p in pars) {
			if (!is.null(fit[[p]])) trace_data[[p]] <- fit[[p]]
		}
	}
	trace_data$iteration <- seq_len(nrow(trace_data))
	####

	# ggplot2
	if (requireNamespace("ggplot2", quietly = TRUE)) {
		df_long <- do.call(rbind, lapply(pars, function(p) {
			if (p %in% names(trace_data)) {
				run_mean <- cumsum(trace_data[[p]]) / seq_along(trace_data[[p]])
				post_mean <- param_df[param_df$parameter == p, "mean"]
				data.frame(
					iteration = trace_data$iteration,
					value = trace_data[[p]],
					running_mean = run_mean,
					posterior_mean = post_mean,
					parameter = p,
					type = "trace"
				)
			}
		}))

		p <- ggplot2::ggplot(df_long, ggplot2::aes(x = iteration)) +
			ggplot2::geom_line(ggplot2::aes(y = value), color = "gray40", alpha = 0.7) +
			ggplot2::geom_line(ggplot2::aes(y = running_mean), color = "red", linewidth = 1) +
			ggplot2::geom_hline(ggplot2::aes(yintercept = posterior_mean),
				color = "blue", linetype = "dashed") +
			ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = ncol) +
			ggplot2::labs(title = "Parameter Trace Plots", x = "Iteration", y = "Value") +
			ggplot2::theme_bw() +
			ggplot2::theme(
				panel.border = ggplot2::element_blank(),
				axis.ticks = ggplot2::element_blank(),
				strip.background = ggplot2::element_rect(fill = "black", color = "black"),
				strip.text.x = ggplot2::element_text(color = "white", hjust = 0, size = 10, face = "bold")
			)
		return(p)
	}
	####

	# base R
	n_pars <- length(pars)
	nrow <- ceiling(n_pars / ncol)
	if (n_pars > 1) {
		oldpar <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))
		on.exit(par(oldpar))
	}
	for (p in pars) {
		if (p %in% names(trace_data)) {
			plot(trace_data[[p]], type = "l", main = p,
				xlab = "Iteration", ylab = "Value", col = "gray40")
			run_mean <- cumsum(trace_data[[p]]) / seq_along(trace_data[[p]])
			lines(run_mean, col = "red", lwd = 2)
			abline(h = param_df[param_df$parameter == p, "mean"], col = "blue", lty = 2)
		}
	}
	invisible(NULL)
	####
}
####

####
#' Plot Theta Heatmap
#'
#' @description Posterior mean of Theta as a heatmap
#' @param fit A dbn model fit object
#' @param time Time point to plot
#' @param rel Relation to plot
#' @param fun Summary function (default: mean)
#' @param ... Additional arguments to theta_summary
#' @return A ggplot object or NULL
#' @seealso \code{\link{theta_summary}}, \code{\link{plot_trace}}, \code{\link{param_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot_theta(fit)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
plot_theta <- function(fit, time = 1, rel = 1, fun = mean, ...) {
	theta_sum <- theta_summary(fit, fun = fun, rel = rel, time = time, ...)
	if (is.null(theta_sum) || nrow(theta_sum) == 0) {
		cli::cli_warn("No Theta values found for specified indices.")
		return(invisible(NULL))
	}

	dims <- fit$meta$dims %||% fit$dims
	n_row <- dims$n_row
	n_col <- dims$n_col

	theta_mat <- matrix(NA, n_row, n_col)
	for (row in seq_len(nrow(theta_sum))) {
		theta_mat[theta_sum$i[row], theta_sum$j[row]] <- theta_sum$value[row]
	}
	if (all(is.na(theta_mat))) return(invisible(NULL))
	theta_mat[!is.finite(theta_mat)] <- NA

	# ggplot2
	if (requireNamespace("ggplot2", quietly = TRUE)) {
		df <- expand.grid(sender = seq_len(n_row), receiver = seq_len(n_col))
		df$value <- as.vector(theta_mat)
		df <- df[!is.na(df$value), ]

		p <- ggplot2::ggplot(df, ggplot2::aes(x = receiver, y = sender, fill = value)) +
			ggplot2::geom_tile() +
			ggplot2::scale_fill_gradient2(
				low = "blue", mid = "white", high = "red",
				midpoint = 0, na.value = "grey90") +
			ggplot2::scale_y_reverse() +
			ggplot2::labs(
				title = paste0("Theta: Relation ", rel, ", Time ", time),
				x = "Receiver (j)", y = "Sender (i)",
				fill = expression(theta[ij])) +
			ggplot2::theme_bw() +
			ggplot2::theme(
				panel.border = ggplot2::element_blank(),
				axis.ticks = ggplot2::element_blank(),
				panel.grid = ggplot2::element_blank()
			)
		if (n_row == n_col) p <- p + ggplot2::coord_equal()
		return(p)
	}
	####

	# base R
	image(seq_len(n_col), seq_len(n_row), t(theta_mat),
		main = paste0("Theta: Relation ", rel, ", Time ", time),
		xlab = "Receiver (j)", ylab = "Sender (i)", col = heat.colors(100))
	invisible(NULL)
	####
}
####

####
#' Plot Regime Probabilities
#'
#' @description Stacked area plot of posterior regime probabilities over time (HMM only)
#' @param fit A dbn_hmm model fit object
#' @return A ggplot object or NULL
#' @seealso \code{\link{regime_probs}}, \code{\link{plot_trace}}, \code{\link{dbn_hmm}}
#' @examples
#' \donttest{
#' sim <- simulate_hmm_dbn(n = 6, time = 10, R = 2, seed = 1)
#' fit <- dbn(sim$Y, model = "hmm", R = 2, nscan = 200, burn = 100, verbose = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot_regime_probs(fit)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
plot_regime_probs <- function(fit) {
	probs <- regime_probs(fit)
	if (is.null(probs)) {
		cli::cli_warn("Not an HMM model or no regime information found.")
		return(invisible(NULL))
	}

	times <- seq_len(nrow(probs))
	R <- ncol(probs)

	# ggplot2
	if (requireNamespace("ggplot2", quietly = TRUE)) {
		df <- do.call(rbind, lapply(1:R, function(r) {
			data.frame(time = times, prob = probs[, r], regime = paste("Regime", r))
		}))

		p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = prob, fill = regime)) +
			ggplot2::geom_area(alpha = 0.7) +
			ggplot2::scale_fill_brewer(palette = "Set2") +
			ggplot2::labs(title = "Posterior Regime Probabilities",
				x = "Time", y = "Probability", fill = NULL) +
			ggplot2::theme_bw() +
			ggplot2::theme(
				panel.border = ggplot2::element_blank(),
				axis.ticks = ggplot2::element_blank(),
				legend.position = "top"
			)
		return(p)
	}
	####

	# base R
	plot(times, probs[, 1], type = "n", ylim = c(0, 1),
		main = "Posterior Regime Probabilities",
		xlab = "Time", ylab = "Probability")
	y_bottom <- rep(0, length(times))
	cols <- rainbow(R)
	for (r in 1:R) {
		y_top <- y_bottom + probs[, r]
		polygon(c(times, rev(times)), c(y_bottom, rev(y_top)),
			col = cols[r], border = NA)
		y_bottom <- y_top
	}
	legend("topright", paste("Regime", 1:R), fill = cols)
	invisible(NULL)
	####
}
####
