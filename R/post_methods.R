####
# post_methods.R
# plot, summary, convergence, and prediction methods for dbn model objects.
# see the S3 dispatch methods for the generic plot.dbn entry point
####

####
#' Plot Static DBN Results
#'
#' @description Creates diagnostic plots for static model results using ggplot2
#' @param results Output from dbn_static()
#' @param alpha Significance level for edge detection (default 0.01)
#' @return A ggplot2 object with multiple panels
#' @keywords internal
plot_static <- function(results, alpha = 0.01) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort(c("Package {.pkg ggplot2} is required for this function.", "i" = "Install with {.code install.packages(\"ggplot2\")}"))
	# suppress R CMD check notes for NSE variables
	iteration <- value <- parameter <- s2 <- t2 <- g2 <- NULL
	from_x <- from_y <- to_x <- to_y <- color <- x <- y <- label <- NULL

	if (is.null(results$params)) {
		cli::cli_abort("results$params is NULL -- nothing to plot")
	}

	# trace plots for scalar parameters
	params_df <- as.data.frame(results$params)
	params_df$iteration <- seq_len(nrow(params_df))

	params_long <- data.frame(
		iteration = rep(params_df$iteration, 3),
		value = c(params_df$s2, params_df$t2, params_df$g2),
		parameter = factor(rep(c("s2", "t2", "g2"), each = nrow(params_df)))
	)

	p_traces <- ggplot2::ggplot(params_long, ggplot2::aes(x = iteration, y = value, color = parameter)) +
		ggplot2::geom_line() +
		ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 1) +
		ggplot2::labs(title = "Parameter Traces", x = "Iteration", y = "Value") +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			legend.position = "none",
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0)
		)

	# posterior histograms
	p_hist_s2 <- ggplot2::ggplot(params_df, ggplot2::aes(x = s2)) +
		ggplot2::geom_histogram(bins = 30) +
		ggplot2::labs(title = "s2 Posterior", x = "s2", y = "Count") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

	p_hist_t2 <- ggplot2::ggplot(params_df, ggplot2::aes(x = t2)) +
		ggplot2::geom_histogram(bins = 30) +
		ggplot2::labs(title = "t2 Posterior", x = "t2", y = "Count") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

	p_hist_g2 <- ggplot2::ggplot(params_df, ggplot2::aes(x = g2)) +
		ggplot2::geom_histogram(bins = 30) +
		ggplot2::labs(title = "g2 Posterior", x = "g2", y = "Count") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

	# network plot of the first B matrix
	if (length(results$B) >= 1 && !is.null(results$B[[1]])) {
		B_post <- results$B[[1]]
		EB <- apply(B_post, c(1, 2), mean)
		diag(EB) <- 0
		SB <- apply(B_post, c(1, 2), sd)

		# find edges with large posterior mean and narrow credible intervals
		BIG <- abs(EB) > quantile(abs(EB), 1 - alpha, na.rm = TRUE)
		SIG <- abs(EB) > qnorm(1 - alpha / 2) * SB
		if (all(SB == 0, na.rm = TRUE)) {
			SIG[] <- FALSE
		}
		BSG <- BIG & SIG

		csig <- (rowSums(BSG) + colSums(BSG)) > 0

		if (sum(csig) > 1) {
			BSG_sub <- BSG[csig, csig]
			EB_sub <- EB[csig, csig]
			nodes_sub <- which(csig)

			# arrange nodes in a circle
			n_nodes <- nrow(BSG_sub)
			angles <- seq(0, 2 * pi, length.out = n_nodes + 1)[1:n_nodes]
			xy <- cbind(x = cos(angles), y = sin(angles))

			# build edge data frame (vectorized)
			idx <- which(BSG_sub != 0, arr.ind = TRUE)
			if (nrow(idx) > 0) {
				edges <- data.frame(
					from_x = xy[idx[, 1], 1], from_y = xy[idx[, 1], 2],
					to_x   = xy[idx[, 2], 1], to_y   = xy[idx[, 2], 2],
					weight = EB_sub[idx],
					color  = ifelse(EB_sub[idx] > 0, "Positive", "Negative")
				)
			} else {
				edges <- data.frame(from_x = numeric(0), from_y = numeric(0),
									to_x = numeric(0), to_y = numeric(0),
									weight = numeric(0), color = character(0))
			}

			nodes <- data.frame(
				x = xy[, 1], y = xy[, 2],
				label = nodes_sub
			)

			p_network <- ggplot2::ggplot() +
				ggplot2::geom_segment(
					data = edges,
					ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y, color = color),
					arrow = grid::arrow(length = grid::unit(0.2, "cm"))
				) +
				ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y), size = 4) +
				ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y, label = label), size = 3) +
				ggplot2::scale_color_manual(values = c("Positive" = "green3", "Negative" = "red3")) +
				ggplot2::labs(title = "Network B[[1]]", color = "Edge Type") +
				ggplot2::theme_bw() +
				ggplot2::theme(
					panel.border = ggplot2::element_blank(),
					axis.ticks = ggplot2::element_blank(),
					axis.text = ggplot2::element_blank(),
					axis.title = ggplot2::element_blank(),
					panel.grid = ggplot2::element_blank()
				) +
				ggplot2::coord_fixed()
		} else {
			p_network <- ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No significant edges for B[[1]]") +
				ggplot2::theme_void()
		}
	} else {
		cli::cli_warn("No B samples available - skipping network panel")
		p_network <- ggplot2::ggplot() +
			ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No B samples available") +
			ggplot2::theme_void()
	}
	####

	# arrange all panels into a single figure
	if (requireNamespace("gridExtra", quietly = TRUE)) {
		tryCatch(
			{
				gridExtra::grid.arrange(p_network, p_traces, p_hist_s2, p_hist_t2, p_hist_g2,
					layout_matrix = rbind(c(1, 1, 2), c(3, 4, 5))
				)
			},
			error = function(e) {
				cli::cli_warn("Could not arrange plots with gridExtra: {e$message}")
				list(
					network = p_network, traces = p_traces,
					hist_s2 = p_hist_s2, hist_t2 = p_hist_t2, hist_g2 = p_hist_g2
				)
			}
		)
	} else {
		list(
			network = p_network, traces = p_traces,
			hist_s2 = p_hist_s2, hist_t2 = p_hist_t2, hist_g2 = p_hist_g2
		)
	}
}
####

####
#' Plot Dynamic DBN Results
#'
#' @description Creates diagnostic plots for dynamic model results using ggplot2
#' @param results Output from dbn_dynamic()
#' @param time_points Which time points to display (default: start, middle, end)
#' @return A list of ggplot2 objects
#' @keywords internal
plot_dynamic <- function(results, time_points = NULL) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort(c("Package {.pkg ggplot2} is required for this function.", "i" = "Install with {.code install.packages(\"ggplot2\")}"))
	# suppress R CMD check notes for NSE variables
	iteration <- value <- parameter <- g2 <- NULL

	dims <- results$meta$dims %||% results$dims
	Tt <- dims$Tt

	if (is.null(time_points)) {
		time_points <- c(1, floor(Tt / 2), Tt)
	}

	# build data frame of scalar parameter traces
	if (!is.null(results$draws$pars)) {
		pars <- results$draws$pars
		trace_df <- data.frame(
			iteration = seq_len(nrow(pars)),
			sigma2 = pars$sigma2,
			tauA2 = pars$tau_A2,
			tauB2 = pars$tau_B2
		)
		if ("g2" %in% names(pars)) trace_df$g2 <- pars$g2
		if ("rho_A" %in% names(pars)) trace_df$rhoA <- pars$rho_A
		if ("rho_B" %in% names(pars)) trace_df$rhoB <- pars$rho_B
	} else {
		trace_df <- data.frame(
			iteration = seq_along(results$sigma2),
			sigma2 = results$sigma2,
			tauA2 = results$tau_A2 %||% results$tauA2,
			tauB2 = results$tau_B2 %||% results$tauB2
		)
		if (!is.null(results$g2)) trace_df$g2 <- results$g2
		if (!is.null(results$rho_A)) trace_df$rhoA <- results$rho_A
		if (!is.null(results$rho_B)) trace_df$rhoB <- results$rho_B
	}
	####

	# reshape to long format for ggplot2
	plist <- c("sigma^2", "tauA^2", "tauB^2")
	trace_long <- data.frame(
		iteration = rep(trace_df$iteration, 3),
		value = c(trace_df$sigma2, trace_df$tauA2, trace_df$tauB2),
		parameter = factor(rep(plist, each = nrow(trace_df)))
	)

	if (!is.null(results$g2)) {
		trace_long <- rbind(
			trace_long,
			data.frame(
				iteration = trace_df$iteration,
				value = trace_df$g2,
				parameter = "g^2"
			)
		)
		plist <- c(plist, "g^2")
	}

	p_traces <- ggplot2::ggplot(trace_long, ggplot2::aes(x = iteration, y = value)) +
		ggplot2::geom_line() +
		ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 1) +
		ggplot2::labs(title = "Parameter Traces", x = "Iteration", y = "Value") +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0)
		)

	# distribution of A matrix entries at selected time points
	if (!is.null(results$draws$misc$A)) {
		last_A <- results$draws$misc$A[[length(results$draws$misc$A)]]
	} else {
		last_A <- results$A[[length(results$A)]]
	}

	time_idx_avail <- seq_len(dim(last_A)[3])
	time_points <- intersect(time_points, time_idx_avail)
	if (length(time_points) == 0) {
		time_points <- time_idx_avail[c(1, floor(length(time_idx_avail) / 2), length(time_idx_avail))]
	}

	hist_list <- vector("list", length(time_points))
	for (i in seq_along(time_points)) {
		t_idx <- time_points[i]
		orig_t <- if (!is.null(results$settings$time_thin) && results$settings$time_thin > 1) {
			(t_idx - 1) * results$settings$time_thin + 1
		} else {
			t_idx
		}
		hist_list[[i]] <- data.frame(
			value = as.vector(last_A[, , t_idx]),
			time = paste("t =", orig_t)
		)
	}
	A_hist_data <- do.call(rbind, hist_list)

	p_A_hist <- ggplot2::ggplot(A_hist_data, ggplot2::aes(x = value)) +
		ggplot2::geom_histogram(bins = 20) +
		ggplot2::facet_wrap(~time, ncol = length(time_points)) +
		ggplot2::labs(title = "Distribution of A Matrix Elements", x = "Value", y = "Count") +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0)
		)
	####

	# AR(1) parameter traces
	if ("rhoA" %in% names(trace_df) && "rhoB" %in% names(trace_df)) {
		rho_df <- data.frame(
			iteration = trace_df$iteration,
			rhoA = trace_df$rhoA,
			rhoB = trace_df$rhoB
		)

		rho_long <- data.frame(
			iteration = rep(rho_df$iteration, 2),
			value = c(rho_df$rhoA, rho_df$rhoB),
			parameter = factor(rep(c("rhoA", "rhoB"), each = nrow(rho_df)))
		)

		p_rho <- ggplot2::ggplot(rho_long, ggplot2::aes(x = iteration, y = value, color = parameter)) +
			ggplot2::geom_line() +
			ggplot2::facet_wrap(~parameter, ncol = 2) +
			ggplot2::labs(title = "AR(1) Parameters", x = "Iteration", y = "Value") +
			ggplot2::theme_bw() +
			ggplot2::theme(
				panel.border = ggplot2::element_blank(),
				axis.ticks = ggplot2::element_blank(),
				legend.position = "none",
				strip.background = ggplot2::element_rect(fill = "black", color = "black"),
				strip.text.x = ggplot2::element_text(color = "white", hjust = 0)
			)
	} else if ("rhoA" %in% names(trace_df)) {
		p_rho <- ggplot2::ggplot(
			data.frame(iteration = trace_df$iteration, rhoA = trace_df$rhoA),
			ggplot2::aes(x = iteration, y = rhoA)
		) +
			ggplot2::geom_line() +
			ggplot2::labs(title = "AR(1) Parameter - rhoA", x = "Iteration", y = "rhoA") +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
	} else if ("rhoB" %in% names(trace_df)) {
		p_rho <- ggplot2::ggplot(
			data.frame(iteration = trace_df$iteration, rhoB = trace_df$rhoB),
			ggplot2::aes(x = iteration, y = rhoB)
		) +
			ggplot2::geom_line() +
			ggplot2::labs(title = "AR(1) Parameter - rhoB", x = "Iteration", y = "rhoB") +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
	} else {
		p_rho <- NULL
	}
	####

	# g2 (mean variance) trace
	if ("g2" %in% names(trace_df)) {
		p_g2 <- ggplot2::ggplot(
			data.frame(iteration = trace_df$iteration, g2 = trace_df$g2),
			ggplot2::aes(x = iteration, y = g2)
		) +
			ggplot2::geom_line() +
			ggplot2::labs(title = "g^2 (tau_mu^2) Trace", x = "Iteration", y = "g^2") +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
	} else {
		p_g2 <- NULL
	}
	####

	# collect panels and arrange
	plots <- list(traces = p_traces, A_hist = p_A_hist)
	if (!is.null(p_rho)) plots$rho <- p_rho
	if (!is.null(p_g2)) plots$g2 <- p_g2

	if (requireNamespace("gridExtra", quietly = TRUE)) {
		do.call(gridExtra::grid.arrange, c(plots, ncol = 1))
	} else {
		plots
	}
}
####

####
#' Summary for static DBN fits
#'
#' Append a plain-language label to a scalar-parameter name
#'
#' @description Maps cryptic parameter names (`s2`, `tau_A2`, ...) used in
#'   `summary()` output to a `name (description)` display string, so a reader
#'   need not memorise the abbreviations.
#' @param par Character scalar parameter name.
#' @return Character display string.
#' @keywords internal
#' @noRd
.dbn_param_disp <- function(par, model = NULL) {
	# `t2` carries no innovations in the static model (A, B are time-invariant)
	# -- relabel it accordingly. Same for tau_A2 / tau_B2 in non-time-varying
	# variants. The dynamic / hmm / lowrank / piecewise variants do have
	# time-varying operators and the "innovation variance" label is correct.
	is_static <- identical(model, "static")
	lbl <- switch(par,
		s2 = , sigma2 = , sigma2_proc = "process variance",
		t2 = , tau2 = if (is_static) "A/B prior variance" else "operator innovation variance",
		tau_A2 = if (is_static) "A prior variance" else "A innovation variance",
		tau_B2 = if (is_static) "B prior variance" else "B innovation variance",
		g2 = "baseline-mean variance",
		sigma2_obs = "observation-noise variance",
		rhoA = , rho_A = "A persistence",
		rhoB = , rho_B = "B persistence",
		NULL
	)
	if (is.null(lbl)) par else paste0(par, " (", lbl, ")")
}

#' Print summary for a static DBN fit
#'
#' Internal helper invoked by `summary.dbn()` when `fit$model == "static"`.
#'
#' @param object Object of class "dbn" with model="static"
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisible object
#' @keywords internal
#' @noRd
summary_static <- function(object, digits = 3, ...) {
	if (object$model != "static") cli::cli_abort("This function requires a static model, but got {.val {object$model}}.")

	cli::cli_h1("Static Bilinear Network Model")
	cli::cli_h3("Dimensions")

	dims <- object$meta$dims %||% object$dims
	cli::cli_bullets(c(
		" " = if (isTRUE(dims$is_bipartite)) "Senders: {dims$n_row}, Receivers: {dims$n_col}" else "Nodes: {dims$n_row}",
		" " = "Relations: {dims$p}",
		" " = "Time points: {dims$Tt}"
	))

	cli::cli_h3("Parameter estimates (mean [95% CI])")

	if (exists("param_summary")) {
		param_df <- param_summary(object, probs = c(0.025, 0.975))
		if (!is.null(param_df)) {
			for (i in seq_len(nrow(param_df))) {
				par <- .dbn_param_disp(param_df$parameter[i], model = "static")
				if ("q2.5" %in% names(param_df)) {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q2.5[i], digits)}, {round(param_df$q97.5[i], digits)}]")
				} else if ("q5" %in% names(param_df)) {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q5[i], digits)}, {round(param_df$q95[i], digits)}]")
				} else if ("q50" %in% names(param_df)) {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q5[i], digits)}, {round(param_df$q95[i], digits)}]")
				} else {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)}")
				}
			}
		}
	} else {
		params <- object$params
		for (par in colnames(params)) {
			vals <- params[, par]
			cli::cli_inform("  {par}: {round(mean(vals, na.rm = TRUE), digits)} [{round(quantile(vals, 0.025, na.rm = TRUE), digits)}, {round(quantile(vals, 0.975, na.rm = TRUE), digits)}]")
		}
	}

	cli::cli_h3("Settings")
	# prefer the inner `settings` block (flat scalars) over the meta wrapper,
	# which contains nested dims/settings lists. plain {x} interpolation on
	# a list collapses to a comma-string.
	settings <- if (!is.null(object$settings)) object$settings
				else if (!is.null(object$meta$settings)) object$meta$settings
				else object$meta
	for (s in names(settings)) {
		v <- settings[[s]]
		if (is.null(v) || s == "dims") next
		if (is.list(v) || length(v) > 1L) {
			cli::cli_inform("  {s}:")
			for (sn in names(v)) {
				if (length(v[[sn]]) == 1L) cli::cli_inform("    {sn}: {v[[sn]]}")
			}
		} else {
			cli::cli_inform("  {s}: {v}")
		}
	}

	invisible(object)
}
####

####
#' Summary for dynamic DBN fits
#'
#' @description Prints summary statistics for dynamic DBN model results
#' @param object Object of class "dbn" with model="dynamic"
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisible object
#' @keywords internal
summary_dynamic <- function(object, digits = 3, ...) {
	if (object$model != "dynamic") cli::cli_abort("This function requires a dynamic model, but got {.val {object$model}}.")

	cli::cli_h1("Dynamic Bilinear Network Model")
	cli::cli_h3("Dimensions")

	dims <- object$meta$dims %||% object$dims
	cli::cli_bullets(c(
		" " = if (isTRUE(dims$is_bipartite)) "Senders: {dims$n_row}, Receivers: {dims$n_col}" else "Nodes: {dims$n_row}",
		" " = "Relations: {dims$p}",
		" " = "Time points: {dims$Tt}"
	))

	# A point-estimate fit (ALS) carries a single pseudo-draw of the scalar
	# parameters, so a "95% CI" on them would collapse to [point, point].
	# bootstrap covers A/B/M but does not produce scalar-parameter draws,
	# so it also reports point estimates here. The bootstrap is advertised
	# separately below.
	su <- object$meta$sampler_used; is_als <- length(su) == 1L && !is.na(su) && su %in% c("als", "als_tv")
	has_unc <- !isFALSE(object$meta$uncertainty_available) && !is_als
	cli::cli_h3(if (has_unc) "Parameter estimates (mean [95% CI])" else "Parameter estimates (point estimate, no credible intervals)")

	if (exists("param_summary")) {
		param_df <- param_summary(object, probs = c(0.025, 0.975))
		if (!is.null(param_df)) {
			for (i in seq_len(nrow(param_df))) {
				par <- .dbn_param_disp(param_df$parameter[i])
				if (!has_unc) {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)}")
				} else if ("q2.5" %in% names(param_df)) {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q2.5[i], digits)}, {round(param_df$q97.5[i], digits)}]")
				} else if ("q5" %in% names(param_df)) {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q5[i], digits)}, {round(param_df$q95[i], digits)}]")
				} else if ("q50" %in% names(param_df)) {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q5[i], digits)}, {round(param_df$q95[i], digits)}]")
				} else {
					cli::cli_inform("  {par}: {round(param_df$mean[i], digits)}")
				}
			}
		}
	} else {
		for (par in c("sigma2", "tau_A2", "tau_B2", "g2")) {
			if (!is.null(object[[par]])) {
				vals <- object[[par]]
				if (!has_unc) {
					cli::cli_inform("  {par}: {round(mean(vals, na.rm = TRUE), digits)}")
				} else {
					cli::cli_inform("  {par}: {round(mean(vals, na.rm = TRUE), digits)} [{round(quantile(vals, 0.025, na.rm = TRUE), digits)}, {round(quantile(vals, 0.975, na.rm = TRUE), digits)}]")
				}
			}
		}

		if (!is.null(object$rho_A)) {
			cli::cli_inform("  rho_A: {round(mean(object$rho_A, na.rm = TRUE), digits)} [{round(quantile(object$rho_A, 0.025, na.rm = TRUE), digits)}, {round(quantile(object$rho_A, 0.975, na.rm = TRUE), digits)}]")
			cli::cli_inform("  rho_B: {round(mean(object$rho_B, na.rm = TRUE), digits)} [{round(quantile(object$rho_B, 0.025, na.rm = TRUE), digits)}, {round(quantile(object$rho_B, 0.975, na.rm = TRUE), digits)}]")
		}
	}

	# bootstrap advertisement: when an ALS fit was run with `bootstrap = N`,
	# point the user at the dbn_boot object that carries A/B/M CIs.
	if (is_als && !is.null(object$bootstrap) && inherits(object$bootstrap, "dbn_boot")) {
		cli::cli_h3("Bootstrap")
		cli::cli_inform("  {object$bootstrap$n_valid}/{object$bootstrap$n_total} valid {object$bootstrap$type}-bootstrap replicates")
		cli::cli_inform("  CIs for A, B, M available via {.code fit$bootstrap} (try {.fun print} or {.fun summary} on it)")
	}

	cli::cli_h3("Settings")
	# prefer the inner `settings` block (flat scalars) over the meta
	# wrapper, which contains nested dims/settings/tv lists. plain {x}
	# interpolation on a list collapses to a comma-string.
	settings <- if (!is.null(object$settings)) object$settings
				else if (!is.null(object$meta$settings)) object$meta$settings
				else object$meta
	for (s in names(settings)) {
		v <- settings[[s]]
		if (is.null(v) || s == "dims") next
		if (is.list(v) || length(v) > 1L) {
			cli::cli_inform("  {s}:")
			for (sn in names(v)) {
				if (length(v[[sn]]) == 1L) cli::cli_inform("    {sn}: {v[[sn]]}")
			}
		} else {
			cli::cli_inform("  {s}: {v}")
		}
	}

	invisible(object)
}
####

####
#' Check MCMC Convergence
#'
#' @description Prints effective sample sizes (ESS) and Geweke diagnostics
#'   for the scalar variance parameters in a fitted DBN model. Always run
#'   this before interpreting results to verify the MCMC sampler has
#'   converged.
#'
#' @details
#' **Effective sample size (ESS):** The number of effectively independent
#' posterior draws. Values above 200-400 are generally adequate. Low ESS
#' means the chain is highly autocorrelated and you should increase
#' `nscan` or `odens`.
#'
#' **Geweke diagnostic:** Tests whether the first and last portions of
#' the chain come from the same distribution. Absolute values above 2
#' suggest the chain has not converged and you should increase `burn`.
#'
#' For visual diagnostics, use [plot_trace()] to inspect trace plots.
#'
#' \strong{Note on Rhat:} `dbn()` runs a single-chain conjugate Gibbs
#' sampler, so there is no between-chain variance to compute Rhat against.
#' This routine instead reports per-parameter Geweke z-scores (comparing
#' the first 10% to the last 50% of the chain) and effective sample size
#' from the chain's autocorrelation function. For multi-chain inference,
#' run `dbn()` several times with different `seed` values and compare the
#' resulting posterior summaries manually.
#'
#' @param results Output from [dbn()]
#' @return Invisibly, a data frame with one row per sampled parameter and
#'   columns `parameter`, `ess`, `geweke_z`, `ess_ok` (ESS >= 200), and
#'   `geweke_ok` (|Geweke z| <= 2). Diagnostics are also printed to the
#'   console, with threshold violations flagged. Returns `NULL` invisibly when
#'   there are no sampled parameters to diagnose.
#' @seealso \code{\link{dbn}}, \code{\link{compare_dbn}},
#'   \code{\link{param_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' check_convergence(fit)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
check_convergence <- function(results) {
	# dispatch onto dbn_multichain by computing cross-chain Rhat/ESS on
	# the scalar variances; fall through to single-chain logic otherwise.
	if (inherits(results, "dbn_multichain")) {
		return(.check_convergence_multichain(results))
	}
	if (!inherits(results, "dbn")) {
		cli::cli_abort(c(
			"{.arg results} must be a {.cls dbn} object returned by {.fun dbn}.",
			"i" = "For multi-chain diagnostics, pass a {.cls dbn_multichain} (from {.fun dbn_multichain}).",
			"x" = "Got {.obj_type_friendly {results}}."
		))
	}
	if (isTRUE(results$meta$sampler_used %in% c("als", "als_tv")) ||
		isFALSE(results$meta$uncertainty_available)) {
		cli::cli_abort(c(
			"Convergence diagnostics require a posterior MCMC chain.",
			"x" = "This is a point-estimate fit from {.fun dbn_als} -- there is no chain to diagnose.",
			"i" = "Refit with {.fun dbn} for full MCMC inference."
		))
	}
	if (!requireNamespace("coda", quietly = TRUE)) {
		cli::cli_inform(c(
			"!" = "Package 'coda' is required for full convergence diagnostics.",
			"i" = "Install it with: {.code install.packages('coda')}",
			" " = "Showing basic diagnostics instead..."
		))

		# basic diagnostics when coda is not installed
		if (results$model %in% c("static", "piecewise")) {
			params <- results$params
		} else {
			params <- cbind(
				sigma2 = results$sigma2 %||% results$sigma2_proc,
				tau_A2 = results$tau_A2 %||% results$tauA2,
				tau_B2 = results$tau_B2 %||% results$tauB2
			)
			if (!is.null(results$g2)) {
				params <- cbind(params, g2 = results$g2)
			}
		}

		cli::cli_h3("Parameter Summary")
		print(summary(params))

		cli::cli_h3("Autocorrelation at lag 1")
		for (j in 1:ncol(params)) {
			ac1 <- cor(params[-nrow(params), j], params[-1, j])
			cli::cli_inform("  {colnames(params)[j]}: {round(ac1, 3)}")
		}

		return(invisible(NULL))
	}

	# coda-based diagnostics (ESS, Geweke, autocorrelation)
	if (results$model %in% c("static", "piecewise")) {
		params_mcmc <- coda::mcmc(results$params)
	} else {
		params_df <- cbind(
			sigma2 = results$sigma2,
			tau_A2 = results$tau_A2 %||% results$tauA2,
			tau_B2 = results$tau_B2 %||% results$tauB2
		)
		if (!is.null(results$g2)) {
			params_df <- cbind(params_df, g2 = results$g2)
		}
		# low-rank fits carry the factor-strength innovation variance
		# tau_alpha2 -- the defining parameter of the model. include it so it
		# is not silently omitted from the diagnostic.
		if (!is.null(results$tau_alpha2) && length(results$tau_alpha2) == nrow(params_df)) {
			params_df <- cbind(params_df, tau_alpha2 = results$tau_alpha2)
		}
		if (!is.null(results$sigma2_obs)) {
			params_df <- cbind(params_df, sigma2_obs = results$sigma2_obs)
		}
		# include AR(1) persistence parameters when present (`ar1 = TRUE`).
		# without these, a user with poor rho mixing would see "all parameters
		# converged" and never know.
		rhoA_draws <- results$rhoA %||% results$rho_A
		rhoB_draws <- results$rhoB %||% results$rho_B
		if (!is.null(rhoA_draws) && length(rhoA_draws) == nrow(params_df)) {
			params_df <- cbind(params_df, rhoA = rhoA_draws)
		}
		if (!is.null(rhoB_draws) && length(rhoB_draws) == nrow(params_df)) {
			params_df <- cbind(params_df, rhoB = rhoB_draws)
		}
		params_mcmc <- coda::mcmc(params_df)
	}

	# drop columns with zero variance (fixed parameters, e.g. ordinal s2)
	col_var <- apply(as.matrix(params_mcmc), 2, var, na.rm = TRUE)
	fixed_cols <- names(col_var)[!is.na(col_var) & col_var < .Machine$double.eps]
	varying_mcmc <- params_mcmc
	if (length(fixed_cols) > 0) {
		cli::cli_alert_info("Fixed parameters (not sampled): {.val {fixed_cols}}")
		keep <- setdiff(colnames(params_mcmc), fixed_cols)
		if (length(keep) > 0) {
			varying_mcmc <- coda::mcmc(as.matrix(params_mcmc)[, keep, drop = FALSE])
		} else {
			cli::cli_alert_warning("All parameters are fixed; no convergence diagnostics to compute.")
			return(invisible(NULL))
		}
	}

	ess <- coda::effectiveSize(varying_mcmc)
	gew <- coda::geweke.diag(varying_mcmc)$z
	gew <- gew[names(ess)]

	cli::cli_h3("Effective Sample Sizes")
	print(round(ess, 1))

	cli::cli_h3("Geweke Diagnostic (z-scores)")
	print(round(gew, 3))

	# flag parameters that fail the documented thresholds (ESS >= 200,
	# |Geweke z| <= 2) instead of leaving the reader to check by eye
	low_ess <- names(ess)[ess < 200]
	bad_gew <- names(gew)[is.finite(gew) & abs(gew) > 2]
	if (length(low_ess) > 0) {
		cli::cli_alert_warning("Low effective sample size (< 200): {.val {low_ess}}. Increase {.arg nscan} (raising {.arg odens} alone reduces the stored draws and will not help).")
	}
	if (length(bad_gew) > 0) {
		cli::cli_alert_warning("Geweke |z| > 2: {.val {bad_gew}}. The chain may not have converged.")
	}
	if (length(low_ess) == 0 && length(bad_gew) == 0) {
		cli::cli_alert_success("All sampled parameters pass ESS >= 200 and |Geweke z| <= 2.")
	}

	n_pars <- ncol(as.matrix(varying_mcmc))
	if (n_pars > 0) {
		par(mfrow = c(min(n_pars, 2), min(n_pars, 2)))
		coda::autocorr.plot(varying_mcmc, auto.layout = FALSE)
	}

	invisible(data.frame(
		parameter = names(ess),
		ess       = as.numeric(ess),
		geweke_z  = as.numeric(gew),
		ess_ok    = as.numeric(ess) >= 200,
		geweke_ok = !(is.finite(gew) & abs(gew) > 2),
		row.names = NULL, stringsAsFactors = FALSE
	))
}
####

####
#' Compare fitted DBN models with information criteria
#'
#' Compare two or more fitted \pkg{dbn} models using WAIC, AIC, and BIC.
#' The returned object is a data frame ordered by the criterion selected in
#' \code{type}. Smaller WAIC/AIC/BIC values are preferred. Delta columns are
#' criterion differences from the best model under that criterion, and weight
#' columns are normalized \code{exp(-0.5 * delta)} weights.
#'
#' @param ... Two or more fitted objects of class \code{"dbn"}, or a single
#'   list of fitted \code{"dbn"} objects.
#' @param type Character string selecting the criterion used to order rows:
#'   \code{"waic"}, \code{"aic"}, or \code{"bic"}. All available criteria are
#'   still returned as columns.
#' @param names Optional character vector of display names. Must have the same
#'   length as the number of fitted models.
#'
#' @return A data frame with class \code{c("dbn_compare", "data.frame")} and
#'   columns \code{model}, \code{elpd_waic}, \code{se_waic}, \code{p_waic},
#'   \code{waic}, \code{aic}, \code{bic}, \code{delta_waic},
#'   \code{delta_aic}, \code{delta_bic}, \code{weight_waic},
#'   \code{weight_aic}, and \code{weight_bic}.
#'
#' @seealso \code{\link{dbn}}, \code{\link{check_convergence}},
#'   \code{\link{compute_waic_dbn}}
#'
#' @examples
#' \dontrun{
#' fit_dyn <- dbn(Y, model = "dynamic", family = "gaussian")
#' fit_static <- dbn(Y, model = "static", family = "gaussian")
#' compare_dbn(dynamic = fit_dyn, static = fit_static)
#' compare_dbn(list(fit_dyn, fit_static), names = c("dynamic", "static"),
#'     type = "aic")
#' }
#'
#' @author Tosin Salau and Shahryar Minhas
#' @export
compare_dbn <- function(..., type = c("waic", "aic", "bic"), names = NULL) {
	type <- match.arg(type)

	args <- list(...)
	if (length(args) == 0L) {
		cli::cli_abort(c(
			"At least two fitted {.cls dbn} objects are required.",
			"i" = "Use {.code compare_dbn(fit_a, fit_b)} or {.code compare_dbn(list(fit_a, fit_b))}."
		))
	}

	dot_names <- base::names(args)
	if (length(args) == 1L && is.list(args[[1L]]) && !inherits(args[[1L]], "dbn")) {
		fits <- args[[1L]]
		if (is.null(names)) {
			dot_names <- base::names(fits)
		}
	} else {
		fits <- args
	}

	if (length(fits) < 2L) {
		cli::cli_abort(c(
			"At least two fitted {.cls dbn} objects are required.",
			"i" = "A single fit cannot be compared to itself."
		))
	}

	bad <- which(!vapply(fits, .dbn_compare_is_dbn, logical(1L)))
	if (length(bad) > 0L) {
		cli::cli_abort(c(
			"Every object supplied to {.fn compare_dbn} must be a fitted {.cls dbn} object.",
			"x" = paste0("Non-dbn object positions: ", paste(bad, collapse = ", "), ".")
		))
	}

	model_names <- .dbn_compare_names(fits, dot_names = dot_names, names = names)

	results <- lapply(seq_along(fits), function(i) {
		.dbn_compare_one_fit(fits[[i]], model_name = model_names[[i]])
	})
	out <- do.call(rbind, results)

	out$delta_waic <- .dbn_compare_delta(out$waic)
	out$delta_aic <- .dbn_compare_delta(out$aic)
	out$delta_bic <- .dbn_compare_delta(out$bic)

	out$weight_waic <- .dbn_compare_weights(out$delta_waic)
	out$weight_aic <- .dbn_compare_weights(out$delta_aic)
	out$weight_bic <- .dbn_compare_weights(out$delta_bic)

	sort_col <- switch(type, waic = "waic", aic = "aic", bic = "bic")

	if (!any(is.finite(out[[sort_col]]))) {
		cli::cli_abort(c(
			"Cannot rank models by {.val {type}} because no finite values were available.",
			"i" = "Check that each fit supports the selected information criterion.",
			"i" = "For ALS or other point-estimate paths without WAIC, try {.code type = \"aic\"} or {.code type = \"bic\"}."
		))
	}

	out <- out[order(out[[sort_col]], na.last = TRUE), , drop = FALSE]
	rownames(out) <- NULL

	class(out) <- c("dbn_compare", "data.frame")
	attr(out, "criterion") <- type
	attr(out, "n_models") <- nrow(out)

	out
}

#' @keywords internal
#' @noRd
.dbn_compare_is_dbn <- function(x) {
	inherits(x, "dbn")
}

#' @keywords internal
#' @noRd
.dbn_compare_names <- function(fits, dot_names = NULL, names = NULL) {
	if (!is.null(names)) {
		if (!is.character(names)) {
			cli::cli_abort(c(
				"{.arg names} must be a character vector.",
				"x" = paste0("Got {.cls ", class(names)[[1L]], "}.")
			))
		}
		if (length(names) != length(fits)) {
			cli::cli_abort(c(
				"{.arg names} must have one entry per fitted model.",
				"x" = paste0("Got ", length(names), " names for ", length(fits), " fits.")
			))
		}
		if (anyNA(names) || any(!nzchar(names))) {
			cli::cli_abort(c(
				"{.arg names} must not contain missing or empty values.",
				"i" = "Use explicit labels such as {.code names = c(\"dynamic\", \"static\")}."
			))
		}
		return(make.unique(names, sep = "_"))
	}

	if (is.null(dot_names)) {
		dot_names <- rep("", length(fits))
	}
	if (length(dot_names) != length(fits)) {
		dot_names <- rep("", length(fits))
	}

	out <- dot_names
	empty <- is.na(out) | !nzchar(out)
	if (any(empty)) {
		out[empty] <- vapply(which(empty), function(i) {
			.dbn_compare_default_name(fits[[i]], i)
		}, character(1L))
	}

	make.unique(out, sep = "_")
}

#' @keywords internal
#' @noRd
.dbn_compare_default_name <- function(fit, i) {
	label <- NULL
	if (!is.null(fit$model)) label <- fit$model
	if (is.null(label) && !is.null(fit$meta) && !is.null(fit$meta$model))
		label <- fit$meta$model
	if (is.null(label) && !is.null(fit$meta) && !is.null(fit$meta$sampler_used))
		label <- fit$meta$sampler_used

	if (is.null(label) || length(label) == 0L || is.na(label[[1L]]) ||
			!nzchar(as.character(label[[1L]]))) {
		return(paste0("dbn_", i))
	}
	paste0(as.character(label[[1L]]), "_", i)
}

#' @keywords internal
#' @noRd
.dbn_compare_one_fit <- function(fit, model_name) {
	waic_vals <- .dbn_compare_waic_values(fit, model_name = model_name)
	data.frame(
		model = model_name,
		elpd_waic = waic_vals[["elpd_waic"]],
		se_waic = waic_vals[["se_waic"]],
		p_waic = waic_vals[["p_waic"]],
		waic = waic_vals[["waic"]],
		aic = .dbn_compare_ic_value(fit, criterion = "AIC", model_name = model_name),
		bic = .dbn_compare_ic_value(fit, criterion = "BIC", model_name = model_name),
		stringsAsFactors = FALSE
	)
}

#' @keywords internal
#' @noRd
.dbn_compare_waic_values <- function(fit, model_name) {
	na_vals <- c(elpd_waic = NA_real_, se_waic = NA_real_,
				 p_waic = NA_real_, waic = NA_real_)

	res <- tryCatch(compute_waic_dbn(fit), error = function(e) e)
	if (inherits(res, "error")) {
		cli::cli_warn(c(
			"Could not compute WAIC for one DBN fit.",
			"x" = paste0("Model ", shQuote(model_name), ": ", conditionMessage(res)),
			"i" = "WAIC columns are set to NA for this model."
		))
		return(na_vals)
	}

	vals <- tryCatch(.dbn_compare_extract_waic(res), error = function(e) e)
	if (inherits(vals, "error")) {
		cli::cli_warn(c(
			"Could not parse the WAIC object returned by {.fn compute_waic_dbn}.",
			"x" = paste0("Model ", shQuote(model_name), ": ", conditionMessage(vals)),
			"i" = "WAIC columns are set to NA for this model."
		))
		return(na_vals)
	}
	vals
}

#' @keywords internal
#' @noRd
.dbn_compare_extract_waic <- function(waic_obj) {
	if (is.list(waic_obj) && !is.null(waic_obj$estimates)) {
		est <- waic_obj$estimates
	} else {
		est <- waic_obj
	}
	if (is.data.frame(est)) est <- as.matrix(est)

	if (!is.matrix(est)) {
		cli::cli_abort(c(
			"WAIC result does not contain an estimates matrix.",
			"i" = "Expected a {.pkg loo}-style object with {.field estimates}."
		))
	}

	required_rows <- c("elpd_waic", "p_waic", "waic")
	missing_rows <- setdiff(required_rows, rownames(est))
	if (length(missing_rows) > 0L) {
		cli::cli_abort(c(
			"WAIC estimates matrix is missing required rows.",
			"x" = paste0("Missing: ", paste(missing_rows, collapse = ", "), ".")
		))
	}

	estimate_col <- if ("Estimate" %in% colnames(est)) "Estimate" else 1L
	se_col <- if ("SE" %in% colnames(est)) "SE" else if ("se" %in% colnames(est)) "se" else NA

	out <- c(
		elpd_waic = as.numeric(est["elpd_waic", estimate_col]),
		se_waic = NA_real_,
		p_waic = as.numeric(est["p_waic", estimate_col]),
		waic = as.numeric(est["waic", estimate_col])
	)
	if (!is.na(se_col)) out[["se_waic"]] <- as.numeric(est["waic", se_col])
	out
}

#' @keywords internal
#' @noRd
.dbn_compare_ic_value <- function(fit, criterion = c("AIC", "BIC"), model_name) {
	criterion <- match.arg(criterion)
	fun <- switch(criterion, AIC = stats::AIC, BIC = stats::BIC)

	res <- tryCatch(fun(fit), error = function(e) e)
	if (inherits(res, "error")) {
		cli::cli_warn(c(
			paste0("Could not compute ", criterion, " for one DBN fit."),
			"x" = paste0("Model ", shQuote(model_name), ": ", conditionMessage(res)),
			"i" = paste0(criterion, " is set to NA for this model.")
		))
		return(NA_real_)
	}

	if (is.data.frame(res)) {
		if (criterion %in% colnames(res)) return(as.numeric(res[[criterion]][[1L]]))
		num_cols <- vapply(res, is.numeric, logical(1L))
		if (any(num_cols)) return(as.numeric(res[[which(num_cols)[[1L]]]][[1L]]))
		return(NA_real_)
	}

	val <- suppressWarnings(as.numeric(res))
	if (length(val) == 0L) return(NA_real_)
	val[[1L]]
}

#' @keywords internal
#' @noRd
.dbn_compare_delta <- function(x) {
	delta <- rep(NA_real_, length(x))
	ok <- is.finite(x)
	if (!any(ok)) return(delta)
	delta[ok] <- x[ok] - min(x[ok])
	delta
}

#' @keywords internal
#' @noRd
.dbn_compare_weights <- function(delta) {
	weights <- rep(NA_real_, length(delta))
	ok <- is.finite(delta)
	if (!any(ok)) return(weights)
	raw <- exp(-0.5 * delta[ok])
	denom <- sum(raw)
	if (!is.finite(denom) || denom <= 0) return(weights)
	weights[ok] <- raw / denom
	weights
}

#' Print a DBN model-comparison table
#'
#' @param x An object returned by \code{\link{compare_dbn}}.
#' @param ... Unused.
#' @param digits Number of significant digits to display.
#'
#' @return Invisibly returns \code{x}.
#'
#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_compare <- function(x, ..., digits = 3) {
	criterion <- attr(x, "criterion")
	if (is.null(criterion) || !nzchar(criterion)) criterion <- "waic"

	cli::cli_h3("DBN model comparison")
	cli::cli_inform(paste0(
		"Ranked by ", toupper(criterion),
		"; smaller information-criterion values are preferred."
	))

	cols <- c("model", "elpd_waic", "se_waic", "p_waic", "waic", "aic", "bic",
			  "delta_waic", "delta_aic", "delta_bic",
			  "weight_waic", "weight_aic", "weight_bic")
	cols <- intersect(cols, colnames(x))
	display <- x[, cols, drop = FALSE]

	num_cols <- vapply(display, is.numeric, logical(1L))
	display[num_cols] <- lapply(display[num_cols], .dbn_compare_format_number,
		digits = digits)

	# strip the dbn_compare class on `display` so this print() dispatches to
	# base R's print.data.frame instead of recursing into print.dbn_compare.
	class(display) <- "data.frame"
	print(display, row.names = FALSE, right = TRUE)

	if (anyNA(x$waic)) {
		cli::cli_inform(c(
			"i" = "At least one WAIC value is NA. This usually means the fit does not expose a pointwise posterior log-likelihood."
		))
	}

	invisible(x)
}

#' @keywords internal
#' @noRd
.dbn_compare_format_number <- function(x, digits = 3) {
	out <- rep(NA_character_, length(x))
	ok <- !is.na(x)
	if (any(ok)) {
		# use scientific notation when the magnitude is < 1e-3 or >= 1e6;
		# otherwise fixed. Avoids 2.7e-21 printing as 20 leading zeros.
		v <- signif(x[ok], digits = digits)
		mag <- abs(v)
		use_sci <- mag > 0 & (mag < 1e-3 | mag >= 1e6)
		fmt <- character(length(v))
		if (any(use_sci)) {
			fmt[use_sci] <- format(v[use_sci], trim = TRUE, scientific = TRUE, digits = digits)
		}
		if (any(!use_sci)) {
			fmt[!use_sci] <- format(v[!use_sci], trim = TRUE, scientific = FALSE)
		}
		out[ok] <- fmt
	}
	out
}
####

####
#' Posterior Predictive ECDF Overlay
#'
#' @description Creates posterior predictive checks using ECDF overlays
#' @param fit Object of class "dbn"
#' @param n_rep Number of replicated datasets to draw (default: 20)
#' @return A ggplot2 object
#' @keywords internal
ppc_ecdf <- function(fit, n_rep = 20) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort(c("Package {.pkg ggplot2} is required for this function.", "i" = "Install with {.code install.packages(\"ggplot2\")}"))
	if (!inherits(fit, "dbn")) {
		cli::cli_abort("{.arg fit} must be of class {.cls dbn}")
	}
	if (!fit$model %in% c("static", "dynamic")) {
		cli::cli_abort("Model must be either 'static' or 'dynamic'")
	}

	if (!exists("tprod")) {
		cli::cli_abort("tprod function not found. Please load the dbn package.")
	}

	y_orig <- fit$Y %||% fit$R
	ecdf_orig <- ecdf(c(y_orig))
	stat_orig <- ecdf_orig(sort(unique(c(y_orig))))

	# draw one replicated dataset from the posterior
	draw_rep <- function() {
		if (fit$model == "static") {
			idx <- sample(seq_len(dim(fit$B[[1]])[3]), 1)
			Bsim <- lapply(fit$B, function(x) x[, , idx])
			Yhat <- tprod(fit$M, Bsim) + array(rnorm(length(fit$M)), dim = dim(fit$M))
		} else {
			idx <- sample(seq_along(fit$A), 1)
			A <- fit$A[[idx]]
			B <- fit$B[[idx]]
			nr <- fit$dims$n_row
			ncl <- fit$dims$n_col
			Theta <- array(0, dim = c(nr, ncl, fit$dims$p, dim(A)[3]))
			Theta[, , , 1] <- 0
			# safe sigma2 index: bootstrap-expanded fits have length-1 sigma2
			sig2_draw <- if (length(fit$sigma2) >= idx) fit$sigma2[idx] else fit$sigma2[1]
			if (!is.finite(sig2_draw)) sig2_draw <- 1
			for (t in 2:dim(Theta)[4]) {
				for (rel in 1:fit$dims$p) {
					Theta[, , rel, t] <- A[, , t] %*% Theta[, , rel, t - 1] %*% t(B[, , t]) +
						sqrt(sig2_draw) * matrix(rnorm(nr * ncl), nr, ncl)
				}
			}
			M_draw <- if (length(dim(fit$M)) == 4) fit$M[, , , idx] else fit$M[[idx]]
			Yhat <- sweep(Theta, 1:3, M_draw, "+")
		}
		# round and clip to observed ordinal range
		vals <- sort(unique(c(y_orig)))
		Yrep <- pmax(min(vals), pmin(max(vals), round(Yhat)))
		ecdf(c(Yrep))(vals)
	}

	reps <- replicate(n_rep, draw_rep())

	df <- data.frame(
		quant = rep(sort(unique(c(y_orig))), n_rep + 1),
		ecdf = c(stat_orig, as.vector(reps)),
		sel = rep(c("Observed", paste0("Rep", 1:n_rep)), each = length(stat_orig))
	)

	ggplot2::ggplot(df, ggplot2::aes(
		x = quant, y = ecdf, group = sel,
		colour = sel == "Observed"
	)) +
		ggplot2::geom_line(alpha = 0.5) +
		ggplot2::scale_color_manual(values = c("grey70", "black"), guide = "none") +
		ggplot2::labs(
			title = "Posterior predictive ECDF overlay",
			x = "Ordinal category", y = "ECDF"
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
}
####

####
#' Plot Dyad Trajectory
#'
#' @description Plots posterior mean and 95% bands for a single dyad through time
#' @param fit Dynamic dbn object
#' @param i Actor i index
#' @param j Actor j index
#' @param rel Relation indices (default: NULL = all relations)
#' @param facet Whether to facet by relation (default: TRUE)
#' @param cred Credible interval quantiles (default: c(0.025, 0.975))
#' @return A ggplot2 object
#' @seealso \code{\link{dbn}}, \code{\link{net_snapshot}},
#'   \code{\link{role_trajectory}}, \code{\link{theta_summary}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' dyad_path(fit, i = 1, j = 2)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
dyad_path <- function(fit, i, j, rel = NULL, facet = TRUE, cred = c(0.025, 0.975)) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort(c("Package {.pkg ggplot2} is required for this function.", "i" = "Install with {.code install.packages(\"ggplot2\")}"))
	if (!fit$model %in% c("dynamic", "lowrank", "hmm")) {
		cli::cli_abort("This function requires a time-varying model (dynamic, lowrank, or hmm)")
	}
	p <- fit$dims$p

	if (is.null(rel)) {
		rel <- seq_len(p)
	}
	# validate `rel` as integer indices into the relation axis; an
	# out-of-range value would otherwise plot relation 1's draws under the
	# wrong facet label.
	if (!all(is.numeric(rel)) || any(rel < 1) || any(rel > p) || any(rel != as.integer(rel))) {
		cli::cli_abort(c(
			"{.arg rel} must be (an integer or integer vector of) relation indices in {.code 1:p}.",
			"x" = "Got {.val {rel}}; fit has {.code p = {p}} relation{?s}."
		))
	}
	rel <- as.integer(rel)
	# read the latent-state draws directly from fit$Theta, which is what
	# the function name and the y-axis label promise (the lag operator
	# (A_t B_t^T)[i,j] would be a different quantity).
	#
	# different model variants store theta differently:
	#  * dynamic / piecewise: fit$Theta (5D top-level cube)
	#  * hmm / lowrank: fit$draws$theta (list of [n_row, n_col, p, T] arrays)
	# fall back to the list-of-arrays form so dyad_path works on hmm fits.
	Theta <- fit$Theta
	if (is.null(Theta) && !is.null(fit$draws$theta) && is.list(fit$draws$theta) &&
		length(fit$draws$theta) > 0L && length(dim(fit$draws$theta[[1]])) == 4L) {
		# stack list -> 5D array [n_row, n_col, p, T, draws]
		d4 <- dim(fit$draws$theta[[1]])
		Theta <- array(unlist(fit$draws$theta),
			dim = c(d4, length(fit$draws$theta)))
	}
	if (is.null(Theta) || length(dim(Theta)) != 5L) {
		cli::cli_abort(c(
			"{.fun dyad_path} requires the latent-state draws.",
			"x" = "Neither {.code fit$Theta} (5D array) nor {.code fit$draws$theta} (list of 4D arrays) is available.",
			"i" = "The fit may have been produced with {.code keep = ...} that dropped Theta. Refit without that exclusion."
		))
	}
	dim_Theta <- dim(Theta)  # [n_row, n_col, p, T, draws]
	n_row <- dim_Theta[1]; n_col <- dim_Theta[2]
	# accept character actor names against dimnames(fit$Y), so substantive
	# users can write dyad_path(fit, "USA", "China") instead of guessing
	# integer indices.
	actor_names_row <- dimnames(fit$Y)[[1L]]
	actor_names_col <- dimnames(fit$Y)[[2L]]
	i_label <- if (is.character(i)) i else NULL
	j_label <- if (is.character(j)) j else NULL
	if (is.character(i)) {
		if (is.null(actor_names_row))
			cli::cli_abort("{.arg i} is a character name but {.code dimnames(fit$Y)[[1]]} is NULL.")
		idx <- match(i, actor_names_row)
		if (length(idx) != 1L || is.na(idx))
			cli::cli_abort(c("{.arg i} actor name not found.", "x" = "Got {.val {i}}; available: {.val {head(actor_names_row, 6)}}"))
		i <- idx
	}
	if (is.character(j)) {
		if (is.null(actor_names_col))
			cli::cli_abort("{.arg j} is a character name but {.code dimnames(fit$Y)[[2]]} is NULL.")
		idx <- match(j, actor_names_col)
		if (length(idx) != 1L || is.na(idx))
			cli::cli_abort(c("{.arg j} actor name not found.", "x" = "Got {.val {j}}; available: {.val {head(actor_names_col, 6)}}"))
		j <- idx
	}
	if (!is.numeric(i) || length(i) != 1L || i < 1 || i > n_row || i != as.integer(i))
		cli::cli_abort("{.arg i} must be a single integer in {.code 1:{n_row}} (or a character name).")
	if (!is.numeric(j) || length(j) != 1L || j < 1 || j > n_col || j != as.integer(j))
		cli::cli_abort("{.arg j} must be a single integer in {.code 1:{n_col}} (or a character name).")
	i <- as.integer(i); j <- as.integer(j)
	# fall back to dimnames if labels weren't supplied as characters
	if (is.null(i_label) && !is.null(actor_names_row)) i_label <- actor_names_row[i]
	if (is.null(j_label) && !is.null(actor_names_col)) j_label <- actor_names_col[j]

	# posterior trajectory for one relation -- read theta_{i,j,r,t,s} directly
	compute_trajectory <- function(r) {
		thetas <- t(Theta[i, j, r, , ])  # [draws, T]
		Tt <- ncol(thetas)

		# convert stored indices to original time scale
		time_vals <- 1:Tt
		if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
			time_vals <- (time_vals - 1) * fit$settings$time_thin + 1
		}

		data.frame(
			relation = paste0("Relation ", r),
			time = time_vals,
			mean = colMeans(thetas),
			lo = apply(thetas, 2, quantile, cred[1]),
			hi = apply(thetas, 2, quantile, cred[2])
		)
	}
	####

	df_all <- do.call(rbind, lapply(rel, compute_trajectory))

	title_str <- if (!is.null(i_label) && !is.null(j_label)) {
		sprintf("Dyad (%s, %s) trajectory", i_label, j_label)
	} else {
		sprintf("Dyad (%d,%d) trajectory", i, j)
	}
	g <- ggplot2::ggplot(df_all, ggplot2::aes(time, mean)) +
		ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80") +
		ggplot2::geom_line(linewidth = 0.7) +
		ggplot2::labs(
			x = "Time", y = expression(theta[list(i, j, t)]),
			title = title_str
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0)
		)

	if (facet && length(rel) > 1) {
		g <- g + ggplot2::facet_wrap(~relation, scales = "free_y")
	}

	g
}
####

####
#' Per-actor singular-vector trajectory
#'
#' @description At each time point, take the posterior-mean operator
#'   (`A_t` for senders, `B_t` for receivers), run its SVD, and read off
#'   the `comp`-th singular vector. This gives one coordinate per actor
#'   per time, tracing how the leading singular direction of the
#'   operator rotates over the panel.
#'
#'   This is a low-dimensional summary of how the operator's row (or
#'   column) space moves; it is not actor influence or coupling. For
#'   per-actor coupling with credible bands, use
#'   [coupling_trajectory()].
#'
#' @param fit Dynamic dbn object
#' @param mat "A" (senders) or "B" (receivers)
#' @param comp Singular vector index to track (default: 1)
#' @param plot Logical. If `TRUE` (default), also draw the trajectory
#'   plot via [ggplot2::ggplot()] and return the result invisibly;
#'   if `FALSE`, return the underlying data frame without plotting.
#' @return A data frame of class `dbn_role_trajectory` with columns
#'   `time`, `actor`, `score`, with attributes `mat` and `comp`. A
#'   companion `plot()` method renders it as a ggplot.
#' @seealso [coupling_trajectory()] for per-actor coupling with bands;
#'   [actor_embedding()] for a single-shot per-actor coordinate.
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' rt <- role_trajectory(fit, mat = "A", comp = 1, plot = FALSE)
#' plot(rt)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
role_trajectory <- function(fit, mat = c("A", "B"), comp = 1, plot = TRUE) {
	mat <- match.arg(mat)
	n_keep <- length(fit[[mat]])
	Tt <- dim(fit[[mat]][[1]])[3]
	m <- if (mat == "A") fit$dims$n_row else fit$dims$n_col
	scores <- matrix(NA_real_, Tt, m)
	for (t in seq_len(Tt)) {
		Mbar <- Reduce(`+`, lapply(fit[[mat]], function(M) M[, , t])) / n_keep
		sv <- svd(Mbar)
		scores[t, ] <- if (mat == "A") sv$u[, comp] else sv$v[, comp]
	}
	actor_names <- if (mat == "A") {
		dimnames(fit$Y)[[1L]] %||% paste0("actor_", seq_len(m))
	} else {
		dimnames(fit$Y)[[2L]] %||% paste0("actor_", seq_len(m))
	}
	out <- data.frame(
		time = rep(seq_len(Tt), times = m),
		actor = factor(rep(actor_names, each = Tt), levels = actor_names),
		score = as.vector(scores),
		stringsAsFactors = FALSE
	)
	attr(out, "mat") <- mat
	attr(out, "comp") <- comp
	class(out) <- c("dbn_role_trajectory", "data.frame")
	if (isTRUE(plot)) {
		p <- plot(out)
		print(p)
		return(invisible(out))
	}
	out
}
####

####
#' Plot method for `dbn_role_trajectory`
#'
#' @param x A `dbn_role_trajectory` object from [role_trajectory()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @author Tosin Salau and Shahryar Minhas
#' @export
plot.dbn_role_trajectory <- function(x, ...) {
	if (!requireNamespace("ggplot2", quietly = TRUE))
		cli::cli_abort("{.pkg ggplot2} is required for plotting role trajectories.")
	mat <- attr(x, "mat") %||% "A"
	comp <- attr(x, "comp") %||% 1
	ggplot2::ggplot(x, ggplot2::aes(x = .data$time, y = .data$score,
	                                colour = .data$actor, group = .data$actor)) +
		ggplot2::geom_line(linewidth = 0.7) +
		ggplot2::labs(
			title = sprintf("%s_t singular vector %d", mat, comp),
			subtitle = "Per-actor coordinate in the leading singular direction over time",
			x = "Time", y = "Singular-vector coordinate", colour = "Actor"
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank())
}
####

####
#' Network Snapshot
#'
#' @description Heat map of Theta at given time (posterior mean)
#' @param fit Dynamic dbn object
#' @param time Time point
#' @param rel Relation index (default: 1)
#' @param sparse Auto-switch to sparse visualization for large networks
#' @param eps Threshold for sparse plotting
#' @param show_significant Logical, whether to show only significant effects (default: FALSE)
#' @param cred_level Credible level for significance (default corresponds to 95% CI)
#' @return A `ggplot` object.
#' @seealso \code{\link{dbn}}, \code{\link{dyad_path}},
#'   \code{\link{theta_summary}}, \code{\link{role_trajectory}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' net_snapshot(fit, time = 5)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
net_snapshot <- function(fit, time, rel = 1, sparse = NULL, eps = 1e-4,
						show_significant = FALSE, cred_level = 0.025) {
	if (!inherits(fit, "dbn")) {
		cli::cli_abort(c(
			"{.fun net_snapshot} expects a fitted {.cls dbn} object.",
			"x" = "Got {.cls {class(fit)[1]}}.",
			"i" = "For raw-data snapshots, use {.code Y[, , rel, t]} directly."
		))
	}
	if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort(c("Package {.pkg ggplot2} is required for this function.", "i" = "Install with {.code install.packages(\"ggplot2\")}"))
	if (!fit$model %in% c("dynamic", "lowrank", "hmm")) {
		cli::cli_abort("This function requires a time-varying model (dynamic, lowrank, or hmm)")
	}
	Tt_fit <- fit$dims$Tt
	if (length(time) != 1L || !is.numeric(time) || !is.finite(time) || time < 1 ||
	    time != round(time) || time > Tt_fit) {
		cli::cli_abort(c(
			"{.arg time} must be a single whole number between 1 and {Tt_fit}.",
			"x" = "Got {.val {time}}."
		))
	}
	n_row <- fit$dims$n_row
	n_col <- fit$dims$n_col
	Th <- matrix(0, n_row, n_col)

	if (show_significant) {
		Th_all <- array(NA, dim = c(n_row, n_col, length(fit$A)))
	}

	# use sparse point plot for large networks
	if (is.null(sparse)) {
		sparse <- max(n_row, n_col) > 150
	}

	# adjust index for time thinning
	t_idx <- time
	if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
		t_idx <- ceiling(time / fit$settings$time_thin)
	}

	# average Theta across posterior draws, model-specific computation
	if (fit$model == "dynamic") {
		# prefer the stored Theta when present (handles bipartite correctly:
		# `A_t %*% t(B_t)` requires n_row == n_col, but Theta_t has the right
		# n_row x n_col shape regardless of topology). fall back to the
		# `A %*% t(B)` reconstruction when Theta is absent.
		if (!is.null(fit$Theta) && length(dim(fit$Theta)) == 5L) {
			d <- dim(fit$Theta)
			n_keep <- d[5]
			if (t_idx > d[4]) {
				cli::cli_abort("Time {time} not available (time_thin={fit$settings$time_thin})")
			}
			for (s in seq_len(n_keep)) {
				Th_s <- fit$Theta[, , rel, t_idx, s]
				if (show_significant) Th_all[, , s] <- Th_s
				Th <- Th + Th_s
			}
			Th <- Th / n_keep
		} else {
			n_keep <- length(fit$A)
			if (t_idx > dim(fit$A[[1]])[3]) {
				cli::cli_abort("Time {time} not available (time_thin={fit$settings$time_thin})")
			}
			if (n_row != n_col) {
				cli::cli_abort(c(
					"{.fun net_snapshot} on a bipartite fit needs {.code fit$Theta} to be stored.",
					"i" = "Refit with default settings ({.code keep} including Theta) to enable this view."
				))
			}
			for (s in seq_len(n_keep)) {
				Th_s <- fit$A[[s]][, , t_idx] %*% t(fit$B[[s]][, , t_idx])
				if (show_significant) Th_all[, , s] <- Th_s
				Th <- Th + Th_s
			}
			Th <- Th / n_keep
		}
	} else if (fit$model == "lowrank") {
		n_keep <- length(fit$U)
		for (s in seq_len(n_keep)) {
			U_s <- fit$U[[s]]
			alpha_s <- fit$alpha[[s]][, t_idx]
			A_s <- U_s %*% diag(alpha_s, nrow = length(alpha_s)) %*% t(U_s)
			B_s <- fit$B[[s]][, , t_idx]
			Th_s <- A_s %*% t(B_s)
			if (show_significant) {
				Th_all[, , s] <- Th_s
			}
			Th <- Th + Th_s
		}
		Th <- Th / n_keep
	} else if (fit$model == "hmm") {
		n_keep <- length(fit$S)
		for (s in seq_len(n_keep)) {
			regime <- fit$S[[s]][t_idx]
			A_s <- fit$A[[s]][, , regime]
			B_s <- fit$B[[s]][, , regime]
			Th_s <- A_s %*% t(B_s)
			if (show_significant) {
				Th_all[, , s] <- Th_s
			}
			Th <- Th + Th_s
		}
		Th <- Th / n_keep
	}
	####

	# zero out cells whose credible interval contains zero
	if (show_significant) {
		alpha <- 1 - cred_level
		lower_q <- alpha / 2
		upper_q <- 1 - alpha / 2

		is_significant <- matrix(TRUE, n_row, n_col)
		for (i in 1:n_row) {
			for (j in 1:n_col) {
				cell_vals <- Th_all[i, j, ]
				ci <- quantile(cell_vals, probs = c(lower_q, upper_q), na.rm = TRUE)
				if (ci[1] <= 0 && ci[2] >= 0) {
					is_significant[i, j] <- FALSE
					Th[i, j] <- NA
				}
			}
		}
	}
	####

	# render the heatmap or sparse point plot
	if (sparse) {
		# threshold small values, keep NAs from significance masking
		if (!show_significant) {
			Th[abs(Th) < eps] <- 0
		} else {
			non_na_mask <- !is.na(Th)
			Th[non_na_mask & abs(Th) < eps] <- 0
		}

		if (show_significant) {
			idx <- which(Th != 0 | is.na(Th), arr.ind = TRUE)
		} else {
			idx <- which(Th != 0, arr.ind = TRUE)
		}
		if (nrow(idx) == 0) {
			cli::cli_warn("No edges above threshold eps = {eps}")
			df <- data.frame(i = 1, j = 1, val = 0)
		} else {
			df <- data.frame(
				i = idx[, 1],
				j = idx[, 2],
				val = Th[idx]
			)
		}

		ggplot2::ggplot(df, ggplot2::aes(i, j, color = val)) +
			ggplot2::geom_point(size = 1) +
			ggplot2::scale_color_gradient2(
				low = "blue", mid = "white", high = "red",
				midpoint = 0,
				name = expression(theta[ij]),
				na.value = "grey80"
			) +
			ggplot2::scale_y_reverse() +
			ggplot2::coord_equal() +
			ggplot2::labs(
				title = paste("Network snapshot at t =", time),
				subtitle = paste("Relation", rel, "- Sparse view"),
				x = "Actor i", y = "Actor j"
			) +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
	} else {
		# tile heatmap for smaller networks
		df <- expand.grid(i = 1:n_row, j = 1:n_col)
		df$val <- c(Th)

		ggplot2::ggplot(df, ggplot2::aes(i, j, fill = val)) +
			ggplot2::geom_tile() +
			ggplot2::scale_fill_gradient2(
				low = "navy", mid = "white", high = "firebrick",
				midpoint = 0,
				name = expression(theta[ij]),
				na.value = "grey80"
			) +
			ggplot2::scale_y_reverse() +
			ggplot2::coord_equal() +
			ggplot2::labs(
				title = paste("Network snapshot at t =", time),
				subtitle = paste("Relation", rel),
				x = "Actor i", y = "Actor j"
			) +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
	}
	####
}
####

####
#' Tidy DBN Summary
#'
#' @description Extract posterior means in tidy format
#' @param fit DBN object
#' @param what Components to extract
#' @param time_subset Time points to include (dynamic model)
#' @return List of posterior mean arrays
#' @seealso [tidy_dbn_lowrank()], [param_summary()], [theta_summary()]
#' @author Tosin Salau and Shahryar Minhas
#' @export
tidy_dbn <- function(fit, what = c("A", "B", "Theta"), time_subset = NULL) {
	what <- match.arg(what, several.ok = TRUE)

	# piecewise stores per-block 2D operators in fit$A_blocks / fit$B_blocks;
	# the downstream reduction below assumes a 3D time-indexed array. return
	# the block-level posterior means directly to sidestep the dim mismatch.
	if (fit$model == "piecewise") {
		out <- list()
		if ("A" %in% what && !is.null(fit$A_blocks)) out$A <- fit$A_blocks
		if ("B" %in% what && !is.null(fit$B_blocks)) out$B <- fit$B_blocks
		if ("Theta" %in% what) {
			if (!is.null(fit$Theta)) {
				out$Theta <- if (length(dim(fit$Theta)) == 5L)
					apply(fit$Theta, 1:4, mean) else fit$Theta
			}
		}
		return(out)
	}

	n_keep <- ifelse(fit$model == "static", dim(fit$B[[1]])[3], length(fit$A))

	if (is.null(time_subset)) {
		time_subset <- if (fit$model == "static") 1 else seq_len(dim(fit$A[[1]])[3])
	}

	out <- list()

	# posterior mean of A
	if ("A" %in% what && fit$model == "dynamic") {
		time_idx <- time_subset
		if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
			time_idx <- unique(ceiling(time_subset / fit$settings$time_thin))
			time_idx <- time_idx[time_idx <= dim(fit$A[[1]])[3]]
		}

		Amean <- Reduce(`+`, lapply(fit$A, function(a) a[, , time_idx, drop = FALSE])) / n_keep
		out$A <- Amean
	}
	####

	# posterior mean of B
	if ("B" %in% what && fit$model == "dynamic") {
		time_idx <- time_subset
		if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
			time_idx <- unique(ceiling(time_subset / fit$settings$time_thin))
			time_idx <- time_idx[time_idx <= dim(fit$B[[1]])[3]]
		}

		Bmean <- Reduce(`+`, lapply(fit$B, function(b) b[, , time_idx, drop = FALSE])) / n_keep
		out$B <- Bmean
	}
	####

	# posterior mean of Theta
	if ("Theta" %in% what) {
		if (fit$model == "static") {
			out$Theta <- fit$M
		} else if (!is.null(fit$Theta) && length(dim(fit$Theta)) == 5L) {
			# prefer the stored latent-state draws when available -- the
			# reconstruction below uses `A %*% t(B)`, which is non-conformable
			# for bipartite (n_row != n_col) and is anyway a different quantity
			# (the lag operator) from the latent state Theta
			out$Theta <- apply(fit$Theta, 1:4, mean)
		} else if (isTRUE(fit$dims$is_bipartite)) {
			cli::cli_abort(c(
				"{.fun tidy_dbn} cannot reconstruct {.code Theta} for bipartite fits without {.code fit$Theta}.",
				"i" = "Refit without dropping {.code \"Theta\"} from {.arg keep}, or read {.code fit$A} / {.code fit$B} directly."
			))
		} else {
			time_idx <- time_subset
			if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
				time_idx <- unique(ceiling(time_subset / fit$settings$time_thin))
				time_idx <- time_idx[time_idx <= dim(fit$A[[1]])[3]]
			}

			nr <- fit$dims$n_row
			ncl <- fit$dims$n_col
			Th <- array(0, dim = c(nr, ncl, fit$dims$p, length(time_idx)))
			for (s in seq_len(n_keep)) {
				for (i in seq_along(time_idx)) {
					for (rel in 1:fit$dims$p) {
						Th[, , rel, i] <- Th[, , rel, i] + fit$A[[s]][, , time_idx[i]] %*% t(fit$B[[s]][, , time_idx[i]])
					}
				}
			}
			out$Theta <- Th / n_keep
		}
	}
	####

	out
}
####

####
#' Plot Group Influence Profile
#'
#' @description Plots posterior group influence over time for dynamic models
#' @param fit A "dbn" object from dbn_dynamic()
#' @param group Integer vector of actor indices
#' @param type "sender" (rows of A_t) or "target" (columns of B_t)
#' @param fun Aggregation across actors: "mean" or "sum"
#' @param measure Per-actor metric: "rowsum" (default), "rowmean", "l2"
#' @param cred Credible band level (0.95 gives 95% bands)
#' @return A ggplot2 object
#' @seealso \code{\link{get_group_influence}},
#'   \code{\link{compare_group_influence}}, \code{\link{dbn}}
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 10, time = 10, seed = 6886)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'     nscan = 200, burn = 100, verbose = FALSE)
#' plot_group_influence(fit, group = c(1, 3, 5), type = "sender")
#' }
plot_group_influence <- function(fit,
								 group,
								 type = c("sender", "target"),
								 fun = c("mean", "sum"),
								 measure = c("rowsum", "rowmean", "l2"),
								 cred = 0.95) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort(c("Package {.pkg ggplot2} is required for this function.", "i" = "Install with {.code install.packages(\"ggplot2\")}"))
	if (fit$model != "dynamic") {
		cli::cli_abort("Group influence is defined for dynamic fit objects only.")
	}

	type <- match.arg(type)
	fun <- match.arg(fun)
	measure <- match.arg(measure)

	# sender groups index rows of A, target groups index columns of B
	n_row <- fit$dims$n_row
	n_col <- fit$dims$n_col
	m_check <- if (type == "sender") n_row else n_col
	if (any(group < 1) || any(group > m_check)) {
		cli::cli_abort("Group indices must be between 1 and {m_check}")
	}

	S <- length(fit$A)
	Tt <- dim(fit$A[[1]])[3]
	influence <- matrix(NA, S, Tt)

	row_fun <- switch(measure,
		rowsum = rowSums,
		rowmean = rowMeans,
		l2 = function(M) sqrt(rowSums(M^2))
	)

	col_fun <- switch(measure,
		rowsum = colSums,
		rowmean = colMeans,
		l2 = function(M) sqrt(colSums(M^2))
	)

	agg_fun <- if (fun == "mean") base::mean else base::sum

	# accumulate influence across MCMC draws
	for (s in seq_len(S)) {
		if (type == "sender") {
			for (t in seq_len(Tt)) {
				a <- fit$A[[s]][, , t]
				influence[s, t] <- agg_fun(row_fun(a[group, , drop = FALSE]))
			}
		} else {
			for (t in seq_len(Tt)) {
				b <- fit$B[[s]][, , t]
				influence[s, t] <- agg_fun(col_fun(b[, group, drop = FALSE]))
			}
		}
	}

	# compute posterior credible bands
	band <- apply(influence, 2, quantile, probs = c((1 - cred) / 2, 0.5, 1 - (1 - cred) / 2))

	# convert to original time scale (undo time thinning)
	time_vals <- 1:Tt
	if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
		time_vals <- (time_vals - 1) * fit$settings$time_thin + 1
	}

	df <- data.frame(
		time = time_vals,
		lo = band[1, ],
		med = band[2, ],
		hi = band[3, ]
	)

	gtitle <- sprintf(
		"%s group influence (%s %s)",
		if (type == "sender") "Sender" else "Target",
		fun, measure
	)

	ggplot2::ggplot(df, ggplot2::aes(time, med)) +
		ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80") +
		ggplot2::geom_line(linewidth = 0.8) +
		ggplot2::labs(
			title = gtitle,
			subtitle = paste("Actors:", paste(group, collapse = ", ")),
			y = sprintf("Posterior median %s %d%% CI", "\u00B1", round(cred * 100)),
			x = "Time"
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
}
####

####
#' Extract Group Influence Trajectories
#'
#' @description Computes group influence trajectories with posterior quantiles
#' @param fit A "dbn" object from dbn_dynamic()
#' @param group Integer vector of actor indices
#' @param type "sender" or "target"
#' @param measure Per-actor metric: "rowsum", "rowmean", "l2"
#' @param fun Aggregation across actors: "mean" or "sum"
#' @param probs Quantile probabilities to compute
#' @return Data frame with time, posterior quantiles, and mean
#' @seealso \code{\link{plot_group_influence}},
#'   \code{\link{compare_group_influence}}, \code{\link{dbn}}
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 10, time = 10, seed = 6886)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'     nscan = 200, burn = 100, verbose = FALSE)
#' inf_data <- get_group_influence(fit, group = c(1, 3, 5), type = "sender")
#' }
get_group_influence <- function(fit,
								group,
								type = c("sender", "target"),
								measure = c("rowsum", "rowmean", "l2"),
								fun = c("mean", "sum"),
								probs = c(0.025, 0.5, 0.975)) {
	if (fit$model != "dynamic") {
		cli::cli_abort("Group influence is defined for dynamic fit objects only.")
	}

	type <- match.arg(type)
	measure <- match.arg(measure)
	fun <- match.arg(fun)

	# sender groups index rows of A, target groups index columns of B
	n_row <- fit$dims$n_row
	n_col <- fit$dims$n_col
	m_check <- if (type == "sender") n_row else n_col
	if (any(group < 1) || any(group > m_check)) {
		cli::cli_abort("Group indices must be between 1 and {m_check}")
	}

	S <- length(fit$A)
	Tt <- dim(fit$A[[1]])[3]
	influence <- matrix(NA, S, Tt)

	row_fun <- switch(measure,
		rowsum = rowSums,
		rowmean = rowMeans,
		l2 = function(M) sqrt(rowSums(M^2))
	)

	col_fun <- switch(measure,
		rowsum = colSums,
		rowmean = colMeans,
		l2 = function(M) sqrt(colSums(M^2))
	)

	agg_fun <- if (fun == "mean") base::mean else base::sum

	# accumulate influence across MCMC draws
	for (s in seq_len(S)) {
		if (type == "sender") {
			for (t in seq_len(Tt)) {
				a <- fit$A[[s]][, , t]
				influence[s, t] <- agg_fun(row_fun(a[group, , drop = FALSE]))
			}
		} else {
			for (t in seq_len(Tt)) {
				b <- fit$B[[s]][, , t]
				influence[s, t] <- agg_fun(col_fun(b[, group, drop = FALSE]))
			}
		}
	}

	# posterior quantiles and mean across draws
	quants <- apply(influence, 2, quantile, probs = probs)
	means <- colMeans(influence)

	# convert to original time scale (undo time thinning)
	time_vals <- 1:Tt
	if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
		time_vals <- (time_vals - 1) * fit$settings$time_thin + 1
	}

	df <- data.frame(time = time_vals, mean = means)
	for (i in seq_along(probs)) {
		df[[paste0("q", probs[i])]] <- quants[i, ]
	}

	attr(df, "group") <- group
	attr(df, "type") <- type
	attr(df, "measure") <- measure
	attr(df, "fun") <- fun

	df
}
####

####
#' Compare Group Influences
#'
#' @description Compares influence trajectories of multiple groups
#' @param fit A "dbn" object from dbn_dynamic()
#' @param groups List of integer vectors, each defining a group
#' @param group_names Optional character vector of group names
#' @param type "sender" or "target"
#' @param measure Per-actor metric: "rowsum", "rowmean", "l2"
#' @param fun Aggregation: "mean" or "sum"
#' @param cred Credible band level
#' @return A ggplot2 object
#' @seealso \code{\link{plot_group_influence}},
#'   \code{\link{get_group_influence}}, \code{\link{dbn}}
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 10, time = 10, seed = 6886)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'     nscan = 200, burn = 100, verbose = FALSE)
#' compare_group_influence(fit,
#'     groups = list(c(1, 3, 5), c(2, 4, 6)),
#'     group_names = c("Group A", "Group B"))
#' }
compare_group_influence <- function(fit,
									groups,
									group_names = NULL,
									type = c("sender", "target"),
									measure = c("rowsum", "rowmean", "l2"),
									fun = c("mean", "sum"),
									cred = 0.95) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort(c("Package {.pkg ggplot2} is required for this function.", "i" = "Install with {.code install.packages(\"ggplot2\")}"))
	type <- match.arg(type)
	measure <- match.arg(measure)
	fun <- match.arg(fun)
	if (!is.list(groups)) {
		cli::cli_abort("groups must be a list of integer vectors")
	}

	n_groups <- length(groups)
	if (is.null(group_names)) {
		group_names <- paste("Group", 1:n_groups)
	} else if (length(group_names) != n_groups) {
		cli::cli_abort("group_names must have same length as groups")
	}

	all_data <- data.frame()

	for (i in seq_len(n_groups)) {
		inf_data <- get_group_influence(fit, groups[[i]], type, measure, fun,
			probs = c((1 - cred) / 2, 0.5, 1 - (1 - cred) / 2)
		)

		df_group <- data.frame(
			time = inf_data$time,
			median = inf_data[[paste0("q", 0.5)]],
			lo = inf_data[[paste0("q", (1 - cred) / 2)]],
			hi = inf_data[[paste0("q", 1 - (1 - cred) / 2)]],
			group = group_names[i]
		)

		all_data <- rbind(all_data, df_group)
	}

	ggplot2::ggplot(all_data, ggplot2::aes(
		x = time, y = median,
		color = group, fill = group
	)) +
		ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.3) +
		ggplot2::geom_line(linewidth = 1) +
		ggplot2::labs(
			title = sprintf(
				"%s group influence comparison (%s %s)",
				if (type == "sender") "Sender" else "Target", fun, measure
			),
			x = "Time",
			y = sprintf("Posterior median %s %d%% CI", "\u00B1", round(cred * 100)),
			color = "Group",
			fill = "Group"
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			legend.position = "top"
		)
}
####

####
#' Simulate data from static DBN model
#'
#' @description Generate posterior predictive samples from a fitted static
#'   DBN model.
#' @param fit A fitted `dbn` object from `dbn_static()`.
#' @param draws Number of posterior samples to generate.
#' @param summary Character string specifying summary type:
#'   - `"none"`: Return full array of simulations (default)
#'   - `"mean"`: Return posterior mean across simulations
#' @return If `summary = "none"`, returns a 5D array with dimensions
#'   `[nodes, nodes, relations, time, draws]`. If `summary = "mean"`,
#'   returns a 4D array with dimensions `[nodes, nodes, relations, time]`.
#' @keywords internal
simulate_static <- function(fit, draws = 100L, summary = "none") {
	n_row <- fit$dims$n_row
	n_col <- fit$dims$n_col
	p <- fit$dims$p
	Tt <- fit$dims$Tt

	if (!exists("tprod")) {
		cli::cli_abort("tprod function not found. Please load the dbn package.")
	}

	out <- array(0, dim = c(n_row, n_col, p, Tt, draws))

	n_saved <- dim(fit$B[[1]])[3]
	idx <- sample(n_saved, draws, replace = TRUE)

	for (s in seq_len(draws)) {
		Bdraw <- lapply(fit$B, function(b) b[, , idx[s]])

		Yrep <- tprod(fit$M, Bdraw)

		s2 <- fit$params[idx[s], "s2"]
		Yrep <- Yrep + array(rnorm(prod(dim(fit$M)), sd = sqrt(s2)),
			dim = dim(fit$M)
		)

		out[, , , , s] <- Yrep
	}

	if (identical(summary, "mean")) {
		apply(out, 1:4, mean)
	} else {
		out
	}
}
####

####
#' Forecast future network states from dynamic DBN model
#'
#' @description Generate H-step ahead forecasts from a fitted dynamic DBN
#'   model.
#' @param fit A fitted `dbn` object from `dbn_dynamic()`.
#' @param H Number of time steps to forecast ahead.
#' @param draws Number of posterior samples to generate.
#' @param summary Character string specifying summary type:
#'   - `"none"`: Return full array of forecasts (default)
#'   - `"mean"`: Return posterior mean forecasts
#' @return If `summary = "none"`, returns a 5D array with dimensions
#'   `[nodes, nodes, relations, horizon, draws]`. If `summary = "mean"`,
#'   returns a 4D array with dimensions `[nodes, nodes, relations, horizon]`.
#' @keywords internal
simulate_dynamic <- function(fit, H, draws = 100L, summary = "none", seed = NULL) {
	if (!is.null(seed)) set.seed(seed)
	n_row <- fit$dims$n_row
	n_col <- fit$dims$n_col
	p <- fit$dims$p
	Tt <- fit$dims$Tt

	if (!missing(H) && is.numeric(H) && length(H) == 1L && is.finite(H) && H > Tt) {
		cli::cli_warn(c(
			"Forecast horizon {.arg H} ({H}) exceeds the observed series length ({Tt}).",
			"i" = "Steps beyond {Tt} are pure extrapolation; credible intervals widen and a non-contractive operator will diverge.",
			"i" = "Check operator stability with {.fun dbn_operator} before trusting long-horizon forecasts."
		))
	}

	Theta_pred <- array(0, c(n_row, n_col, p, H, draws))

	n_saved <- length(fit$A)
	idx <- sample(n_saved, draws, replace = TRUE)

	# the forecast must propagate forward from the end of the observed
	# series, so seed the latent state at the last estimated Theta (not 0).
	# the last stored time index handles time-thinned fits correctly.
	have_theta <- !is.null(fit$Theta) && is.array(fit$Theta) &&
		length(dim(fit$Theta)) == 5
	t_last_theta <- if (have_theta) dim(fit$Theta)[4] else NA_integer_
	if (!have_theta) {
		cli::cli_warn(c(
			"Fitted {.code Theta} is not stored; seeding the forecast at zero.",
			"i" = "Refit with {.code store_theta = TRUE} to forecast from the observed end-of-series state."
		))
	}

	# propagate the operator forward via RW(1) for h > 1:
	# A_{T+h} = A_T + sum_{k=1..h} epsilon_k with epsilon ~ N(0, tau_A^2).
	# this gives the right h-step CIs under the RW prior. for ALS / TV-ALS
	# fits tau_A2 / tau_B2 may be NA; in that case we hold the operator
	# constant and warn that the forecast uncertainty is conditional on
	# the terminal operator value.
	tau_A2_avail <- !all(is.na(fit$tau_A2))
	tau_B2_avail <- !all(is.na(fit$tau_B2))
	propagate_operator <- tau_A2_avail && tau_B2_avail
	if (!propagate_operator && H > 1L &&
	    !isTRUE(getOption("dbn.terminal_operator_inform_fired", FALSE))) {
		cli::cli_inform(c(
			"i" = "Forecast holds the terminal operator (A_T, B_T) constant for h > 1.",
			"i" = "This fit has no RW innovation variance (typical for ALS); intervals are conditional on the terminal operator value.",
			"i" = "This message fires only once per session; clear with {.code options(dbn.terminal_operator_inform_fired = FALSE)}."
		))
		options(dbn.terminal_operator_inform_fired = TRUE)
	}

	# fit$sigma2 / tau_A2 / tau_B2 may be length-1 point estimates (ALS
	# bootstrap-expanded fits) or full per-draw chains (MCMC fits). use
	# safe indexing so a length-1 vector is reused across draws instead
	# of giving NA for s > 1.
	idx_scalar <- function(vec, m) {
		if (is.null(vec) || length(vec) == 0L) return(NA_real_)
		if (length(vec) == 1L) return(vec[1L])
		vec[m]
	}
	for (s in seq_len(draws)) {
		t_last_A <- dim(fit$A[[idx[s]]])[3]
		t_last_B <- dim(fit$B[[idx[s]]])[3]
		A_curr_op <- fit$A[[idx[s]]][, , t_last_A]
		B_curr_op <- fit$B[[idx[s]]][, , t_last_B]
		sigma2 <- idx_scalar(fit$sigma2, idx[s])
		if (!is.finite(sigma2)) sigma2 <- 1
		if (propagate_operator) {
			# scalar RW innovation sd per draw (uniform across entries)
			tauA <- sqrt(max(idx_scalar(fit$tau_A2, idx[s]), 0, na.rm = TRUE))
			tauB <- sqrt(max(idx_scalar(fit$tau_B2, idx[s]), 0, na.rm = TRUE))
		}
		# AR(1) coefficient for the operator state. when the fit was AR(1)
		# the prior is A_t = rho * A_{t-1} + eps, eps ~ N(0, tauA2 I); we
		# must honor that during forecast projection instead of defaulting
		# to pure random walk (rho = 1)
		rhoA <- if (is.numeric(fit$rhoA)) fit$rhoA[idx[s]] else 1
		rhoB <- if (is.numeric(fit$rhoB)) fit$rhoB[idx[s]] else 1
		if (!is.finite(rhoA)) rhoA <- 1
		if (!is.finite(rhoB)) rhoB <- 1

		if (have_theta) {
			Theta_curr <- fit$Theta[, , , t_last_theta, idx[s], drop = FALSE]
			dim(Theta_curr) <- c(n_row, n_col, p)
		} else {
			Theta_curr <- array(0, c(n_row, n_col, p))
		}
		# per-draw baseline mean for the model's centering form
		# Theta_t = M + A_t (Theta_{t-1} - M) B_t' + eps. fit$M is
		# [n_row, n_col, p, draws]; pull this draw's slice (M may be near
		# zero in the dynamic model, where the level sits in Theta).
		M_s <- array(0, c(n_row, n_col, p))
		if (is.array(fit$M) && length(dim(fit$M)) == 4L && dim(fit$M)[4] >= idx[s]) {
			M_s <- fit$M[, , , idx[s], drop = FALSE]
			dim(M_s) <- c(n_row, n_col, p)
		}

		for (h in seq_len(H)) {
			# propagate operator: A_{T+h} = rhoA * A_{T+h-1} + eps ~ N(0, tauA2 I)
			# RW case (default) has rhoA = 1; AR(1) fits give rhoA < 1
			if (propagate_operator && h >= 1L) {
				if (tauA > 0)
					A_curr_op <- rhoA * A_curr_op +
						matrix(rnorm(n_row * n_row, 0, tauA), n_row, n_row)
				else if (rhoA != 1)
					A_curr_op <- rhoA * A_curr_op
				if (tauB > 0)
					B_curr_op <- rhoB * B_curr_op +
						matrix(rnorm(n_col * n_col, 0, tauB), n_col, n_col)
				else if (rhoB != 1)
					B_curr_op <- rhoB * B_curr_op
			}
			Theta_new <- array(0, c(n_row, n_col, p))

			for (rel in seq_len(p)) {
				dev <- Theta_curr[, , rel] - M_s[, , rel]
				Theta_new[, , rel] <- M_s[, , rel] +
					A_curr_op %*% dev %*% t(B_curr_op) +
					sqrt(sigma2) * matrix(rnorm(n_row * n_col), n_row, n_col)
			}

			Theta_pred[, , , h, s] <- Theta_new
			Theta_curr <- Theta_new
		}
	}

	# map the latent forecast to the observation scale of the fit's family:
	# binary -> probability via probit (pnorm); ordinal -> empirical-cut
	# discretisation; gaussian -> identity.
	fam_name <- if (is.list(fit$family)) fit$family$name else fit$family
	if (identical(fam_name, "binary")) {
		Theta_pred[] <- stats::pnorm(Theta_pred)
	} else if (identical(fam_name, "ordinal")) {
		# discretise via the empirical cut points of the observed Y if available
		yvals <- sort(unique(stats::na.omit(as.numeric(fit$Y))))
		if (length(yvals) >= 2L) {
			cuts <- stats::quantile(as.numeric(fit$Theta), na.rm = TRUE,
				probs = seq(0, 1, length.out = length(yvals) + 1L))
			cuts[1] <- -Inf; cuts[length(cuts)] <- Inf
			tmp <- as.numeric(Theta_pred)
			disc <- cut(tmp, breaks = cuts, labels = FALSE, include.lowest = TRUE)
			Theta_pred[] <- yvals[disc]
		}
	}
	if (identical(summary, "mean")) {
		apply(Theta_pred, 1:4, mean)
	} else {
		Theta_pred
	}
}
####

####
#' Posterior Predictive Ordinal Data
#'
#' @description Generate ordinal data from posterior predictive distribution
#' @param fit A "dbn" object
#' @param draws Number of draws
#' @param H Forecast horizon (dynamic models)
#' @return Array of ordinal predictions
#' @keywords internal
predict_ordinal <- function(fit, draws = 100, H = NULL) {
	if (fit$model == "static") {
		Z_pred <- predict(fit, draws = draws, summary = "none")

		for (s in seq_len(draws)) {
			Z_pred[, , , , s] <- sweep(Z_pred[, , , , s], 1:3, fit$M, "+")
		}

		# map latent scores to ordinal categories via empirical CDF
		vals <- sort(unique(c(fit$Y %||% fit$R)))
		Z_vec <- c(Z_pred)

		R_pred <- vals[findInterval(
			Z_vec,
			quantile(c(fit$Y %||% fit$R),
				probs = seq(0, 1, length = length(vals) + 1)[-c(1, length(vals) + 1)],
				na.rm = TRUE
			)
		)]

		array(R_pred, dim = dim(Z_pred))
	} else {
		if (is.null(H)) H <- 1

		Theta_pred <- predict(fit, H = H, draws = draws, summary = "none")

		n_row <- fit$dims$n_row
		n_col <- fit$dims$n_col
		p <- fit$dims$p
		R_pred <- array(NA, dim = c(n_row, n_col, p, H, draws))

		vals <- sort(unique(c(fit$Y %||% fit$R)))

		for (s in seq_len(draws)) {
			M_sample <- fit$M[[sample(length(fit$M), 1)]]

			for (h in seq_len(H)) {
				Z_h <- sweep(Theta_pred[, , , h, s], 1:3, M_sample, "+")

				Z_vec <- c(Z_h)
				R_vec <- vals[findInterval(
					Z_vec,
					quantile(c(fit$Y %||% fit$R),
						probs = seq(0, 1, length = length(vals) + 1)[-c(1, length(vals) + 1)],
						na.rm = TRUE
					)
				)]

				R_pred[, , , h, s] <- array(R_vec, dim = c(n_row, n_col, p))
			}
		}

		R_pred
	}
}
####

####
#' Compare Posterior Means Across Different Samplers
#'
#' @description For two fits to the same data using different samplers
#'   (e.g., exact PCG vs approximate FFBS), compute the Frobenius-norm difference
#'   in posterior mean A and B matrices. Useful for assessing how much the
#'   approximation differs from the exact method.
#'
#' @param fit1 First fitted `dbn` object (typically exact PCG)
#' @param fit2 Second fitted `dbn` object (typically approximate FFBS)
#' @param verbose Logical. If TRUE, print summary.
#'
#' @return Invisibly, a list with:
#'   \item{A_diff_frob}{Frobenius norm of difference in posterior mean A matrices}
#'   \item{B_diff_frob}{Frobenius norm of difference in posterior mean B matrices}
#'   \item{sampler1}{Sampler used in fit1}
#'   \item{sampler2}{Sampler used in fit2}
#'
#' @author Tosin Salau and Shahryar Minhas
#' @export
compare_samplers <- function(fit1, fit2, verbose = TRUE) {
	if (!inherits(fit1, "dbn") || !inherits(fit2, "dbn")) {
		cli::cli_abort("Both arguments must be fitted {.cls dbn} objects.")
	}
	# A/B comparison only makes sense between two fits that both expose a
	# top-level A draws list. Static fits store A inside fit$B[[1]] (so $A is
	# NULL) and piecewise stores per-block matrices, so cross-comparisons
	# need the same model layout.
	if (!is.list(fit1$A) || !is.list(fit1$B) ||
	    !is.list(fit2$A) || !is.list(fit2$B)) {
		cli::cli_abort(c(
			"{.fun compare_samplers} requires two fits with top-level {.code $A} and {.code $B} draw lists (dynamic, hmm, or lowrank).",
			"x" = "Got models {.val {fit1$model}} and {.val {fit2$model}}.",
			"i" = "Static fits store A under {.code fit$B[[1]]}; piecewise stores per-block operators. These layouts are not directly comparable."
		))
	}

	sampler1 <- fit1$meta$sampler_used %||% "unknown"
	sampler2 <- fit2$meta$sampler_used %||% "unknown"

	# compute posterior means for A and B
	A1_mean <- apply(simplify2array(fit1$A), 1:2, mean)
	A2_mean <- apply(simplify2array(fit2$A), 1:2, mean)
	B1_mean <- apply(simplify2array(fit1$B), 1:2, mean)
	B2_mean <- apply(simplify2array(fit2$B), 1:2, mean)

	# frobenius norms of differences
	A_diff_frob <- norm(A1_mean - A2_mean, type = "F")
	B_diff_frob <- norm(B1_mean - B2_mean, type = "F")

	if (verbose) {
		cli::cli_h3("Sampler Comparison")
		cli::cli_inform("Sampler 1: {sampler1}")
		cli::cli_inform("Sampler 2: {sampler2}")
		cli::cli_inform("")
		cli::cli_text("Frobenius norm of difference in posterior mean A: {sprintf('%.4f', A_diff_frob)}")
		cli::cli_text("Frobenius norm of difference in posterior mean B: {sprintf('%.4f', B_diff_frob)}")
		if (A_diff_frob < 0.1 && B_diff_frob < 0.1) {
			cli::cli_alert_success("Samplers agree closely (low difference).")
		} else if (A_diff_frob < 0.5 || B_diff_frob < 0.5) {
			cli::cli_alert_info("Samplers show moderate agreement.")
		} else {
			cli::cli_alert_warning("Samplers differ substantially (check convergence).")
		}
	}

	invisible(list(
		A_diff_frob = A_diff_frob,
		B_diff_frob = B_diff_frob,
		sampler1 = sampler1,
		sampler2 = sampler2
	))
}
####
