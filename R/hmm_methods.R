####
#' HMM-DBN Posterior Analysis Methods
#'
#' @description Posterior analysis utilities for HMM-DBN models
#' @name hmm-methods
#' @keywords internal
NULL
####

####
#' Plot HMM-DBN Results
#'
#' @description Diagnostic plots for HMM-DBN model fits
#' @param x A dbn object with model="hmm"
#' @return List of ggplot objects or combined plot
#' @keywords internal
plot_hmm <- function(x) {
	if (x$model != "hmm") cli::cli_abort("Not an HMM fit.")
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_abort(c(
			"package {.pkg ggplot2} is required for this function.",
			"i" = "install with {.code install.packages(\"ggplot2\")}"
		))
	}

	# regime probability heatmap
	S_mat <- do.call(cbind, lapply(x$S, identity))
	R <- x$settings$R
	Tt_ <- nrow(S_mat)
	probs <- sapply(1:R, function(r) rowMeans(S_mat == r))

	df_reg <- data.frame()
	for (t in 1:Tt_) {
		for (r in 1:R) {
			df_reg <- rbind(df_reg, data.frame(time = t, regime = r, prob = probs[t, r]))
		}
	}
	df_reg$regime <- factor(df_reg$regime)

	p_reg <- ggplot2::ggplot(df_reg, ggplot2::aes(time, regime, fill = prob)) +
		ggplot2::geom_tile() +
		ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
		ggplot2::labs(title = "Posterior regime probabilities", x = "Time", y = "Regime") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank())
	####

	# transition matrix
	Pi_bar <- Reduce(`+`, x$Pi) / length(x$Pi)

	df_Pi <- data.frame()
	for (i in 1:nrow(Pi_bar)) {
		for (j in 1:ncol(Pi_bar)) {
			df_Pi <- rbind(df_Pi, data.frame(from = i, to = j, prob = Pi_bar[i, j]))
		}
	}
	df_Pi$from <- factor(df_Pi$from)
	df_Pi$to <- factor(df_Pi$to)

	p_Pi <- ggplot2::ggplot(df_Pi, ggplot2::aes(from, to, fill = prob)) +
		ggplot2::geom_tile() +
		ggplot2::scale_fill_gradient(low = "white", high = "firebrick") +
		ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", prob)), size = 3) +
		ggplot2::labs(title = "Posterior mean transition matrix Pi", x = "From", y = "To") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank())
	####

	# parameter trace plots
	trace <- data.frame(
		iter = seq_along(x$sigma2),
		sigma2 = x$sigma2,
		tau_A2 = x$tau_A2,
		tau_B2 = x$tau_B2,
		g2 = x$g2
	)

	tl <- data.frame()
	for (var in c("sigma2", "tau_A2", "tau_B2", "g2")) {
		tl <- rbind(tl, data.frame(iter = trace$iter, variable = var, value = trace[[var]]))
	}

	p_trace <- ggplot2::ggplot(tl, ggplot2::aes(iter, value)) +
		ggplot2::geom_line(colour = "darkgreen") +
		ggplot2::facet_wrap(~variable, scales = "free_y", ncol = 1) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0)
		) +
		ggplot2::labs(title = "MCMC traces", x = "Iteration", y = NULL)
	####

	if (requireNamespace("gridExtra", quietly = TRUE)) {
		gridExtra::grid.arrange(p_reg, p_Pi, p_trace, ncol = 1)
	} else {
		list(regime = p_reg, Pi = p_Pi, trace = p_trace)
	}
}
####

####
#' Summary for HMM-DBN Fits
#'
#' @description Posterior summaries for HMM-DBN models
#' @param object A dbn object with model="hmm"
#' @param digits Number of digits to display
#' @param ... Ignored
#' @return Invisibly returns the input object
#' @keywords internal
summary_hmm <- function(object, digits = 3, ...) {
	if (object$model != "hmm") cli::cli_abort("Not an HMM fit.")

	cat("Regime-switching (HMM) DBN model\n")
	nr <- object$dims$n_row
	ncl <- object$dims$n_col
	cat(
		if (isTRUE(object$dims$is_bipartite)) paste0("  senders   : ", nr, "\n  receivers : ", ncl) else paste0("  nodes     : ", nr),
		"\n  relations :", object$dims$p,
		"\n  time pts  :", object$dims$Tt,
		"\n  regimes   :", object$settings$R, "\n\n"
	)

	ss <- function(x) {
		formatC(c(mean(x, na.rm = TRUE), quantile(x, c(.025, .975), na.rm = TRUE)),
			digits = digits, format = "f")
	}
	sm <- rbind(sigma2 = ss(object$sigma2), tau_A2 = ss(object$tau_A2),
		tau_B2 = ss(object$tau_B2), g2 = ss(object$g2))
	colnames(sm) <- c("mean", "2.5%", "97.5%")
	print(sm, quote = FALSE)

	cat("\nPosterior mean transition matrix Pi:\n")
	print(round(Reduce(`+`, object$Pi) / length(object$Pi), 3))
	invisible(object)
}
####

####
#' Predict from HMM-DBN Fit
#'
#' @description Generates H-step-ahead forecasts from an HMM-DBN model
#' @param object A dbn object with model="hmm"
#' @param H Number of forecast steps
#' @param draws Number of posterior draws
#' @param summary "mean" for posterior mean, "none" for all draws
#' @return Prediction array
#' @keywords internal
predict_hmm <- function(object, H = 1, draws = 100,
						summary = c("mean", "none")) {
	if (object$model != "hmm") cli::cli_abort("Not an HMM fit.")
	summary <- match.arg(summary)

	n_row <- object$dims$n_row
	n_col <- object$dims$n_col
	p <- object$dims$p
	R <- object$settings$R
	S_saved <- length(object$A)
	pick <- sample(S_saved, draws, TRUE)

	Theta_pred <- array(0, c(n_row, n_col, p, H, draws))

	for (d in seq_len(draws)) {
		s <- pick[d]
		A_list <- object$A[[s]]
		B_list <- object$B[[s]]
		Pi <- object$Pi[[s]]
		sigma2 <- object$sigma2[s]

		regime <- sample(R, 1, prob = colMeans(Pi))
		Theta_now <- array(0, c(n_row, n_col, p))

		for (h in seq_len(H)) {
			regime <- sample(R, 1, prob = Pi[regime, ])
			A_t <- A_list[, , regime]
			B_t <- B_list[, , regime]
			for (rel in seq_len(p)) {
				Theta_now[, , rel] <- A_t %*% Theta_now[, , rel] %*% t(B_t) +
					matrix(rnorm(n_row * n_col, 0, sqrt(sigma2)), n_row, n_col)
			}
			Theta_pred[, , , h, d] <- Theta_now
		}
	}

	if (summary == "mean") apply(Theta_pred, 1:4, mean) else Theta_pred
}
####
