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
		ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
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
		ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
	####

	# parameter trace plots -- collect whichever scalar chains are present.
	# the HMM model stores process variance as `sigma2_proc` (not `sigma2`),
	# and only the gaussian family carries `sigma2_obs`; build the trace from
	# the chains that actually exist so the plot works for every family.
	cand <- list(
		sigma2_proc = x$sigma2_proc, sigma2_obs = x$sigma2_obs,
		tau_A2 = x$tau_A2, tau_B2 = x$tau_B2, g2 = x$g2
	)
	niter <- max(lengths(cand), 0L)
	cand <- cand[vapply(cand, function(v) length(v) == niter && niter > 0L, logical(1))]
	tl <- do.call(rbind, lapply(names(cand), function(v)
		data.frame(iter = seq_len(niter), variable = v, value = cand[[v]])))

	p_trace <- ggplot2::ggplot(tl, ggplot2::aes(iter, value)) +
		ggplot2::geom_line() +
		ggplot2::facet_wrap(~variable, scales = "free_y", ncol = 1) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
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
	# the HMM model stores process variance as `sigma2_proc`, not `sigma2`
	sm <- rbind(sigma2 = ss(object$sigma2 %||% object$sigma2_proc),
		tau_A2 = ss(object$tau_A2),
		tau_B2 = ss(object$tau_B2), g2 = ss(object$g2))
	if (!is.null(object$sigma2_obs)) {
		sm <- rbind(sm, sigma2_obs = ss(object$sigma2_obs))
	}
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
						summary = c("mean", "none"), seed = NULL) {
	if (object$model != "hmm") cli::cli_abort("Not an HMM fit.")
	if (!is.null(seed)) set.seed(seed)
	summary <- match.arg(summary)

	n_row <- object$dims$n_row
	n_col <- object$dims$n_col
	p <- object$dims$p
	R <- object$settings$R
	S_saved <- length(object$A)
	pick <- sample(S_saved, draws, TRUE)

	Theta_pred <- array(0, c(n_row, n_col, p, H, draws))

	# seed the forecast from the last estimated latent state, not from zero
	# (a zero seed makes the h = 1 step A %*% 0 %*% B' = pure noise). HMM
	# stores per-draw latent draws in $draws$theta and per-draw baselines in
	# $draws$misc$M; fall back to the observed last slice / zero baseline.
	theta_draws <- object$draws$theta
	M_draws <- object$draws$misc$M
	Tt <- object$dims$Tt
	Y_last <- if (!is.null(object$Y)) object$Y[, , , Tt, drop = FALSE] else NULL

	for (d in seq_len(draws)) {
		s <- pick[d]
		A_list <- object$A[[s]]
		B_list <- object$B[[s]]
		Pi <- object$Pi[[s]]
		sigma2 <- (object$sigma2 %||% object$sigma2_proc)[s]

		# per-draw baseline mean
		M_s <- array(0, c(n_row, n_col, p))
		if (is.list(M_draws) && length(M_draws) >= s && !is.null(M_draws[[s]])) {
			M_s <- array(M_draws[[s]], dim = c(n_row, n_col, p))
		}
		# lag-0 latent state: last stored Theta draw, else last observed slice
		Theta_now <- array(0, c(n_row, n_col, p))
		if (is.list(theta_draws) && length(theta_draws) >= s &&
		    !is.null(theta_draws[[s]]) && length(dim(theta_draws[[s]])) == 4L) {
			Theta_now <- array(theta_draws[[s]][, , , Tt], dim = c(n_row, n_col, p))
		} else if (!is.null(Y_last)) {
			yl <- Y_last; yl[!is.finite(yl)] <- 0
			Theta_now <- array(yl, dim = c(n_row, n_col, p))
		}

		regime <- sample(R, 1, prob = colMeans(Pi))

		for (h in seq_len(H)) {
			regime <- sample(R, 1, prob = Pi[regime, ])
			A_t <- A_list[, , regime]
			B_t <- B_list[, , regime]
			Theta_next <- array(0, c(n_row, n_col, p))
			for (rel in seq_len(p)) {
				dev <- Theta_now[, , rel] - M_s[, , rel]
				Theta_next[, , rel] <- M_s[, , rel] +
					A_t %*% dev %*% t(B_t) +
					matrix(rnorm(n_row * n_col, 0, sqrt(sigma2)), n_row, n_col)
			}
			Theta_pred[, , , h, d] <- Theta_next
			Theta_now <- Theta_next
		}
	}

	if (identical(summary, "mean")) apply(Theta_pred, 1:4, mean) else Theta_pred
}
####
