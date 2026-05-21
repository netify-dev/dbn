####
#' Low-Rank DBN Posterior Analysis Methods
#'
#' @description Posterior analysis utilities for low-rank DBN models
#' @name lowrank-methods
#' @keywords internal
NULL
####

####
#' Plot Low-Rank DBN Results
#'
#' @description Diagnostic plots for low-rank DBN model fits
#' @param x A dbn object with model="lowrank"
#' @param factors_show Number of factors to display
#' @param time_points Optional subset of time points
#' @return List of ggplot objects or combined plot
#' @keywords internal
plot_lowrank <- function(x,
						 factors_show = min(3, x$settings$r),
						 time_points = NULL) {
	if (x$model != "lowrank") cli::cli_abort("Input is not a low-rank fit.")
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_abort(c(
			"package {.pkg ggplot2} is required for this function.",
			"i" = "install with {.code install.packages(\"ggplot2\")}"
		))
	}
	# avoid R CMD check notes
	alpha <- NULL

	# parameter trace plots
	df_trace <- data.frame(
		iter = seq_along(x$sigma2),
		sigma2 = x$sigma2,
		tau_alpha2 = x$tau_alpha2,
		tau_B2 = x$tau_B2,
		g2 = x$g2
	)

	trace_long <- do.call(rbind, lapply(c("sigma2", "tau_alpha2", "tau_B2", "g2"), function(par) {
		data.frame(iter = df_trace$iter, par = par, val = df_trace[[par]])
	}))

	p_trace <- ggplot2::ggplot(trace_long, ggplot2::aes(iter, val)) +
		ggplot2::geom_line() +
		ggplot2::facet_wrap(~par, scales = "free_y", ncol = 1) +
		ggplot2::labs(title = "MCMC traces", x = "Iteration", y = NULL) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0)
		)
	####

	# factor paths
	r <- x$settings$r
	S <- length(x$alpha)
	Tt_ <- ncol(x$alpha[[1]])

	time_vals <- if (isTRUE(x$settings$time_thin > 1)) {
		(seq_len(Tt_) - 1) * x$settings$time_thin + 1
	} else {
		seq_len(Tt_)
	}

	arr <- array(NA_real_, c(r, Tt_, S))
	for (s in seq_len(S)) arr[, , s] <- x$alpha[[s]]

	qs <- apply(arr, 1:2, quantile, probs = c(.025, .5, .975))

	alpha_df <- do.call(rbind, lapply(seq_len(factors_show), function(k) {
		data.frame(time = time_vals, lo = qs[1, k, ], med = qs[2, k, ],
			hi = qs[3, k, ], k = factor(k))
	}))

	p_alpha <- ggplot2::ggplot(alpha_df, ggplot2::aes(time, med)) +
		ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80") +
		ggplot2::geom_line() +
		ggplot2::facet_wrap(~k, scales = "free_y", ncol = 1,
			labeller = ggplot2::label_bquote(alpha[.(k)])) +
		ggplot2::labs(title = "Latent factor trajectories", x = "Time", y = expression(alpha)) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0)
		)
	####

	# node loading heatmap
	U_bar <- Reduce(`+`, x$U) / S

	grid <- expand.grid(actor = seq_len(nrow(U_bar)), factor = seq_len(ncol(U_bar)))
	df_U <- data.frame(actor = grid$actor, factor = grid$factor,
						loading = as.vector(U_bar))

	p_U <- ggplot2::ggplot(df_U, ggplot2::aes(factor, actor, fill = loading)) +
		ggplot2::geom_tile() +
		ggplot2::scale_fill_gradient2(low = "navy", mid = "white", high = "darkred", midpoint = 0) +
		ggplot2::coord_equal() +
		ggplot2::labs(title = "Posterior mean of U", x = "Factor k", y = "Actor i") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank())
	####

	if (requireNamespace("gridExtra", quietly = TRUE)) {
		gridExtra::grid.arrange(p_trace, p_alpha, p_U, ncol = 1)
	} else {
		list(trace = p_trace, alpha = p_alpha, U = p_U)
	}
}
####

####
#' Summary for Low-Rank Fits
#'
#' @description Posterior summaries for low-rank DBN models
#' @param object A dbn object with model="lowrank"
#' @param digits Number of digits to display
#' @param ... Ignored
#' @return Invisibly returns the input object
#' @keywords internal
summary_lowrank <- function(object, digits = 3, ...) {
	if (object$model != "lowrank") cli::cli_abort("Not a low-rank fit.")

	cat("Low-rank Dynamic Bilinear Network model\n")
	nr <- object$dims$n_row
	ncl <- object$dims$n_col
	cat(
		if (isTRUE(object$dims$is_bipartite)) paste0("  senders   : ", nr, "\n  receivers : ", ncl) else paste0("  nodes     : ", nr),
		"\n  relations :", object$dims$p,
		"\n  time pts  :", object$dims$Tt,
		"\n  rank      :", object$settings$r, "\n\n"
	)

	ss <- function(x) {
		formatC(c(mean(x, na.rm = TRUE), quantile(x, c(.025, .975), na.rm = TRUE)),
			digits = digits, format = "f")
	}
	sm <- rbind(sigma2 = ss(object$sigma2), tau_alpha2 = ss(object$tau_alpha2),
		tau_B2 = ss(object$tau_B2), g2 = ss(object$g2))
	colnames(sm) <- c("mean", "2.5%", "97.5%")
	print(sm, quote = FALSE)
	invisible(object)
}
####

####
#' Predict from Low-Rank Fit
#'
#' @description H-step-ahead forecasts using random-walk extension of alpha
#' @param object A dbn object with model="lowrank"
#' @param H Number of forecast steps
#' @param draws Number of posterior draws
#' @param summary "mean" for posterior mean, "none" for all draws
#' @return Prediction array
#' @keywords internal
predict_lowrank <- function(object, H = 1, draws = 100,
							summary = c("mean", "none")) {
	if (object$model != "lowrank") cli::cli_abort("Not a low-rank fit.")
	summary <- match.arg(summary)

	n_row <- object$dims$n_row
	n_col <- object$dims$n_col
	p <- object$dims$p
	r <- object$settings$r
	S_saved <- length(object$U)
	pick <- sample(S_saved, draws, TRUE)

	Theta_pred <- array(0, c(n_row, n_col, p, H, draws))

	for (d in seq_len(draws)) {
		s <- pick[d]
		sigma2 <- object$sigma2[s]
		U <- object$U[[s]]
		B_now <- object$B[[s]][, , dim(object$B[[s]])[3]]
		alpha_path <- object$alpha[[s]]

		for (h in seq_len(H)) {
			a_t <- alpha_path[, ncol(alpha_path)] +
				rnorm(r, 0, sqrt(object$tau_alpha2[s]))
			A_t <- U %*% diag(a_t) %*% t(U)

			for (rel in seq_len(p)) {
				Theta_pred[, , rel, h, d] <- A_t %*% Theta_pred[, , rel, max(h - 1, 1), d] %*% t(B_now) +
					matrix(rnorm(n_row * n_col, 0, sqrt(sigma2)), n_row, n_col)
			}
			alpha_path <- cbind(alpha_path, a_t)
		}
	}

	if (summary == "mean") apply(Theta_pred, 1:4, mean) else Theta_pred
}
####

####
#' Tidy Extractor for Low-Rank Factor Trajectories
#'
#' @description Extracts the time-varying factor strengths \eqn{\alpha_t}
#'   from a fitted low-rank model in tidy data-frame format, ready for
#'   plotting with ggplot2. Each factor captures how strongly one latent
#'   dimension of the influence structure drives the dynamics at each time
#'   point. The returned data frame includes posterior means and 95\%
#'   credible intervals.
#' @param fit A fitted `dbn` object with `model = "lowrank"`.
#' @param factors Integer vector specifying which factors to extract
#'   (default: all `r` factors). For example, `factors = 1:2` extracts the
#'   first two factors.
#' @return Data frame with columns: `time`, `mean` (posterior mean of
#'   \eqn{\alpha_k(t)}), `lo` (2.5th percentile), `hi` (97.5th
#'   percentile), and `factor` (factor index).
#' @seealso [dbn()] with `model = "lowrank"`, [plot.dbn()] for built-in
#'   diagnostic plots
#' @examples
#' \donttest{
#' sim <- simulate_lowrank_dbn(n = 8, time = 8, r = 2, seed = 1)
#' fit <- dbn_lowrank(sim$Y, r = 2, nscan = 200, burn = 100, verbose = FALSE)
#' df_alpha <- tidy_dbn_lowrank(fit)
#' head(df_alpha)
#' }
#' @export
tidy_dbn_lowrank <- function(fit, factors = NULL) {
	if (is.null(factors)) factors <- 1:fit$settings$r
	if (fit$model != "lowrank") cli::cli_abort("Not a low-rank fit.")
	S <- length(fit$alpha)
	Tt_ <- ncol(fit$alpha[[1]])

	time_vals <- if (isTRUE(fit$settings$time_thin > 1)) {
		(seq_len(Tt_) - 1) * fit$settings$time_thin + 1
	} else {
		seq_len(Tt_)
	}

	out <- lapply(factors, function(k) {
		mat <- sapply(fit$alpha, function(a) a[k, ])
		if (!is.matrix(mat)) mat <- matrix(mat, nrow = Tt_, ncol = S)
		data.frame(
			time = time_vals,
			mean = rowMeans(mat),
			lo = apply(mat, 1, quantile, .025),
			hi = apply(mat, 1, quantile, .975),
			factor = k
		)
	})
	do.call(rbind, out)
}
####
