#' Effective SNR Audit for a Dynamic Bilinear AR Fit
#'
#' @description Computes the per-period effective signal-to-noise ratio for a
#'   fitted symmetric dynamic bilinear AR model. The signal is the off-diagonal
#'   Frobenius norm of the bilinear projection \eqn{A_t \Theta_{t-1} A_t^\top};
#'   the noise is the off-diagonal Frobenius norm of the residual
#'   \eqn{\Theta_t - A_t \Theta_{t-1} A_t^\top}. Computed per posterior draw and
#'   summarized across draws.
#'
#' @param fit A \code{dbn} object returned by \code{\link{dbn_dynamic}} fit with
#'   \code{symmetric = TRUE}.
#' @param windows Optional named list of integer time indices defining windows
#'   over which to average. If NULL, only the full-sample average is reported.
#'
#' @return List with components:
#'   \describe{
#'     \item{snr_per_draw_per_t}{Matrix of size \code{n_draws x (Tt-1)}.}
#'     \item{full_sample}{Named vector with mean, sd, q025, q975 across draws of
#'       the across-time average.}
#'     \item{windows}{Optional list, one per window, each with the same fields.}
#'   }
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 8, p = 1, time = 12, seed = 1)
#' for (t in seq_len(dim(sim$Y)[4])) {
#'   Yt <- sim$Y[, , 1, t]; Yt <- 0.5 * (Yt + t(Yt)); diag(Yt) <- 0
#'   sim$Y[, , 1, t] <- Yt
#' }
#' fit <- dbn_dynamic(sim$Y, family = "gaussian", symmetric = TRUE,
#'                    lambda_diag = 0.5,
#'                    nscan = 300, burn = 100, verbose = FALSE)
#' snr <- dbn_compute_snr(fit, windows = list(early = 1:6, late = 7:12))
#' snr$full_sample
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_compute_snr <- function(fit, windows = NULL) {
	if (!inherits(fit, "dbn")) {
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object returned by {.fun dbn_dynamic} or {.fun dbn}.")
	}
	# ALS+bootstrap fits have per-replicate A and Theta draws after expansion;
	# SNR is conceptually well-defined under fixed-design bootstrap. ALS without
	# bootstrap has only a single point estimate -- refuse with a hint.
	su <- fit$meta$sampler_used
	is_als <- !is.null(su) && su %in% c("als", "als_tv", "als_piecewise")
	if (is_als) {
		if (is.null(fit$bootstrap)) {
			cli::cli_abort(c(
				"{.fun dbn_compute_snr} needs replicate draws of {.code A} and {.code Theta}.",
				"x" = "This is a point-estimate fit from {.fun dbn_als} with no bootstrap.",
				"i" = "Refit with {.code dbn_als(..., bootstrap = N)} to get bootstrap replicate SNR, or use {.fun dbn} for posterior draws."
			))
		}
		if (!isTRUE(fit$meta$bootstrap_expanded)) {
			fit <- .bootstrap_expand(fit)
		}
	}
	if (is.null(fit$A) || is.null(fit$Theta)) {
		cli::cli_abort(c(
			"{.fun dbn_compute_snr} requires both {.code fit$A} and {.code fit$Theta}.",
			"x" = "{.code fit$A}: {ifelse(is.null(fit$A), \"NULL\", \"OK\")} ; {.code fit$Theta}: {ifelse(is.null(fit$Theta), \"NULL\", \"OK\")}",
			"i" = "Refit the model and ensure {.code keep} does not drop {.code \"Theta\"}."
		))
	}
	# the signal term below is A_t Theta A_t^T, which assumes the symmetric
	# specification (B = A); on an asymmetric fit it is the wrong estimand
	if (isFALSE(fit$dims$is_symmetric)) {
		cli::cli_abort(c(
			"{.fun dbn_compute_snr} requires a symmetric fit ({.code symmetric = TRUE}).",
			"i" = "The signal-to-residual decomposition assumes B = A; on an asymmetric fit it would compute the wrong quantity."
		))
	}

	n_draws <- length(fit$A)
	Theta_arr <- fit$Theta
	# theta is n_row x n_col x p x Tt x n_draws; collapse p (must be 1)
	if (length(dim(Theta_arr)) == 5 && dim(Theta_arr)[3] == 1) {
		Theta4 <- array(Theta_arr, dim = c(dim(Theta_arr)[1], dim(Theta_arr)[2],
			dim(Theta_arr)[4], dim(Theta_arr)[5]))
	} else {
		stop("dbn_compute_snr() requires a unipartite (p = 1) symmetric fit.")
	}
	n <- dim(Theta4)[1]
	Tt <- dim(Theta4)[3]
	off_mask <- upper.tri(matrix(0, n, n))

	snr <- matrix(NA_real_, n_draws, Tt - 1)
	for (m in seq_len(n_draws)) {
		A_m <- fit$A[[m]]
		# A draws may be n x n x Tt or n x n x time_keep; assume aligned
		for (t in 2:Tt) {
			A_t <- A_m[, , min(t, dim(A_m)[3])]
			Th_prev <- Theta4[, , t - 1, m]
			Th_now <- Theta4[, , t, m]
			signal <- A_t %*% Th_prev %*% t(A_t)
			resid <- Th_now - signal
			s_f <- sqrt(sum(signal[off_mask]^2))
			n_f <- sqrt(sum(resid[off_mask]^2))
			snr[m, t - 1] <- s_f / max(n_f, 1e-12)
		}
	}

	summarize <- function(x) {
		c(mean = mean(x), sd = sd(x),
		  q025 = unname(quantile(x, 0.025)),
		  q975 = unname(quantile(x, 0.975)))
	}

	full_avg <- rowMeans(snr)
	out <- list(
		snr_per_draw_per_t = snr,
		full_sample = summarize(full_avg)
	)

	if (!is.null(windows)) {
		# auto-label unnamed windows so a list like list(c(1,6), c(7,12)) is
		# iterable by name; without this the names(windows) loop below
		# silently produces an empty result.
		if (is.null(names(windows))) {
			names(windows) <- paste0("window_", seq_along(windows))
		} else {
			blanks <- !nzchar(names(windows))
			if (any(blanks)) names(windows)[blanks] <- paste0("window_", which(blanks))
		}
		win_results <- list()
		for (wnm in names(windows)) {
			w_idx <- windows[[wnm]]
			w_idx <- w_idx[w_idx >= 2 & w_idx <= Tt] - 1L
			if (length(w_idx) > 0) {
				w_avg <- rowMeans(snr[, w_idx, drop = FALSE])
				win_results[[wnm]] <- summarize(w_avg)
			}
		}
		out$windows <- win_results
	}

	class(out) <- c("dbn_snr", "list")
	out
}


#' Posterior Coupling Rank Probabilities
#'
#' @description For a symmetric dynamic bilinear AR fit, computes the
#'   posterior probabilities that each actor's time-averaged coupling
#'   \eqn{\bar C_i(\mathcal{T}) = T^{-1}\sum_{t} \|W_{t,i\cdot}\|_2} occupies
#'   the top-K rank, and pairwise dominance probabilities. Implements the
#'   rank-probability reporting framework recommended in Section 5.5 of the
#'   methods paper.
#'
#' @param fit A \code{dbn} object with \code{symmetric = TRUE}.
#' @param windows Optional named list of integer time indices defining windows.
#'   If NULL, the full sample is used.
#' @param K_grid Integer vector of top-K thresholds at which to report
#'   \eqn{p_{i,K}}.
#' @param actor_names Optional character vector of actor names; otherwise
#'   uses \code{actor_1, actor_2, ...}.
#'
#' @return List, one entry per window plus \code{full_sample}, with
#'   \code{post_mean}, \code{post_sd}, \code{post_q025}, \code{post_q975},
#'   \code{p_topK} (named list keyed by K), and \code{pi_dom} (n x n
#'   pairwise dominance matrix).
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 8, p = 1, time = 12, seed = 1)
#' for (t in seq_len(dim(sim$Y)[4])) {
#'   Yt <- sim$Y[, , 1, t]; Yt <- 0.5 * (Yt + t(Yt)); diag(Yt) <- 0
#'   sim$Y[, , 1, t] <- Yt
#' }
#' fit <- dbn_dynamic(sim$Y, family = "gaussian", symmetric = TRUE,
#'                    nscan = 300, burn = 100, verbose = FALSE)
#' rp <- dbn_coupling_rank_probs(fit, K_grid = c(1, 3, 5))
#' rp$full_sample$p_topK$p_top3
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_coupling_rank_probs <- function(fit, windows = NULL, K_grid = c(1, 3, 5),
	actor_names = NULL) {
	if (!inherits(fit, "dbn")) {
		stop("fit must be a dbn object returned by dbn_dynamic().")
	}
	if (is.null(fit$A)) {
		stop("fit must contain $A draws.")
	}
	# bootstrap-aware path: an ALS fit with `fit$bootstrap` carries R A draws
	# in vectorised form. Expand them so the per-draw rank loop below sees R
	# A's and the resulting probabilities are honest bootstrap rates.
	su <- fit$meta$sampler_used
	is_als <- length(su) == 1L && !is.na(su) && su %in% c("als", "als_tv")
	if (is_als && !is.null(fit$bootstrap) && !isTRUE(fit$meta$bootstrap_expanded)) {
		fit <- .bootstrap_expand(fit)
	}
	n_draws <- length(fit$A)
	A_first <- fit$A[[1]]
	n <- dim(A_first)[1]
	Tt <- dim(A_first)[3]
	if (is.null(actor_names)) actor_names <- paste0("actor_", seq_len(n))

	# coupling below is built from W_t = A_t A_t^T, i.e. the symmetric
	# (B = A) specification; on an asymmetric fit B is ignored and the result
	# is not the A_t B_t^T coupling
	if (isFALSE(fit$dims$is_symmetric)) {
		cli::cli_abort(c(
			"{.fun dbn_coupling_rank_probs} requires a symmetric fit.",
			"x" = "This fit is asymmetric; coupling built from {.code A A^T} ignores B and is the wrong quantity.",
			"i" = "Refit with {.code symmetric = TRUE} for a valid coupling diagnostic."
		))
	}
	# warn only when we have neither MCMC draws nor a bootstrap expansion.
	is_bootstrap <- isTRUE(fit$meta$bootstrap_expanded)
	if (n_draws < 2L || (is_als && !is_bootstrap) ||
		(isFALSE(fit$meta$uncertainty_available) && !is_bootstrap)) {
		cli::cli_warn(c(
			"Rank probabilities from a point-estimate fit are degenerate.",
			"i" = "{.code p_topK} and {.code pi_dom} are 0/1 indicators from a single estimate, not posterior probabilities.",
			"i" = "Refit with {.code bootstrap = N} or use a full {.fun dbn} MCMC fit."
		))
	}

	# per-draw, per-actor, per-time coupling: C_{m,i,t} = ||W_t row i||_2
	coupling <- array(0, dim = c(n_draws, n, Tt),
		dimnames = list(draw = NULL, actor = actor_names, t = NULL))
	for (m in seq_len(n_draws)) {
		A_m <- fit$A[[m]]
		for (t in seq_len(Tt)) {
			A_t <- A_m[, , t]
			W_t <- A_t %*% t(A_t)
			coupling[m, , t] <- sqrt(rowSums(W_t^2))
		}
	}

	if (is.null(windows)) windows <- list(full_sample = seq_len(Tt))
	# auto-label unnamed windows -- matches the behavior of dbn_compute_snr()
	if (is.null(names(windows))) {
		names(windows) <- paste0("window_", seq_along(windows))
	} else {
		blanks <- !nzchar(names(windows))
		if (any(blanks)) names(windows)[blanks] <- paste0("window_", which(blanks))
	}

	results <- list()
	for (wnm in names(windows)) {
		w_idx <- windows[[wnm]]
		Cbar <- apply(coupling[, , w_idx, drop = FALSE], c(1, 2), mean)
		ranks_draws <- t(apply(-Cbar, 1, rank, ties.method = "min"))
		p_topK <- lapply(K_grid, function(K) colMeans(ranks_draws <= K))
		names(p_topK) <- paste0("p_top", K_grid)

		pi_dom <- matrix(0, n, n, dimnames = list(actor_names, actor_names))
		for (i in seq_len(n)) for (j in seq_len(n)) {
			if (i != j) pi_dom[i, j] <- mean(Cbar[, i] > Cbar[, j])
		}

		results[[wnm]] <- list(
			window = wnm, t_indices = w_idx,
			Cbar = Cbar,
			post_mean = colMeans(Cbar),
			post_sd = apply(Cbar, 2, sd),
			post_q025 = apply(Cbar, 2, quantile, 0.025),
			post_q975 = apply(Cbar, 2, quantile, 0.975),
			p_topK = p_topK,
			pi_dom = pi_dom
		)
	}

	class(results) <- c("dbn_rank_probs", "list")
	results
}

#' Print method for DBN coupling rank-probability objects
#'
#' @description Prints a compact, per-window summary of an object returned by
#'   [dbn_coupling_rank_probs()] instead of dumping the full nested list
#'   (which includes the large per-draw `Cbar` matrix).
#' @param x A `dbn_rank_probs` object.
#' @param digits Number of digits for coupling values (default 3).
#' @param top Number of top-coupled actors to show per window (default 6).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method print dbn_rank_probs
print.dbn_rank_probs <- function(x, digits = 3, top = 6, ...) {
	cli::cli_h2("DBN actor coupling rank probabilities")
	for (wnm in names(x)) {
		w <- x[[wnm]]
		n <- length(w$post_mean)
		anames <- names(w$post_mean)
		if (is.null(anames)) anames <- paste0("actor_", seq_len(n))
		ord <- order(w$post_mean, decreasing = TRUE)
		cli::cli_h3("Window {.field {wnm}} ({n} actor{?s})")
		for (i in ord[seq_len(min(top, n))]) {
			cli::cli_text("  {anames[i]}: coupling {round(w$post_mean[i], digits)} [{round(w$post_q025[i], digits)}, {round(w$post_q975[i], digits)}]")
		}
		if (n > top) cli::cli_text("  {cli::symbol$ellipsis} {n - top} more actor{?s}")
	}
	cli::cli_text("")
	cli::cli_text("Top-K probabilities in {.code $<window>$p_topK}; pairwise dominance in {.code $<window>$pi_dom}.")
	invisible(x)
}

#' Print method for DBN signal-to-residual audit objects
#'
#' @description Prints a compact summary of an object returned by
#'   [dbn_compute_snr()] instead of dumping the full per-draw matrix.
#' @param x A `dbn_snr` object.
#' @param digits Number of digits (default 3).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method print dbn_snr
print.dbn_snr <- function(x, digits = 3, ...) {
	cli::cli_h2("DBN signal-to-residual audit")
	fs <- x$full_sample
	cli::cli_text("Full sample: {round(fs[['mean']], digits)} [{round(fs[['q025']], digits)}, {round(fs[['q975']], digits)}]")
	if (!is.null(x$windows) && length(x$windows)) {
		cli::cli_h3("By window")
		for (wnm in names(x$windows)) {
			w <- x$windows[[wnm]]
			cli::cli_text("  {wnm}: {round(w[['mean']], digits)} [{round(w[['q025']], digits)}, {round(w[['q975']], digits)}]")
		}
	}
	d <- dim(x$snr_per_draw_per_t)
	cli::cli_text("")
	cli::cli_text("Per-draw, per-transition values in {.code $snr_per_draw_per_t} ({d[1]} x {d[2]}).")
	invisible(x)
}

#' Identifiability diagnostics for a symmetric DBN fit
#'
#' @description Computes three numerical diagnostics that bear on the
#'   identification analysis of the time-varying bilinear AR model:
#'
#'   1. \strong{Lagged-state design rank.} The rank of the matrix whose
#'      rows are the posterior-mean off-diagonal entries of $\\Theta_{t-1}$
#'      for $t = 2, \\ldots, T$. Maximum is $\\min(T-1, n(n-1)/2)$.
#'      Compared to the input-space dimension $n(n-1)/2$, this reports
#'      whether classical operator identification (the population-kernel
#'      / blockwise-static version) is empirically available.
#'
#'   2. \strong{Symmetric Jacobian rank.} The numerical rank of the
#'      stacked Jacobian $H \\mapsto P_\\mathrm{off}(H \\Theta_{t-1} A^\\top +
#'      A \\Theta_{t-1} H^\\top)$ over symmetric perturbations $H$ and
#'      lagged states $t = 2, \\ldots, T$, evaluated at the posterior-mean
#'      $A$. Full rank ($n(n+1)/2$) indicates local identifiability of
#'      $A_t$ at the posterior mean.
#'
#'   3. \strong{Augmented Jacobian rank.} Same as 2 but with diagonal
#'      entries of the Jacobian output included (diagonal-penalty case).
#'      Tests whether the diagonal penalty closes a local rank deficit.
#'
#' @param fit A `dbn` object from `dbn_dynamic(..., symmetric = TRUE)`.
#' @param relative_tol Relative singular-value tolerance for the
#'   effective-rank diagnostic (default: 1e-6).
#' @return A list with components `design`, `jacobian_sym`,
#'   `jacobian_aug`, each containing dimensions, numerical and effective
#'   ranks, singular values, and condition numbers. Class
#'   `"dbn_identifiability_diag"`.
#' @details Local identifiability at the posterior mean is the Bayesian
#'   counterpart of the classical Jacobian rank condition for symmetric
#'   operator identification (Proposition 3 of the PA methods paper
#'   supplement). When the design rank is below the input-space
#'   dimension, classical per-period identification is unavailable and
#'   the posterior on $W_t$ is the natural estimand under the random-
#'   walk-smoothed model; full Jacobian rank indicates the symmetric
#'   posterior at the fitted $A$ is locally well-posed.
#' @seealso [dbn_compute_snr()], [dbn_coupling_rank_probs()],
#'   [dbn_dynamic()]
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 8, p = 1, time = 12, symmetric = TRUE, seed = 1)
#' Y <- sim$Z[, , 1, , drop = FALSE]
#' fit <- dbn_dynamic(Y, family = "gaussian", symmetric = TRUE,
#'   nscan = 200, burn = 100, verbose = FALSE)
#' diag_out <- dbn_identifiability_diagnostics(fit)
#' diag_out$design$numerical_rank
#' diag_out$jacobian_sym$smallest_sv
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_identifiability_diagnostics <- function(fit, relative_tol = 1e-6) {
	if (!inherits(fit, "dbn")) {
		stop("fit must be a dbn object returned by dbn_dynamic().")
	}
	if (is.null(fit$A) || is.null(fit$Theta)) {
		stop("fit must contain $A and $Theta. Refit with default time_thin ",
			"and store_z settings to retain theta draws.")
	}

	# detect symmetric specification
	is_sym <- isTRUE(fit$dims$is_symmetric) ||
		isTRUE(identical(fit$A[[1]], fit$B[[1]]))
	if (!is_sym) {
		cli::cli_abort(c(
			"{.fun dbn_identifiability_diagnostics} requires a symmetric fit.",
			"x" = "This fit is asymmetric; the Jacobian rank condition is defined for the symmetric (B = A) specification.",
			"i" = "Refit with {.code symmetric = TRUE} for a valid identifiability diagnostic."
		))
	}

	n_draws <- length(fit$A)
	Theta_arr <- fit$Theta
	if (length(dim(Theta_arr)) == 5 && dim(Theta_arr)[3] == 1) {
		Theta_post <- apply(Theta_arr, c(1, 2, 4), mean)
	} else {
		stop("dbn_identifiability_diagnostics() requires a single-relation (p = 1) ",
			"fit; got Theta with dim ",
			paste(dim(Theta_arr), collapse = "x"), ".")
	}
	n <- dim(Theta_post)[1]
	Tt <- dim(Theta_post)[3]

	# posterior mean A
	A_arr <- array(0, c(n, n, Tt))
	for (m in seq_len(n_draws)) A_arr <- A_arr + fit$A[[m]]
	A_arr <- A_arr / n_draws
	A_mean <- apply(A_arr, c(1, 2), mean)
	A_mean <- 0.5 * (A_mean + t(A_mean))

	off_mask <- upper.tri(matrix(0, n, n))
	m_up <- sum(off_mask)
	upper_vec <- function(M) M[off_mask]
	upper_diag_vec <- function(M) c(M[off_mask], diag(M))

	# diagnostic 1: lagged-state design rank
	X_design <- matrix(0, Tt - 1, m_up)
	for (t in 2:Tt) X_design[t - 1, ] <- upper_vec(Theta_post[, , t - 1])
	sv_design <- svd(X_design, nu = 0, nv = 0)$d
	num_tol_design <- max(dim(X_design)) * .Machine$double.eps * max(sv_design)
	rel_tol_design <- max(sv_design) * relative_tol
	design <- list(
		dim = dim(X_design),
		input_space_dim = m_up,
		max_possible_rank = min(Tt - 1, m_up),
		numerical_rank = sum(sv_design > num_tol_design),
		effective_rank = sum(sv_design > rel_tol_design),
		singular_values = sv_design,
		smallest_sv = min(sv_design),
		condition = sv_design[1] / max(sv_design[length(sv_design)], 1e-12)
	)

	# diagnostics 2 and 3: symmetric and augmented Jacobian
	sym_pairs <- which(upper.tri(matrix(0, n, n), diag = TRUE), arr.ind = TRUE)
	n_sym <- nrow(sym_pairs)
	jacobian_at_t <- function(A_t, X, augment = FALSE) {
		out_len <- if (augment) m_up + n else m_up
		J <- matrix(0, out_len, n_sym)
		A_t_T <- t(A_t)
		for (b in seq_len(n_sym)) {
			i <- sym_pairs[b, 1]; j <- sym_pairs[b, 2]
			H <- matrix(0, n, n); H[i, j] <- 1
			if (i != j) H[j, i] <- 1
			M <- H %*% X %*% A_t_T + A_t %*% X %*% t(H)
			J[, b] <- if (augment) upper_diag_vec(M) else upper_vec(M)
		}
		J
	}

	# stacked symmetric Jacobian
	J_stack <- matrix(0, (Tt - 1) * m_up, n_sym)
	for (t in 2:Tt) {
		J_stack[((t - 2) * m_up + 1):((t - 1) * m_up), ] <-
			jacobian_at_t(A_mean, Theta_post[, , t - 1], augment = FALSE)
	}
	sv_jac <- svd(J_stack, nu = 0, nv = 0)$d
	jacobian_sym <- list(
		dim = dim(J_stack),
		max_possible_rank = n_sym,
		numerical_rank = sum(sv_jac >
			max(dim(J_stack)) * .Machine$double.eps * max(sv_jac)),
		effective_rank = sum(sv_jac > max(sv_jac) * relative_tol),
		singular_values = sv_jac,
		smallest_sv = min(sv_jac),
		condition = sv_jac[1] / max(sv_jac[length(sv_jac)], 1e-12)
	)

	# stacked augmented Jacobian (diagonal penalty case)
	m_full <- m_up + n
	J_stack_aug <- matrix(0, (Tt - 1) * m_full, n_sym)
	for (t in 2:Tt) {
		J_stack_aug[((t - 2) * m_full + 1):((t - 1) * m_full), ] <-
			jacobian_at_t(A_mean, Theta_post[, , t - 1], augment = TRUE)
	}
	sv_jac_aug <- svd(J_stack_aug, nu = 0, nv = 0)$d
	jacobian_aug <- list(
		dim = dim(J_stack_aug),
		max_possible_rank = n_sym,
		numerical_rank = sum(sv_jac_aug >
			max(dim(J_stack_aug)) * .Machine$double.eps * max(sv_jac_aug)),
		effective_rank = sum(sv_jac_aug > max(sv_jac_aug) * relative_tol),
		singular_values = sv_jac_aug,
		smallest_sv = min(sv_jac_aug),
		condition = sv_jac_aug[1] / max(sv_jac_aug[length(sv_jac_aug)], 1e-12)
	)

	out <- list(
		n = n, Tt = Tt, n_draws = n_draws,
		design = design,
		jacobian_sym = jacobian_sym,
		jacobian_aug = jacobian_aug
	)
	class(out) <- c("dbn_identifiability_diag", "list")
	out
}

#' Print method for DBN identifiability diagnostics
#'
#' @description Prints a compact verdict from [dbn_identifiability_diagnostics()]
#'   instead of dumping the full nested list of singular-value blocks.
#' @param x A `dbn_identifiability_diag` object.
#' @param digits Number of digits (default 3).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method print dbn_identifiability_diag
print.dbn_identifiability_diag <- function(x, digits = 3, ...) {
	cli::cli_h2("DBN identifiability diagnostics")
	cli::cli_text("Network: {x$n} actors, {x$Tt} time points, {x$n_draws} draw{?s}.")
	report <- function(label, j) {
		full <- j$numerical_rank >= j$max_possible_rank
		cli::cli_h3(label)
		cli::cli_text("  numerical rank {j$numerical_rank} / {j$max_possible_rank}{if (full) ' (full rank)' else ' (rank-deficient)'}")
		cli::cli_text("  smallest singular value {format(signif(j$smallest_sv, digits))}, condition {format(signif(j$condition, digits))}")
		if (full) {
			cli::cli_alert_success("Jacobian has full rank: the symmetric operator is locally identified.")
		} else {
			cli::cli_alert_warning("Jacobian is rank-deficient: the operator is not locally identified from these data.")
		}
	}
	if (!is.null(x$design)) {
		d <- x$design
		cli::cli_h3("Design rank (classical per-period identification)")
		cli::cli_text("  numerical rank {d$numerical_rank} / {d$input_space_dim} input dimensions")
		if (d$numerical_rank >= d$input_space_dim) {
			cli::cli_alert_success("The lagged-state design has full column rank.")
		} else {
			cli::cli_alert_warning("The lagged-state design is rank-deficient ({d$numerical_rank} of {d$input_space_dim}): classical per-period identification is unavailable, and the random-walk prior is doing the identifying work.")
		}
	}
	if (!is.null(x$jacobian_sym)) report("Symmetric Jacobian", x$jacobian_sym)
	if (!is.null(x$jacobian_aug)) report("Diagonal-augmented Jacobian", x$jacobian_aug)
	cli::cli_text("")
	cli::cli_text("Full singular-value blocks in {.code $design}, {.code $jacobian_sym}, and {.code $jacobian_aug}.")
	invisible(x)
}

#' Extract the time-varying operator and assess its stability
#'
#' @description
#' Returns the posterior-mean time-varying lag operator
#' \eqn{W_t = A_t B_t^\top} (for a symmetric fit, \eqn{B_t = A_t}) together
#' with a per-period stability summary. The latent state evolves as
#' \eqn{\mathrm{vec}(\Theta_t) = (B_t \otimes A_t)\,\mathrm{vec}(\Theta_{t-1})},
#' which is locally stable when the gain \eqn{\rho(A_t)\,\rho(B_t)} is below 1.
#'
#' Because the operator drifts as a random walk, individual periods and
#' individual posterior draws can be locally explosive even when the
#' posterior-mean operator looks contractive (averaging operators whose
#' eigenvalues point in different directions cancels them). This function
#' therefore reports the gain per draw, so an average is not mistaken for
#' every draw.
#'
#' @param fit A fitted dynamic `dbn` object.
#' @return A list of class `dbn_operator` with `W` (an `[n, n, T]`
#'   posterior-mean operator array) and `stability` (a data frame, one row per
#'   time point, with the posterior mean and 95% interval of the gain
#'   \eqn{\rho(A_t)\rho(B_t)} and `frac_explosive`, the fraction of draws whose
#'   gain exceeds 1).
#' @seealso [dbn_coupling_rank_probs()], [compute_irf()]
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 8, time = 12, seed = 1)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'            nscan = 200, burn = 100, verbose = FALSE)
#' op <- dbn_operator(fit)
#' op$stability
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_operator <- function(fit) {
	if (!inherits(fit, "dbn")) cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	if (fit$model != "dynamic" || is.null(fit$A) || !is.list(fit$A) ||
	    is.null(fit$B) || !is.list(fit$B)) {
		cli::cli_abort(c(
			"{.fun dbn_operator} requires a dynamic fit with time-varying A/B draws.",
			"i" = "Got model type {.val {fit$model}}."
		))
	}
	# the interaction summary W_t = A_t B_t^T is well-defined only for
	# square (unipartite) networks; for bipartite (n_row != n_col), A and
	# B have different sizes and the product is non-conformable.
	if (isTRUE(fit$dims$is_bipartite)) {
		cli::cli_abort(c(
			"{.fun dbn_operator} requires a unipartite (square) fit; got a bipartite fit.",
			"x" = "{.code A_t} is {.val {fit$dims$n_row}} x {.val {fit$dims$n_row}} but {.code B_t} is {.val {fit$dims$n_col}} x {.val {fit$dims$n_col}}.",
			"i" = "The interaction summary {.code W_t = A_t B_t^T} is defined only for square networks. Inspect {.code fit$A} and {.code fit$B} separately, or use {.fun dbn_leverage}, {.fun compute_irf} for bipartite analysis."
		))
	}
	n_draws <- length(fit$A)
	d <- dim(fit$A[[1]])
	n <- d[1]
	Tt <- d[3]
	# signed mean, sign-flip-invariant absolute mean, and per-draw frobenius
	W_sum <- array(0, dim = c(n, n, Tt))
	W_abs_sum <- array(0, dim = c(n, n, Tt))
	W_frob <- matrix(NA_real_, n_draws, Tt)
	gain <- matrix(NA_real_, n_draws, Tt)
	for (m in seq_len(n_draws)) {
		Am <- fit$A[[m]]
		Bm <- fit$B[[m]]
		for (t in seq_len(Tt)) {
			At <- Am[, , t]
			Bt <- Bm[, , t]
			Wt <- At %*% t(Bt)
			W_sum[, , t] <- W_sum[, , t] + Wt
			W_abs_sum[, , t] <- W_abs_sum[, , t] + abs(Wt)
			W_frob[m, t] <- sqrt(sum(Wt^2))
			ra <- max(Mod(eigen(At, only.values = TRUE)$values))
			rb <- max(Mod(eigen(Bt, only.values = TRUE)$values))
			gain[m, t] <- ra * rb
		}
	}
	# flag t = 1 as the anchor: in the asymmetric MCMC path, A_1 / B_1 are
	# pinned to identity (src/dynamic.cpp:436, 505, 737, 793) rather than
	# sampled. The per-iteration scale rebalancer (R/models.R:1698-1706) then
	# scales the identity by some constant c per draw, so A_1 looks like c*I
	# rather than exactly I. Detect "scaled identity with zero off-diagonal
	# variance across draws" -- that's the anchor signature, regardless of
	# rebalance scale.
	is_anchor <- logical(Tt)
	if (Tt >= 1L && n_draws >= 1L) {
		# check that A_1 across all draws is diagonal and B_1 likewise
		A1_diag_only <- TRUE
		B1_diag_only <- TRUE
		for (m in seq_len(min(n_draws, 5L))) {
			A1 <- fit$A[[m]][, , 1]
			B1 <- fit$B[[m]][, , 1]
			# off-diagonal mass relative to total mass
			off_A <- sum(A1[row(A1) != col(A1)]^2)
			tot_A <- sum(A1^2)
			if (tot_A > 0 && off_A / tot_A > 1e-8) { A1_diag_only <- FALSE; break }
			off_B <- sum(B1[row(B1) != col(B1)]^2)
			tot_B <- sum(B1^2)
			if (tot_B > 0 && off_B / tot_B > 1e-8) { B1_diag_only <- FALSE; break }
		}
		if (A1_diag_only && B1_diag_only) is_anchor[1] <- TRUE
	}
	stability <- data.frame(
		t = seq_len(Tt),
		gain_mean = colMeans(gain, na.rm = TRUE),
		gain_q025 = apply(gain, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
		gain_q975 = apply(gain, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
		frac_explosive = colMeans(gain > 1, na.rm = TRUE),
		is_anchor = is_anchor,
		row.names = NULL
	)
	W_frob_summary <- data.frame(
		t        = seq_len(Tt),
		frob_mean = colMeans(W_frob, na.rm = TRUE),
		frob_q025 = apply(W_frob, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
		frob_q975 = apply(W_frob, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
		row.names = NULL
	)
	out <- list(
		W           = W_sum / n_draws,
		W_abs_mean  = W_abs_sum / n_draws,
		W_frob      = W_frob_summary,
		stability   = stability,
		n_draws     = n_draws,
		has_t1_anchor = any(is_anchor)
	)
	class(out) <- c("dbn_operator", "list")
	out
}

#' Print method for the DBN time-varying operator
#'
#' @param x A `dbn_operator` object.
#' @param digits Number of digits (default 3).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method print dbn_operator
print.dbn_operator <- function(x, digits = 3, ...) {
	cli::cli_h2("DBN time-varying interaction summary W_t = A_t B_t^T")
	d <- dim(x$W)
	cli::cli_text("Posterior-mean interaction matrix in {.code $W}: {d[1]} x {d[2]} x {d[3]} (actors x actors x time).")
	cli::cli_text("{.strong Note:} {.code $W} is the {.emph interaction summary} A_t B_t^T -- not the state-transition operator. The relevant operator for stability is the Kronecker product B_t {.strong tensor} A_t whose spectral radius equals rho(A_t) * rho(B_t); the {.code $stability$gain_mean} column below reports this gain, while {.code eigen($W)} does not.")
	cli::cli_text("{.strong Sign cancellation:} the rank-r factorization admits per-latent sign flips across draws, so {.code $W} (signed mean) can be cancelled toward zero entry-by-entry while individual draws are far from zero. Use {.code $W_abs_mean} for entry-wise coupling and {.code $W_frob} for per-time scale.")
	s <- x$stability
	if (isTRUE(x$has_t1_anchor)) {
		cli::cli_alert_info("t = 1 is a structural identity anchor (A_1 = B_1 = I); the gain = 1.000 entry there is NOT an estimate. Summaries below exclude t = 1.")
		s_summary <- s[!s$is_anchor, , drop = FALSE]
	} else {
		s_summary <- s
	}
	cli::cli_text("Operator gain rho(A_t) rho(B_t) over time: mean {round(mean(s_summary$gain_mean), digits)}, range [{round(min(s_summary$gain_mean), digits)}, {round(max(s_summary$gain_mean), digits)}].")
	n_expl <- sum(s_summary$frac_explosive > 0.5)
	if (n_expl > 0) {
		cli::cli_alert_warning("Most posterior draws are locally explosive (gain > 1) in {n_expl} of {nrow(s)} period{?s}; impulse responses and forecasts will not decay there.")
	} else if (mean(s$frac_explosive) > 0.05) {
		cli::cli_alert_info("Some draws are locally explosive; see {.code $stability$frac_explosive}.")
	} else {
		cli::cli_alert_success("The operator is locally stable across draws and time.")
	}
	invisible(x)
}

#' Plot a dbn_rank_probs object
#'
#' @description Renders the per-actor coupling-rank probabilities as a
#'   point-range or tile plot, depending on `kind`. The default
#'   summarises the per-actor coupling CI (`Cbar` posterior with
#'   `q025`/`q975`) over a chosen window.
#'
#' @param x A `dbn_rank_probs` object from [dbn_coupling_rank_probs()].
#' @param kind Character. One of `"coupling"` (point-range per actor) or
#'   `"top1"` (per-actor probability of being top-1).
#' @param ... Unused.
#' @return Invisibly, the ggplot object (or `NULL` if ggplot2 not loaded).
#' @author Tosin Salau and Shahryar Minhas
#' @export
plot.dbn_rank_probs <- function(x, kind = c("coupling", "top1"), ...) {
	if (!inherits(x, "dbn_rank_probs")) {
		cli::cli_abort(c(
			"{.fn plot.dbn_rank_probs} requires a {.cls dbn_rank_probs} object.",
			"x" = "Got {.cls {class(x)[1]}}.",
			"i" = "Produce one with {.fn dbn_coupling_rank_probs} on a symmetric fit."
		))
	}
	kind <- match.arg(kind)
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_warn("Install {.pkg ggplot2} for the dbn_rank_probs plot.")
		return(invisible(NULL))
	}
	full <- x$full_sample %||% x[[1]]
	if (is.null(full)) {
		cli::cli_abort("Could not locate a sample slot on the {.cls dbn_rank_probs} object.")
	}
	if (kind == "coupling") {
		df <- data.frame(
			actor = full$actor %||% paste0("a", seq_along(full$post_mean)),
			mean = full$post_mean,
			lo = full$post_q025 %||% full$q025,
			hi = full$post_q975 %||% full$q975
		)
		df <- df[order(df$mean), , drop = FALSE]
		df$actor <- factor(df$actor, levels = df$actor)
		p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$actor, y = .data$mean)) +
			ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$lo, ymax = .data$hi)) +
			ggplot2::coord_flip() +
			ggplot2::labs(x = NULL, y = "Posterior coupling (mean and 95% CI)") +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.border = ggplot2::element_blank(),
				axis.ticks = ggplot2::element_blank(),
				legend.position = "top")
	} else {
		p_top <- full$p_topK
		p_top1 <- if (is.list(p_top)) p_top$p_top1 else p_top[, 1]
		df <- data.frame(actor = names(p_top1) %||% paste0("a", seq_along(p_top1)),
			p_top1 = as.numeric(p_top1))
		df <- df[order(df$p_top1), , drop = FALSE]
		df$actor <- factor(df$actor, levels = df$actor)
		p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$actor, y = .data$p_top1)) +
			ggplot2::geom_col() +
			ggplot2::coord_flip() +
			ggplot2::labs(x = NULL, y = "P(actor is top-1)") +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.border = ggplot2::element_blank(),
				axis.ticks = ggplot2::element_blank(),
				legend.position = "top")
	}
	print(p)
	invisible(p)
}

#' Summary method for dbn_rank_probs
#'
#' Returns a data frame of per-actor coupling means / 95% intervals /
#' top-1 / top-3 probabilities, sorted by posterior mean coupling
#' (descending). Useful for compact reporting and downstream tables.
#' @param object A `dbn_rank_probs` object.
#' @param ... Unused.
#' @return A data frame, one row per actor.
#' @author Tosin Salau and Shahryar Minhas
#' @export
summary.dbn_rank_probs <- function(object, ...) {
	if (!inherits(object, "dbn_rank_probs")) {
		cli::cli_abort(c(
			"{.fn summary.dbn_rank_probs} requires a {.cls dbn_rank_probs} object.",
			"x" = "Got {.cls {class(object)[1]}}.",
			"i" = "Produce one with {.fn dbn_coupling_rank_probs} on a symmetric fit."
		))
	}
	full <- object$full_sample %||% object[[1]]
	if (is.null(full)) {
		cli::cli_abort("Could not locate a sample slot on the {.cls dbn_rank_probs} object.")
	}
	actor <- full$actor %||% paste0("a", seq_along(full$post_mean))
	p_top <- full$p_topK
	p_top1 <- if (is.list(p_top)) p_top$p_top1 else if (is.matrix(p_top)) p_top[, 1] else rep(NA_real_, length(actor))
	p_top3 <- if (is.list(p_top)) p_top$p_top3 else if (is.matrix(p_top) && ncol(p_top) >= 2L) p_top[, 2] else rep(NA_real_, length(actor))
	out <- data.frame(
		actor = actor,
		mean = full$post_mean,
		lower = full$post_q025 %||% full$q025,
		upper = full$post_q975 %||% full$q975,
		p_top1 = as.numeric(p_top1),
		p_top3 = as.numeric(p_top3),
		stringsAsFactors = FALSE
	)
	out[order(-out$mean), , drop = FALSE]
}

#' autoplot alias for plot.dbn_rank_probs
#'
#' Lets `ggplot2::autoplot(x)` dispatch to [plot.dbn_rank_probs()]
#' when ggplot2 is loaded. Registered conditionally on `ggplot2::autoplot`
#' at package load so ggplot2 stays in `Suggests`.
#' @param object A `dbn_rank_probs` object.
#' @param ... Forwarded to [plot.dbn_rank_probs()].
#' @return A ggplot object (invisibly, via the underlying plot method).
#' @keywords internal
autoplot.dbn_rank_probs <- function(object, ...) {
	if (!inherits(object, "dbn_rank_probs")) {
		cli::cli_abort(c(
			"{.fn autoplot.dbn_rank_probs} requires a {.cls dbn_rank_probs} object.",
			"x" = "Got {.cls {class(object)[1]}}.",
			"i" = "Produce one with {.fn dbn_coupling_rank_probs} on a symmetric fit."
		))
	}
	plot.dbn_rank_probs(object, ...)
}
