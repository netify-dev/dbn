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
#' @export
dbn_compute_snr <- function(fit, windows = NULL) {
	if (!inherits(fit, "dbn")) {
		stop("fit must be a dbn object returned by dbn_dynamic().")
	}
	if (is.null(fit$A) || is.null(fit$Theta)) {
		stop("fit must contain $A and $Theta. Refit with store_z = TRUE or use ",
			"a fit that retains theta draws.")
	}

	n_draws <- length(fit$A)
	Theta_arr <- fit$Theta
	# Theta is n_row x n_col x p x Tt x n_draws; collapse p (must be 1)
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
#' @export
dbn_coupling_rank_probs <- function(fit, windows = NULL, K_grid = c(1, 3, 5),
	actor_names = NULL) {
	if (!inherits(fit, "dbn")) {
		stop("fit must be a dbn object returned by dbn_dynamic().")
	}
	if (is.null(fit$A)) {
		stop("fit must contain $A draws.")
	}
	n_draws <- length(fit$A)
	A_first <- fit$A[[1]]
	n <- dim(A_first)[1]
	Tt <- dim(A_first)[3]
	if (is.null(actor_names)) actor_names <- paste0("actor_", seq_len(n))

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
#' sim <- simulate_dynamic_dbn(n = 8, time = 12, symmetric = TRUE, seed = 1)
#' Y <- sim$Z[, , 1, , drop = FALSE]
#' fit <- dbn_dynamic(Y, family = "gaussian", symmetric = TRUE,
#'   nscan = 200, burn = 100, verbose = FALSE)
#' diag_out <- dbn_identifiability_diagnostics(fit)
#' diag_out$design$numerical_rank
#' diag_out$jacobian_sym$smallest_sv
#' }
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
		warning("dbn_identifiability_diagnostics() is designed for the ",
			"symmetric specification; results on asymmetric fits may not ",
			"reflect the intended diagnostic quantities.")
	}

	n_draws <- length(fit$A)
	Theta_arr <- fit$Theta
	if (length(dim(Theta_arr)) == 5 && dim(Theta_arr)[3] == 1) {
		Theta_post <- apply(Theta_arr, c(1, 2, 4), mean)
	} else {
		stop("dbn_identifiability_diagnostics() requires a unipartite (p = 1) ",
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

	# Diagnostic 1: lagged-state design rank
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

	# Diagnostics 2 and 3: symmetric and augmented Jacobian
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
