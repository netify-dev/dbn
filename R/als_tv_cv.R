####
# lambda-path with continuation + rolling-origin cv.
#   lambda = "path" -> fit a full grid; return all
#   lambda = "cv"   -> fit grid, pick by rolling-origin cv with 1-se rule
#   lambda = "eb"   -> single fit at empirical-bayes plug-in lambda
# the continuation strategy fits from large lambda to small, warm-starting
# each fit from the previous. large-lambda fits sit near the static
# solution and relaxing lambda gradually walks the trajectory through
# identifiable basins.
####

#' Empirical-Bayes lambda pilot with measurement-noise bias correction
#'
#' Fit static ALS on non-overlapping windows. Estimate sigma2 from one-step
#' residuals and raw tau_A2 / tau_B2 from window-to-window operator
#' differences. Subtract the expected sampling-variance contribution of the
#' windowed estimator so the corrected tau2 reflects true operator drift
#' instead of pilot-fit noise. Return both raw and corrected quantities plus
#' separate lambda_A, lambda_B, and a precision-weighted lambda_common.
#'
#' @keywords internal
#' @noRd
.tv_als_pilot_lambda <- function(Y, family, symmetric, n_row, n_col, p, Tt,
                                  window = NULL) {
	# longer windows give a more reliable sigma2 estimate; cap at T/2 so we still
	# get at least 2 non-overlapping windows
	if (is.null(window)) window <- max(4L, min(floor(Tt / 2), ceiling(Tt / 3)))
	floor_lambda <- 0.01
	floor_tau2 <- 1e-8

	fallback <- function() list(
		lambda_pilot = 1.0, lambda_A = 1.0, lambda_B = 1.0,
		lambda_common = 1.0, lambda_raw = 1.0,
		sigma2 = 0.5,
		tau_A2_raw = 0.5, tau_B2_raw = 0.5,
		tau_A2_corrected = 0.5, tau_B2_corrected = 0.5,
		tau_A2_meas_trace = 0, tau_B2_meas_trace = 0,
		correction_factor = 1.0,
		lambda_boundary_hit = FALSE
	)

	if (Tt < window) return(fallback())

	# non-overlapping windows so successive measurement errors are independent
	starts <- seq(1L, Tt - window + 1L, by = window)
	A_wins <- list()
	B_wins <- list()
	XtX_inv_traces <- numeric(0)
	BBt_inv_traces <- numeric(0)
	AtA_inv_traces <- numeric(0)
	win_failures <- character(0)

	for (start in starts) {
		end <- start + window - 1L
		Y_win <- Y[, , , start:end, drop = FALSE]
		fit_win <- tryCatch(
			dbn_als(Y_win, family = family, symmetric = symmetric,
				bootstrap = 0, verbose = FALSE),
			error = function(e) {
				win_failures[length(win_failures) + 1L] <<- conditionMessage(e)
				NULL
			}
		)
		if (is.null(fit_win)) next
		A_k <- fit_win$A[[1]][, , 1]
		B_k <- fit_win$B[[1]][, , 1]
		A_wins[[length(A_wins) + 1L]] <- A_k
		B_wins[[length(B_wins) + 1L]] <- B_k
		# window-specific design statistics used for hessian-based V_W trace
		XX <- matrix(0, n_row, n_row)
		YY <- matrix(0, n_col, n_col)
		for (r in 1:p) for (t in start:(end - 1L)) {
			Z <- Y[, , r, t]
			Z[is.na(Z)] <- 0
			XX <- XX + Z %*% t(Z)
			YY <- YY + t(Z) %*% Z
		}
		XtX_inv_traces <- c(XtX_inv_traces, .safe_trace_inv(XX))
		BBt_inv_traces <- c(BBt_inv_traces, .safe_trace_inv(B_k %*% t(B_k)))
		AtA_inv_traces <- c(AtA_inv_traces, .safe_trace_inv(t(A_k) %*% A_k))
	}

	if (length(A_wins) < 2L) {
		if (length(win_failures) > 0L) {
			cli::cli_warn(c(
				"Empirical-Bayes pilot fit too few windowed static models.",
				"x" = "First failure: {substr(win_failures[1], 1, 100)}",
				"i" = "Falling back to lambda_pilot = 1.0; pass a numeric {.arg lambda} for control."
			))
		}
		return(fallback())
	}

	# sigma2 from one-step prediction residuals
	resid_sq <- 0
	n_resid <- 0
	for (k in seq_along(A_wins)) {
		start <- (k - 1L) * window + 1L
		end <- min(Tt, start + window - 1L)
		for (r in 1:p) for (t in (start + 1L):end) {
			pred <- A_wins[[k]] %*% Y[, , r, t - 1L] %*% t(B_wins[[k]])
			obs <- Y[, , r, t]
			obs[is.na(obs)] <- 0
			pred[is.na(pred)] <- 0
			resid_sq <- resid_sq + sum((obs - pred)^2)
			n_resid <- n_resid + length(obs)
		}
	}
	sigma2_hat <- max(resid_sq / max(n_resid, 1), 1e-12)

	# raw tau2 estimates from successive window operators
	n_diffs <- length(A_wins) - 1L
	dA_sq <- 0
	dB_sq <- 0
	for (k in 2:length(A_wins)) {
		dA_sq <- dA_sq + sum((A_wins[[k]] - A_wins[[k - 1L]])^2)
		dB_sq <- dB_sq + sum((B_wins[[k]] - B_wins[[k - 1L]])^2)
	}
	tau_A2_raw <- max(dA_sq / (n_diffs * n_row * n_row), floor_tau2)
	tau_B2_raw <- max(dB_sq / (n_diffs * n_col * n_col), floor_tau2)

	# measurement variance of windowed vec(A) = sigma2 * tr((XX')^{-1}) * tr((BB')^{-1})
	mean_trace_VA <- sigma2_hat * mean(XtX_inv_traces) * mean(BBt_inv_traces)
	mean_trace_VB <- sigma2_hat * mean(XtX_inv_traces) * mean(AtA_inv_traces)
	meas_A <- 2 * mean_trace_VA / (n_row * n_row)
	meas_B <- 2 * mean_trace_VB / (n_col * n_col)

	# safety cap: the correction nudges lambda, not blows it up; cap so
	# tau2_corrected >= tau2_raw / 3, lambda bumped by at most factor 3
	meas_A <- min(meas_A, (2 / 3) * tau_A2_raw)
	meas_B <- min(meas_B, (2 / 3) * tau_B2_raw)

	tau_A2_corrected <- max(tau_A2_raw - meas_A, floor_tau2)
	tau_B2_corrected <- max(tau_B2_raw - meas_B, floor_tau2)

	# symmetric: B = A; force tau_B = tau_A
	if (isTRUE(symmetric)) {
		tau_A2_corrected <- (tau_A2_corrected + tau_B2_corrected) / 2
		tau_B2_corrected <- tau_A2_corrected
		tau_A2_raw <- (tau_A2_raw + tau_B2_raw) / 2
		tau_B2_raw <- tau_A2_raw
	}

	lambda_A <- max(sigma2_hat / tau_A2_corrected, floor_lambda)
	lambda_B <- max(sigma2_hat / tau_B2_corrected, floor_lambda)

	# precision-weighted common lambda
	p_A <- n_row * n_row
	p_B <- n_col * n_col
	lambda_common <- (p_A * lambda_A + p_B * lambda_B) / (p_A + p_B)

	# raw lambda for reporting (uses raw tau2)
	lambda_raw <- max(sigma2_hat / max(tau_A2_raw, floor_tau2), floor_lambda)

	correction_factor <- lambda_common / lambda_raw

	# warn if the two side-specific lambdas differ by more than factor of e
	if (max(lambda_A, lambda_B) / min(lambda_A, lambda_B) > exp(1)) {
		cli::cli_alert_info(
			"EB pilot suggests asymmetric smoothing: lambda_A = {sprintf('%.3g', lambda_A)}, lambda_B = {sprintf('%.3g', lambda_B)}. Using precision-weighted common lambda = {sprintf('%.3g', lambda_common)}."
		)
	}

	list(
		lambda_pilot = lambda_common,
		lambda_A = lambda_A,
		lambda_B = lambda_B,
		lambda_common = lambda_common,
		lambda_raw = lambda_raw,
		sigma2 = sigma2_hat,
		tau_A2_raw = tau_A2_raw,
		tau_B2_raw = tau_B2_raw,
		tau_A2_corrected = tau_A2_corrected,
		tau_B2_corrected = tau_B2_corrected,
		tau_A2_meas_trace = meas_A,
		tau_B2_meas_trace = meas_B,
		correction_factor = correction_factor,
		lambda_boundary_hit = lambda_common <= floor_lambda
	)
}

#' Trace of inverse with regularization for ill-conditioned matrices
#'
#' If eigen() fails (rare, e.g. non-finite entries) we substitute the
#' largest plausible trace (`n / 1e-6`), but only after warning once
#' per session so a corrupted upstream M doesn't silently bias the EB
#' pilot toward over-shrinkage.
#' @keywords internal
#' @noRd
.safe_trace_inv <- function(M) {
	n <- nrow(M)
	if (is.null(n) || n == 0L) return(0)
	ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
		error = function(e) NULL)
	if (is.null(ev) || any(!is.finite(ev))) {
		if (!isTRUE(getOption("dbn.safe_trace_inv_warned", FALSE))) {
			cli::cli_warn(c(
				"{.fun .safe_trace_inv}: eigendecomposition failed or produced non-finite eigenvalues.",
				"i" = "Substituting a large-trace sentinel; EB pilot estimate may be biased toward over-shrinkage.",
				"i" = "Inspect upstream {.code M} for non-finite entries. Suppress via {.code options(dbn.safe_trace_inv_warned = TRUE)}."
			))
			options(dbn.safe_trace_inv_warned = TRUE)
		}
		return(n / 1e-6)
	}
	floor_ev <- max(1e-8, .Machine$double.eps * max(abs(ev), 1))
	sum(1 / pmax(ev, floor_ev))
}

#' Fit the lambda-path via continuation
#'
#' @param Y data array
#' @param family,symmetric,mu,eta,gauge passed through
#' @param lambda_grid numeric vector of lambda values (sorted descending recommended)
#' @param tv_max_iter,tv_tol_obj,tv_tol_par convergence settings
#' @return list with element `fits` (list of TV-ALS fits per lambda), `lambda_grid`
#' @keywords internal
#' @noRd
.tv_als_lambda_path <- function(Y, family, symmetric, lambda_grid, mu, eta,
                                  gauge, tv_max_iter, tv_tol_obj, tv_tol_par,
                                  init_static_fit = NULL, verbose = FALSE) {
	dims <- dim(Y); n_row <- dims[1]; n_col <- dims[2]; p <- dims[3]; Tt <- dims[4]

	# compute static warm-start once
	if (is.null(init_static_fit)) {
		init_static_fit <- dbn_als(Y, family = family, symmetric = symmetric,
			bootstrap = 0, verbose = FALSE)
	}
	A_static <- init_static_fit$A[[1]][, , 1]
	B_static <- init_static_fit$B[[1]][, , 1]
	M_static <- array(NA, dim = c(n_row, n_col, p))
	for (r in 1:p) M_static[, , r] <- init_static_fit$M[, , r, 1]

	# pre-compute mu from data scale if not given
	if (is.null(mu)) {
		Z_for_M <- Y; Z_for_M[is.na(Z_for_M)] <- 0
		# h_bar via static B_static
		h_bar <- mean(sapply(1:Tt, function(t) {
			M_t <- Z_for_M[, , 1, t] %*% t(B_static)
			sum(M_t^2) / n_row
		}))
		mu <- eta * h_bar
	}

	# sort lambda_grid descending so continuation goes large -> small
	lambda_grid <- sort(lambda_grid, decreasing = TRUE)

	# initial: broadcast static across t
	A_curr <- array(rep(c(A_static), Tt), dim = c(n_row, n_row, Tt))
	B_curr <- array(rep(c(B_static), Tt), dim = c(n_col, n_col, Tt))

	fits <- vector("list", length(lambda_grid))
	for (i in seq_along(lambda_grid)) {
		lam <- lambda_grid[i]
		if (verbose) cli::cli_inform("path iter {i}/{length(lambda_grid)}: lambda = {sprintf('%.4g', lam)}")
		fit_lam <- .dbn_als_tv_fixed_lambda(
			Y, family = family, symmetric = symmetric,
			lambda = lam, mu = mu,
			max_iter = tv_max_iter, tol_obj = tv_tol_obj, tol_par = tv_tol_par,
			gauge = gauge,
			init = list(A_init = A_curr, B_init = B_curr, M_init = M_static),
			verbose = FALSE
		)
		fits[[i]] <- fit_lam
		# warm-start next lambda from this fit
		A_curr <- fit_lam$A
		B_curr <- fit_lam$B
	}

	list(fits = fits, lambda_grid = lambda_grid, mu = mu)
}

#' Compute rolling-origin CV loss for a fixed lambda
#'
#' Trains on transitions t in [2, t*], predicts transition t*+1, accumulates
#' loss across t*.
#'
#' @keywords internal
#' @noRd
.tv_als_cv_loss <- function(Y, family, symmetric, lambda, mu, eta, gauge,
                              tv_max_iter, tv_tol_obj, tv_tol_par,
                              t_min = NULL, verbose = FALSE) {
	dims <- dim(Y); n_row <- dims[1]; n_col <- dims[2]; p <- dims[3]; Tt <- dims[4]
	# min training size: at least 5 transitions before first CV evaluation
	if (is.null(t_min)) t_min <- max(5L, ceiling(Tt / 3))
	if (t_min >= Tt - 1L) {
		# not enough data for rolling CV; return Inf
		return(list(loss = Inf, n_folds = 0L, per_fold = numeric(0)))
	}

	per_fold <- numeric(0)
	for (t_star in t_min:(Tt - 1L)) {
		# train on Y[, , , 1:t_star]
		Y_train <- Y[, , , 1:t_star, drop = FALSE]
		dim_train <- dim(Y_train)
		# fit TV-ALS and collect any error message so a per-fold failure
		# is reported rather than silently producing Inf scores.
		fit_train <- tryCatch(
			.dbn_als_tv_fixed_lambda(Y_train, family = family, symmetric = symmetric,
				lambda = lambda, mu = mu,
				max_iter = tv_max_iter, tol_obj = tv_tol_obj, tol_par = tv_tol_par,
				gauge = gauge, init = NULL, verbose = FALSE),
			error = function(e) structure(list(error = conditionMessage(e)),
				class = "cv_fold_error"))
		if (inherits(fit_train, "cv_fold_error") || is.null(fit_train)) {
			if (inherits(fit_train, "cv_fold_error")) {
				cli::cli_warn(c(
					"CV fold (t* = {t_star}, lambda = {sprintf('%.4g', lambda)}) refit failed.",
					"x" = "{fit_train$error}",
					"i" = "Marking fold loss as Inf; CV selection will treat this fold as invalid."
				))
			}
			per_fold <- c(per_fold, Inf)
			next
		}
		# one-step prediction at t_star + 1: extend last operator
		A_last <- fit_train$A[, , dim(fit_train$A)[3]]
		B_last <- fit_train$B[, , dim(fit_train$B)[3]]
		# state at t_star (centered)
		Phi_last <- fit_train$Phi[, , dim(fit_train$Phi)[3]]
		Omega_last <- fit_train$Omega[, , dim(fit_train$Omega)[3]]
		# prediction (use family-specific later; for Gaussian, mean prediction)
		pred <- A_last %*% Phi_last %*% t(B_last)
		# observed at t_star + 1, centered
		# for Gaussian: Phi_{t*+1} = Y_{t*+1} - M
		if (family == "gaussian") {
			Y_next <- Y[, , 1, t_star + 1L]  # assume p=1 for simplicity
			M_hat <- fit_train$M[, , 1]
			Phi_next <- Y_next - M_hat
			# mask
			obs_mask <- !is.na(Phi_next)
			if (n_row == n_col) diag(obs_mask) <- FALSE
			Phi_next[!obs_mask] <- 0
			pred[!obs_mask] <- 0
			loss_t <- sum((Phi_next - pred)^2) / max(sum(obs_mask), 1)
			per_fold <- c(per_fold, loss_t)
		} else if (family == "binary") {
			# held-out negative log-likelihood under the probit link.
			# theta_pred = M + A_t Phi_{t-1} B_t^T (one-step prediction on the
			# latent scale). For each observed Y_{t+1} in {0, 1}:
			#   loss_ij = -log P(Y_ij | theta_ij)
			#           = -log Phi(theta) if Y=1, -log Phi(-theta) if Y=0
			Y_next <- Y[, , 1, t_star + 1L]
			M_hat <- fit_train$M[, , 1]
			theta_pred <- M_hat + pred
			obs_mask <- !is.na(Y_next)
			if (n_row == n_col) diag(obs_mask) <- FALSE
			# clamp theta to avoid log(0)
			th <- pmin(pmax(theta_pred, -8), 8)
			ll <- ifelse(Y_next == 1, log(pmax(pnorm(th), 1e-12)),
			                          log(pmax(pnorm(-th), 1e-12)))
			ll[!obs_mask] <- 0
			loss_t <- -sum(ll) / max(sum(obs_mask), 1)
			per_fold <- c(per_fold, loss_t)
		} else if (family == "ordinal") {
			# held-out predictive ordinal probability via per-cell truncated-
			# normal: P(Y = c | theta) = Phi(gamma_c - theta) - Phi(gamma_{c-1} - theta)
			# loss = -mean log P over observed cells.
			Y_next <- Y[, , 1, t_star + 1L]
			M_hat <- fit_train$M[, , 1]
			theta_pred <- M_hat + pred
			obs_mask <- !is.na(Y_next)
			if (n_row == n_col) diag(obs_mask) <- FALSE
			# construct thresholds from observed Y across training (heuristic)
			Y_int <- as.integer(round(Y[, , 1, 1:t_star]))
			cats <- sort(unique(Y_int[!is.na(Y_int)]))
			K <- length(cats)
			if (K < 2L) {
				per_fold <- c(per_fold, NA_real_)
				next
			}
			probs <- seq(1 / K, 1 - 1 / K, length.out = K - 1L)
			thr <- quantile(theta_pred[obs_mask], probs, na.rm = TRUE)
			gammas <- c(-Inf, thr, Inf)
			# per-observed-cell log-prob
			ll_sum <- 0; n_obs <- 0
			Y_next_int <- as.integer(round(Y_next))
			for (k in seq_along(cats)) {
				c <- cats[k]
				idx <- which(obs_mask & Y_next_int == c)
				if (!length(idx)) next
				th_c <- pmin(pmax(theta_pred[idx], -8), 8)
				lo <- gammas[k] - th_c
				hi <- gammas[k + 1L] - th_c
				p_c <- pmax(pnorm(hi) - pnorm(lo), 1e-12)
				ll_sum <- ll_sum + sum(log(p_c))
				n_obs <- n_obs + length(idx)
			}
			loss_t <- if (n_obs > 0) -ll_sum / n_obs else NA_real_
			per_fold <- c(per_fold, loss_t)
		} else {
			per_fold <- c(per_fold, NA_real_)
		}
	}
	per_fold_finite <- per_fold[is.finite(per_fold)]
	if (length(per_fold_finite) == 0) return(list(loss = Inf, n_folds = 0L, per_fold = per_fold))
	list(loss = mean(per_fold_finite), n_folds = length(per_fold_finite),
	     per_fold = per_fold,
	     loss_sd = sd(per_fold_finite))
}

#' Select lambda by rolling-origin CV with 1-SE rule
#'
#' @keywords internal
#' @noRd
.tv_als_select_lambda_cv <- function(Y, family, symmetric, lambda_grid,
                                       mu, eta, gauge,
                                       tv_max_iter, tv_tol_obj, tv_tol_par,
                                       verbose = FALSE) {
	cv_results <- vector("list", length(lambda_grid))
	# always show a progress bar for CV: this loop fits a model per lambda
	# value (and per fold inside), which can take minutes on moderate
	# data. silent waiting is the worst-of-both UX.
	show_progress <- length(lambda_grid) >= 2L &&
		!isTRUE(getOption("dbn.suppress_cv_progress", FALSE))
	if (show_progress) {
		cli::cli_progress_bar(
			name = "TV-ALS CV",
			total = length(lambda_grid),
			format = "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
			format_done = "{cli::pb_name} done ({cli::pb_total} lambda values).",
			.envir = parent.frame()
		)
	}
	for (i in seq_along(lambda_grid)) {
		lam <- lambda_grid[i]
		if (verbose) cli::cli_inform("CV at lambda = {sprintf('%.4g', lam)}")
		cv_results[[i]] <- .tv_als_cv_loss(Y, family, symmetric, lam, mu, eta,
			gauge, tv_max_iter, tv_tol_obj, tv_tol_par, verbose = FALSE)
		if (show_progress) {
			cli::cli_progress_update(.envir = parent.frame())
		}
	}
	if (show_progress) cli::cli_progress_done(.envir = parent.frame())
	losses <- sapply(cv_results, function(r) r$loss)
	loss_sds <- sapply(cv_results, function(r) if (is.null(r$loss_sd)) NA_real_ else r$loss_sd)
	n_folds <- sapply(cv_results, function(r) r$n_folds)
	se_losses <- loss_sds / sqrt(pmax(n_folds, 1))

	# 1-SE rule: pick the LARGEST lambda whose loss is within 1 SE of the minimum
	min_loss <- min(losses, na.rm = TRUE)
	idx_min <- which.min(losses)
	se_min <- if (idx_min <= length(se_losses)) se_losses[idx_min] else 0
	threshold <- min_loss + (if (!is.na(se_min)) se_min else 0)
	# lambda_grid is sorted descending, so "largest lambda" is the first satisfying
	candidates <- which(losses <= threshold & !is.na(losses))
	if (length(candidates) == 0) candidates <- idx_min
	# pick the smallest index (= largest lambda) among candidates
	idx_sel <- min(candidates)

	# boundary warning: if CV picks the smoothest or roughest end of the grid,
	# the grid likely didn't extend far enough -- warn so users can extend it.
	# lambda_grid is sorted DESCENDING: idx = 1 is the largest (smoothest) lambda;
	# idx = length(grid) is the smallest (roughest) lambda.
	if (length(lambda_grid) >= 2L) {
		if (idx_sel == 1L) {
			cli::cli_warn(c(
				"CV picked the {.emph largest} {.arg lambda} in the grid ({.val {signif(lambda_grid[1], 4)}}).",
				"i" = "The smoothest end may not be smooth enough; consider extending the grid upward via {.arg lambda_grid = c(lambda * 10, lambda * 100, ...)}.",
				"i" = "An even-larger lambda might yield a better-fitting model."
			))
		} else if (idx_sel == length(lambda_grid)) {
			cli::cli_warn(c(
				"CV picked the {.emph smallest} {.arg lambda} in the grid ({.val {signif(lambda_grid[length(lambda_grid)], 4)}}).",
				"i" = "The roughest end may not be rough enough; consider extending the grid downward.",
				"i" = "An even-smaller lambda might fit the data better -- check whether the model is under-regularized."
			))
		}
	}

	list(
		lambda_grid = lambda_grid,
		losses = losses,
		loss_sds = loss_sds,
		n_folds = n_folds,
		idx_min = idx_min,
		idx_1se = idx_sel,
		lambda_min = lambda_grid[idx_min],
		lambda_1se = lambda_grid[idx_sel],
		boundary_hit = (idx_sel == 1L) || (idx_sel == length(lambda_grid)),
		cv_results = cv_results
	)
}

#' three-point local CV refinement around a pilot lambda
#'
#' @keywords internal
#' @noRd
.tv_als_select_lambda_local <- function(Y, family, symmetric, lambda_candidates,
                                          mu, eta, gauge,
                                          tv_max_iter, tv_tol_obj, tv_tol_par,
                                          verbose = FALSE) {
	losses <- numeric(length(lambda_candidates))
	for (i in seq_along(lambda_candidates)) {
		lam <- lambda_candidates[i]
		res <- .tv_als_cv_loss(Y, family, symmetric, lam, mu, eta,
			gauge, tv_max_iter, tv_tol_obj, tv_tol_par, verbose = FALSE)
		losses[i] <- res$loss
		if (verbose) cli::cli_inform("local CV at lambda = {sprintf('%.4g', lam)} -> loss = {sprintf('%.4g', losses[i])}")
	}
	if (all(is.na(losses))) return(lambda_candidates[2])
	lambda_candidates[which.min(losses)]
}

#' Build the default lambda grid centered on the pilot
#'
#' @keywords internal
#' @noRd
.tv_als_make_lambda_grid <- function(lambda_pilot, n_lambdas = 12L,
                                       min_ratio = 1e-2, max_ratio = 1e2) {
	# log-spaced from lambda_pilot * min_ratio to lambda_pilot * max_ratio
	exp(seq(log(lambda_pilot * max_ratio), log(lambda_pilot * min_ratio),
		length.out = n_lambdas))
}
