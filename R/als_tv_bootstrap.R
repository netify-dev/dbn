####
# bootstrap for tv-als with stability auto-switch.
# default fixed-design parametric bootstrap keeps the observed lags fixed
# and resamples only the errors:
#   Z_t^{(b)} = A_hat_t * Phi_{t-1}^state * B_hat_t^T + epsilon_t^{(b)}
# the recursive parametric variant auto-switches to fixed-design when
# max_t rho(A_hat_t) * rho(B_hat_t) >= 1 - delta.
# a moving-block residual bootstrap is also available.
####

#' Estimate M from a Y array as the per-cell time-and-mask average
#'
#' Matches the convention used by `.tv_als_warm_start_static()` so each
#' bootstrap replicate gets a per-rep M estimate (rather than reusing the
#' original fit's M). Used by the TV-ALS bootstrap to inject the correct
#' level-term variance into the resampling distribution; without this,
#' empirical CI coverage on M collapses to near zero.
#'
#' @keywords internal
#' @noRd
.tv_als_estimate_M_from_Y <- function(Y, n_row, n_col, p, Tt) {
	M <- array(0, dim = c(n_row, n_col, p))
	for (r in seq_len(p)) {
		num <- matrix(0, n_row, n_col)
		den <- matrix(0, n_row, n_col)
		for (t in seq_len(Tt)) {
			y_t <- Y[, , r, t]
			ok <- !is.na(y_t)
			if (n_row == n_col) diag(ok) <- FALSE
			num[ok] <- num[ok] + y_t[ok]
			den[ok] <- den[ok] + 1
		}
		# avoid 0/0 -- structural-diagonal cells get M = 0 (their predictions
		# are masked anyway)
		den_safe <- den
		den_safe[den_safe == 0] <- 1
		M_r <- num / den_safe
		M_r[den == 0] <- 0
		M[, , r] <- M_r
	}
	M
}

#' Compute max spectral radius product across time for a TV-ALS fit
#'
#' @keywords internal
#' @noRd
.tv_als_max_spectral_product <- function(A_arr, B_arr) {
	Tt <- dim(A_arr)[3]
	if (Tt < 2L) return(0)
	rhos <- numeric(Tt - 1L)
	for (t in 2:Tt) {
		# eigen() may legitimately fail on slices with non-finite entries
		# (e.g. divergent draws). NA propagates honestly through `max`.
		ra <- tryCatch(max(Mod(eigen(A_arr[, , t], only.values = TRUE)$values)),
			error = function(e) NA_real_)
		rb <- tryCatch(max(Mod(eigen(B_arr[, , t], only.values = TRUE)$values)),
			error = function(e) NA_real_)
		rhos[t - 1L] <- if (is.finite(ra) && is.finite(rb)) ra * rb else NA_real_
	}
	max(rhos, na.rm = TRUE)
}

#' Bootstrap for TV-ALS fits
#'
#' @param fit a dbn TV-ALS fit object (meta$sampler_used == "als_tv")
#' @param R number of bootstrap replicates
#' @param type "fixed" (default), "recursive", or "residual_block"
#' @param block_length block length for residual block bootstrap (default ceil(T^(1/3)))
#' @param stability_delta safety margin for auto-switch (default 0.05)
#' @param seed optional integer seed
#' @param verbose progress
#' @param parallel logical. If TRUE (default) and `future.apply` is loaded
#'   with a non-sequential plan active, the per-rep refits run in parallel.
#'   Otherwise falls back to a sequential for-loop.
#' @return dbn_boot object
#' @keywords internal
#' @noRd
.dbn_als_tv_bootstrap <- function(fit, R = 200L,
                                    type = c("fixed", "recursive", "residual_block"),
                                    block_length = NULL,
                                    stability_delta = 0.05,
                                    seed = NULL,
                                    verbose = TRUE,
                                    parallel = TRUE) {
	type <- match.arg(type)
	if (!inherits(fit, "dbn") ||
		!(identical(fit$meta$sampler_used, "als_tv") ||
		  identical(fit$meta$sampler_used, "als_piecewise"))) {
		cli::cli_abort(c(
			"{.arg fit} must have time-varying operators ({.code meta$sampler_used} in {.val als_tv} or {.val als_piecewise}).",
			"i" = "Got {.val {fit$meta$sampler_used}}."
		))
	}
	if (!is.null(seed)) set.seed(seed)

	# suppress the per-rep diagonal warning: the top-level dbn_als() call
	# already fired once if needed; firing it 200x more inside the loop just
	# trains users to ignore warnings.
	.dbn_prev_in_boot <- getOption("dbn.in_bootstrap", FALSE)
	options(dbn.in_bootstrap = TRUE)
	on.exit(options(dbn.in_bootstrap = .dbn_prev_in_boot), add = TRUE)

	# pull point estimates
	A_arr <- fit$A[[1]]  # [n_row, n_row, Tt]
	B_arr <- fit$B[[1]]  # [n_col, n_col, Tt]
	M     <- fit$M[, , , 1, drop = FALSE]
	dim(M) <- c(dim(A_arr)[1], dim(B_arr)[1], dim(fit$M)[3])
	dims <- fit$dims
	n_row <- dims$n_row; n_col <- dims$n_col; p <- dims$p; Tt <- dims$Tt
	family <- fit$family
	# lambda_hat / mu_hat are only used for TV refits; piecewise refits use
	# its own per-block static solve and ignores these.
	lambda_hat <- fit$meta$tv$lambda_used %||% NA_real_
	mu_hat <- fit$meta$tv$mu %||% NA_real_
	# honor the outer fit's tv_max_iter for the per-rep inner refits so a
	# user who raises it sees that change propagate; default to 100 when
	# the outer fit didn't store the setting.
	tv_max_iter_inner <- fit$settings$tv_max_iter %||%
		fit$meta$settings$tv_max_iter %||% 100L

	Y <- fit$Y
	# state: Phi_t^state used as predictor
	# compute mask
	Omega <- array(1, dim = c(n_row, n_col, Tt))
	if (n_row == n_col) for (t in 1:Tt) diag(Omega[, , t]) <- 0
	for (r in 1:p) for (t in 1:Tt) {
		Omega[, , t][is.na(Y[, , r, t])] <- 0
	}
	Y_filled <- Y; Y_filled[is.na(Y_filled)] <- 0
	Phi_state <- array(0, dim = c(n_row, n_col, Tt))
	for (t in 1:Tt) {
		Phi_state[, , t] <- Y_filled[, , 1, t] - M[, , 1]
		Phi_state[, , t][Omega[, , t] == 0] <- 0
	}

	# compute residuals on the state scale
	resid_arr <- array(0, dim = c(n_row, n_col, Tt))
	resid_var <- 0; n_resid <- 0
	for (t in 2:Tt) {
		pred <- A_arr[, , t] %*% Phi_state[, , t - 1L] %*% t(B_arr[, , t])
		resid_arr[, , t] <- Omega[, , t] * (Phi_state[, , t] - pred)
		resid_var <- resid_var + sum(resid_arr[, , t]^2)
		n_resid <- n_resid + sum(Omega[, , t])
	}
	sigma_hat <- sqrt(resid_var / max(n_resid, 1))

	# stability check + auto-switch
	max_rho <- .tv_als_max_spectral_product(A_arr, B_arr)
	switched <- FALSE
	if (type == "recursive" && is.finite(max_rho) && max_rho >= 1 - stability_delta) {
		cli::cli_warn(c(
			"Recursive bootstrap requested but max_t rho(A_t)*rho(B_t) = {sprintf('%.3f', max_rho)} >= 1 - delta.",
			"i" = "Auto-switching to {.val fixed} design to avoid explosive bootstrap replicates.",
			"i" = "Override with {.code bootstrap_design = \"recursive\"} explicitly to ignore."
		))
		type <- "fixed"
		switched <- TRUE
	}

	if (verbose) {
		cli::cli_h3("TV-ALS bootstrap: type = {.val {type}}, R = {R}")
		cli::cli_inform("max rho(A_t) * rho(B_t) = {sprintf('%.3f', max_rho)}")
	}

	# prepare storage: per-time A and B replicates
	# use t = 2..Tt slices (the meaningful operator slices)
	S <- Tt - 1L
	coefs_A <- array(NA, dim = c(R, n_row, n_row, S))
	coefs_B <- array(NA, dim = c(R, n_col, n_col, S))
	coefs_M <- array(NA, dim = c(R, n_row, n_col, p))
	# also store sigma2 (constant) and per-slice norms for diagnostics
	rep_max_rho <- numeric(R)
	rep_failed <- logical(R)
	# track which replicates actually converged versus those that returned
	# but hit max-iter, so the user can judge effective sample size honestly.
	rep_converged <- logical(R)
	# collect bootstrap per-rep error messages so failures don't disappear
	# silently. The first N unique messages are warned at the end (so the
	# user sees what went wrong, not just a low n_valid count).
	rep_errors <- character(R)

	# per-rep closure: builds Y_b, refits, returns a list of {A, B, M, max_rho,
	# converged, failed, err}. defined inside the function so it inherits all
	# the locals captured from the outer scope without explicit passing.
	# stripped of <<- side effects so future_lapply can copy the closure to
	# workers cleanly.
	do_rep <- function(b) {
		if (verbose && b %% max(1L, R %/% 10L) == 0L)
			cli::cli_inform("  rep {b}/{R}")

		# build Y_b based on type.
		# for gaussian: simulate on the latent-state (centered Phi) scale
		#   then add back M to get Y_b.
		# for binary/ordinal: simulate latent Z_b at the Y-level (Bernoulli /
		#   ordinal categorical draws from the implied probit / rank model),
		#   so the bootstrap replicate has the right OBSERVED-DATA distribution
		#   and the ECM Z-update inside the refit will re-impute correctly.
		Y_b <- Y
		.observe_from_theta <- function(theta_t, fam) {
			# theta_t is the predicted mean on the latent scale, n_row x n_col.
			# return an observed-Y matrix of the same shape, respecting Omega.
			if (fam == "gaussian") {
				eps <- matrix(rnorm(n_row * n_col, 0, sigma_hat), n_row, n_col)
				M[, , 1] + theta_t + eps
			} else if (fam == "binary") {
				th <- pmin(pmax(theta_t, -8), 8)
				p_pos <- pnorm(th)
				matrix(rbinom(n_row * n_col, 1, p_pos), n_row, n_col)
			} else if (fam == "ordinal") {
				# generate categorical Y_b from the implied rank-likelihood:
				#   draw Z_b ~ N(theta_t, 1), assign to category by empirical
				#   cuts of the observed Y_t (matching the original family's
				#   discretization). For each time t we use the GLOBAL category
				#   set of the observed Y for stability.
				cats <- sort(unique(as.integer(round(Y[!is.na(Y)]))))
				K <- length(cats)
				if (K < 2L) {
					eps <- matrix(rnorm(n_row * n_col, 0, sigma_hat), n_row, n_col)
					return(theta_t + eps)
				}
				probs <- seq(1 / K, 1 - 1 / K, length.out = K - 1L)
				cuts <- quantile(theta_t, probs, na.rm = TRUE)
				Z_b <- theta_t + matrix(rnorm(n_row * n_col), n_row, n_col)
				out <- matrix(NA_real_, n_row, n_col)
				out[Z_b <= cuts[1]] <- cats[1]
				for (k in seq_along(cats)[-1]) {
					if (k <= length(cuts)) {
						out[Z_b > cuts[k - 1] & Z_b <= cuts[k]] <- cats[k]
					} else {
						out[Z_b > cuts[k - 1]] <- cats[k]
					}
				}
				out
			}
		}

		if (type == "fixed") {
			# fixed-design: use observed Phi_{t-1}^state as lag
			Phi_b <- array(0, dim = c(n_row, n_col, Tt))
			Phi_b[, , 1] <- Phi_state[, , 1]  # init same
			for (t in 2:Tt) {
				pred <- A_arr[, , t] %*% Phi_state[, , t - 1L] %*% t(B_arr[, , t])
				if (family == "gaussian") {
					eps <- matrix(rnorm(n_row * n_col, 0, sigma_hat), n_row, n_col)
					Phi_b[, , t] <- pred + eps
					Phi_b[, , t][Omega[, , t] == 0] <- 0
					Y_b[, , 1, t] <- M[, , 1] + Phi_b[, , t]
				} else {
					# theta = M + pred on the latent scale
					theta_t <- M[, , 1] + pred
					Y_b[, , 1, t] <- .observe_from_theta(theta_t, family)
				}
			}
			if (family == "gaussian") {
				Y_b[, , 1, 1] <- M[, , 1] + Phi_state[, , 1]
			} else {
				# t = 1: use observed Y (no transition prediction for it)
				Y_b[, , 1, 1] <- Y[, , 1, 1]
			}
			# restore NA where Omega = 0
			for (t in 1:Tt) Y_b[, , 1, t][Omega[, , t] == 0] <- NA_real_
		} else if (type == "recursive") {
			# recursive: Phi_t^{(b)} depends on Phi_{t-1}^{(b)}.
			# for non-gaussian, generate Y_b at the Y-level from theta_pred
			# computed from the simulated Phi_b.
			Phi_b <- array(0, dim = c(n_row, n_col, Tt))
			Phi_b[, , 1] <- Phi_state[, , 1]
			for (t in 2:Tt) {
				pred <- A_arr[, , t] %*% Phi_b[, , t - 1L] %*% t(B_arr[, , t])
				if (family == "gaussian") {
					eps <- matrix(rnorm(n_row * n_col, 0, sigma_hat), n_row, n_col)
					Phi_b[, , t] <- pred + eps
				} else {
					# need a latent draw to evolve the recursion; use the
					# generated Y as the new "observed" for the next step
					theta_t <- M[, , 1] + pred
					Y_b[, , 1, t] <- .observe_from_theta(theta_t, family)
					# update Phi_b on the latent scale by re-imputing Z from
					# Y_b given theta_t (deterministic mean)
					Phi_b[, , t] <- .expected_Z(Y_b[, , 1, t], theta_t, family) - M[, , 1]
				}
				Phi_b[, , t][Omega[, , t] == 0] <- 0
			}
			if (family == "gaussian") {
				for (t in 1:Tt) Y_b[, , 1, t] <- M[, , 1] + Phi_b[, , t]
			}
			for (t in 1:Tt) Y_b[, , 1, t][Omega[, , t] == 0] <- NA_real_
		} else if (type == "residual_block") {
			# moving-block residual: resample blocks of residual matrices
			if (is.null(block_length)) block_length <- max(2L, ceiling(Tt^(1 / 3)))
			# get residual slices from t=2..Tt
			resid_times <- 2:Tt
			# construct sampled block starts
			n_blocks_needed <- ceiling(length(resid_times) / block_length)
			Phi_b <- array(0, dim = c(n_row, n_col, Tt))
			Phi_b[, , 1] <- Phi_state[, , 1]
			# resample block starts
			possible_starts <- 2:(Tt - block_length + 1L)
			if (length(possible_starts) < 1L) possible_starts <- 2:Tt
			block_starts <- sample(possible_starts, n_blocks_needed, replace = TRUE)
			# build a contiguous resid sequence
			resid_seq <- list()
			for (bs in block_starts) {
				for (offset in 0:(block_length - 1L)) {
					tt <- bs + offset
					if (tt <= Tt) resid_seq[[length(resid_seq) + 1L]] <- resid_arr[, , tt]
				}
			}
			# truncate to S
			resid_seq <- resid_seq[1:S]
			for (t in 2:Tt) {
				pred <- A_arr[, , t] %*% Phi_state[, , t - 1L] %*% t(B_arr[, , t])
				Phi_b[, , t] <- pred + resid_seq[[t - 1L]]
				Phi_b[, , t][Omega[, , t] == 0] <- 0
			}
			if (family == "gaussian") {
				for (t in 1:Tt) Y_b[, , 1, t] <- M[, , 1] + Phi_b[, , t]
			} else {
				# for non-gaussian, the residual-block scheme produces a
				# replicate latent trajectory; observe Y_b from theta = M + Phi_b
				for (t in 2:Tt) {
					theta_t <- M[, , 1] + Phi_b[, , t]
					Y_b[, , 1, t] <- .observe_from_theta(theta_t, family)
				}
				Y_b[, , 1, 1] <- Y[, , 1, 1]
				for (t in 1:Tt) Y_b[, , 1, t][Omega[, , t] == 0] <- NA_real_
			}
		}

		# re-fit at the inherited lambda (or per-block static for piecewise).
		# route by (sampler_used, symmetric). Binary/ordinal go through the
		# directed solver with Y_obs (for ECM Z-update) or symmetric L-BFGS
		# with Y_obs.
		rep_err <- ""
		fit_b <- tryCatch({
			is_sym <- isTRUE(fit$dims$is_symmetric)
			is_piecewise <- identical(fit$meta$sampler_used, "als_piecewise")
			if (is_piecewise) {
				# refit the same breakpoint structure on Y_b.
				# per-rep seed: NULL inherits the outer RNG state (advanced per
				# rep by the .observe_from_theta calls); each rep gets distinct
				# RNG draws naturally.
				breaks_b <- fit$meta$breaks
				pf <- .dbn_als_piecewise(Y_b, family = family, symmetric = is_sym,
					breaks = breaks_b, ridge = 3e-3, max_iter = tv_max_iter_inner,
					verbose = FALSE, seed = NULL)
				# repackage to match the loop's expectations
				list(A = pf$A[[1]], B = pf$B[[1]])
			} else if (!is_sym) {
				# directed: handles gaussian, binary, ordinal uniformly via
				# the ECM Z-update inside .dbn_als_tv_fixed_lambda when
				# Y_obs is supplied.
				A_init_b <- A_arr; B_init_b <- B_arr
				# the inner TV-ALS solver holds M fixed at M_init (the
				# closed-form M conditional is degenerate under Phi = Z - M),
				# so we have to inject the per-rep level-term variance from
				# outside. for gaussian Y_b (already on the latent scale) we
				# re-estimate M_init_b from the per-cell time-and-mask
				# average of Y_b. for non-gaussian families Y_b is on the
				# observation scale (not the M scale), so we drop the init
				# entirely and let the inner warm-start recompute M directly.
				if (family == "gaussian") {
					M_init_b <- .tv_als_estimate_M_from_Y(Y_b, n_row, n_col, p, Tt)
					init_b <- list(A_init = A_init_b, B_init = B_init_b, M_init = M_init_b)
					Y_obs_b <- NULL
				} else {
					init_b <- NULL  # warm-start re-estimates M from Y_b
					Y_obs_b <- Y_b
				}
				.dbn_als_tv_fixed_lambda(
					Y_b, family = family, symmetric = FALSE,
					lambda = lambda_hat, mu = mu_hat,
					max_iter = tv_max_iter_inner, tol_obj = 1e-4, tol_par = 1e-3,
					gauge = "penalty_min",
					init = init_b,
					Y_obs = Y_obs_b,
					verbose = FALSE
				)
			} else {
				# symmetric: L-BFGS solver with ECM Z-update for non-gaussian.
				# re-estimate M from Y_b for gaussian; drop init for non-
				# gaussian (see asymmetric branch above).
				if (family == "gaussian") {
					M_init_b <- .tv_als_estimate_M_from_Y(Y_b, n_row, n_col, p, Tt)
					init_b <- list(A_init = A_arr, M_init = M_init_b)
					Y_obs_b <- NULL
				} else {
					init_b <- NULL
					Y_obs_b <- Y_b
				}
				.dbn_als_tv_symmetric_lbfgs(
					Y_b, lambda = lambda_hat, mu = mu_hat,
					max_iter = tv_max_iter_inner, tol_obj = 1e-4,
					init = init_b,
					family = family,
					Y_obs = Y_obs_b,
					verbose = FALSE
				)
			}
		}, error = function(e) {
			rep_err <<- conditionMessage(e)
			NULL
		})

		if (is.null(fit_b)) {
			return(list(failed = TRUE, err = rep_err, A = NULL, B = NULL,
				M = NULL, converged = FALSE, max_rho = NA_real_))
		}
		Ab <- fit_b$A
		Bb <- if (!is.null(fit_b$B)) fit_b$B else Ab  # symmetric
		Mb <- fit_b$M
		list(failed = FALSE, err = "",
			A = Ab, B = Bb, M = Mb,
			converged = isTRUE(fit_b$converged),
			max_rho = .tv_als_max_spectral_product(Ab, Bb))
	}

	# decide between sequential and future_lapply based on parallel + plan.
	# we only parallelise when future.apply is available AND a non-sequential
	# plan is active; otherwise lapply (sequential) is identical in semantics
	# and avoids future.apply's overhead.
	use_parallel_boot <- isTRUE(parallel) &&
		requireNamespace("future.apply", quietly = TRUE) &&
		requireNamespace("future", quietly = TRUE)
	if (use_parallel_boot) {
		plan_cls <- class(future::plan())
		if ("sequential" %in% plan_cls || "uniprocess" %in% plan_cls) {
			use_parallel_boot <- FALSE
		}
	}
	# devtools::load_all caveat: workers spawned by future::plan(multisession)
	# call library(dbn) (via future.packages = "dbn"), which loads the
	# INSTALLED dbn from .libPaths(), NOT the parent's load_all'd development
	# state. If the dev tree has internals (like .tv_als_estimate_M_from_Y)
	# that the installed version lacks, every replicate errors silently and
	# n_valid = 0. Detect and refuse cleanly with a directive.
	if (use_parallel_boot &&
	    requireNamespace("pkgload", quietly = TRUE) &&
	    pkgload::is_dev_package("dbn")) {
		cli::cli_warn(c(
			"Detected {.fn devtools::load_all} state for {.pkg dbn}.",
			"i" = "Parallel bootstrap workers would load the INSTALLED {.pkg dbn} from {.path .libPaths()}, not the dev tree, and any dev-only internals would be missing.",
			"i" = "Falling back to a sequential bootstrap so the replicates run against the same code you're auditing.",
			"i" = "For real parallel speedup, install the package first: {.code R CMD INSTALL <pkg>} (or {.code devtools::install()})."
		))
		use_parallel_boot <- FALSE
	}
	# future.packages = "dbn" makes the workers load the package before
	# evaluating do_rep, so internal helpers like .dbn_als_tv_fixed_lambda
	# and .tv_als_estimate_M_from_Y resolve correctly. without this the
	# closure references unbound symbols on workers spawned from a fresh
	# session (which is what `future::plan(future::multisession)` does).
	reps_out <- if (use_parallel_boot) {
		future.apply::future_lapply(seq_len(R), do_rep,
			future.seed = if (!is.null(seed)) seed else TRUE,
			future.packages = "dbn")
	} else {
		lapply(seq_len(R), do_rep)
	}

	for (b in seq_len(R)) {
		ro <- reps_out[[b]]
		if (isTRUE(ro$failed)) {
			rep_failed[b] <- TRUE
			rep_errors[b] <- ro$err
			next
		}
		Ab <- ro$A; Bb <- ro$B; Mb <- ro$M
		rep_converged[b] <- ro$converged
		rep_max_rho[b] <- ro$max_rho
		# store slices t=2..Tt
		coefs_A[b, , , ] <- Ab[, , 2:Tt]
		coefs_B[b, , , ] <- Bb[, , 2:Tt]
		if (!is.null(Mb)) {
			# Mb may be [n_row, n_col, p] or [n_row, n_col] (p=1 with drop)
			if (length(dim(Mb)) == 2L) {
				coefs_M[b, , , 1] <- Mb
			} else {
				coefs_M[b, , , ] <- Mb
			}
		}
	}

	# surface per-rep errors loudly so a low n_valid is explained by the
	# top unique error classes rather than left to the caller to diagnose.
	if (any(rep_failed)) {
		err_msgs <- rep_errors[rep_failed]
		err_table <- sort(table(err_msgs), decreasing = TRUE)
		# show up to 3 unique error classes
		n_shown <- min(3L, length(err_table))
		cli::cli_warn(c(
			"{sum(rep_failed)} of {R} bootstrap replicates failed to refit.",
			"i" = "Top unique error{?s} (showing {n_shown}):",
			setNames(paste0("(", as.integer(err_table[1:n_shown]), "x) ",
				substr(names(err_table[1:n_shown]), 1, 80)),
				rep("x", n_shown))
		))
	}

	# aggregate: pointwise quantiles, SE, means
	valid_idx <- which(!rep_failed)
	n_valid <- length(valid_idx)
	n_converged <- sum(rep_converged & !rep_failed)
	if (n_valid < 2L) {
		cli::cli_warn("Bootstrap had fewer than 2 valid replicates ({n_valid}).")
	}
	# surface non-convergence rate: a rep that returned-but-hit-tv_max_iter
	# contributes biased quantiles. Warn when many reps fail to converge.
	n_nonconverged <- n_valid - n_converged
	if (n_valid >= 2L && n_nonconverged >= max(1L, n_valid %/% 5L)) {
		pct <- round(100 * n_nonconverged / n_valid)
		tv_max_iter_hint <- fit$meta$tv$tv_max_iter
		if (is.null(tv_max_iter_hint)) tv_max_iter_hint <- 50L
		cli::cli_warn(c(
			"{n_nonconverged} of {n_valid} bootstrap replicates did not converge ({pct}%).",
			"i" = "These replicates hit the per-rep iteration cap ({tv_max_iter_hint}) without satisfying the tolerance.",
			"i" = "Intervals based on these reps may be biased; raise the outer fit's {.arg tv_max_iter} or loosen {.arg tv_tol_obj}.",
			"i" = "Inspect {.code fit$bootstrap$n_converged} vs {.code fit$bootstrap$n_valid}."
		))
	}

	if (n_valid >= 2L) {
		ci_A_lo <- apply(coefs_A[valid_idx, , , , drop = FALSE], c(2, 3, 4), quantile, 0.025, na.rm = TRUE)
		ci_A_hi <- apply(coefs_A[valid_idx, , , , drop = FALSE], c(2, 3, 4), quantile, 0.975, na.rm = TRUE)
		se_A    <- apply(coefs_A[valid_idx, , , , drop = FALSE], c(2, 3, 4), sd, na.rm = TRUE)
		ci_B_lo <- apply(coefs_B[valid_idx, , , , drop = FALSE], c(2, 3, 4), quantile, 0.025, na.rm = TRUE)
		ci_B_hi <- apply(coefs_B[valid_idx, , , , drop = FALSE], c(2, 3, 4), quantile, 0.975, na.rm = TRUE)
		se_B    <- apply(coefs_B[valid_idx, , , , drop = FALSE], c(2, 3, 4), sd, na.rm = TRUE)
		ci_M_lo <- apply(coefs_M[valid_idx, , , , drop = FALSE], c(2, 3, 4), quantile, 0.025, na.rm = TRUE)
		ci_M_hi <- apply(coefs_M[valid_idx, , , , drop = FALSE], c(2, 3, 4), quantile, 0.975, na.rm = TRUE)
		se_M    <- apply(coefs_M[valid_idx, , , , drop = FALSE], c(2, 3, 4), sd, na.rm = TRUE)
	} else {
		ci_A_lo <- ci_A_hi <- se_A <- array(NA, dim = c(n_row, n_row, S))
		ci_B_lo <- ci_B_hi <- se_B <- array(NA, dim = c(n_col, n_col, S))
		ci_M_lo <- ci_M_hi <- se_M <- array(NA, dim = c(n_row, n_col, p))
	}

	result <- list(
		coefs_A = coefs_A,
		coefs_B = coefs_B,
		coefs_M = coefs_M,
		se_A = se_A, se_B = se_B, se_M = se_M,
		ci_A_lo = ci_A_lo, ci_A_hi = ci_A_hi,
		ci_B_lo = ci_B_lo, ci_B_hi = ci_B_hi,
		ci_M_lo = ci_M_lo, ci_M_hi = ci_M_hi,
		point_est_A = A_arr[, , 2:Tt],
		point_est_B = B_arr[, , 2:Tt],
		point_est_M = M,
		n_converged = n_converged,
		rep_converged = rep_converged,
		n_valid = n_valid,
		n_total = R,
		type = type,
		switched_from_recursive = switched,
		max_rho_orig = max_rho,
		rep_max_rho = rep_max_rho,
		family = family,
		symmetric = isTRUE(fit$dims$is_symmetric),
		dims = list(n_row = n_row, n_col = n_col, p = p, Tt = Tt),
		lambda = lambda_hat,
		mu = mu_hat
	)
	class(result) <- c("dbn_tv_boot", "dbn_boot")
	result
}

#' Print method for TV-ALS bootstrap
#' @export
#' @keywords internal
print.dbn_tv_boot <- function(x, ...) {
	cli::cli_h2("TV-ALS bootstrap ({x$type})")
	cli::cli_inform("  Replicates: {x$n_valid} valid / {x$n_total} total")
	n_conv <- x$n_converged %||% NA_integer_
	if (!is.na(n_conv)) {
		cli::cli_inform("  Converged (under tv_max_iter): {n_conv} / {x$n_valid}")
		if (n_conv == 0L) {
			cli::cli_alert_warning("No bootstrap replicate converged within {.code tv_max_iter}; CIs are based on iterate-at-cap point estimates and may understate uncertainty. Raise {.arg tv_max_iter} or loosen {.arg tv_tol_obj} / {.arg tv_tol_par}.")
		} else if (n_conv < x$n_valid) {
			cli::cli_alert_warning("{x$n_valid - n_conv} of {x$n_valid} bootstrap replicates hit {.code tv_max_iter} without converging.")
		}
	}
	cli::cli_inform("  Network: {x$dims$n_row} x {x$dims$n_col}, {x$dims$Tt} time points")
	cli::cli_inform("  lambda = {sprintf('%.4g', x$lambda)}, mu = {sprintf('%.4g', x$mu)}")
	cli::cli_inform("  Max rho(A)*rho(B) at point estimate: {sprintf('%.3f', x$max_rho_orig)}")
	if (x$switched_from_recursive)
		cli::cli_inform("  (Auto-switched from recursive to fixed-design due to stability)")
	cli::cli_inform("  Pointwise CIs available at {.code $ci_A_lo}, {.code $ci_A_hi}, etc.")
	invisible(x)
}
