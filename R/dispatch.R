####
# s3 method dispatch for dbn objects
####

#' S3 Methods for DBN Objects
#'
#' @description These methods dispatch to model-specific implementations
#' @name dbn-methods
#' @details
#' The following S3 methods are available for dbn objects:
#' \itemize{
#'   \item \code{print}: Display a concise summary of the dbn model
#'   \item \code{plot}: Create model-specific diagnostic plots
#'   \item \code{summary}: Generate detailed model summaries
#'   \item \code{predict}: Generate predictions or simulations
#' }
NULL

#' Plot Diagnostics for a Fitted DBN Model
#'
#' @description Creates model-specific diagnostic plots for a fitted DBN
#'   object. The type of plot depends on the model variant:
#'   \itemize{
#'     \item **Static**: trace plots, posterior histograms, and a network
#'       summary of the estimated B matrix.
#'     \item **Dynamic**: trace plots for variance parameters and, if
#'       available, time-varying A/B summaries.
#'     \item **Lowrank**: trace plots, estimated factor trajectories
#'       (\eqn{\alpha_t}), and the posterior mean node-loading matrix U.
#'     \item **HMM**: regime probabilities over time, the estimated
#'       transition matrix, and MCMC trace plots.
#'     \item **Piecewise**: trace plots for each regime block.
#'   }
#'
#' @param x A fitted `dbn` object returned by [dbn()].
#' @param ... Additional arguments passed to model-specific plot functions
#'   (e.g., `alpha` for edge significance in the static model).
#' @return A ggplot2 object (or arranged multi-panel plot) is printed and
#'   returned invisibly.
#' @seealso [dbn()], [plot_trace()], [check_convergence()], [summary.dbn()]
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' plot(fit)
#' }
#' @export
#' @method plot dbn
plot.dbn <- function(x, ...) {
	if (!inherits(x, "dbn")) cli::cli_abort("{.arg x} must be a {.cls dbn} object.")
	plot_fun <- switch(x$model,
		static = "plot_static",
		dynamic = "plot_dynamic",
		lowrank = "plot_lowrank",
		hmm = "plot_hmm",
		piecewise = "plot_piecewise",
		cli::cli_abort("Unknown model: {.val {x$model}}.")
	)
	do.call(plot_fun, list(x, ...))
}

#' Summarize a Fitted DBN Model
#'
#' @description Prints a structured summary of a fitted DBN model including
#'   data dimensions, posterior means and 95\% credible intervals for scalar
#'   parameters, and model-specific details (e.g., transition matrix for HMM,
#'   block structure for piecewise, DIC for Gaussian models).
#'
#' @param object A fitted `dbn` object returned by [dbn()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `object`. The summary is printed to the console.
#' @seealso [dbn()], [param_summary()], [check_convergence()], [plot.dbn()]
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
#' fit <- dbn(sim$Y, model = "static", nscan = 200, burn = 100, verbose = FALSE)
#' summary(fit)
#' }
#' @export
#' @method summary dbn
summary.dbn <- function(object, ...) {
	if (!inherits(object, "dbn")) cli::cli_abort("{.arg object} must be a {.cls dbn} object.")
	summary_fun <- switch(object$model,
		static = "summary_static",
		dynamic = "summary_dynamic",
		lowrank = "summary_lowrank",
		hmm = "summary_hmm",
		piecewise = "summary_piecewise",
		cli::cli_abort("Unknown model: {.val {object$model}}.")
	)
	do.call(summary_fun, list(object, ...))
}

#' Predict from a Fitted DBN Model
#'
#' @description Generates predictions from a fitted DBN model. With no
#'   forecasting arguments, returns posterior predictive replications of the
#'   observed data for model checking. With `H` (and optionally `draws` and
#'   `summary`), propagates the estimated dynamics forward to produce
#'   multi-step-ahead forecasts.
#'
#' @details
#' **Posterior predictive distribution (default):** Simulates new datasets
#' from the fitted model to compare against the observed data. Use
#' [plot_ppc_ecdf()] or [plot_ppc_density()] for visual checks.
#'
#' **Forecasting:** Pass `H` (horizon), `draws` (number of posterior draws
#' to use), and optionally `summary = "mean"` to average across draws. The
#' forecast propagates \eqn{A_t}, \eqn{B_t}, and process noise forward from
#' the final observed time point.
#'
#' @param object A fitted `dbn` object returned by [dbn()].
#' @param ... Model-specific arguments:
#'   \describe{
#'     \item{H}{Forecast horizon (number of steps ahead). Triggers
#'       forecasting mode.}
#'     \item{draws}{Number of posterior draws to use (for both forecasting
#'       and posterior predictive checks). Default depends on the entry
#'       point used downstream.}
#'     \item{summary}{If `"mean"`, return the pointwise mean across draws.
#'       If `NULL` (default), return the full array of simulated
#'       trajectories.}
#'     \item{seed}{Random seed for reproducibility.}
#'   }
#' @return For forecasting (`H` specified): an array of dimensions
#'   `[n_row, n_col, p, H]` (if `summary = "mean"`) or
#'   `[n_row, n_col, p, H, draws]`.
#'   For posterior predictive checks: a list of class `"dbn_ppd"` containing
#'   replicated datasets.
#' @seealso [dbn()], [posterior_predict_dbn()], [plot_ppc_ecdf()],
#'   [plot_ppc_density()]
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 6886)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'            nscan = 200, burn = 100, verbose = FALSE)
#'
#' # forecasting: 5 steps ahead
#' fc <- predict(fit, H = 5, draws = 50, summary = "mean")
#' dim(fc)
#'
#' # posterior predictive check
#' ppd <- predict(fit, draws = 10)
#' }
#' @export
#' @method predict dbn
predict.dbn <- function(object, ...) {
	if (!inherits(object, "dbn")) cli::cli_abort("{.arg object} must be a {.cls dbn} object.")
	if (object$model == "hmm") return(predict_hmm(object, ...))
	if (object$model == "lowrank") return(predict_lowrank(object, ...))
	dots <- list(...)
	is_forecast <- any(c("H", "summary") %in% names(dots))
	if (!is_forecast && length(dots) > 0) {
		known_ppd <- c("draws", "seed")
		unknown <- setdiff(names(dots), known_ppd)
		if (length(unknown) > 0) {
			cli::cli_warn(c(
				"Unrecognized argument{?s}: {.arg {unknown}}.",
				"i" = "For forecasting use {.arg H}, {.arg draws}, or {.arg summary}.",
				"i" = "Falling back to posterior predictive distribution."
			))
		}
	}
	if (is_forecast) {
		predict_fun <- switch(object$model,
			static = "simulate_static",
			dynamic = "simulate_dynamic",
			piecewise = "simulate_piecewise",
			cli::cli_abort("Unknown model: {.val {object$model}}.")
		)
		return(do.call(predict_fun, list(object, ...)))
	}
	ppd_args <- dots[intersect(names(dots), c("draws", "seed"))]
	fam <- tryCatch(get_family(object), error = function(e) NULL)
	if (!is.null(fam)) {
		return(do.call(posterior_predict_dbn, c(list(object), ppd_args)))
	}
	predict_fun <- switch(object$model,
		static = "simulate_static",
		dynamic = "simulate_dynamic",
		piecewise = "simulate_piecewise",
		cli::cli_abort("Unknown model: {.val {object$model}}.")
	)
	do.call(predict_fun, c(list(object), ppd_args))
}

#' Plot method for DBN objects (internal wrapper)
#'
#' @description Internal wrapper that calls [plot.dbn()]. Use `plot(fit)`
#'   directly; this bare wrapper is retained for internal package use.
#' @param x A `dbn` object
#' @param ... Additional arguments passed to the S3 method
#' @return Invisibly returns the plot object
#' @keywords internal
#' @noRd
plot_dbn <- function(x, ...) plot.dbn(x, ...)

#' Summary method for DBN objects (internal wrapper)
#'
#' @description Internal wrapper that calls [summary.dbn()]. Use
#'   `summary(fit)` directly; this bare wrapper is retained for internal
#'   package use.
#' @param object A `dbn` object
#' @param ... Additional arguments passed to the S3 method
#' @return A summary object specific to the model type
#' @keywords internal
#' @noRd
summary_dbn <- function(object, ...) summary.dbn(object, ...)

####
# print method for dbn objects
####

#' @export
#' @method print dbn
print.dbn <- function(x, ...) {
	if (!inherits(x, "dbn")) cli::cli_abort("{.arg x} must be a {.cls dbn} object.")

	# header
	cat("Dynamic Bilinear Network (DBN) Model\n")
	cat(rep("-", 40), "\n", sep = "")

	# model type
	cat("Model Type: ", toupper(x$model), "\n", sep = "")

	# distribution family
	if (!is.null(x$family)) {
		fam_name <- if (is.list(x$family) && !is.null(x$family$name)) x$family$name
					else if (is.character(x$family)) x$family
					else "unknown"
		cat("Family: ", fam_name, "\n", sep = "")
	}

	# dimensions
	if (!is.null(x$dims)) {
		cat("\nData Dimensions:\n")
		if (isTRUE(x$dims$is_bipartite)) {
			cat("  Senders: ", x$dims$n_row, "\n", sep = "")
			cat("  Receivers: ", x$dims$n_col, "\n", sep = "")
		} else {
			cat("  Nodes: ", x$dims$n_row, "\n", sep = "")
		}
		cat("  Relations: ", x$dims$p, "\n", sep = "")
		if (!is.null(x$dims$Tt)) {
			cat("  Time points: ", x$dims$Tt, "\n", sep = "")
		}
	}

	# MCMC settings
	if (!is.null(x$settings)) {
		cat("\nMCMC Settings:\n")
		cat("  Iterations: ", x$settings$nscan, "\n", sep = "")
		cat("  Burn-in: ", x$settings$burn, "\n", sep = "")
		cat("  Thinning: ", x$settings$odens, "\n", sep = "")
		cat("  Saved draws: ", x$settings$draws, "\n", sep = "")
	}

	# model-specific details
	if (x$model == "hmm" && !is.null(x$R)) {
		cat("\nRegimes: ", x$R, "\n", sep = "")
	} else if (x$model == "lowrank" && !is.null(x$rank)) {
		cat("\nRank: ", x$rank, "\n", sep = "")
	} else if (x$model == "piecewise" && !is.null(x$blocks)) {
		cat("\nBlocks: ", x$blocks$K, "\n", sep = "")
		cat("Boundaries: ", paste(x$blocks$boundaries, collapse = ", "), "\n", sep = "")
	}

	# stored components
	cat("\nModel Components:\n")
	component_names <- names(x)
	exclude <- c("model", "family", "dims", "settings", "meta", "diagnostics")
	components <- setdiff(component_names, exclude)

	for (comp in components) {
		if (!is.null(x[[comp]])) {
			cat("  ", comp, ": ", sep = "")
			if (is.list(x[[comp]])) {
				cat("list of length ", length(x[[comp]]), "\n", sep = "")
			} else if (is.array(x[[comp]]) || is.matrix(x[[comp]])) {
				dims <- dim(x[[comp]])
				cat("[", paste(dims, collapse = " x "), "]\n", sep = "")
			} else if (is.vector(x[[comp]])) {
				cat("vector of length ", length(x[[comp]]), "\n", sep = "")
			} else {
				cat(class(x[[comp]])[1], "\n", sep = "")
			}
		}
	}

	# diagnostics
	if (!is.null(x$diagnostics) && !is.null(x$diagnostics$dic)) {
		cat("\nModel Diagnostics:\n")
		cat("  DIC: ", round(x$diagnostics$dic, 2), "\n", sep = "")
	}

	invisible(x)
}
####
#' S3 Generic for Sampler Description
#'
#' @description Describe a fitted sampler (ALS, MCMC, etc.) and its uncertainty properties.
#'
#' @param object A fitted model object (dbn, dbn_boot, etc.)
#' @param verbose Logical. If TRUE, print full description.
#' @param ... Additional arguments for methods.
#'
#' @return Invisibly, information about the sampler.
#'
#' @keywords internal
#' @export
sampler_describe <- function(object, verbose = TRUE) {
	UseMethod("sampler_describe", object)
}

#' @export
sampler_describe.dbn <- function(object, verbose = TRUE) {
	sampler <- object$meta$sampler_used %||% "unknown"
	cli::cli_h3("Sampler Information")
	cli::cli_text("Type: {.val {sampler}}")
	if (!is.na(object$meta$uncertainty_available)) {
		cli::cli_text("Uncertainty available: {.val {object$meta$uncertainty_available}}")
	}
	if (!is.null(object$meta$approximation_note)) {
		cli::cli_text("Note: {object$meta$approximation_note}")
	}
	invisible(object)
}

#' @export
sampler_describe.dbn_boot <- function(object, verbose = TRUE) {
	cli::cli_h3("Bootstrap Sampler")
	cli::cli_text("Type: {.field {object$type}}")
	cli::cli_text("Valid replicates: {.val {object$n_valid}}/{.val {object$n_total}}")
	cli::cli_text("Standard errors available for A, B, M")
	invisible(object)
}

####
#' Fast Exploratory ALS Fitting for Dynamic Bilinear Networks
#'
#' @description
#' Fits a dynamic bilinear network model via alternating least squares (ALS),
#' providing a fast point estimate without MCMC. Useful for exploring data,
#' diagnosing model fit, and choosing initial values for full Bayesian inference.
#' No credible intervals are produced.
#'
#' @param data Numeric array or file path. Same format as [dbn()]:
#'   3-dimensional `[actors, actors, time]` (single relation) or
#'   4-dimensional `[actors, actors, relations, time]` (multiple relations).
#'   Diagonal entries should be `NA`.
#' @param family Character string: `"gaussian"` (default), `"ordinal"`, or `"binary"`.
#'   For ordinal/binary, ALS is applied to latent continuous scores from [shared_preprocess()].
#' @param symmetric Logical. If TRUE, enforce B = A (symmetric/undirected networks).
#' @param max_iter Maximum ALS iterations (default: 50).
#' @param tol Convergence tolerance: relative change in Frobenius norm of A (default: 1e-4).
#' @param verbose Logical. If TRUE, print iteration progress (default: TRUE).
#' @param seed Random seed for reproducibility (default: 6886).
#'
#' @return
#' A list of class `"dbn"` with structure matching a full MCMC fit, but with:
#' - Single pseudo-draw (n_keep = 1)
#' - No MCMC variance parameters (tau_A2, tau_B2 are NA)
#' - `meta$sampler_used = "als"`
#' - `meta$uncertainty_available = FALSE`
#'
#' Compatible with downstream tools: [compute_irf()], [dbn_compute_snr()],
#' [dbn_coupling_rank_probs()], etc.
#'
#' @seealso [dbn()], [compute_irf()]
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 8, time = 10, seed = 1)
#' fit_als <- dbn_explore(sim$Y, family = "gaussian")
#' sampler_describe(fit_als)
#' }
#'
#' @export
dbn_explore <- function(data,
                        family = c("gaussian", "ordinal", "binary"),
                        symmetric = FALSE,
                        max_iter = 50,
                        tol = 1e-4,
                        verbose = TRUE,
                        seed = 6886) {

	family <- match.arg(family)
	if (!is.null(seed)) set.seed(seed)

	# load data from file path or use the array directly
	if (is.character(data)) {
		cli::cli_inform("Loading data from: {.path {data}}")
		env <- new.env()
		load(data, envir = env)
		Y <- env$Y
		if (is.null(Y)) cli::cli_abort("Data file must contain object {.var Y}")
	} else {
		Y <- data
	}

	# convert 3D to 4D if needed
	if (length(dim(Y)) == 3) {
		dim_orig <- dim(Y)
		Y <- array(Y, dim = c(dim_orig[1], dim_orig[2], 1, dim_orig[3]))
	} else if (length(dim(Y)) != 4) {
		cli::cli_abort("Data must be a 3D or 4D array.")
	}

	# preprocess data
	pre <- shared_preprocess(Y, family = family)
	Z <- pre$Z  # latent continuous scores
	M <- pre$M  # baseline
	dims <- pre$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	p <- dims$p
	Tt <- dims$Tt

	if (symmetric && n_row != n_col) {
		cli::cli_abort("Symmetric model requires square network (n_row == n_col).")
	}

	if (verbose) {
		cli::cli_h3("ALS Exploration")
		cli::cli_text("Network: {n_row} × {n_col}, {p} relation{?s}, {Tt} time point{?s}")
		cli::cli_text("Tolerance: {tol}, Max iterations: {max_iter}")
	}

	# Initialize A and B: identity or small random
	A_est <- diag(n_row) + 0.01 * matrix(rnorm(n_row^2), n_row, n_row)
	if (symmetric) {
		A_est <- 0.5 * (A_est + t(A_est))  # symmetrize
		B_est <- A_est
	} else {
		B_est <- diag(n_col) + 0.01 * matrix(rnorm(n_col^2), n_col, n_col)
	}

	# Extract latent Z: from 4D [n_row, n_col, p, Tt] to list of T matrices
	# Z[, , 1, t] is the first relation at time t; use all relations averaged
	Z_list <- list()
	Z_prev_list <- list()
	for (t in 1:Tt) {
		# Average latent scores across relations (alternative: Z[, , 1, t] for first relation only)
		Z_t_all <- Z[, , , t]  # n_row x n_col x p
		Z_list[[t]] <- apply(Z_t_all, c(1, 2), mean)  # Average to get n_row x n_col matrix
		if (t > 1) {
			Z_prev_list[[t]] <- Z_list[[t-1]]
		}
	}

	# ALS iteration
	converged <- FALSE
	for (iter in seq_len(max_iter)) {

		# Update A (fix B, solve for each row of A)
		A_old <- A_est
		for (i in seq_len(n_row)) {
			# Stack equations: Z[i, j, t] = (A @ Z_{t-1})_{i,j} + noise for all t, j
			# For fixed i: Z[i, :, t] = a_i^T @ Z_{t-1} @ B^T + noise^T
			# Transpose: Z[i, :, t]^T = B @ Z_{t-1}^T @ a_i

			X_lhs <- list()
			y_rhs <- c()
			for (t in 2:Tt) {
				# (Z_{t-1} @ B^T)^T = B @ Z_{t-1}^T
				X_t_T <- t(Z_prev_list[[t]] %*% t(B_est))  # n_col x n_col
				# Stack by rows: append all n_col rows
				X_lhs[[t-1]] <- X_t_T
				# Stack corresponding RHS values (Z[i, :, t] as column)
				y_rhs <- c(y_rhs, Z_list[[t]][i, ])  # n_col elements
			}
			X_design <- do.call(rbind, X_lhs)  # (Tt-1)*n_col x n_col

			# Solve: a_i = (X^T X)^{-1} X^T y
			if (ncol(X_design) > 0 && nrow(X_design) > 0) {
				XtX <- crossprod(X_design)
				Xty <- crossprod(X_design, y_rhs)
				a_i <- tryCatch({
					solve(XtX, Xty)
				}, error = function(e) {
					# If solve fails (singular matrix), use pseudoinverse
					MASS::ginv(XtX) %*% Xty
				})
				A_est[i, ] <- as.numeric(a_i)
			}
		}

		# Update B (fix A, solve for each column of B)
		B_old <- B_est
		for (j in seq_len(n_col)) {
			# Stack equations: Z[i, j, t] = (A @ Z_{t-1} @ B^T)_{i,j} + noise for all t, i
			# For fixed j: Z[:, j, t] = A @ Z_{t-1} @ b_j + noise

			X_lhs <- list()
			y_rhs <- c()
			for (t in 2:Tt) {
				# A @ Z_{t-1} is n_row x n_col
				X_t <- A_est %*% Z_prev_list[[t]]  # n_row x n_col
				# Stack by rows: append all n_row rows
				X_lhs[[t-1]] <- X_t
				# Stack corresponding RHS values (Z[:, j, t] as column)
				y_rhs <- c(y_rhs, Z_list[[t]][, j])  # n_row elements
			}
			X_design <- do.call(rbind, X_lhs)  # (Tt-1)*n_row x n_col

			# Solve: b_j = (X^T X)^{-1} X^T y
			if (ncol(X_design) > 0 && nrow(X_design) > 0) {
				XtX <- crossprod(X_design)
				Xty <- crossprod(X_design, y_rhs)
				b_j <- tryCatch({
					solve(XtX, Xty)
				}, error = function(e) {
					# If solve fails (singular matrix), use pseudoinverse
					MASS::ginv(XtX) %*% Xty
				})
				B_est[, j] <- as.numeric(b_j)
			}
		}

		# Enforce symmetry if needed
		if (symmetric) {
			A_sym <- 0.5 * (A_est + t(A_est))
			A_est <- A_sym
			B_est <- A_sym
		}

		# Check convergence
		rel_change_A <- norm(A_est - A_old, "F") / (norm(A_old, "F") + 1e-10)
		rel_change_B <- norm(B_est - B_old, "F") / (norm(B_old, "F") + 1e-10)
		rel_change <- max(rel_change_A, rel_change_B)

		if (verbose && (iter %% 5 == 0 || iter == 1)) {
			cli::cli_text("Iter {iter}: rel_change = {sprintf('%.6f', rel_change)}")
		}

		if (rel_change < tol) {
			converged <- TRUE
			if (verbose) cli::cli_text("Converged at iteration {iter}.")
			break
		}
	}

	if (!converged && verbose) {
		cli::cli_warn("ALS did not converge after {max_iter} iterations.")
	}

	# Estimate Theta_t from current A, B, Z
	Theta_est <- array(NA, dim = c(n_row, n_col, p, Tt))
	for (t in 1:Tt) {
		Theta_t <- A_est %*% Z_list[[t]] %*% t(B_est)
		for (r in 1:p) {
			Theta_est[, , r, t] <- Theta_t
		}
	}

	# Estimate M from residuals
	M_est <- M  # keep as-is from preprocessing

	# Estimate sigma2 from residuals
	resid_sq <- 0
	n_resid <- 0
	for (t in 2:Tt) {
		pred_t <- A_est %*% Z_prev_list[[t]] %*% t(B_est)
		resid_t <- Z_list[[t]] - pred_t
		resid_sq <- resid_sq + sum(resid_t^2)
		n_resid <- n_resid + length(resid_t)
	}
	sigma2_est <- resid_sq / max(n_resid, 1)
	if (sigma2_est < 1e-6) sigma2_est <- 1.0

	# Build output object (single pseudo-draw, n_keep = 1)
	n_keep <- 1L
	time_keep <- seq_len(Tt)

	# Reshape for output: each element in list is [n_row, n_row, n_time_keep]
	A_store <- list(array(NA, dim = c(n_row, n_row, length(time_keep))))
	A_store[[1]][, , ] <- A_est

	B_store <- list(array(NA, dim = c(n_col, n_col, length(time_keep))))
	B_store[[1]][, , ] <- B_est

	M_store <- array(NA, dim = c(n_row, n_col, p, n_keep))
	M_store[, , , 1] <- M_est

	Theta_store <- array(NA, dim = c(n_row, n_col, p, length(time_keep), n_keep))
	Theta_store[, , , , 1] <- Theta_est

	sigma2_store <- c(sigma2_est)
	tau_A2_store <- c(NA_real_)
	tau_B2_store <- c(NA_real_)
	g2_store <- c(NA_real_)

	# Build output list (matching dbn_dynamic structure)
	out <- list(
		model = "dynamic",
		family = family,
		Y = Y,
		dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
					is_bipartite = (n_row != n_col), is_symmetric = symmetric),
		settings = list(
			nscan = 0,
			burn = 0,
			odens = 1,
			time_thin = 1,
			store_z = FALSE,
			draws = n_keep
		),
		meta = list(
			family = family,
			dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
						is_bipartite = (n_row != n_col), is_symmetric = symmetric),
			draws = n_keep,
			settings = list(nscan = 0, burn = 0, odens = 1, time_thin = 1),
			sampler_requested = "explore",
			sampler_used = "als",
			uncertainty_available = FALSE,
			approximation_note = "Point estimate only; no credible intervals"
		),
		params = data.frame(sigma2 = sigma2_store),
		draws = list(
			pars = data.frame(sigma2 = sigma2_store),
			misc = list(A = A_store, B = B_store, M = M_store)
		),
		time_kept = time_keep,
		A = A_store,
		B = B_store,
		M = M_store,
		Theta = Theta_store,
		sigma2 = sigma2_store,
		tau_A2 = tau_A2_store,
		tau_B2 = tau_B2_store,
		g2 = g2_store,
		n_iter = max_iter,
		burn = 0,
		thin = 1,
		time_thin = 1
	)

	class(out) <- "dbn"
	out
}

####
#' Internal ALS Refitting with Warm Start
#'
#' @description
#' ALS fitting with optional warm initialization from a prior fit.
#' Used by bootstrap to stabilize convergence on resampled data.
#'
#' @param Z Latent scores matrix (n_row x n_col)
#' @param Z_prev List of previous Z matrices for AR structure
#' @param A_init Optional initial A matrix (n_row x n_row)
#' @param B_init Optional initial B matrix (n_col x n_col)
#' @param symmetric Logical, enforce B = A
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#'
#' @return List with fitted (A, B) matrices
#' @keywords internal
als_refit_warmstart <- function(Z_list, Z_prev_list, A_init, B_init,
                                 symmetric, max_iter, tol, n_row, n_col, Tt) {
	# Use provided initializations or defaults
	if (!is.null(A_init)) {
		A_est <- A_init
	} else {
		A_est <- diag(n_row)
	}

	if (!is.null(B_init)) {
		B_est <- B_init
	} else {
		B_est <- diag(n_col)
	}

	if (symmetric) {
		B_est <- A_est
	}

	# Run ALS with early stopping (fewer iterations for warm start)
	converged <- FALSE
	for (iter in seq_len(max_iter)) {
		A_old <- A_est

		# Update A
		for (i in seq_len(n_row)) {
			X_lhs <- list()
			y_rhs <- c()
			for (t in 2:Tt) {
				X_t_T <- t(Z_prev_list[[t]] %*% t(B_est))
				X_lhs[[t-1]] <- X_t_T
				y_rhs <- c(y_rhs, Z_list[[t]][i, ])
			}
			X_design <- do.call(rbind, X_lhs)
			if (ncol(X_design) > 0 && nrow(X_design) > 0) {
				XtX <- crossprod(X_design)
				Xty <- crossprod(X_design, y_rhs)
				a_i <- tryCatch({
					solve(XtX, Xty)
				}, error = function(e) {
					# If solve fails (singular matrix), use pseudoinverse
					MASS::ginv(XtX) %*% Xty
				})
				A_est[i, ] <- as.numeric(a_i)
			}
		}

		# Update B
		B_old <- B_est
		for (j in seq_len(n_col)) {
			X_lhs <- list()
			y_rhs <- c()
			for (t in 2:Tt) {
				X_t <- A_est %*% Z_prev_list[[t]]
				X_lhs[[t-1]] <- X_t
				y_rhs <- c(y_rhs, Z_list[[t]][, j])
			}
			X_design <- do.call(rbind, X_lhs)
			if (ncol(X_design) > 0 && nrow(X_design) > 0) {
				XtX <- crossprod(X_design)
				Xty <- crossprod(X_design, y_rhs)
				b_j <- tryCatch({
					solve(XtX, Xty)
				}, error = function(e) {
					# If solve fails (singular matrix), use pseudoinverse
					MASS::ginv(XtX) %*% Xty
				})
				B_est[, j] <- as.numeric(b_j)
			}
		}

		# Enforce symmetry
		if (symmetric) {
			A_sym <- 0.5 * (A_est + t(A_est))
			A_est <- A_sym
			B_est <- A_sym
		}

		# Check convergence
		rel_change <- norm(A_est - A_old, "F") / (norm(A_old, "F") + 1e-10)
		if (rel_change < tol) {
			converged <- TRUE
			break
		}
	}

	list(A = A_est, B = B_est, converged = converged)
}

####
#' Bootstrap Inference for ALS Fits
#'
#' @description
#' Computes bootstrap standard errors and confidence intervals for ALS point
#' estimates from [dbn_explore()]. Two bootstrap strategies are available:
#'
#' - **Block bootstrap**: Resamples time periods (t indices) with replacement.
#'   Preserves within-period dependence structure. Recommended when Tt >= 10.
#' - **Parametric bootstrap**: Simulates new outcome arrays from the fitted model
#'   using estimated parameters and the specified family distribution. Better
#'   when Tt is small but the model is well-specified.
#'
#' Each replicate refits the full ALS model. Replicates that fail to converge
#' are dropped and reported. Standard errors are column standard deviations of
#' the successful replicates. Confidence intervals use the percentile method.
#'
#' @param als_fit A fitted dbn object from [dbn_explore()].
#' @param R Integer. Number of bootstrap replicates. Default is 200. Increase
#'   to 500-1000 for publication-quality intervals.
#' @param type Character. Bootstrap type: `"block"` (default) resamples time
#'   periods with replacement; `"parametric"` simulates new outcomes from the
#'   fitted model.
#' @param seed Optional integer for reproducibility. Sets the random seed
#'   before resampling.
#' @param verbose Logical. If TRUE, prints progress every 10 replicates and
#'   warnings about convergence issues.
#'
#' @return An object of class `"dbn_boot"` with components:
#' - `coefs_A`: R x (n_row * n_row) matrix of vectorized bootstrap A estimates
#' - `coefs_B`: R x (n_col * n_col) matrix of vectorized bootstrap B estimates
#' - `coefs_M`: R x (n_row * n_col * p) matrix of vectorized bootstrap M estimates
#' - `se_A`, `se_B`, `se_M`: Standard errors (column sds of valid replicates)
#' - `ci_A_lo`, `ci_A_hi`: 2.5% and 97.5% percentile bounds for A
#' - `ci_B_lo`, `ci_B_hi`: Percentile bounds for B
#' - `ci_M_lo`, `ci_M_hi`: Percentile bounds for M
#' - `point_est_A`, `point_est_B`, `point_est_M`: Original ALS estimates
#' - `n_valid`: Number of successful bootstrap replicates
#' - `n_total`: Total number of replicates attempted
#' - `type`: The bootstrap type used (`"block"` or `"parametric"`)
#' - `family`: The distribution family
#' - `symmetric`: Whether the symmetric constraint was used
#' - `dims`: Network dimensions (n_row, n_col, p, Tt)
#'
#' @seealso [dbn_explore()], [print.dbn_boot()]
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dynamic_dbn(n = 8, time = 10, seed = 1)
#' fit_als <- dbn_explore(sim$Y, family = "gaussian", verbose = FALSE)
#'
#' # Block bootstrap with 200 replicates
#' boot_result <- dbn_explore_bootstrap(fit_als, R = 200, seed = 42)
#' print(boot_result)
#' summary(boot_result)
#'
#' # Parametric bootstrap
#' boot_par <- dbn_explore_bootstrap(fit_als, R = 200, type = "parametric")
#' }
#'
#' @export
dbn_explore_bootstrap <- function(als_fit, R = 200,
                                  type = c("block", "parametric"),
                                  seed = NULL, verbose = TRUE) {

	if (!inherits(als_fit, "dbn")) {
		cli::cli_abort("{.arg als_fit} must be a dbn object from {.fun dbn_explore}.")
	}
	if (als_fit$meta$sampler_used != "als") {
		cli::cli_abort("{.arg als_fit} must be fitted with ALS (dbn_explore).")
	}

	type <- match.arg(type)
	if (!is.null(seed)) set.seed(seed)

	# Note: parametric bootstrap not recommended for asymmetric models due to scale ambiguity
	if (type == "parametric" && !als_fit$dims$is_symmetric && verbose) {
		cli::cli_warn("Parametric bootstrap for asymmetric models: scale ambiguity may cause inflated SE for B. Consider block bootstrap instead.")
	}

	Y <- als_fit$Y
	family <- als_fit$family
	symmetric <- als_fit$dims$is_symmetric
	n_row <- als_fit$dims$n_row
	n_col <- als_fit$dims$n_col
	p <- als_fit$dims$p
	Tt <- als_fit$dims$Tt

	# Original point estimates
	A_orig <- als_fit$A[[1]][, , 1]  # n_row x n_row
	B_orig <- als_fit$B[[1]][, , 1]  # n_col x n_col
	M_orig <- als_fit$M[, , , 1]     # n_row x n_col x p

	n_params_A <- n_row * n_row
	n_params_B <- n_col * n_col
	n_params_M <- n_row * n_col * p

	boot_coefs_A <- matrix(NA, R, n_params_A)
	boot_coefs_B <- matrix(NA, R, n_params_B)
	boot_coefs_M <- matrix(NA, R, n_params_M)

	if (verbose) {
		cli::cli_h3("Bootstrap for ALS Fit")
		cli::cli_text("Type: {.field {type}} | Replicates: {.val {R}}")
		cli::cli_text("Network: {n_row} × {n_col}, {p} relation{?s}, {Tt} time point{?s}")
	}

	for (b in 1:R) {
		if (verbose && b %% 10 == 0) {
			cli::cli_inform("Bootstrap {.val {b}}/{.val {R}}")
		}

		if (type == "block") {
			# Block bootstrap: resample time indices with replacement
			t_idx <- sample(1:Tt, Tt, replace = TRUE)
			Y_b <- Y[, , , t_idx, drop = FALSE]
		} else {
			# Parametric bootstrap: simulate from fitted model
			# First, extract Z for prediction
			pre <- shared_preprocess(Y, family = family)
			Z <- pre$Z
			M <- pre$M

			# Compute predicted latent state for all times
			Theta_pred <- array(NA, dim = c(n_row, n_col, p, Tt))
			Z_list <- list()
			for (t in 1:Tt) {
				Z_t_all <- Z[, , , t]
				Z_list[[t]] <- apply(Z_t_all, c(1, 2), mean)
			}
			for (t in 1:Tt) {
				Theta_t <- A_orig %*% Z_list[[t]] %*% t(B_orig)
				for (r in 1:p) {
					Theta_pred[, , r, t] <- Theta_t
				}
			}

			# Simulate new Y from fitted model
			if (family == "gaussian") {
				sigma <- sqrt(als_fit$sigma2[1])
				Y_b <- array(rnorm(length(Theta_pred), mean = Theta_pred, sd = sigma),
							 dim = dim(Y))
			} else if (family == "ordinal") {
				# Ordinal: simulate latent Z, then discretize
				Y_b <- array(NA, dim = dim(Y))
				for (t in 1:Tt) {
					for (r in 1:p) {
						Z_sim <- rnorm(n_row * n_col, mean = c(Theta_pred[, , r, t]), sd = sqrt(als_fit$sigma2[1]))
						Y_b[, , r, t] <- matrix(pmin(pmax(round(Z_sim), 1), 5), n_row, n_col)
					}
				}
			} else if (family == "binary") {
				# Binary: Bernoulli from logit link
				Y_b <- array(NA, dim = dim(Y))
				for (t in 1:Tt) {
					for (r in 1:p) {
						prob <- 1 / (1 + exp(-Theta_pred[, , r, t]))
						Y_b[, , r, t] <- matrix(rbinom(n_row * n_col, 1, prob), n_row, n_col)
					}
				}
			}
			# Preserve diagonal NA structure
			for (tt in 1:Tt) {
				for (rr in 1:p) {
					diag(Y_b[, , rr, tt]) <- NA
				}
			}
		}

		# Refit ALS on bootstrap sample using warm-start from original solution
		# This prevents converging to different modes/local optima
		tryCatch({
			# Preprocess bootstrap data
			pre_b <- shared_preprocess(Y_b, family = family)
			Z_b <- pre_b$Z
			M_b <- pre_b$M

			# Extract latent Z
			Z_list_b <- list()
			Z_prev_list_b <- list()
			for (t in 1:Tt) {
				Z_t_all <- Z_b[, , , t]
				Z_list_b[[t]] <- apply(Z_t_all, c(1, 2), mean)
				if (t > 1) {
					Z_prev_list_b[[t]] <- Z_list_b[[t-1]]
				}
			}

			# Warm-start ALS refinement from original estimates
			als_result <- als_refit_warmstart(
				Z_list_b, Z_prev_list_b,
				A_init = A_orig, B_init = B_orig,  # Use original as starting point
				symmetric = symmetric,
				max_iter = 30,  # Fewer iterations needed with warm start
				tol = 1e-5,
				n_row = n_row, n_col = n_col, Tt = Tt
			)

			A_boot <- als_result$A
			B_boot <- als_result$B
			M_boot <- M_b

			# Vectorize and store
			boot_coefs_A[b, ] <- c(A_boot)
			boot_coefs_B[b, ] <- c(B_boot)
			boot_coefs_M[b, ] <- c(M_boot)

		}, error = function(e) {
			# Leave as NA for failed replicates
		})
	}

	# Compute statistics from successful replicates
	valid_A <- apply(boot_coefs_A, 1, function(r) !any(is.na(r)))
	valid_B <- apply(boot_coefs_B, 1, function(r) !any(is.na(r)))
	valid_M <- apply(boot_coefs_M, 1, function(r) !any(is.na(r)))

	# Use intersection of all valid replicates
	valid <- valid_A & valid_B & valid_M
	n_valid <- sum(valid)

	# Quality control: flag replicates that are outliers (likely different local optima)
	# by checking Frobenius norm of A relative to others
	if (n_valid > 3) {
		A_norms <- sqrt(rowSums(boot_coefs_A[valid, ]^2))
		A_norm_median <- median(A_norms)
		A_norm_mad <- mad(A_norms)
		# Flag replicates with very different norms (>3 MAD from median)
		outlier_threshold <- A_norm_median + 3 * A_norm_mad
		outliers <- which(A_norms > outlier_threshold)

		if (length(outliers) > 0 && verbose) {
			n_outliers <- length(outliers)
			cli::cli_warn("Bootstrap: {n_outliers} replicate{?s} with large coefficient norms (possible local optima). Consider larger R or MCMC.")
		}
	}

	if (n_valid < 10 && verbose) {
		cli::cli_warn("Only {.val {n_valid}} valid bootstrap replicates (of {.val {R}}). Results unreliable.")
	}

	# Compute standard errors and CIs
	boot_se_A <- apply(boot_coefs_A[valid, , drop = FALSE], 2, sd)
	boot_se_B <- apply(boot_coefs_B[valid, , drop = FALSE], 2, sd)
	boot_se_M <- apply(boot_coefs_M[valid, , drop = FALSE], 2, sd)

	boot_ci_A <- apply(boot_coefs_A[valid, , drop = FALSE], 2,
					   quantile, probs = c(0.025, 0.975))
	boot_ci_B <- apply(boot_coefs_B[valid, , drop = FALSE], 2,
					   quantile, probs = c(0.025, 0.975))
	boot_ci_M <- apply(boot_coefs_M[valid, , drop = FALSE], 2,
					   quantile, probs = c(0.025, 0.975))

	result <- list(
		coefs_A = boot_coefs_A,
		coefs_B = boot_coefs_B,
		coefs_M = boot_coefs_M,
		se_A = boot_se_A,
		se_B = boot_se_B,
		se_M = boot_se_M,
		ci_A_lo = boot_ci_A[1, ],
		ci_A_hi = boot_ci_A[2, ],
		ci_B_lo = boot_ci_B[1, ],
		ci_B_hi = boot_ci_B[2, ],
		ci_M_lo = boot_ci_M[1, ],
		ci_M_hi = boot_ci_M[2, ],
		point_est_A = c(A_orig),
		point_est_B = c(B_orig),
		point_est_M = c(M_orig),
		n_valid = n_valid,
		n_total = R,
		type = type,
		family = family,
		symmetric = symmetric,
		dims = list(n_row = n_row, n_col = n_col, p = p, Tt = Tt)
	)

	class(result) <- "dbn_boot"
	result
}

####
#' Print Bootstrap Results for ALS
#'
#' @param x A `dbn_boot` object from [dbn_explore_bootstrap()].
#' @param digits Number of significant digits.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the `dbn_boot` object.
#' @export
print.dbn_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cli::cli_text("\n")
	cli::cli_h3("Bootstrap Results for ALS Fit")
	cli::cli_text("Type: {.field {x$type}} | Replicates: {.val {x$n_valid}}/{.val {x$n_total}} valid")
	cli::cli_text("Network: {x$dims$n_row} × {x$dims$n_col}, {x$dims$p} relation{?s}, {x$dims$Tt} time point{?s}")
	cli::cli_text("Family: {.field {x$family}} {ifelse(x$symmetric, '(symmetric)', '')}")

	cli::cli_rule()

	# Summarize A matrix
	cli::cli_text("{.strong Operator A ({x$dims$n_row} × {x$dims$n_row})}")
	A_mean_se <- mean(x$se_A)
	A_mean_est <- mean(x$point_est_A)
	cli::cli_text("  Mean estimate: {sprintf('%.4f', A_mean_est)}, Mean SE: {sprintf('%.4f', A_mean_se)}")
	cli::cli_text("  Diagonal mean: {sprintf('%.4f', mean(matrix(x$point_est_A, x$dims$n_row, x$dims$n_row)[cbind(1:x$dims$n_row, 1:x$dims$n_row)]))}")

	# Summarize B matrix
	cli::cli_text("{.strong Operator B ({x$dims$n_col} × {x$dims$n_col})}")
	B_mean_se <- mean(x$se_B)
	B_mean_est <- mean(x$point_est_B)
	cli::cli_text("  Mean estimate: {sprintf('%.4f', B_mean_est)}, Mean SE: {sprintf('%.4f', B_mean_se)}")
	cli::cli_text("  Diagonal mean: {sprintf('%.4f', mean(matrix(x$point_est_B, x$dims$n_col, x$dims$n_col)[cbind(1:x$dims$n_col, 1:x$dims$n_col)]))}")

	# Summarize M baseline
	cli::cli_text("{.strong Baseline M ({x$dims$n_row} × {x$dims$n_col} × {x$dims$p})}")
	M_mean_se <- mean(x$se_M)
	M_mean_est <- mean(x$point_est_M)
	cli::cli_text("  Mean estimate: {sprintf('%.4f', M_mean_est)}, Mean SE: {sprintf('%.4f', M_mean_se)}")

	cli::cli_text("")
	invisible(x)
}

####
#' Summary of Bootstrap Results for ALS
#'
#' @param object A `dbn_boot` object from [dbn_explore_bootstrap()].
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the `dbn_boot` object.
#' @export
summary.dbn_boot <- function(object, ...) {
	cli::cli_h1("Bootstrap Results for ALS Fit")
	cli::cli_text("Type: {.field {object$type}} | Family: {.field {object$family}}")
	cli::cli_text("Replicates: {.val {object$n_valid}} valid of {.val {object$n_total}} total ({.val {round(100 * object$n_valid / object$n_total, 1)}}%)")
	cli::cli_text("Network: {object$dims$n_row} × {object$dims$n_col}, {object$dims$p} relation{?s}, {object$dims$Tt} time point{?s}")

	cli::cli_rule()

	# A matrix summary
	cli::cli_text("{.strong Operator A}")
	A_mat <- matrix(object$point_est_A, object$dims$n_row, object$dims$n_row)
	A_se_mat <- matrix(object$se_A, object$dims$n_row, object$dims$n_row)
	cli::cli_text("Point estimate (diagonal): {paste(sprintf('%.4f', diag(A_mat)), collapse=', ')}")
	cli::cli_text("Bootstrap SE (diagonal): {paste(sprintf('%.4f', diag(A_se_mat)), collapse=', ')}")

	# B matrix summary
	cli::cli_text("{.strong Operator B}")
	B_mat <- matrix(object$point_est_B, object$dims$n_col, object$dims$n_col)
	B_se_mat <- matrix(object$se_B, object$dims$n_col, object$dims$n_col)
	cli::cli_text("Point estimate (diagonal): {paste(sprintf('%.4f', diag(B_mat)), collapse=', ')}")
	cli::cli_text("Bootstrap SE (diagonal): {paste(sprintf('%.4f', diag(B_se_mat)), collapse=', ')}")

	# Coverage check: what percentage of params have CIs that include the point est
	A_covers <- (object$ci_A_lo <= object$point_est_A) & (object$point_est_A <= object$ci_A_hi)
	B_covers <- (object$ci_B_lo <= object$point_est_B) & (object$point_est_B <= object$ci_B_hi)
	M_covers <- (object$ci_M_lo <= object$point_est_M) & (object$point_est_M <= object$ci_M_hi)

	cli::cli_rule()
	cli::cli_text("{.strong CI Coverage}")
	cli::cli_text("A: {sum(A_covers)}/{length(A_covers)} params have point est in 95% CI")
	cli::cli_text("B: {sum(B_covers)}/{length(B_covers)} params have point est in 95% CI")
	cli::cli_text("M: {sum(M_covers)}/{length(M_covers)} params have point est in 95% CI")

	invisible(object)
}
