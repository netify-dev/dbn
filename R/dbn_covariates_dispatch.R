####
# covariate dispatch for non-dynamic models
#
# the dynamic path has its own joint Gibbs sampler in dbn_with_covariates.R.
# for static / piecewise / hmm / lowrank we implement the simpler shift-trick
# block-Gibbs: at each outer iteration, fix beta and run the base sampler
# on the shifted observation Z - L, then update beta from the observation
# residual Z - Theta (post-fit) using the conjugate Gaussian update.
# the result is a `dbn_covariates_fit` carrying beta + variance + operator
# draws.
####

#' @keywords internal
#' @noRd
.dbn_with_covariates_static <- function(Y, covariates, family = "gaussian",
                                          nscan = 2000L, burn = 500L,
                                          odens = 1L, symmetric = FALSE,
                                          prior_beta_scale = 2.5,
                                          prior_beta_mean = NULL,
                                          actor_effects = "none",
                                          seed = NULL, verbose = TRUE, ...) {
	.dbn_with_covariates_generic(Y, covariates,
		family = family, model = "static",
		nscan = nscan, burn = burn, odens = odens,
		symmetric = symmetric,
		prior_beta_scale = prior_beta_scale, prior_beta_mean = prior_beta_mean,
		actor_effects = actor_effects, seed = seed, verbose = verbose, ...)
}

#' @keywords internal
#' @noRd
.dbn_with_covariates_piecewise <- function(Y, covariates, family = "gaussian",
                                             blocks = NULL,
                                             nscan = 2000L, burn = 500L,
                                             odens = 1L, symmetric = FALSE,
                                             prior_beta_scale = 2.5,
                                             prior_beta_mean = NULL,
                                             actor_effects = "none",
                                             seed = NULL, verbose = TRUE, ...) {
	if (is.null(blocks))
		cli::cli_abort("{.arg blocks} is required for piecewise covariate fits.")
	.dbn_with_covariates_generic(Y, covariates,
		family = family, model = "piecewise",
		blocks = blocks,
		nscan = nscan, burn = burn, odens = odens, symmetric = symmetric,
		prior_beta_scale = prior_beta_scale, prior_beta_mean = prior_beta_mean,
		actor_effects = actor_effects, seed = seed, verbose = verbose, ...)
}

#' @keywords internal
#' @noRd
.dbn_with_covariates_hmm <- function(Y, covariates, family = "gaussian",
                                       R = 2L,
                                       nscan = 2000L, burn = 500L,
                                       odens = 1L, symmetric = FALSE,
                                       prior_beta_scale = 2.5,
                                       prior_beta_mean = NULL,
                                       actor_effects = "none",
                                       seed = NULL, verbose = TRUE, ...) {
	# the block-Gibbs driver needs `fitted(fit_inner)` to back out the
	# residual state R for the beta update. HMM fits don't implement
	# fitted() (per-regime operators applied conditional on the inferred
	# regime path); falling back to a single initial OLS beta would
	# produce a fake posterior chain. abort directly until a proper
	# Theta-reconstruction path lands.
	cli::cli_abort(c(
		"Covariate fits on {.code model = \"hmm\"} are not yet supported.",
		"x" = "HMM lacks a {.fn fitted} method, which the block-Gibbs covariate driver needs to update beta.",
		"i" = "Use {.code model = \"dynamic\"} or {.code model = \"piecewise\"} for time-varying operator covariate fits."
	))
}

#' @keywords internal
#' @noRd
.dbn_with_covariates_lowrank <- function(Y, covariates, family = "gaussian",
                                           r = 2L,
                                           nscan = 2000L, burn = 500L,
                                           odens = 1L, symmetric = FALSE,
                                           prior_beta_scale = 2.5,
                                           prior_beta_mean = NULL,
                                           actor_effects = "none",
                                           seed = NULL, verbose = TRUE, ...) {
	cli::cli_abort(c(
		"Covariate fits on {.code model = \"lowrank\"} are not yet supported.",
		"x" = "Lowrank does not expose a {.fn fitted} method that the block-Gibbs covariate driver can use to update beta.",
		"i" = "Use {.code model = \"dynamic\"} for time-varying operator covariate fits with the joint Gibbs sampler."
	))
}

####
# generic block-Gibbs shift-trick driver for non-dynamic covariate fits.
# runs the base sampler on Z - L at the current beta; updates beta from
# obs residual Z - Theta; iterates. each "outer pass" is an MCMC chain
# of length n_inner = nscan / n_outer_passes. defaults: 5 outer passes
# of 400 inner iter each = nscan = 2000.
.dbn_with_covariates_generic <- function(Y, covariates,
                                            family, model,
                                            nscan, burn, odens, symmetric,
                                            prior_beta_scale, prior_beta_mean,
                                            actor_effects = "none",
                                            seed = NULL, verbose = TRUE,
                                            n_outer = NULL,
                                            outer_burn = NULL, ...) {
	if (!is.null(seed)) set.seed(seed)
	# refuse actor_effects on the non-dynamic block-Gibbs driver: the inner
	# samplers (dbn_static, dbn_piecewise, etc.) don't accept actor-effect
	# arguments, so we can't honestly thread them through. abort directly
	# rather than silently dropping the user's argument
	if (!identical(actor_effects, "none"))
		cli::cli_abort(c(
			"{.arg actor_effects = \"{actor_effects}\"} is currently supported only for {.code model = \"dynamic\"}.",
			"i" = "On {.val {model}} the non-dynamic block-Gibbs driver does not yet thread actor effects through the inner sampler.",
			"i" = "Either refit with {.code model = \"dynamic\"}, or drop {.arg actor_effects}."
		))
	# validate Y
	dY <- dim(Y)
	if (length(dY) != 4L)
		cli::cli_abort("{.arg Y} must be a 4D array.")
	n_row <- dY[1L]; n_col <- dY[2L]; p <- dY[3L]; Tt <- dY[4L]
	if (p != 1L)
		cli::cli_abort("Covariate support currently requires {.code p = 1}.")
	# validate covariates (covariates may be NULL for actor-effects-only)
	if (is.null(covariates)) {
		covars <- list(K = 0L, n_row = n_row, n_col = n_col, validated = TRUE,
			entries = list(),
			terms = data.frame(name = character(0), stringsAsFactors = FALSE))
		class(covars) <- c("dbn_covariates", "list")
		K <- 0L
	} else {
		covars <- .dbn_validate_covariates(covariates, n_row, n_col, p, Tt,
			symmetric = symmetric)
		K <- covars$K
	}

	# initial beta from per-time pooled OLS on Z (post-NA-fill)
	pre <- shared_preprocess(Y, family = family)
	Z <- pre$Z
	if (K > 0L) {
		HtH <- matrix(0, K, K); Htc <- numeric(K)
		off_diag <- if (n_row == n_col) as.logical(diag(n_row) == 0) else
			rep(TRUE, n_row * n_col)
		for (t in seq_len(Tt)) {
			D_t <- .dbn_build_design_t(covars, t)
			z_t <- as.numeric(Z[, , 1, t])
			HtH <- HtH + crossprod(D_t[off_diag, , drop = FALSE])
			Htc <- Htc + crossprod(D_t[off_diag, , drop = FALSE], z_t[off_diag])
		}
		prior_prec <- 1 / prior_beta_scale ^ 2
		beta <- as.numeric(solve(diag(prior_prec, K) + HtH, Htc))
	} else {
		beta <- numeric(0)
	}

	build_L <- function(beta) {
		L <- array(0, dim = c(n_row, n_col, Tt))
		if (K == 0L) return(L)
		for (t in seq_len(Tt))
			L[, , t] <- .dbn_build_linear_predictor_t(covars, beta, t)
		L
	}

	# outer block-Gibbs: alternate between base sampler on Y - L and beta
	# update. each outer pass produces ONE beta draw from its full
	# conditional given the current Theta posterior. those draws form the
	# actual posterior chain for beta, with the first outer_burn passes
	# discarded as burn-in.
	# defaults: total beta chain length = nscan / odens, plus a 20% burn-in
	if (is.null(n_outer)) n_outer <- max(50L, as.integer(nscan / odens))
	if (is.null(outer_burn)) outer_burn <- max(10L, as.integer(0.2 * n_outer))
	n_outer    <- as.integer(n_outer)
	outer_burn <- as.integer(outer_burn)
	if (outer_burn >= n_outer)
		cli::cli_abort("{.arg outer_burn} must be less than {.arg n_outer}.")
	# each outer pass runs a SHORT inner chain to mix R | beta. inner_nscan
	# trades wall-clock for mixing of the residual state; small is fine
	# because R is conjugate-Gibbs sampled
	inner_nscan <- 50L
	inner_burn  <- 10L
	# storage: one beta row per outer pass; final fit's $beta is the
	# post-burn portion of this chain (NOT a single replicated value)
	beta_path <- if (K > 0L) matrix(NA_real_, n_outer, K) else NULL

	# choose the base fitter
	base_fit_fn <- switch(model,
		static    = "dbn_static",
		piecewise = "dbn_piecewise",
		hmm       = "dbn_hmm",
		lowrank   = "dbn_lowrank",
		cli::cli_abort("Unsupported model {.val {model}} in covariate generic dispatch.")
	)

	dots <- list(...)
	# strip seed from dots: the inner sampler's `.use_seed_locally` resets
	# the RNG to the same value at the start of every outer pass and
	# restores it on exit, freezing the chain. instead, generate a fresh
	# per-iter seed from the outer RNG so the inner sampler advances state
	# independently each pass while remaining reproducible overall
	dots$seed <- NULL
	final_fit <- NULL
	for (op in seq_len(n_outer)) {
		L <- build_L(beta)
		# build the shifted Y in the right family scale; for gaussian we
		# subtract directly; for binary / ordinal the standard approach is
		# to shift the *latent* Z (handled inside the sampler). this
		# generic driver only fully supports the gaussian path; non-Gaussian
		# falls back to a "subtract before fitting" approximation
		Y_shift <- Y
		for (t in seq_len(Tt)) Y_shift[, , 1, t] <- Y[, , 1, t] - L[, , t]
		# run base sampler on shifted data with a fresh seed per outer pass
		# (drawn from the outer RNG state so the overall fit is still
		# reproducible given the user's outer seed)
		inner_seed <- as.integer(.Machine$integer.max *
			stats::runif(1, -1, 1))
		base_args <- c(list(Y = Y_shift, family = family,
			nscan = inner_nscan, burn = inner_burn, odens = odens,
			verbose = FALSE, symmetric = symmetric,
			seed = inner_seed), dots)
		fit_inner <- suppressMessages(suppressWarnings(
			do.call(base_fit_fn, base_args)
		))
		# extract R-hat from the inner fit via fitted(). this returns the
		# model's posterior-mean prediction E[Y_shift | params] which is
		# exactly the residual state R (since Y_shift = R + obs_noise in
		# the shift-trick frame). use fitted() because the per-model
		# inner samplers don't all expose a top-level $Theta slot
		if (K > 0L) {
			R_post_arr <- tryCatch(
				suppressMessages(suppressWarnings(stats::fitted(fit_inner))),
				error = function(e) NULL
			)
			if (!is.null(R_post_arr) && length(dim(R_post_arr)) == 4L) {
				Z_3d <- array(Z[, , 1, ], dim = c(n_row, n_col, Tt))
				R_post <- R_post_arr[, , 1, , drop = TRUE]
				dim(R_post) <- c(n_row, n_col, Tt)
				sig2_obs <- if (!is.null(fit_inner$sigma2_obs))
					mean(fit_inner$sigma2_obs)
				else if (!is.null(fit_inner$sigma2))
					mean(fit_inner$sigma2)
				else 1
				# guard against numerically-zero residual variance
				# (which would collapse the posterior to a point mass)
				sig2_obs <- max(sig2_obs, 1e-4)
				bres <- .dbn_update_beta_obs(covars, Z_3d, R_post,
					sigma2_obs = sig2_obs,
					prior_mean = prior_beta_mean,
					prior_sd   = prior_beta_scale)
				beta <- bres$beta
			} else if (op == 1L) {
				cli::cli_warn(c(
					"Could not extract fitted values from inner {.val {model}} fit.",
					"i" = "Falling back to the initial OLS estimate; beta chain will not mix."
				))
			}
		}
		if (K > 0L) beta_path[op, ] <- beta
		final_fit <- fit_inner
	}

	# attach covariate info to the final fit
	final_fit$covariates <- covars
	# real posterior chain: discard outer_burn rows, keep the rest. each
	# row is a draw from p(beta | Theta_s, sigma2_obs_s, prior) where
	# Theta_s is the inner sampler's posterior mean at outer iteration s.
	# this is NOT marginal p(beta | data) but a block-Gibbs draw; honest
	# uncertainty quantification given the inner-sampler resolution.
	final_fit$beta <- if (K > 0L) {
		kept <- beta_path[(outer_burn + 1L):n_outer, , drop = FALSE]
		colnames(kept) <- covars$terms$name
		kept
	} else NULL
	final_fit$beta_outer_path <- beta_path
	final_fit$actor_effects <- actor_effects
	final_fit$time_varying_beta <- FALSE
	final_fit$settings <- c(final_fit$settings,
		list(outer_n = n_outer, outer_burn = outer_burn, inner_nscan = inner_nscan))
	# preserve the new class so accessors pick up the beta slot
	class(final_fit) <- c("dbn_covariates_fit", class(final_fit))
	final_fit
}
