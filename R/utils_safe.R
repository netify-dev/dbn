####
#' Safe Statistical Utilities
#'
#' @description Numerically stable implementations of common statistical operations
#' @name utils_safe
#' @keywords internal
NULL
####

####
#' Expand a bootstrap ALS fit into per-draw A/B/M lists
#'
#' @description An ALS fit produced by `dbn_als(..., bootstrap = R)` carries
#'   one pseudo-draw in `fit$A`, `fit$B`, `fit$M` (the point estimate), and
#'   the bootstrap A/B/M coefficients in `fit$bootstrap$coefs_*`. Downstream
#'   accessors (`compute_irf`, `dbn_compute_snr`, `dbn_coupling_rank_probs`,
#'   forecast paths) iterate over draws via `fit$A[[m]]`, etc. This helper
#'   reshapes the bootstrap coefficients into per-draw arrays so those
#'   functions can treat a bootstrap fit as if it had genuine posterior draws.
#'
#'   Static (2D bootstrap) fits tile each bootstrap A_b / B_b across the
#'   time dimension; TV-ALS (4D bootstrap) keeps per-time A_t / B_t draws
#'   directly. Bootstrap does not resample residual or innovation variances:
#'   `sigma2` is kept as the length-1 point estimate and `tau_A2`, `tau_B2`,
#'   `g2` are set to `NA_real_` to flag that no replicate-level variation
#'   exists for them. Per-draw consumers must index `fit$sigma2[1]` (not
#'   `fit$sigma2[m]`) and treat `tau_A2` / `tau_B2` / `g2` as absent. M is
#'   resampled in the static path but not in TV-ALS, so the static-boot
#'   `M_arr` has trailing dimension `R_valid` while the TV-boot `M_arr` has
#'   trailing dimension 1.
#'
#'   The returned object has `meta$uncertainty_available = TRUE`,
#'   `meta$uncertainty_source = "bootstrap"`, `meta$n_bootstrap_draws`
#'   recording the number of valid replicates used, and
#'   `meta$bootstrap_expanded = TRUE` so downstream code can branch on it.
#'   The original `fit$bootstrap` is preserved on the result.
#'
#' @param fit A `dbn` object with `meta$sampler_used = "als"` and a
#'   non-null `fit$bootstrap` of class `dbn_boot`.
#' @return A `dbn` object whose `A`, `B`, `M` are lists/arrays of length
#'   `n_bootstrap_draws` (M is length-1 for TV-boot), `sigma2` is the
#'   length-1 point estimate, and `tau_A2` / `tau_B2` / `g2` are
#'   `NA_real_`. If the input is not a bootstrap fit, returns it unchanged.
#' @keywords internal
#' @noRd
.bootstrap_expand <- function(fit) {
	if (is.null(fit$bootstrap) || !(fit$meta$sampler_used %in% c("als", "als_tv")) ||
		!inherits(fit$bootstrap, "dbn_boot")) {
		return(fit)
	}
	boot <- fit$bootstrap
	n_row <- fit$dims$n_row %||% boot$dims$n_row
	n_col <- fit$dims$n_col %||% boot$dims$n_col
	p     <- fit$dims$p     %||% boot$dims$p
	Tt    <- fit$dims$Tt    %||% boot$dims$Tt

	# coefs_A may be:
	#   2D [R, n_row*n_row]: static ALS bootstrap (each rep has time-invariant A)
	#   4D [R, n_row, n_row, S=Tt-1]: TV-ALS bootstrap (per-time draws of A_t)
	is_tv_boot <- length(dim(boot$coefs_A)) == 4L

	if (is_tv_boot) {
		# 4D: per-time A draws, no need to tile
		Rtotal <- dim(boot$coefs_A)[1]
		valid <- vapply(seq_len(Rtotal), function(r)
			!any(is.na(boot$coefs_A[r, , , ])) &&
			!any(is.na(boot$coefs_B[r, , , ])),
			logical(1))
		# coefs_M may not exist in TV bootstrap; fall back to point estimate
		idx <- which(valid)
		R_valid <- length(idx)
		if (R_valid < 2L) cli::cli_warn(c(
			"Bootstrap-expanded fit has fewer than 2 valid replicates ({R_valid}).",
			"i" = "Raise {.arg bootstrap}."))

		# TV bootstrap doesn't resample M -- only A_t and B_t. keep M as a
		# length-1 point estimate (not tiled across replicates) so the
		# downstream uncertainty over M is honestly zero rather than
		# silently-zero. M_arr below is shape [n_row, n_col, p, 1].
		M_point <- if (is.array(fit$M) && length(dim(fit$M)) >= 3L)
			fit$M[, , , 1, drop = FALSE] else fit$M
		dim(M_point) <- c(n_row, n_col, p)
		M_arr <- array(M_point, dim = c(n_row, n_col, p, 1L))

		A_list <- vector("list", R_valid)
		B_list <- vector("list", R_valid)
		for (i in seq_len(R_valid)) {
			b <- idx[i]
			# reconstruct full Tt slices: slice 1 = identity, slices 2..Tt from bootstrap
			A_full <- array(0, dim = c(n_row, n_row, Tt))
			A_full[, , 1] <- diag(n_row)
			A_full[, , 2:Tt] <- boot$coefs_A[b, , , ]
			B_full <- array(0, dim = c(n_col, n_col, Tt))
			B_full[, , 1] <- diag(n_col)
			B_full[, , 2:Tt] <- boot$coefs_B[b, , , ]
			A_list[[i]] <- A_full
			B_list[[i]] <- B_full
		}
	} else {
		# 2D static-bootstrap format
		valid <- !apply(boot$coefs_A, 1, function(r) any(is.na(r))) &
		         !apply(boot$coefs_B, 1, function(r) any(is.na(r))) &
		         !apply(boot$coefs_M, 1, function(r) any(is.na(r)))
		idx <- which(valid)
		R_valid <- length(idx)
		if (R_valid < 2L) cli::cli_warn(c(
			"Bootstrap-expanded fit has fewer than 2 valid replicates ({R_valid}).",
			"i" = "Raise {.arg bootstrap}."))

		A_list <- vector("list", R_valid)
		B_list <- vector("list", R_valid)
		M_arr  <- array(0, dim = c(n_row, n_col, p, R_valid))
		for (i in seq_len(R_valid)) {
			b <- idx[i]
			A_b <- matrix(boot$coefs_A[b, ], n_row, n_row)
			B_b <- matrix(boot$coefs_B[b, ], n_col, n_col)
			M_b <- array(boot$coefs_M[b, ], dim = c(n_row, n_col, p))
			A_list[[i]] <- array(rep(c(A_b), Tt), dim = c(n_row, n_row, Tt))
			B_list[[i]] <- array(rep(c(B_b), Tt), dim = c(n_col, n_col, Tt))
			M_arr[, , , i] <- M_b
		}
	}

	fit$A <- A_list
	fit$B <- B_list
	fit$M <- M_arr
	# bootstrap gives draws over (A, B, M) but NOT over the scalar variance
	# parameters. keep sigma2 as the length-1 point estimate (not a fake
	# replicated chain); downstream consumers that index `fit$sigma2[m]`
	# must use `fit$sigma2[1]` instead. tau_A2 / tau_B2 / g2 are likewise
	# NA -- the bootstrap has no opinion on them.
	if (length(fit$sigma2) > 1L) fit$sigma2 <- fit$sigma2[1]
	fit$tau_A2 <- NA_real_
	fit$tau_B2 <- NA_real_
	fit$g2     <- NA_real_
	# recompute Theta per replicate using the per-rep (A, B, M). this is
	# the fixed-design behaviour: the lag state Phi_{t-1} is held at the
	# latent point estimate (or observed Y for gaussian), and uncertainty
	# comes from variation in (A_t, B_t, M). for binary/ordinal families
	# fit$Y is on the {0,1}/integer scale but the bilinear dynamics live
	# on the latent Z scale, so the non-gaussian path uses fit$Theta as
	# the lag state.
	if (!is.null(fit$Y) && !is.null(fit$Theta) && length(dim(fit$Theta)) == 5L) {
		family_name <- if (is.list(fit$family)) fit$family$name else fit$family
		Theta_pt <- fit$Theta[, , , , 1, drop = FALSE]
		dim(Theta_pt) <- c(n_row, n_col, p, Tt)
		# gaussian keeps the Y-based lag (Y == Z up to noise); non-gaussian
		# uses the latent Theta from the point fit.
		use_theta_lag <- !identical(family_name, "gaussian")
		Y_filled <- fit$Y; Y_filled[is.na(Y_filled)] <- 0
		Theta_arr <- array(NA_real_, dim = c(n_row, n_col, p, Tt, R_valid))
		for (i in seq_len(R_valid)) {
			A_i <- A_list[[i]]
			B_i <- B_list[[i]]
			# M_arr may be replicate-indexed (static-boot path) or length-1
			# (TV-boot path, where M is not resampled). pick the right slice
			M_idx <- if (dim(M_arr)[4] >= i) i else 1L
			for (r in seq_len(p)) {
				M_ir <- M_arr[, , r, M_idx]
				if (use_theta_lag) {
					# latent-scale lag from the ALS point-estimate Theta
					Theta_arr[, , r, 1, i] <- Theta_pt[, , r, 1]
					for (t in 2:Tt) {
						Phi_prev <- Theta_pt[, , r, t - 1L] - M_ir
						Theta_arr[, , r, t, i] <- M_ir +
							A_i[, , t] %*% Phi_prev %*% t(B_i[, , t])
					}
				} else {
					Theta_arr[, , r, 1, i] <- Y_filled[, , r, 1]
					for (t in 2:Tt) {
						Phi_prev <- Y_filled[, , r, t - 1L] - M_ir
						Theta_arr[, , r, t, i] <- M_ir +
							A_i[, , t] %*% Phi_prev %*% t(B_i[, , t])
					}
				}
			}
		}
		fit$Theta <- Theta_arr
	}

	fit$meta$uncertainty_available <- TRUE
	fit$meta$uncertainty_source <- "bootstrap"
	fit$meta$n_bootstrap_draws <- R_valid
	fit$meta$draws <- R_valid
	# advertise the expansion in fit metadata so downstream code can branch
	# on it (e.g. plotting routines that want to label intervals "bootstrap"
	# rather than "credible")
	fit$meta$bootstrap_expanded <- TRUE
	fit
}
####

####
#' Convert a piecewise-DBN fit into a dynamic-shaped fit for accessors
#'
#' Builds per-draw `[n, n, Tt]` operator cubes by tiling each block's
#' single `A_k` / `B_k` matrix across the time indices that belong to
#' block k. The result has `model = "dynamic"` plus everything else the
#' dynamic-path IRF / forecast code expects, so accessors that work on
#' dynamic fits also work on piecewise via this expansion.
#'
#' @param fit A `dbn` fit with `model = "piecewise"`.
#' @return A `dbn` object with `model = "dynamic"`, ready to feed into
#'   accessors that work for dynamic fits.
#' @keywords internal
#' @noRd
.piecewise_expand_to_dynamic <- function(fit) {
	if (!identical(fit$model, "piecewise")) return(fit)
	dims <- fit$dims
	n_row <- dims$n_row; n_col <- dims$n_col; p <- dims$p %||% 1L; Tt <- dims$Tt
	# t -> k mapping from block_info
	block_indices <- fit$blocks$block_indices
	K <- length(block_indices)
	t2k <- integer(Tt)
	for (k in seq_len(K)) t2k[block_indices[[k]]] <- k
	# per-draw A / B lists in piecewise: fit$draws$misc$A is a list (per draw)
	# of K-element lists of [n_row, n_row] matrices. Tile each into [n,n,Tt].
	A_draws <- fit$draws$misc$A
	B_draws <- fit$draws$misc$B
	if (is.null(A_draws) || is.null(B_draws)) {
		cli::cli_abort(c(
			"Piecewise fit is missing per-draw operator storage.",
			"i" = "Refit without {.code keep = ...} stripping {.field A} / {.field B}."
		))
	}
	n_draws <- length(A_draws)
	A_list <- vector("list", n_draws)
	B_list <- vector("list", n_draws)
	for (s in seq_len(n_draws)) {
		A_cube <- array(0, dim = c(n_row, n_row, Tt))
		B_cube <- array(0, dim = c(n_col, n_col, Tt))
		for (t in seq_len(Tt)) {
			k <- t2k[t]
			A_cube[, , t] <- A_draws[[s]][[k]]
			B_cube[, , t] <- B_draws[[s]][[k]]
		}
		A_list[[s]] <- A_cube
		B_list[[s]] <- B_cube
	}
	fit$A <- A_list
	fit$B <- B_list
	# tile M across draws so downstream code that indexes fit$M[, , , s]
	# (e.g. compute_irf's extract_baseline) sees the right shape. piecewise
	# stores per-draw M under fit$draws$misc$M; assemble it into the
	# dynamic-layout 4D array [n_row, n_col, p, n_draws].
	M_draws <- fit$draws$misc$M
	if (is.list(M_draws) && length(M_draws) == n_draws) {
		M_arr <- array(NA_real_, dim = c(n_row, n_col, p, n_draws))
		for (s in seq_len(n_draws)) {
			M_s <- M_draws[[s]]
			# M_s may be [n_row, n_col, p] or [n_row, n_col] (p = 1 with drop)
			if (length(dim(M_s)) == 2L) {
				M_arr[, , 1L, s] <- M_s
			} else {
				M_arr[, , , s] <- M_s
			}
		}
		fit$M <- M_arr
	} else if (is.array(fit$M) && length(dim(fit$M)) == 3L) {
		# fall back: replicate the point estimate across draws
		M_arr <- array(NA_real_, dim = c(n_row, n_col, p, n_draws))
		for (s in seq_len(n_draws)) M_arr[, , , s] <- fit$M
		fit$M <- M_arr
	}
	# tile Theta across draws if a single point-estimate Theta is stored
	if (!is.null(fit$Theta) && is.array(fit$Theta)) {
		d_th <- dim(fit$Theta)
		if (length(d_th) == 5L && d_th[5] == 1L && n_draws > 1L) {
			Theta_tiled <- array(NA_real_, dim = c(d_th[1:4], n_draws))
			for (s in seq_len(n_draws)) Theta_tiled[, , , , s] <- fit$Theta[, , , , 1L]
			fit$Theta <- Theta_tiled
		}
	}
	# variance scalars: piecewise stores them under draws$pars (s2, t2, g2).
	# map onto the canonical sigma2 name and synthesise tau_A2/tau_B2 from
	# t2 so the downstream code that indexes fit$sigma2[s] sees a vector.
	if (!is.null(fit$draws$pars)) {
		pars <- fit$draws$pars
		if (!is.null(pars$s2)) fit$sigma2 <- as.numeric(pars$s2)
		if (!is.null(pars$t2)) {
			fit$tau_A2 <- as.numeric(pars$t2)
			fit$tau_B2 <- as.numeric(pars$t2)
		}
		if (!is.null(pars$g2)) fit$g2 <- as.numeric(pars$g2)
	}
	fit$model <- "dynamic"
	fit$meta$piecewise_expanded <- TRUE
	fit$meta$piecewise_K <- K
	fit$meta$piecewise_t2k <- t2k
	fit
}
####

####
#' Safe Inverse-Gamma Sampling
#'
#' @description Draws from inverse-Gamma on log-scale to avoid under/overflow
#' @param shape Shape parameter
#' @param rate Rate parameter
#' @param floor Minimum return value
#' @param ceiling Maximum return value
#' @return Sample from inverse-gamma distribution
#' @keywords internal
#' Seed the RNG locally without leaking the global state
#'
#' @description `set.seed(seed)` resets the global `.Random.seed` and leaves
#'   it changed when the function returns -- the caller's RNG stream is
#'   silently corrupted. This helper snapshots `.Random.seed`, calls
#'   `set.seed(seed)`, and returns a thunk the caller invokes via `on.exit()`
#'   to restore the original state.
#' @param seed Integer or NULL. If NULL, `.Random.seed` is left untouched
#'   (and the returned thunk is a no-op).
#' @return A zero-argument function suitable for `on.exit(<>(), add = TRUE)`.
#' @keywords internal
#' @noRd
.use_seed_locally <- function(seed) {
	# seed = NULL: inherit the caller's RNG state and advance it as a side
	# effect, so successive calls without an explicit seed produce
	# monte-Carlo variability.
	if (is.null(seed)) return(function() invisible(NULL))
	old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
		get(".Random.seed", envir = globalenv())
	} else NULL
	set.seed(seed)
	function() {
		if (is.null(old)) {
			suppressWarnings(rm(".Random.seed", envir = globalenv()))
		} else {
			assign(".Random.seed", old, envir = globalenv())
		}
	}
}
####

####
safe_rinv_gamma <- function(shape, rate, floor = 1e-8, ceiling = 1e8) {
	shape <- max(shape, 1e-10)
	rate <- max(rate, 1e-10)

	log_x <- (log(rate) - log(stats::rgamma(1, shape = shape, rate = 1)))
	x <- exp(log_x)
	x <- min(max(x, floor), ceiling)
	return(x)
}
####

####
#' Warn when posterior variance chains are pinned at the safe-IG ceiling
#'
#' @description `safe_rinv_gamma()` silently clamps draws to `[floor, ceiling]`
#'   to keep the Gibbs sweep numerically alive. When a variance posterior
#'   actually wants a value beyond the ceiling -- typically because the data
#'   were never rescaled toward unit variance, or the operator is genuinely
#'   explosive -- every draw is clamped and the chain becomes a degenerate
#'   spike. The reported posterior summary then looks like a valid (huge)
#'   estimate with zero uncertainty. This helper inspects the relevant
#'   variance traces of a finished fit and warns if any have saturated.
#'
#' @param fit A fitted `dbn` object.
#' @param ceiling The same ceiling used by `safe_rinv_gamma()` (1e8).
#' @return Invisibly, a character vector of saturated parameter names (or
#'   `character(0)` if none).
#' @keywords internal
#' @noRd
.warn_variance_saturation <- function(fit, ceiling = 1e8) {
	candidates <- c("sigma2", "sigma2_obs", "sigma2_proc",
		"tau_A2", "tau_B2", "tauA2", "tauB2", "tau_alpha2", "g2")
	# classify a single variance trace: "full" (every draw at the ceiling),
	# "partial" (a non-trivial fraction at the ceiling, so the mean / upper CI
	# is dominated by clamped draws), or "" (fine). partial saturation is the
	# insidious case -- a handful of clamped draws push the posterior mean into
	# the millions while the median looks sane.
	classify <- function(v) {
		if (is.null(v) || !is.numeric(v) || length(v) < 2L) return("")
		if (!all(is.finite(v))) return("")
		at_ceiling <- v >= 0.99 * ceiling
		frac <- mean(at_ceiling)
		if (frac >= 0.999 && sd(v) < 1e-6 * ceiling) return("full")
		med <- stats::median(v)
		# partial: >= 5% of draws clamped, or a heavy right tail where the max
		# dwarfs the median (max/median > 1e4 with a finite positive median)
		heavy_tail <- is.finite(med) && med > 0 && (max(v) / med) > 1e4
		if (frac >= 0.05 || (any(at_ceiling) && heavy_tail)) return("partial")
		""
	}
	full <- character(0); partial <- character(0)
	scan_one <- function(nm, v) {
		cls <- classify(v)
		if (cls == "full") full <<- c(full, nm)
		else if (cls == "partial") partial <<- c(partial, nm)
	}
	for (nm in candidates) scan_one(nm, fit[[nm]])
	# static/piecewise hide their scalar traces under params (data.frame or matrix)
	if (!is.null(fit$params) && (is.data.frame(fit$params) || is.matrix(fit$params))) {
		seen <- c(full, partial)
		for (nm in setdiff(colnames(fit$params), seen)) {
			v <- if (is.data.frame(fit$params)) fit$params[[nm]] else fit$params[, nm]
			scan_one(nm, v)
		}
	}
	if (length(full) > 0) {
		label <- if (length(full) == 1L) "parameter" else "parameters"
		cli::cli_warn(c(
			"Variance {label} {.val {full}} saturated the sampler ceiling for every kept draw.",
			"i" = "The posterior chain is a degenerate spike, not a genuine summary -- the reported credible interval has zero width but is meaningless.",
			"i" = "Rescale the input data toward unit variance, check for explosive dynamics with {.fun dbn_operator}, or shorten the chain to inspect early divergence."
		))
	}
	if (length(partial) > 0) {
		label <- if (length(partial) == 1L) "parameter" else "parameters"
		cli::cli_warn(c(
			"Variance {label} {.val {partial}} hit the sampler ceiling on a subset of draws.",
			"i" = "The posterior mean and upper credible bound are dominated by clamped draws and are not trustworthy (the median may still look reasonable).",
			"i" = "Rescale the input data toward unit variance, or check for explosive dynamics with {.fun dbn_operator}."
		))
	}
	invisible(c(full, partial))
}
####

####
#' Auto-warn when the posterior-mean lag operator is locally explosive
#'
#' @description The random-walk prior does not enforce contractivity, so the
#'   per-period spectral radius rho(A_t) * rho(B_t) can exceed 1 in posterior
#'   draws. Users routinely miss this until they look at an exploded IRF. This
#'   helper inspects the posterior-mean operators on a finished dynamic fit
#'   and warns if a substantial fraction of periods are locally explosive,
#'   pointing the user at `dbn_operator()` for the full profile.
#' @param fit A fitted `dbn` object.
#' @param threshold Fraction of explosive periods that triggers the warning
#'   (default 0.5).
#' @return Invisibly, the per-period gain vector (or NULL if not applicable).
#' @keywords internal
#' @noRd
.warn_operator_stability <- function(fit, threshold = 0.5) {
	if (!identical(fit$model, "dynamic")) return(invisible(NULL))
	if (!is.list(fit$A) || !is.list(fit$B)) return(invisible(NULL))
	if (length(fit$A) < 2L) return(invisible(NULL))  # ALS / point estimate
	if (isTRUE(fit$dims$is_bipartite)) return(invisible(NULL))  # Kronecker case
	d <- dim(fit$A[[1]])
	if (length(d) != 3L || d[3] < 2L) return(invisible(NULL))
	Tt <- d[3]
	# posterior-mean operators (cheap; spectral radius of a mean is a sensible
	# first-pass diagnostic even though it does not equal the mean radius)
	Amean <- Reduce(`+`, fit$A) / length(fit$A)
	Bmean <- Reduce(`+`, fit$B) / length(fit$B)
	gains <- numeric(Tt)
	n_failed <- 0L
	for (t in seq_len(Tt)) {
		# `eigen()` is allowed to fail (e.g. when an operator slice has
		# non-finite entries because of a divergent draw). track failures
		# explicitly so the downstream warning logic isn't silently NA-poisoned
		ra <- tryCatch(max(Mod(eigen(Amean[, , t], only.values = TRUE)$values)),
			error = function(e) NA_real_)
		rb <- tryCatch(max(Mod(eigen(Bmean[, , t], only.values = TRUE)$values)),
			error = function(e) NA_real_)
		if (is.na(ra) || is.na(rb)) n_failed <- n_failed + 1L
		gains[t] <- ra * rb
	}
	# if a substantial fraction of periods failed eigen, warn loudly so the
	# stability diagnostic itself doesn't silently underreport (NA-poisoning
	# would otherwise make `n_expl` drop to zero)
	if (n_failed > 0.2 * Tt) {
		cli::cli_warn(c(
			"Spectral-radius computation failed for {n_failed} of {Tt} period{?s}.",
			"i" = "Likely divergent operator entries in some posterior-mean slices; the explosive-operator diagnostic is unreliable on this fit.",
			"i" = "Inspect {.code coef(fit, \"A\")} for non-finite entries and consider rescaling the data."
		))
	}
	# two complementary triggers: many marginally-explosive periods, or any
	# clearly-explosive period. The zero-mean prior pulls operators toward 0,
	# so a single period with gain > 1.5 is already a strong signal even when
	# most periods are contractive.
	n_expl <- sum(gains > 1.1, na.rm = TRUE)
	max_gain <- max(gains, na.rm = TRUE)
	if (n_expl > threshold * Tt || (is.finite(max_gain) && max_gain > 1.5)) {
		cli::cli_warn(c(
			"Posterior-mean lag operator is locally explosive in {n_expl} of {Tt} period{?s} (max gain = {sprintf('%.2f', max_gain)}).",
			"i" = "Impulse responses and forecasts may diverge.",
			"i" = "Tighten the contractivity prior with a smaller {.arg shrink_rho} (it defaults to 0.9), or pool transitions with {.code model = \"piecewise\"}.",
			"i" = "Use {.fun dbn_operator} for the full per-period stability profile."
		))
	}
	invisible(gains)
}
####

####
#' Sparse-Aware Matrix Multiply
#'
#' @description Multiplies matrices that may be sparse
#' @param A First matrix
#' @param B Second matrix
#' @return Product A %*% B
#' @keywords internal
sparse_mult <- function(A, B) {
	is_sparse_A <- inherits(A, "Matrix") || inherits(A, "sparseMatrix")
	is_sparse_B <- inherits(B, "Matrix") || inherits(B, "sparseMatrix")

	if (is_sparse_A || is_sparse_B) {
		requireNamespace("Matrix", quietly = TRUE)
		return(A %*% B)
	} else {
		return(A %*% B)
	}
}
####

####
#' Sparse-Aware Quadratic Form
#'
#' @description Computes t(A) %*% B %*% A efficiently for diagonal B
#' @param A Matrix
#' @param B Matrix (often diagonal/sparse)
#' @return Quadratic form result
#' @keywords internal
sparse_quad <- function(A, B) {
	if (inherits(B, "ddiMatrix") || inherits(B, "diagonalMatrix") ||
		(is.matrix(B) && all(B[row(B) != col(B)] == 0))) {
		diag_B <- if (inherits(B, "ddiMatrix") || inherits(B, "diagonalMatrix")) {
			Matrix::diag(B)
		} else {
			diag(B)
		}
		return(t(A * sqrt(diag_B)) %*% (A * sqrt(diag_B)))
	} else {
		return(sparse_mult(t(A), sparse_mult(B, A)))
	}
}
####

####
#' Sparse Diagonal Matrix
#'
#' @description Creates diagonal matrix using sparse format for large dimensions
#' @param x Diagonal values or dimension
#' @param sparse_threshold Size threshold for sparse format
#' @return Diagonal matrix
#' @keywords internal
sparse_diag <- function(x, sparse_threshold = 100) {
	n <- if (length(x) == 1) x else length(x)

	if (n >= sparse_threshold && requireNamespace("Matrix", quietly = TRUE)) {
		if (length(x) == 1) {
			return(Matrix::Diagonal(n))
		} else {
			return(Matrix::Diagonal(x = x))
		}
	} else {
		return(diag(x))
	}
}
####

####
#' Simulate Low-Rank Factors
#'
#' @description Generates orthogonal U and V factors for low-rank decomposition
#' @param m Dimension
#' @param r Rank
#' @param sparse Use sparse matrices if beneficial
#' @return List with U and V matrices
#' @keywords internal
simulate_lowrank_factors <- function(m, r, sparse = FALSE) {
	if (r > m) cli::cli_abort("Rank cannot exceed dimension.")

	U_raw <- matrix(rnorm(m * r), m, r)
	V_raw <- matrix(rnorm(m * r), m, r)
	U <- qr.Q(qr(U_raw))
	V <- qr.Q(qr(V_raw))

	if (sparse && m >= 100 && r < m / 2 && requireNamespace("Matrix", quietly = TRUE)) {
		U <- Matrix::Matrix(U, sparse = TRUE)
		V <- Matrix::Matrix(V, sparse = TRUE)
	}

	list(U = U, V = V)
}
####

####
#' Initialize HMM States
#'
#' @description Generates initial hidden state sequence for HMM model
#' @param Tt Number of time points
#' @param K Number of states
#' @param init_probs Initial state probabilities
#' @param trans_probs Transition probability matrix
#' @return Integer vector of state assignments
#' @keywords internal
initialize_hmm_states <- function(Tt, K, init_probs = NULL, trans_probs = NULL) {
	if (is.null(init_probs)) {
		init_probs <- rep(1 / K, K)
	}
	if (is.null(trans_probs)) {
		trans_probs <- matrix(0.1 / (K - 1), K, K)
		diag(trans_probs) <- 0.9
	}

	states <- integer(Tt)
	states[1] <- sample(K, 1, prob = init_probs)
	for (t in 2:Tt) {
		states[t] <- sample(K, 1, prob = trans_probs[states[t - 1], ])
	}
	states
}
####

####
#' Sparse Matrix Multiply Operator
#'
#' @description Infix operator that keeps sparse objects sparse
#' @param A First matrix
#' @param B Second matrix
#' @return Matrix product
#' @keywords internal
`%sp%` <- function(A, B) {
	if (inherits(A, "sparseMatrix") || inherits(B, "sparseMatrix")) {
		return(Matrix::crossprod(Matrix::t(A), B))
	}
	A %*% B
}
####

####
#' Bilinear Multiplication Step
#'
#' @description Computes A %*% Theta %*% t(B)
#' @param A Left matrix
#' @param Theta_prev Middle matrix
#' @param B Right matrix
#' @return Bilinear product
#' @keywords internal
bilinear_step <- function(A, Theta_prev, B) {
	A %sp% Theta_prev %sp% Matrix::t(B)
}
####

####
#' Remove Diagonal from Matrix
#'
#' @description Sets diagonal elements to zero, handling sparse matrices
#' @param mat Sparse or regular matrix
#' @return Matrix with zero diagonal
#' @keywords internal
drop_diag <- function(mat) {
	if (inherits(mat, "sparseMatrix")) {
		diag(mat) <- 0
		return(Matrix::drop0(mat))
	} else {
		diag(mat) <- 0
		return(mat)
	}
}
####

####
#' Regularize Matrix for Cholesky
#'
#' @description Adds small ridge to diagonal for numerical stability
#' @param mat Matrix to regularize
#' @param eps Regularization strength
#' @return Regularized matrix
#' @keywords internal
regularize_matrix <- function(mat, eps = NULL) {
	if (is.null(eps)) {
		eps <- 5e-8 * mean(diag(mat), na.rm = TRUE)
	}
	mat + eps * diag(nrow(mat))
}
####

####
# attach actor / time dimnames from fit$Y onto operator arrays and Theta.
# called once at the end of dbn() so downstream accessors return labeled
# arrays. silently no-ops if fit$Y has no dimnames.
.attach_actor_names <- function(fit) {
	Y <- fit$Y
	if (is.null(Y) || is.null(dimnames(Y))) return(fit)
	dn <- dimnames(Y)
	row_nm <- dn[[1]]; col_nm <- dn[[2]]
	rel_nm <- if (length(dn) >= 3L) dn[[3]] else NULL
	time_nm <- if (length(dn) >= 4L) dn[[4]] else NULL

	# Theta: 5D [n_row, n_col, p, Tt, n_keep] (dynamic / piecewise / hmm) or 4D
	if (is.array(fit$Theta)) {
		d <- dim(fit$Theta)
		if (length(d) == 5L) {
			dimnames(fit$Theta) <- list(row_nm, col_nm, rel_nm, time_nm, NULL)
		} else if (length(d) == 4L) {
			dimnames(fit$Theta) <- list(row_nm, col_nm, rel_nm, time_nm)
		}
	}

	# A / B per-draw lists. for dynamic / piecewise the third dim is time
	# (length == Tt or Tt - 1 with anchor drop); for HMM it is regime
	# (length == R). use regime labels when the third dim isn't a time
	# match — otherwise we attach the wrong vector and break downstream
	# accessors
	attach_3d <- function(arr, rn, cn, tn_default) {
		if (!is.array(arr) || length(dim(arr)) != 3L) return(arr)
		third_n <- dim(arr)[3]
		tn <- if (!is.null(tn_default) && length(tn_default) == third_n)
			tn_default
		else if (!is.null(tn_default) && length(tn_default) - 1L == third_n)
			tn_default[-1L]    # piecewise / anchor-drop convention
		else
			paste0("regime_", seq_len(third_n))  # HMM regimes (or fallback)
		dimnames(arr) <- list(rn, cn, tn)
		arr
	}
	if (is.list(fit$A)) {
		fit$A <- lapply(fit$A, attach_3d, rn = row_nm, cn = row_nm,
			tn_default = time_nm)
	}
	if (is.list(fit$B)) {
		fit$B <- lapply(fit$B, attach_3d, rn = col_nm, cn = col_nm,
			tn_default = time_nm)
	}

	# M baseline mean (3D or 4D)
	if (is.array(fit$M)) {
		d <- dim(fit$M)
		if (length(d) == 3L) {
			dimnames(fit$M) <- list(row_nm, col_nm, rel_nm)
		} else if (length(d) == 4L) {
			dimnames(fit$M) <- list(row_nm, col_nm, rel_nm, NULL)
		}
	}
	fit
}
####
