####
# piecewise ALS: fit a separate static ALS on each segment defined by
# user-supplied breakpoints, then assemble into a time-varying cube
# where A_t is constant within each block. fast operator-per-block
# estimator for users with known regime breaks (e.g. before/after a
# crisis), without TV-ALS smoothing.
####

#' Piecewise-constant ALS by user-supplied breakpoints
#'
#' @param Y The 4D data array.
#' @param family Family string ("gaussian", "binary", "ordinal").
#' @param symmetric Logical.
#' @param breaks Numeric vector of breakpoint time indices in (1, Tt-1).
#'   Blocks are [1, breaks[1]], (breaks[1], breaks[2]], ..., (breaks[k], Tt].
#' @param ridge Ridge regularizer for the inner static fit.
#' @param max_iter Inner static max_iter.
#' @param verbose Logical.
#' @return A `dbn` fit with `meta$sampler_used = "als_piecewise"`, A and B
#'   as length-1 lists of `[n, n, T]` cubes (constant within blocks).
#' @keywords internal
#' @noRd
.dbn_als_piecewise <- function(Y, family, symmetric, breaks,
                                 ridge, max_iter, verbose, seed,
                                 bootstrap = 0L,
                                 bootstrap_type = "block",
                                 parallel = TRUE) {
	dims <- dim(Y); n_row <- dims[1]; n_col <- dims[2]; p <- dims[3]; Tt <- dims[4]

	# validate breaks
	if (!is.numeric(breaks) || any(!is.finite(breaks))) {
		cli::cli_abort("{.arg breaks} must be a numeric vector of breakpoint time indices.")
	}
	breaks <- sort(unique(as.integer(round(breaks))))
	if (any(breaks < 2L) || any(breaks >= Tt)) {
		cli::cli_abort(c(
			"All {.arg breaks} must lie in {.code [2, Tt - 1]} (got {.val {breaks}}, Tt = {Tt}).",
			"i" = "A breakpoint of {.val k} separates block [1..k] from block [k+1..Tt]."
		))
	}

	# construct block ranges: list of integer vectors (each is a contiguous
	# range of time indices).
	block_starts <- c(1L, breaks + 1L)
	block_ends   <- c(breaks, Tt)
	n_blocks <- length(block_starts)

	# each block needs at least 2 time points to fit static ALS
	short <- which((block_ends - block_starts + 1L) < 2L)
	if (length(short) > 0L) {
		cli::cli_abort(c(
			"Block(s) {.val {short}} have fewer than 2 time points -- ALS needs at least one transition.",
			"i" = "Drop the offending breakpoint or widen the block."
		))
	}

	if (verbose) {
		cli::cli_h3("Piecewise ALS")
		cli::cli_text("Network: {n_row} x {n_col}, {p} relation{?s}, {Tt} time point{?s}")
		cli::cli_text("Blocks: {n_blocks} (boundaries at {.val {breaks}})")
	}

	# fit a static ALS on each block
	block_fits <- vector("list", n_blocks)
	for (k in seq_len(n_blocks)) {
		t_idx <- block_starts[k]:block_ends[k]
		Y_blk <- Y[, , , t_idx, drop = FALSE]
		if (verbose) cli::cli_inform("  Block {k}/{n_blocks}: t = {block_starts[k]}-{block_ends[k]}")
		# per-block seed: if outer `seed` is NULL (default), pass NULL so each
		# block inherits global RNG state. If outer is an integer, offset by `k`
		# so blocks don't share the same seed.
		seed_k <- if (is.null(seed)) NULL else seed + k
		block_fits[[k]] <- dbn_als(Y_blk, family = family, symmetric = symmetric,
			lambda = "static", ridge = ridge, max_iter = max_iter,
			bootstrap = 0L, verbose = FALSE, seed = seed_k)
	}

	# assemble TV cube: A_t = A_block(t) for each t.
	A_cube <- array(0, dim = c(n_row, n_row, Tt))
	B_cube <- array(0, dim = c(n_col, n_col, Tt))
	M_arr  <- array(0, dim = c(n_row, n_col, p))
	# pool M across blocks (weighted by block length)
	M_total_weight <- 0
	for (k in seq_len(n_blocks)) {
		t_idx <- block_starts[k]:block_ends[k]
		A_k <- block_fits[[k]]$A[[1]][, , 1]
		B_k <- block_fits[[k]]$B[[1]][, , 1]
		for (t in t_idx) {
			A_cube[, , t] <- A_k
			B_cube[, , t] <- B_k
		}
		w_k <- length(t_idx)
		M_total_weight <- M_total_weight + w_k
		# block_fits[[k]]$M is [n_row, n_col, p, n_keep=1]; extract the
		# n_keep-1 slice as a [n_row, n_col, p] array before accumulating.
		M_k <- block_fits[[k]]$M[, , , 1, drop = FALSE]
		dim(M_k) <- c(n_row, n_col, p)
		M_arr <- M_arr + w_k * M_k
	}
	M_arr <- M_arr / M_total_weight

	# compute residual variance (sigma2) from the assembled fit
	resid_sq <- 0; n_resid <- 0
	for (t in 2:Tt) {
		for (r in 1:p) {
			Phi_prev <- Y[, , r, t - 1L] - M_arr[, , r]
			Phi_prev[is.na(Phi_prev)] <- 0
			Phi_curr <- Y[, , r, t] - M_arr[, , r]
			obs <- !is.na(Phi_curr)
			if (n_row == n_col) diag(obs) <- FALSE
			Phi_curr[!obs] <- 0
			pred <- A_cube[, , t] %*% Phi_prev %*% t(B_cube[, , t])
			pred[!obs] <- 0
			resid_sq <- resid_sq + sum((Phi_curr - pred)^2)
			n_resid <- n_resid + sum(obs)
		}
	}
	sigma2_est <- if (n_resid > 0) resid_sq / n_resid else 1

	# build a dbn fit object matching the TV-ALS layout
	A_store <- list(A_cube)
	B_store <- list(B_cube)
	M_store <- array(NA, dim = c(n_row, n_col, p, 1L))
	M_store[, , , 1] <- M_arr
	# phi (state) for downstream tools
	Theta_store <- array(0, dim = c(n_row, n_col, p, Tt, 1L))
	for (r in 1:p) {
		Theta_store[, , r, 1, 1] <- M_arr[, , r] + (Y[, , r, 1] - M_arr[, , r])
		Theta_store[, , r, 1, 1][is.na(Theta_store[, , r, 1, 1])] <- 0
		for (t in 2:Tt) {
			Phi_prev <- Y[, , r, t - 1L] - M_arr[, , r]
			Phi_prev[is.na(Phi_prev)] <- 0
			Theta_store[, , r, t, 1] <- M_arr[, , r] +
				A_cube[, , t] %*% Phi_prev %*% t(B_cube[, , t])
		}
	}

	out <- list(
		model = "dynamic",
		symmetric = isTRUE(symmetric),
		family = family,
		Y = Y,
		dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
					is_bipartite = (n_row != n_col), is_symmetric = symmetric),
		settings = list(nscan = 0, burn = 0, odens = 1, time_thin = 1,
			draws = 1L, ridge = ridge, breaks = breaks),
		meta = list(
			family = family,
			dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
						is_bipartite = (n_row != n_col), is_symmetric = symmetric),
			draws = 1L,
			sampler_requested = "als",
			sampler_used = "als_piecewise",
			uncertainty_available = FALSE,
			approximation_note = sprintf("Piecewise-constant ALS with %d block%s (boundaries at %s); no credible intervals",
				n_blocks, if (n_blocks == 1L) "" else "s",
				paste(breaks, collapse = ", ")),
			breaks = breaks,
			block_starts = block_starts,
			block_ends = block_ends,
			n_blocks = n_blocks
		),
		params = data.frame(sigma2 = sigma2_est),
		draws = list(pars = data.frame(sigma2 = sigma2_est),
			misc = list(A = A_store, B = B_store, M = M_store)),
		time_kept = seq_len(Tt),
		A = A_store,
		B = B_store,
		M = M_store,
		Theta = Theta_store,
		sigma2 = sigma2_est,
		tau_A2 = NA_real_,
		tau_B2 = NA_real_,
		g2 = NA_real_,
		n_iter = max_iter,
		burn = 0, thin = 1, time_thin = 1
	)
	class(out) <- "dbn"

	# bootstrap support: route through the TV-ALS bootstrap machinery, which
	# already handles per-time A_t, B_t replicates and family-aware refits.
	# each bootstrap replicate refits the piecewise model on simulated data.
	if (bootstrap > 0L) {
		if (verbose) cli::cli_inform("Bootstrapping piecewise fit ({bootstrap} replicates) ...")
		boot <- .dbn_als_tv_bootstrap(out, R = bootstrap,
			type = switch(bootstrap_type,
				"block" = "residual_block",
				"parametric" = "fixed",
				bootstrap_type),
			verbose = isTRUE(verbose), parallel = parallel)
		# route through the shared attacher: downgrades to point-only when
		# every replicate failed instead of silently shipping NA CIs.
		out <- .dbn_attach_bootstrap(out, boot,
			label = sprintf("Piecewise-ALS point estimate plus %d-replicate %s bootstrap CIs.",
				boot$n_valid %||% 0L,
				switch(bootstrap_type, "block" = "residual_block",
					"parametric" = "fixed", bootstrap_type)),
			source = "bootstrap")
	}

	out
}
