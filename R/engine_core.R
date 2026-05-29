####
# preprocessing and model registration
####

#' Generic MCMC Engine for DBN Models
#'
#' @description Core Gibbs sampler that drives all DBN model variants
#' @name engine_core
#' @keywords internal
NULL

#' Warn on a non-NA actor diagonal
#'
#' @description Warns (without failing) when a square network's actor-by-actor
#'   diagonal contains non-`NA` values; self-ties are not modelled and the
#'   diagonal is expected to be `NA`. Shared by `dbn()` and `dbn_als()`.
#'   Suppressed inside bootstrap refits (where the top-level call already
#'   warned and the inner replicates would otherwise log it once per rep).
#' @param Y A 4D data array.
#' @return Invisibly `NULL`.
#' @keywords internal
#' @noRd
warn_filled_diagonal <- function(Y) {
	# skip when called from inside a bootstrap loop. The flag is set by the
	# bootstrap helper before per-rep refits and cleared on exit.
	if (isTRUE(getOption("dbn.in_bootstrap", FALSE))) return(invisible(NULL))
	d <- dim(Y)
	if (length(d) != 4L || d[1] != d[2]) return(invisible(NULL))
	diag_idx <- which(diag(TRUE, d[1]))
	for (rr in seq_len(d[3])) {
		for (tt in seq_len(d[4])) {
			if (any(!is.na(Y[, , rr, tt][diag_idx]))) {
				cli::cli_warn(c(
					"The actor-by-actor diagonal contains non-{.val NA} values.",
					"i" = "Self-ties are not modelled; the diagonal is expected to be {.val NA}.",
					"i" = "Diagonal entries will be treated as observed data unless set to {.val NA}."
				))
				return(invisible(NULL))
			}
		}
	}
	invisible(NULL)
}

#' Shared Preprocessing
#'
#' @description Common preprocessing for all DBN models
#' @param Y Data array
#' @param family Distribution family
#' @return List with preprocessed components
#' @keywords internal
shared_preprocess <- function(Y, family = "ordinal") {
	# handle sparse input
	if (inherits(Y, "dgCMatrix")) {
		Y <- as(Y, "TsparseMatrix")
	}

	n_row <- dim(Y)[1]
	n_col <- dim(Y)[2]
	is_bipartite <- (n_row != n_col)

	dims <- list(
		m = n_row,
		n_row = n_row,
		n_col = n_col,
		p = dim(Y)[3],
		Tt = dim(Y)[4],
		is_bipartite = is_bipartite
	)

	# validate data content. NA is supported (missing data); but the data
	# must contain at least one usable observation, and NaN/Inf are anomalies
	# (NA is the proper missing marker) that would otherwise be coerced
	# silently to zero below.
	if (!is.numeric(Y)) {
		cli::cli_abort(c(
			"Data array must be numeric.",
			"x" = "Got a {.cls {typeof(Y)}} array.",
			"i" = "Convert factors/characters to numeric codes before fitting."
		))
	}
	if (!any(is.finite(Y))) {
		cli::cli_abort(c(
			"Data array contains no usable (finite) observations.",
			"i" = "Every entry is {.val NA}, {.val NaN}, or infinite -- there is nothing to fit."
		))
	}
	n_nan <- sum(is.nan(Y))
	n_inf <- sum(is.infinite(Y))
	if (n_nan > 0 || n_inf > 0) {
		cli::cli_warn(c(
			"Data array contains {n_nan} {.val NaN} and {n_inf} infinite value{?s}.",
			"i" = "These will be treated as missing. Use {.val NA} to mark missing data explicitly."
		))
	}

	# validate data for the chosen family
	if (family == "binary") {
		Y_vals <- Y[!is.na(Y)]
		if (length(Y_vals) > 0 && !all(Y_vals %in% c(0, 1))) {
			cli::cli_abort(c(
				"Binary family requires 0/1 data.",
				"x" = "Found values outside {{0, 1}}.",
				"i" = "Use {.code family = \"gaussian\"} for continuous data or {.code family = \"ordinal\"} for integer data."
			))
		}
	} else if (family == "ordinal") {
		Y_vals <- Y[!is.na(Y)]
		if (length(Y_vals) > 0 && !all(Y_vals == floor(Y_vals))) {
			cli::cli_warn(c(
				"Ordinal family expects integer-valued data.",
				"i" = "Non-integer values detected; data will be treated as continuous ranks."
			))
		}
	} else if (family == "gaussian") {
		# warn when the data look discrete -- a common novice mistake is to
		# pass ordinal/binary-coded data to the gaussian family
		Y_vals <- Y[is.finite(Y)]
		if (length(Y_vals) > 0 && all(Y_vals == floor(Y_vals))) {
			uv <- unique(as.vector(Y_vals))
			if (all(uv %in% c(0, 1))) {
				cli::cli_warn(c(
					"Data are all 0/1 but {.code family = \"gaussian\"} was selected.",
					"i" = "If these are presence/absence ties, use {.code family = \"binary\"}."
				))
			} else if (length(uv) <= 20L) {
				cli::cli_warn(c(
					"Data are integer-valued with only {length(uv)} distinct value{?s}, but {.code family = \"gaussian\"} was selected.",
					"i" = "If these are ordered categories, use {.code family = \"ordinal\"}."
				))
			}
		}
	}

	# rank structure for ordinal sampling
	IR <- precompute_ranks(Y)

	# build the observation mask Omega (1 = observed, 0 = missing or
	# structural diagonal-NA). off-diagonal NA cells represent actor
	# entry/exit, dropouts, structurally absent dyads, etc., and must
	# not be silently coerced to "observed = 0". the MCMC samplers in
	# models.R treat the diagonal as observed-zero (the zero-diagonal
	# convention); the TV-ALS path uses Omega directly to skip masked
	# cells. we surface Omega so downstream consumers (.bootstrap_expand,
	# posterior_predict_dbn, fitted.dbn) can also skip them.
	Omega <- array(1L, dim = dim(Y))
	off_diag_na <- 0L
	for (rr in seq_len(dims$p)) {
		for (tt in seq_len(dims$Tt)) {
			Y_slice <- Y[, , rr, tt, drop = FALSE]
			dim(Y_slice) <- c(n_row, n_col)
			na_mask <- is.na(Y_slice)
			Omega[, , rr, tt] <- as.integer(!na_mask)
			# count off-diagonal NA cells (unipartite case)
			if (n_row == n_col && n_row >= 2L) {
				diag(na_mask) <- FALSE
			}
			off_diag_na <- off_diag_na + sum(na_mask)
		}
	}
	if (off_diag_na > 0L) {
		cli::cli_warn(c(
			"{off_diag_na} off-diagonal cell{?s} of {.arg Y} are NA (actor entry/exit, dropouts, or structurally absent dyads).",
			"i" = "The MCMC sampler treats these as observed zero (binary) or z=0 (ordinal/gaussian) at initialization and during sampling. This is not the right behaviour for cohort or panel data with actor entry/exit.",
			"i" = "For honest missingness, use {.fun dbn_als} (the TV-ALS path applies an observation mask correctly), or impute NAs explicitly before calling {.fun dbn}.",
			"i" = "The observation mask is stored as {.code fit$meta$Omega} for downstream use."
		))
	}

	# initial latent Z based on family
	Z <- Y
	if (family == "ordinal") {
		means <- numeric(dims$p)
		sds <- numeric(dims$p)
		for (j in 1:dims$p) {
			Y_j <- Y[, , j, ]
			means[j] <- mean(Y_j, na.rm = TRUE)
			sds[j] <- sd(Y_j, na.rm = TRUE)
		}
		Y_flat <- array(Y, dim = c(n_row, n_col, dims$p * dims$Tt))
		Z_flat <- compute_zscores_batch(Y_flat, means, sds, n_row, n_col, dims$p, dims$Tt)
		Z <- array(Z_flat, dim = c(n_row, n_col, dims$p, dims$Tt))
		Z[is.na(Z)] <- 0
	} else {
		Z <- Y
		if (family == "gaussian") {
			Z[!is.finite(Z)] <- 0
		} else if (family == "binary") {
			Z[is.na(Z)] <- 0
		}
	}

	# initial mean M (NaN arises from all-NA diagonals in unipartite networks)
	M <- array(apply(Z, c(1, 2, 3), mean, na.rm = TRUE), dim = c(n_row, n_col, dims$p))
	M[is.nan(M)] <- 0

	# initial Theta (centered residual + noise)
	Theta <- sweep(Z, c(1, 2, 3), M, "-") + rsan(dim(Z))

	list(
		Y = Y,
		Z = Z,
		M = M,
		Theta = Theta,
		R = Y,
		IR = IR,
		Omega = Omega,
		dims = dims
	)
}

#' Get package environment
#'
#' `.pkg_env$models` is allocated for diagnostic uses (and so `.onLoad`
#' has a binding to assign into); `dbn()` validates `model` via
#' `match.arg` and does not consult the registry.
#'
#' @keywords internal
#' @noRd
get_pkg_env <- function() {
	ns <- getNamespace("dbn")
	if (!exists(".pkg_env", envir = ns)) {
		assign(".pkg_env", new.env(parent = emptyenv()), envir = ns)
		assign("models", list(), envir = get(".pkg_env", envir = ns))
	}
	get(".pkg_env", envir = ns)
}
