####
# multi-chain wrapper around dbn(): K independent single-chain fits
# bundled into a dbn_multichain that downstream posterior::, bayesplot::,
# and loo:: tools consume.
####

#' Run multiple independent dbn chains
#'
#' Wrapper that calls [dbn()] K times with K different seeds and bundles
#' the results into a `dbn_multichain` object. Each fit retains its full
#' draws; convergence summaries (Rhat, ESS) can then be computed across
#' chains via [check_convergence()] or via the standard `posterior`
#' package after [as_draws.dbn_multichain()].
#'
#' Parallelism follows the `future` package's plan if installed and a
#' multisession plan is active; otherwise chains run sequentially. The
#' default `seeds = NULL` chooses `K` values deterministically from
#' `1:.Machine$integer.max` based on the system seed at call time.
#'
#' @param Y A 4D array `[n_row, n_col, p, Tt]` (passed to [dbn()]).
#' @param chains Integer >= 1. Number of independent samplers to run.
#' @param seeds Integer vector of length `chains`, or `NULL` (default)
#'   to choose deterministic seeds derived from `set.seed(NULL)`.
#' @param ... All other arguments passed through to [dbn()] (e.g.
#'   `model`, `family`, `nscan`, `burn`, `verbose`).
#' @param parallel Logical. If `TRUE` (default) and `future.apply` is
#'   installed and a multisession plan is active, chains run in
#'   parallel; otherwise sequentially.
#' @return A list of class `dbn_multichain` with elements:
#'   - `chains`: list of K `dbn` fit objects.
#'   - `seeds`: integer vector of seeds actually used.
#'   - `n_chains`: K.
#'   - `args`: the captured `...` for reproducibility.
#' @seealso [dbn()], [check_convergence()], [as_draws.dbn_multichain()]
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
#' mc <- dbn_multichain(sim$Z, chains = 2, seeds = c(1, 2),
#'                      model = "dynamic", family = "gaussian",
#'                      nscan = 150, burn = 50, verbose = FALSE)
#' length(mc$chains)
#' }
dbn_multichain <- function(Y, chains = 2L, seeds = NULL, ...,
                            parallel = TRUE) {
	# catch the `seed = ...` (singular) typo before R partial-matches it
	# onto `seeds` and length validation gives a confusing error.
	raw_call <- sys.call()
	raw_names <- names(raw_call) %||% character(0L)
	called_seed <- any(raw_names == "seed")
	if (called_seed) {
		cli::cli_warn(c(
			"{.fun dbn_multichain} takes {.arg seeds} (plural), not {.arg seed}.",
			"i" = "R partial-matched {.code seed = ...} to {.arg seeds}; auto-generating {.arg seeds} of length {chains} and ignoring the single value supplied.",
			"i" = "Pass {.code seeds = c(s1, s2, ...)} with one value per chain to silence."
		))
		seeds <- NULL
	}
	if (length(chains) != 1L || !is.numeric(chains) || !is.finite(chains) ||
	    chains < 1 || chains != round(chains)) {
		cli::cli_abort("{.arg chains} must be a single positive integer.")
	}
	chains <- as.integer(chains)
	if (is.null(seeds)) {
		# deterministic but distinct seeds
		seeds <- as.integer(sample.int(.Machine$integer.max, chains))
	}
	args <- list(...)
	# strip a seed= that survived because seeds= was also passed
	if ("seed" %in% names(args)) args$seed <- NULL
	if (length(seeds) != chains) {
		cli::cli_abort(c(
			"{.arg seeds} must have length {.arg chains} = {chains}.",
			"x" = "Got length {length(seeds)}.",
			"i" = "Pass {.code seeds = c(s1, s2, ...)} with one value per chain, or leave NULL to auto-generate."
		))
	}
	if (!is.numeric(seeds) || any(!is.finite(seeds))) {
		cli::cli_abort("{.arg seeds} must be a finite numeric vector.")
	}
	seeds <- as.integer(seeds)
	# strip chains= if the user accidentally passed it
	args$chains <- NULL
	# multichain is most often run to diagnose tau_A^2 / tau_B^2, which
	# only exist in the dynamic model; default to dynamic here.
	if (is.null(args$model)) args$model <- "dynamic"

	fit_one <- function(s) {
		do.call(dbn, c(list(data = Y, seed = s), args))
	}
	# future.seed = NULL lets each fit_one's inner set.seed(s) control
	# reproducibility; the L'Ecuyer streams from future.seed=TRUE would
	# override it and de-sync parallel from sequential runs.
	use_parallel <- isTRUE(parallel) &&
	                requireNamespace("future.apply", quietly = TRUE)
	# parallel = TRUE without an active multisession plan silently falls
	# back to sequential; tell the user so they can opt in to real parallelism.
	# future::plan() returns an object whose class vector begins with
	# "FutureStrategy" and includes "sequential" / "multisession" / etc., so
	# we check inclusion rather than the first element.
	if (use_parallel && requireNamespace("future", quietly = TRUE)) {
		plan_cls <- class(future::plan())
		if ("sequential" %in% plan_cls || "uniprocess" %in% plan_cls) {
			cli::cli_inform(c(
				"i" = "{.code parallel = TRUE} but no multisession plan is active; chains will run sequentially.",
				"i" = "Set up parallelism with {.code future::plan(future::multisession, workers = {chains})} before calling."
			))
		}
	}
	chain_fits <- if (use_parallel) {
		future.apply::future_lapply(seeds, fit_one, future.seed = NULL)
	} else {
		lapply(seeds, fit_one)
	}
	out <- list(
		chains = chain_fits,
		seeds = seeds,
		n_chains = chains,
		args = args
	)
	class(out) <- c("dbn_multichain", "list")
	out
}

#' Print method for dbn_multichain
#' @param x A `dbn_multichain` object.
#' @param ... Unused.
#' @export
#' @keywords internal
print.dbn_multichain <- function(x, ...) {
	cli::cli_h2("dbn_multichain")
	cli::cli_text("{x$n_chains} chains")
	cli::cli_text("seeds: {paste(x$seeds, collapse = ', ')}")
	if (length(x$chains) >= 1L) {
		f1 <- x$chains[[1]]
		# compute the values before cli interpolation: cli >= 3.4.0
		# parses any {.name(...)} token as a style class.
		fam_str <- .fam_name(f1$family)
		mdl_str <- f1$model %||% NA
		cli::cli_text("model = {.val {mdl_str}}, family = {.val {fam_str}}")
		dims <- f1$dims
		if (!is.null(dims)) {
			cli::cli_text("dims: n_row = {dims$n_row}, n_col = {dims$n_col}, p = {dims$p}, Tt = {dims$Tt}")
		}
	}
	invisible(x)
}

#' @keywords internal
.fam_name <- function(f) {
	if (is.character(f)) f else if (is.list(f) && !is.null(f$name)) f$name else NA_character_
}

#' Convert dbn_multichain to a posterior::draws_array
#'
#' Builds an `[iterations, chains, variables]` array of scalar
#' parameters (sigma^2, sigma^2_obs, tau_A^2, tau_B^2) for downstream
#' use with the `posterior` package: rank-normalised split-Rhat, ESS,
#' summary tables, and bayesplot.
#'
#' Operator parameters (`A_t`, `B_t`, `M`) are not included by default
#' because the `posterior` array layout flattens variables into a
#' character-named axis; including thousands of operator entries
#' explodes the size. Use `include_operators = TRUE` to opt in.
#'
#' Requires the `posterior` package.
#'
#' @param x A `dbn_multichain` object.
#' @param include_operators Logical (default `FALSE`). If `TRUE`,
#'   include per-cell `A_t`, `B_t`, and `M` entries (variable names of
#'   the form `A[i,j,t]`, `B[i,j,t]`, `M[i,j]`).
#' @param ... Unused.
#' @return A `posterior::draws_array` object.
#' @author Tosin Salau and Shahryar Minhas
#' @export
as_draws.dbn_multichain <- function(x, include_operators = FALSE, ...) {
	if (!requireNamespace("posterior", quietly = TRUE)) {
		cli::cli_abort(c(
			"{.fun as_draws.dbn_multichain} requires the {.pkg posterior} package.",
			"i" = "Install with {.code install.packages(\"posterior\")}."
		))
	}
	# nudge: dynamic / piecewise fits have time-varying operators A_t, B_t
	# that don't appear in the default scalar-only output. Users routinely
	# expect them and silently get only sigma2/tau_*2/g2; surface the option.
	# fires once per session via getOption() to avoid log spam.
	mdl1 <- if (length(x$chains) >= 1L) x$chains[[1]]$model %||% NA_character_ else NA_character_
	if (!isTRUE(include_operators) && !is.na(mdl1) &&
	    mdl1 %in% c("dynamic", "piecewise", "hmm", "lowrank") &&
	    !isTRUE(getOption("dbn.silence_as_draws_nudge", FALSE)) &&
	    isTRUE(interactive())) {
		cli::cli_inform(c(
			"i" = "Returning scalar variances only. Pass {.code include_operators = TRUE} to expand {.var A_t}/{.var B_t}/{.var M} per cell.",
			"i" = "Silence with {.code options(dbn.silence_as_draws_nudge = TRUE)}."
		))
	}
	K <- x$n_chains
	# pull scalar traces from each chain and rename to canonical sigma2 / tau_*2
	extract_scalars <- function(f) {
		mdl <- f$model %||% ""
		if (identical(mdl, "piecewise") || identical(mdl, "static")) {
			pars <- f$draws$pars
			if (is.null(pars) || !is.data.frame(pars)) pars <- f$params
			if (is.null(pars) || !is.data.frame(pars)) return(list())
			out <- lapply(pars, as.numeric)
			ren <- c(s2 = "sigma2", t2 = "tau_A2", g2 = "g2")
			nms <- names(out)
			nms[nms %in% names(ren)] <- ren[nms[nms %in% names(ren)]]
			names(out) <- nms
			# the static model has a single variance; `sigma2_obs` is an alias
			# of `sigma2` (s2), not a separately-sampled observation variance.
			# only expose it when it actually differs, so as_draws() doesn't
			# show a duplicated parameter.
			if (is.numeric(f$sigma2_obs) && !"sigma2_obs" %in% names(out) &&
			    !(("sigma2" %in% names(out)) &&
			      isTRUE(all.equal(as.numeric(f$sigma2_obs), out[["sigma2"]],
			                       check.attributes = FALSE)))) {
				out$sigma2_obs <- as.numeric(f$sigma2_obs)
			}
			return(out)
		}
		if (identical(mdl, "hmm")) {
			# HMM uses sigma2_proc internally; expose as sigma2.
			# pull every numeric column from $draws$pars (renaming proc -> sigma2)
			# plus any of tau_A2/tau_B2/g2/sigma2_obs at top level so the user
			# gets the same scalar diagnostics as the dynamic path
			out <- list()
			pars <- f$draws$pars
			if (is.data.frame(pars)) {
				for (nm in names(pars)) {
					v <- pars[[nm]]
					if (!is.numeric(v)) next
					name <- if (nm == "sigma2_proc") "sigma2" else nm
					out[[name]] <- as.numeric(v)
				}
			}
			# augment with top-level vectors when present
			for (nm in c("tau_A2", "tau_B2", "g2", "sigma2_obs")) {
				v <- f[[nm]]
				if (is.numeric(v) && !nm %in% names(out)) out[[nm]] <- as.numeric(v)
			}
			return(out)
		}
		# dynamic / lowrank: scalars stored at top level
		scalar_names <- c(
			"sigma2", "sigma2_obs", "tau_A2", "tau_B2",
			"g2", "rhoA", "rhoB", "tau_alpha2"
		)
		out <- lapply(scalar_names, function(nm) {
			v <- f[[nm]]
			if (is.null(v) || !is.numeric(v)) return(NULL)
			as.numeric(v)
		})
		names(out) <- scalar_names
		out[!vapply(out, is.null, logical(1L))]
	}
	chain_scalars <- lapply(x$chains, extract_scalars)
	# align chains on the minimum scalar iteration count (per chain, then min across chains)
	scalar_n_iter <- min(vapply(chain_scalars, function(s) {
		if (length(s) == 0L) 0L else min(vapply(s, length, integer(1L)))
	}, integer(1L)))

	# when include_operators = TRUE, append per-cell A[i,j,t], B[i,j,t],
	# M[i,j] entries as separate posterior variables so tidybayes and
	# posterior::summarise_draws can address them.
	op_arrays <- list()
	if (isTRUE(include_operators)) {
		extract_op <- function(f, nm) {
			# piecewise: 4d block draws at f$draws$A_blocks / B_blocks
			# (n_draws x K x n_row x n_col); expose with A_block / B_block labels
			pw_name <- paste0(nm, "_blocks")
			pw <- f$draws[[pw_name]]
			if (is.array(pw) && length(dim(pw)) == 4L) return(pw)
			# dynamic / lowrank / hmm: list of per-draw cubes
			a <- f[[nm]]
			if (is.list(a) && length(a) > 0L) {
				d1 <- dim(a[[1]])
				if (is.null(d1) || length(d1) == 0L) return(NULL)
				flat <- do.call(rbind, lapply(a, as.numeric))
				out <- array(flat, dim = c(length(a), d1))
				return(out)
			}
			if (is.array(a)) return(NULL)
			NULL
		}
		op_per_chain <- vector("list", K)
		for (k in seq_len(K)) {
			op_per_chain[[k]] <- list()
			for (nm in c("A", "B")) {
				oa <- extract_op(x$chains[[k]], nm)
				if (!is.null(oa)) {
					# piecewise stores at $draws$A_blocks; label vars A_block / B_block
					label <- if (is.array(x$chains[[k]]$draws[[paste0(nm, "_blocks")]])) {
						paste0(nm, "_block")
					} else {
						nm
					}
					op_per_chain[[k]][[label]] <- oa
				}
			}
			Mk <- x$chains[[k]]$M
			if (is.array(Mk) && length(dim(Mk)) == 4L) {
				op_per_chain[[k]][["M"]] <- aperm(Mk, c(4, 1, 2, 3))
			}
		}
		# operator draws may be longer than scalar chains (e.g. ALS+bootstrap
		# fits have length-1 scalars and R-rep operator draws). use the min
		# operator chain length across chains as the operator n_iter.
		op_n_iter_per_chain <- vapply(op_per_chain, function(opk) {
			if (length(opk) == 0L) 0L else min(vapply(opk, function(a) as.integer(dim(a)[1]), integer(1L)))
		}, integer(1L))
		op_n_iter <- if (length(op_n_iter_per_chain) > 0L) min(op_n_iter_per_chain) else 0L
		ch1 <- op_per_chain[[1]]
		# resolve n_iter: prefer scalar_n_iter when chains agree, but allow
		# op_n_iter to win when scalars are length-1 (bootstrap-expand path).
		n_iter <- if (scalar_n_iter >= 2L) scalar_n_iter else op_n_iter
		# drop t=1 only for dynamic random-walk fits (the t=1 identity anchor).
		# piecewise block 1 is estimated freely and must NOT be dropped.
		ch1_model <- x$chains[[1]]$model %||% NA_character_
		drop_anchor <- identical(ch1_model, "dynamic")
		for (nm in names(ch1)) {
			ar1 <- ch1[[nm]]
			d_op <- dim(ar1)
			cell_dims <- d_op[-1]
			n_cells <- prod(cell_dims)
			per_chain_mats <- vector("list", K)
			for (k in seq_len(K)) {
				ar_k <- op_per_chain[[k]][[nm]]
				if (is.null(ar_k)) {
					per_chain_mats[[k]] <- matrix(NA_real_, nrow = n_iter, ncol = n_cells)
					next
				}
				iters_k <- min(n_iter, dim(ar_k)[1])
				mat_k <- matrix(as.numeric(ar_k), nrow = dim(ar_k)[1])
				mat_k <- mat_k[seq_len(iters_k), , drop = FALSE]
				pad <- matrix(NA_real_, nrow = n_iter - iters_k, ncol = n_cells)
				per_chain_mats[[k]] <- rbind(mat_k, pad)
			}
			# label each cell by its multi-index; skip t=1 only on dynamic A / B
			idx_mat <- arrayInd(seq_len(n_cells), cell_dims)
			has_time_axis <- nm %in% c("A", "B") && length(cell_dims) == 3L
			for (cell_i in seq_len(n_cells)) {
				idx <- as.integer(idx_mat[cell_i, ])
				if (isTRUE(drop_anchor) && has_time_axis && idx[3] == 1L) next
				var_label <- sprintf("%s[%s]", nm, paste(idx, collapse = ","))
				mat <- vapply(per_chain_mats, function(m) m[, cell_i], numeric(n_iter))
				op_arrays[[var_label]] <- mat
			}
		}
	} else {
		# include_operators = FALSE: fall back to scalar-only chain length
		n_iter <- scalar_n_iter
	}
	if (is.null(n_iter) || n_iter < 2L) {
		cli::cli_abort(c(
			"Could not extract aligned chains of length >= 2.",
			"i" = "Scalar chains have length {scalar_n_iter}; pass {.code include_operators = TRUE} for ALS+bootstrap fits.",
			"i" = "Otherwise ensure each chain ran for the same {.arg nscan} / {.arg odens}."
		))
	}

	# only keep scalars that have at least n_iter values per chain. for
	# bootstrap-expanded fits, length-1 scalars are dropped (the operator
	# chain is the real source of variation; a length-1 sigma2 padded to
	# n_iter would be a fake chain).
	var_names <- Filter(function(nm) {
		all(vapply(chain_scalars, function(s) length(s[[nm]]) >= n_iter, logical(1L)))
	}, names(chain_scalars[[1]]))

	full_var_names <- c(var_names, names(op_arrays))
	arr <- array(NA_real_,
	             dim = c(n_iter, K, length(full_var_names)),
	             dimnames = list(iteration = NULL, chain = NULL, variable = full_var_names))
	for (k in seq_len(K)) {
		for (v_idx in seq_along(var_names)) {
			vec <- chain_scalars[[k]][[var_names[v_idx]]]
			arr[, k, v_idx] <- vec[seq_len(n_iter)]
		}
	}
	if (length(op_arrays) > 0L) {
		offset <- length(var_names)
		for (i in seq_along(op_arrays)) {
			arr[, , offset + i] <- op_arrays[[i]]
		}
	}
	posterior::as_draws_array(arr)
}

#' Convert dbn fit to a posterior::draws_array
#'
#' Single-chain analog of [as_draws.dbn_multichain()]. Useful for
#' calling `posterior::summarise_draws()` / `bayesplot` on a one-chain
#' fit. Multi-chain Rhat is meaningless on a single chain but
#' rank-normalised ESS still works.
#'
#' @param x A `dbn` fit.
#' @param include_operators Logical (default `FALSE`). If `TRUE`, expand
#'   per-cell entries (variable names of the form `A[i,j,t]`,
#'   `B[i,j,t]`, `M[i,j]`) as separate posterior variables.
#' @param ... Unused.
#' @return A `posterior::draws_array` (chains = 1).
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'            nscan = 200, burn = 100, verbose = FALSE)
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   da <- as_draws.dbn(fit)
#'   dim(da)
#' }
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
as_draws.dbn <- function(x, include_operators = FALSE, ...) {
	if (!requireNamespace("posterior", quietly = TRUE)) {
		cli::cli_abort(c(
			"{.fun as_draws.dbn} requires the {.pkg posterior} package.",
			"i" = "Install with {.code install.packages(\"posterior\")}."
		))
	}
	# ALS+bootstrap fits have per-rep (A, B, M) but no chain over the
	# scalar variances; expand the bootstrap into a per-rep operator
	# trace and treat the reps as a single chain. scalar variance Rhat
	# is meaningless in this layout so we force operator-only output.
	su <- x$meta$sampler_used %||% ""
	is_als <- su %in% c("als", "als_tv")
	if (is_als) {
		if (is.null(x$bootstrap)) {
			cli::cli_abort(c(
				"{.fun as_draws.dbn} on an ALS fit requires bootstrap draws.",
				"x" = "{.code fit$bootstrap} is NULL.",
				"i" = "Refit with {.code dbn_als(..., bootstrap = N)} for N >= 30."
			))
		}
		if (!isTRUE(x$meta$bootstrap_expanded)) {
			x <- .bootstrap_expand(x)
		}
		if (!isTRUE(include_operators)) {
			cli::cli_inform(c(
				"!" = "ALS+bootstrap has no chain over scalar variances; {.fun as_draws.dbn} returns operator entries only.",
				"i" = "Forcing {.code include_operators = TRUE}."
			))
			include_operators <- TRUE
		}
	}
	mc <- list(chains = list(x), seeds = NA_integer_, n_chains = 1L, args = list())
	class(mc) <- c("dbn_multichain", "list")
	as_draws.dbn_multichain(mc, include_operators = include_operators)
}

#' Cross-chain convergence diagnostics for a dbn_multichain
#'
#' Computes rank-normalised split-Rhat and ESS-bulk on the scalar
#' variance parameters across chains, plus a one-line headline of the
#' worst-mixing parameter. Used by [check_convergence()] when called on
#' a `dbn_multichain` object.
#'
#' @keywords internal
#' @noRd
.check_convergence_multichain <- function(mc) {
	# cross-chain Rhat / ESS using posterior if installed, else a manual
	# between/within estimator
	K <- mc$n_chains
	if (K < 2L) {
		cli::cli_inform("Single-chain {.cls dbn_multichain}; falling back to per-chain {.fn check_convergence}.")
		return(check_convergence(mc$chains[[1]]))
	}
	scalar_names <- c("sigma2", "sigma2_obs", "tau_A2", "tau_B2")
	rows <- list()
	for (nm in scalar_names) {
		mats <- lapply(mc$chains, function(f) {
			v <- f[[nm]]
			if (is.null(v) || !is.numeric(v)) return(NULL)
			as.numeric(v)
		})
		mats <- mats[!vapply(mats, is.null, logical(1L))]
		if (length(mats) < 2L) next
		n_min <- min(vapply(mats, length, integer(1L)))
		if (n_min < 10L) next
		mat <- do.call(cbind, lapply(mats, function(v) v[seq_len(n_min)]))
		# rhat: posterior::rhat if available, else simple between/within.
		if (requireNamespace("posterior", quietly = TRUE)) {
			rh <- tryCatch(as.numeric(posterior::rhat(mat)), error = function(e) NA_real_)
			ess_b <- tryCatch(as.numeric(posterior::ess_bulk(mat)), error = function(e) NA_real_)
		} else {
			chain_means <- colMeans(mat)
			chain_vars <- apply(mat, 2, var)
			B <- n_min * var(chain_means)
			W <- mean(chain_vars)
			var_hat <- ((n_min - 1) / n_min) * W + B / n_min
			rh <- if (W > 0) sqrt(var_hat / W) else NA_real_
			# crude ESS: n_min * K / (1 + 2*ac1)
			ac1 <- mean(apply(mat, 2, function(v) cor(v[-1], v[-length(v)])), na.rm = TRUE)
			ess_b <- if (is.finite(ac1) && ac1 < 1) n_min * K * (1 - ac1) / (1 + ac1) else NA_real_
		}
		rows[[length(rows) + 1L]] <- data.frame(
			parameter = nm, rhat = rh, ess_bulk = ess_b,
			stringsAsFactors = FALSE
		)
	}
	out <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
	rownames(out) <- NULL
	cli::cli_h3("Multi-chain convergence ({K} chains)")
	if (nrow(out) == 0L) {
		cli::cli_alert_warning("No scalar parameters available across all chains.")
		return(invisible(out))
	}
	worst_rhat <- max(out$rhat, na.rm = TRUE)
	min_ess <- min(out$ess_bulk, na.rm = TRUE)
	headline <- if (worst_rhat > 1.05 || min_ess < 100) "x" else if (worst_rhat > 1.01 || min_ess < 400) "!" else "v"
	if (headline == "v") {
		cli::cli_alert_success("All scalars: Rhat <= 1.01 and ESS_bulk >= 400")
	} else if (headline == "!") {
		cli::cli_alert_warning("Worst Rhat = {sprintf('%.3f', worst_rhat)}, min ESS_bulk = {round(min_ess)} -- borderline (target Rhat < 1.01, ESS >= 400).")
	} else {
		cli::cli_alert_danger("Worst Rhat = {sprintf('%.3f', worst_rhat)}, min ESS_bulk = {round(min_ess)} -- NOT converged. Increase {.arg nscan} or check the model.")
	}
	print(out, row.names = FALSE)
	invisible(out)
}

# dbn_multichain accessors ----------------------------------------------
# every standard dbn accessor (tidy, coef, glance, predict, gof) is
# routed through combine_chains() so the multichain object behaves
# identically to a single-chain fit with K * n_draws draws concatenated.

#' Combine the chains of a `dbn_multichain` into a single `dbn` fit
#'
#' Concatenates per-draw operator lists (`A`, `B`, `M`) and scalar
#' parameter vectors (`sigma2`, `tau_A2`, `tau_B2`) across chains. The
#' result is a single `dbn` object with all `K * n_draws` draws stacked,
#' suitable for every standard accessor (`tidy`, `coef`, `predict`,
#' `gof`, `compute_irf`, ...).
#'
#' This is the brms-style "merged-chains" view. Use it when you want a
#' single posterior summary that respects all chains; use the underlying
#' `mc$chains` when you need per-chain access.
#'
#' @param mc A `dbn_multichain` object.
#' @return A single `dbn` fit object with all draws concatenated.
#' @seealso [dbn_multichain()], [check_convergence()] (multichain).
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
#' mc <- dbn_multichain(sim$Z, chains = 2, seeds = c(1, 2),
#'                      model = "dynamic", family = "gaussian",
#'                      nscan = 150, burn = 50, verbose = FALSE)
#' fit_all <- combine_chains(mc)
#' length(fit_all$sigma2)  # 2 * nscan
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
combine_chains <- function(mc) {
	if (!inherits(mc, "dbn_multichain")) {
		cli::cli_abort("{.arg mc} must be a {.cls dbn_multichain} object.")
	}
	if (mc$n_chains < 1L) cli::cli_abort("Empty multichain.")
	if (mc$n_chains == 1L) return(mc$chains[[1L]])
	# use chain 1 as the template and concatenate per-draw structures
	base <- mc$chains[[1L]]
	# A and B are lists of cubes (one per draw)
	if (is.list(base$A)) {
		base$A <- do.call(c, lapply(mc$chains, function(f) f$A))
	}
	if (is.list(base$B)) {
		base$B <- do.call(c, lapply(mc$chains, function(f) f$B))
	}
	# M is [n_row, n_col, p, n_draws] -> concatenate along the 4th axis
	if (is.array(base$M) && length(dim(base$M)) == 4L) {
		M_list <- lapply(mc$chains, function(f) f$M)
		base$M <- do.call(abind::abind, c(M_list, list(along = 4L)))
	}
	# theta is [n_row, n_col, p, Tt, n_draws] -> concatenate along 5th
	if (is.array(base$Theta) && length(dim(base$Theta)) == 5L) {
		Theta_list <- lapply(mc$chains, function(f) f$Theta)
		base$Theta <- do.call(abind::abind, c(Theta_list, list(along = 5L)))
	}
	# scalar parameter vectors
	for (nm in c("sigma2", "sigma2_obs", "tau_A2", "tau_B2", "g2",
	             "rho_A", "rho_B", "rhoA", "rhoB")) {
		if (is.numeric(base[[nm]])) {
			base[[nm]] <- do.call(c, lapply(mc$chains, function(f) f[[nm]]))
		}
	}
	# update kept-draw bookkeeping
	new_n <- if (is.list(base$A)) length(base$A) else NA_integer_
	if (!is.null(base$settings$draws)) base$settings$draws <- new_n
	if (!is.null(base$meta$draws))     base$meta$draws <- new_n
	# record where the draws came from
	base$meta$combined_chains <- mc$n_chains
	base$meta$combined_chain_seeds <- mc$seeds
	base
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
`[[.dbn_multichain` <- function(x, i, ...) {
	# numeric index returns the i-th chain so the multichain behaves
	# like an array of chains; named index falls through to list semantics.
	if (is.numeric(i)) {
		if (i < 1L || i > x$n_chains) {
			cli::cli_abort("Chain index {.val {i}} out of range (1..{x$n_chains}).")
		}
		return(x$chains[[as.integer(i)]])
	}
	# allow named access to the list slots (chains, seeds, n_chains, args)
	# via the original list semantics.
	NextMethod()
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
length.dbn_multichain <- function(x) x$n_chains

#' @author Tosin Salau and Shahryar Minhas
#' @export
`[.dbn_multichain` <- function(x, i, ...) {
	# subset by integer index preserves the dbn_multichain class with the
	# selected chains.
	if (is.numeric(i)) {
		i <- as.integer(i)
		if (any(i < 1L) || any(i > x$n_chains)) {
			cli::cli_abort("Chain index out of range (1..{x$n_chains}).")
		}
		out <- list(
			chains = x$chains[i],
			seeds = x$seeds[i],
			n_chains = length(i),
			args = x$args
		)
		class(out) <- c("dbn_multichain", "list")
		return(out)
	}
	NextMethod()
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
tidy.dbn_multichain <- function(x, level = 0.95, ...) {
	tidy(combine_chains(x), level = level, ...)
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
glance.dbn_multichain <- function(x, ...) {
	g <- glance(combine_chains(x), ...)
	g$n_chains <- x$n_chains
	g
}

#' Methods for `dbn_multichain` objects
#'
#' S3 methods that combine the chains (via [combine_chains()]) and
#' delegate to the corresponding `dbn` method.
#'
#' @param object,fit A `dbn_multichain` object.
#' @param what For `coef`, which operator slot to return.
#' @param ... Forwarded to the underlying single-chain method.
#' @return Whatever the underlying single-chain method returns.
#' @name dbn_multichain-methods
NULL

#' @rdname dbn_multichain-methods
#' @author Tosin Salau and Shahryar Minhas
#' @export
coef.dbn_multichain <- function(object, what = c("all", "A", "B", "M"), ...) {
	coef(combine_chains(object), what = what, ...)
}

#' @rdname dbn_multichain-methods
#' @author Tosin Salau and Shahryar Minhas
#' @export
predict.dbn_multichain <- function(object, ...) {
	predict(combine_chains(object), ...)
}

#' @rdname dbn_multichain-methods
#' @author Tosin Salau and Shahryar Minhas
#' @export
gof.dbn_multichain <- function(fit, ...) {
	gof(combine_chains(fit), ...)
}

####
#' tidybayes-compatible tidy_draws methods for dbn fits
#'
#' @description Returns a tibble of posterior draws shaped for
#'   [tidybayes::spread_draws()] / [tidybayes::gather_draws()]: one row
#'   per `.draw`, with named scalar variables and (optionally)
#'   index-named operator entries like `A[i,j,t]`. This avoids the
#'   `tidy_draws.default(as_draws(fit))` dance and keeps the column
#'   semantics stable across model variants.
#'
#' @param model A `dbn` or `dbn_multichain` fit.
#' @param ... Forwarded to [as_draws.dbn()] /
#'   [as_draws.dbn_multichain()] (e.g. `include_operators = TRUE`).
#' @return A tibble with `.chain`, `.iteration`, `.draw` columns and
#'   one column per posterior variable.
#' @keywords internal
#' @export
tidy_draws.dbn <- function(model, ...) {
	if (!requireNamespace("posterior", quietly = TRUE))
		cli::cli_abort("Package {.pkg posterior} required.")
	posterior::as_draws_df(as_draws.dbn(model, ...))
}

#' @rdname tidy_draws.dbn
#' @author Tosin Salau and Shahryar Minhas
#' @export
tidy_draws.dbn_multichain <- function(model, ...) {
	if (!requireNamespace("posterior", quietly = TRUE))
		cli::cli_abort("Package {.pkg posterior} required.")
	posterior::as_draws_df(as_draws.dbn_multichain(model, ...))
}
