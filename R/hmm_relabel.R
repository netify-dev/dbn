####
# label-switching correction for hmm fits
#
# the gibbs sampler assigns arbitrary labels to regimes at each iteration.
# to compute a meaningful posterior summary of the transition matrix Pi
# we relabel each draw so the regimes are aligned with a reference.
# for R <= 5 we exhaustively search the R! permutations per draw and
# pick the one minimising a Frobenius cost on the sender / receiver
# operators (or a user-supplied alignment statistic). R > 5 returns an
# error directing the user to choose a reference draw manually.
####

####
#' Permutation-aligned posterior summary of the HMM transition matrix
#'
#' @description Corrects label switching in an HMM fit by searching over
#'   the R! permutations of the regime labels and picking, for each
#'   posterior draw, the permutation that minimises the Frobenius
#'   distance between this draw's (A, B) operators and a reference set
#'   (the posterior-mode draw, by default). The returned object reports
#'   the relabeled posterior mean and credible intervals for Pi.
#'
#' @details For \eqn{R} regimes the search is exhaustive over \eqn{R!}
#'   permutations per draw. This is exact and tractable for \eqn{R \le 5}
#'   (\eqn{120} permutations per draw). For larger \eqn{R} the function
#'   refuses to run; manual reference selection or a sequential alignment
#'   (e.g. Stephens, 2000) is required.
#'
#' @param fit An HMM fit produced by \code{dbn(..., model = "hmm")}.
#' @param ref Reference draw selection: \code{"mode"} (the highest
#'   log-likelihood draw, default), \code{"first"} (draw 1), or an
#'   integer index in \code{seq_len(n_draws)}.
#' @param level Coverage for posterior credible intervals (default 0.95).
#' @param verbose Print progress when \code{TRUE} (default \code{FALSE}).
#' @return A list with class \code{dbn_hmm_relabel} containing:
#'   \describe{
#'     \item{\code{Pi_mean}}{R x R posterior mean transition matrix.}
#'     \item{\code{Pi_lo}, \code{Pi_hi}}{R x R credible interval bounds.}
#'     \item{\code{Pi_draws}}{n_draws x R x R relabeled draws of Pi.}
#'     \item{\code{S_relabeled}}{n_draws x Tt regime assignment matrix
#'       after relabeling (each row is the relabeled sequence for one
#'       draw).}
#'     \item{\code{perms}}{n_draws x R permutation indices applied per
#'       draw (\code{perms[d, r]} is the original label that draw d's
#'       regime r was renamed to).}
#'     \item{\code{ref_idx}}{Index of the reference draw used.}
#'   }
#' @author Tosin Salau and Shahryar Minhas
#' @export
hmm_relabel <- function(fit, ref = "mode", level = 0.95, verbose = FALSE) {
	if (!inherits(fit, "dbn") || !identical(fit$model, "hmm"))
		cli::cli_abort("{.arg fit} must be a {.cls dbn} HMM fit ({.code model = \"hmm\"}).")
	R <- fit$R
	if (is.null(R))
		cli::cli_abort("Fit does not record the number of regimes {.code R}.")
	if (R > 5L)
		cli::cli_abort(c(
			"Exhaustive permutation-search relabeling is only supported for {.code R <= 5} ({R}! = {factorial(R)} permutations).",
			"i" = "For {.code R > 5}, choose a sequential alignment (e.g. Stephens 2000) or fix the labels by ordering regimes by a marker statistic before fitting."
		))
	if (level <= 0 || level >= 1)
		cli::cli_abort("{.arg level} must be in (0, 1).")

	n_draws <- length(fit$Pi)
	if (n_draws < 2L)
		cli::cli_abort("Need at least 2 posterior draws to relabel.")

	# enumerate all permutations of {1, ..., R}
	all_perms <- .all_permutations(R)
	n_perms <- nrow(all_perms)

	# choose reference draw
	if (is.numeric(ref)) {
		ref_idx <- as.integer(ref)
		if (ref_idx < 1L || ref_idx > n_draws)
			cli::cli_abort("{.arg ref} index out of range.")
	} else {
		ref_choice <- match.arg(ref, c("mode", "first"))
		if (ref_choice == "first") {
			ref_idx <- 1L
		} else {
			# posterior-mode reference: pick the draw with the largest
			# average log-density of the observed regimes under its Pi
			ref_idx <- .hmm_pick_mode_draw(fit)
		}
	}

	# extract per-draw A / B / Pi / S as arrays for fast alignment
	A_draws <- fit$A
	B_draws <- fit$B
	Pi_draws_raw <- fit$Pi
	S_draws <- fit$S  # n_draws x Tt (regime assignment per time)

	A_ref <- A_draws[[ref_idx]]
	B_ref <- B_draws[[ref_idx]]
	# A and B per draw are 3D arrays (n_row x n_row x R): operator per regime
	# stacked along the third dim
	if (!is.array(A_ref) || length(dim(A_ref)) != 3L || dim(A_ref)[3] != R)
		cli::cli_abort(c(
			"HMM fit's {.code A[[d]]} should be a 3D array {.code [n_row, n_row, R]}.",
			"i" = "Got dims {.val {paste(dim(A_ref), collapse = ' x ')}} with R = {R}."
		))

	# storage for relabeled outputs
	S_len <- length(S_draws[[1L]])
	Pi_draws_rel <- array(NA_real_, dim = c(n_draws, R, R))
	S_rel <- matrix(NA_integer_, nrow = n_draws, ncol = S_len)
	perms_applied <- matrix(NA_integer_, nrow = n_draws, ncol = R)

	# Frobenius cost between operator sets under a permutation sigma:
	# cost(sigma) = sum_r ||A^(d)[, , sigma(r)] - A_ref[, , r]||_F^2 +
	#               ||B^(d)[, , sigma(r)] - B_ref[, , r]||_F^2
	for (d in seq_len(n_draws)) {
		A_d <- A_draws[[d]]
		B_d <- B_draws[[d]]
		Pi_d <- Pi_draws_raw[[d]]
		best_cost <- Inf
		best_perm <- seq_len(R)
		for (k in seq_len(n_perms)) {
			sig <- all_perms[k, ]
			cost <- 0
			for (r in seq_len(R)) {
				cost <- cost + sum((A_d[, , sig[r]] - A_ref[, , r])^2)
				cost <- cost + sum((B_d[, , sig[r]] - B_ref[, , r])^2)
			}
			if (cost < best_cost) {
				best_cost <- cost
				best_perm <- sig
			}
		}
		# relabel: new regime r corresponds to old regime best_perm[r]
		# Pi_new[r1, r2] = Pi_old[best_perm[r1], best_perm[r2]]
		Pi_rel <- Pi_d[best_perm, best_perm, drop = FALSE]
		Pi_draws_rel[d, , ] <- Pi_rel
		# relabel S: if S_old[t] == best_perm[r], then S_new[t] = r
		inv_perm <- order(best_perm)
		S_rel[d, ] <- inv_perm[S_draws[[d]]]
		perms_applied[d, ] <- best_perm
		if (verbose && d %% max(1L, floor(n_draws / 20)) == 0L)
			cli::cli_alert_info("Relabeled draw {d}/{n_draws}; cost = {round(best_cost, 4)}")
	}

	alpha <- (1 - level) / 2
	Pi_mean <- apply(Pi_draws_rel, c(2, 3), mean)
	Pi_lo   <- apply(Pi_draws_rel, c(2, 3), stats::quantile, probs = alpha)
	Pi_hi   <- apply(Pi_draws_rel, c(2, 3), stats::quantile, probs = 1 - alpha)

	rn <- paste0("regime_", seq_len(R))
	dimnames(Pi_mean) <- list(rn, rn)
	dimnames(Pi_lo)   <- list(rn, rn)
	dimnames(Pi_hi)   <- list(rn, rn)

	out <- list(
		Pi_mean = Pi_mean,
		Pi_lo   = Pi_lo,
		Pi_hi   = Pi_hi,
		Pi_draws = Pi_draws_rel,
		S_relabeled = S_rel,
		perms = perms_applied,
		ref_idx = ref_idx,
		R = R,
		level = level
	)
	class(out) <- c("dbn_hmm_relabel", "list")
	out
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_hmm_relabel <- function(x, digits = 3L, ...) {
	cli::cli_h2("HMM transition matrix (label-switching corrected)")
	cli::cli_text("Regimes: {.val {x$R}}; reference draw: {.val {x$ref_idx}}")
	cli::cli_text("Posterior mean Pi:")
	print(round(x$Pi_mean, digits))
	cli::cli_text("{format(100 * x$level, nsmall = 0)}% CI bounds:")
	cat("Lower:\n"); print(round(x$Pi_lo, digits))
	cat("Upper:\n"); print(round(x$Pi_hi, digits))
	invisible(x)
}

####
# enumerate all permutations of seq_len(R)
.all_permutations <- function(R) {
	if (R == 1L) return(matrix(1L, 1, 1))
	perms <- matrix(0L, factorial(R), R)
	# recursive enumeration
	build <- function(remaining, prefix) {
		if (length(remaining) == 0L) {
			return(matrix(prefix, nrow = 1L))
		}
		out <- NULL
		for (k in seq_along(remaining)) {
			sub <- build(remaining[-k], c(prefix, remaining[k]))
			out <- if (is.null(out)) sub else rbind(out, sub)
		}
		out
	}
	build(seq_len(R), integer(0))
}

####
# pick the posterior-mode draw: simple heuristic based on the trace of Pi
# (regimes with stronger self-persistence drive a more "modal" decomposition)
.hmm_pick_mode_draw <- function(fit) {
	n_draws <- length(fit$Pi)
	scores <- numeric(n_draws)
	for (d in seq_len(n_draws)) {
		scores[d] <- sum(diag(fit$Pi[[d]]))
	}
	which.max(scores)
}
