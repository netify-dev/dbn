####
#' Held-out link-prediction metrics for a dbn fit
#'
#' @description Computes the standard binary-classification metrics
#'   used in the graph-ML literature (AUC, average precision,
#'   precision@K, accuracy at a given threshold) against held-out dyads
#'   or a held-out time slice. Designed for users who want to compare
#'   dbn to GNN baselines on the same evaluation protocol.
#'
#' @param fit A fitted `dbn` object (typically `family = "binary"` for
#'   AUC/AP semantics; other families produce continuous predictions
#'   and the metrics are computed against the binarized truth at
#'   `threshold_truth`).
#' @param Y_truth Either a 4D `[n_row, n_col, p, T]` array of truth
#'   (same shape as `fit$Y`) or `NULL`. If `NULL`, `fit$Y` is used and
#'   the user MUST provide `held_out_mask` or `held_out_time` so we
#'   don't score against the training data.
#' @param held_out_mask Logical 4D array marking held-out cells (TRUE
#'   = held out). Mutually exclusive with `held_out_time`.
#' @param held_out_time Integer vector of held-out time indices.
#' @param threshold_truth Numeric. For non-binary families, dichotomize
#'   `Y_truth` at this threshold to define positive ties. Ignored for
#'   `family = "binary"`.
#' @param K Integer or integer vector of K values for precision@K.
#'   Default `c(10, 50)`.
#' @return A list with components `auc`, `ap` (average precision),
#'   `accuracy`, `precision_at_K` (named numeric vector), `n_pos`,
#'   `n_neg`, and `held_out` (the mask actually used). Returns NA
#'   metrics with a directive when there are fewer than 2 positives.
#' @author Tosin Salau and Shahryar Minhas
#' @export
link_prediction_metrics <- function(fit, Y_truth = NULL,
                                    held_out_mask = NULL,
                                    held_out_time = NULL,
                                    threshold_truth = 0,
                                    K = c(10L, 50L)) {
	if (!inherits(fit, "dbn"))
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	if (is.null(Y_truth)) Y_truth <- fit$Y
	if (is.null(Y_truth) || length(dim(Y_truth)) != 4L)
		cli::cli_abort("{.arg Y_truth} (or {.code fit$Y}) must be a 4D array of shape [n_row, n_col, p, T].")
	# require a held-out specification
	if (is.null(held_out_mask) && is.null(held_out_time))
		cli::cli_abort(c(
			"You must specify either {.arg held_out_mask} or {.arg held_out_time}.",
			"i" = "Scoring against training data inflates AUC; pass a held-out specification."
		))
	dY <- dim(Y_truth)
	mask <- array(FALSE, dim = dY)
	if (!is.null(held_out_mask)) {
		if (!identical(dim(held_out_mask), dY))
			cli::cli_abort("{.arg held_out_mask} must have the same shape as {.arg Y_truth}.")
		mask <- as.logical(held_out_mask)
		dim(mask) <- dY
	} else {
		# mark all cells in held-out time slices
		t_idx <- as.integer(held_out_time)
		if (any(t_idx < 1L) || any(t_idx > dY[4L]))
			cli::cli_abort("Some {.arg held_out_time} indices are out of range.")
		for (t in t_idx) mask[, , , t] <- TRUE
	}
	# exclude diagonal cells and NA truth from scoring
	if (dY[1L] == dY[2L])
		for (r in seq_len(dY[3L])) for (t in seq_len(dY[4L])) diag(mask[, , r, t]) <- FALSE
	mask[is.na(Y_truth)] <- FALSE
	if (sum(mask) == 0L)
		cli::cli_abort("No scorable cells after applying mask and NA filter.")
	# predicted scores: forecast on held-out future, in-sample posterior mean otherwise
	fam <- if (is.list(fit$family)) fit$family$name else fit$family
	score <- array(NA_real_, dim = dY)
	if (!is.null(held_out_time) && all(held_out_time > max(setdiff(seq_len(dY[4L]), held_out_time)))) {
		# pure forward forecast
		H <- length(held_out_time)
		pred <- stats::predict(
			fit,
			H = H,
			draws = min(100L, length(fit$A) %||% 50L),
			summary = "mean"
		)
		for (i in seq_along(held_out_time)) score[, , , held_out_time[i]] <- pred[, , , i]
	} else {
		# in-sample posterior mean
		Th <- if (is.array(fit$Theta) && length(dim(fit$Theta)) == 5L) apply(fit$Theta, 1:4, mean) else fit$Theta
		Mm <- if (is.array(fit$M) && length(dim(fit$M)) == 4L) apply(fit$M, 1:3, mean) else fit$M
		for (r in seq_len(dY[3L])) for (t in seq_len(dY[4L])) {
			score[, , r, t] <- (Th[, , r, t] %||% 0) + (Mm[, , r] %||% 0)
		}
	}
	# binary labels
	truth <- if (identical(fam, "binary")) Y_truth > 0.5 else Y_truth > threshold_truth
	y <- as.integer(truth[mask])
	s <- as.numeric(score[mask])
	keep <- !is.na(s) & !is.na(y)
	y <- y[keep]
	s <- s[keep]
	n_pos <- sum(y == 1L)
	n_neg <- sum(y == 0L)
	if (n_pos < 2L || n_neg < 2L) {
		cli::cli_warn(c(
			"Fewer than 2 positives or negatives after masking.",
			"i" = "n_pos = {n_pos}, n_neg = {n_neg}. Returning NA metrics."
		))
		return(list(
			auc = NA_real_,
			ap = NA_real_,
			accuracy = NA_real_,
			precision_at_K = stats::setNames(rep(NA_real_, length(K)), paste0("P@", K)),
			n_pos = n_pos,
			n_neg = n_neg,
			held_out = mask
		))
	}
	# auc via mann-whitney u normalization
	rank_y <- rank(s, ties.method = "average")
	auc <- (sum(rank_y[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
	# average precision
	o <- order(s, decreasing = TRUE)
	cum_tp <- cumsum(y[o] == 1L)
	prec <- cum_tp / seq_along(cum_tp)
	rec <- cum_tp / n_pos
	delta <- diff(c(0, rec))
	ap <- sum(prec * delta)
	# accuracy at threshold
	pred_pos <- if (identical(fam, "binary")) s > 0.5 else s > threshold_truth
	acc <- mean(pred_pos == (y == 1L))
	# precision at k
	p_at_k <- vapply(K, function(k) {
		k_eff <- min(k, length(o))
		mean(y[o[seq_len(k_eff)]] == 1L)
	}, numeric(1L))
	names(p_at_k) <- paste0("P@", K)
	list(
		auc = auc,
		ap = ap,
		accuracy = acc,
		precision_at_K = p_at_k,
		n_pos = n_pos,
		n_neg = n_neg,
		held_out = mask
	)
}

####
#' Coupling trajectory with credible bands
#'
#' @description Returns a tidy data frame of per-actor coupling
#'   `C_i(t) = ||W_t[i, ]||_2` over time with posterior credible
#'   intervals. Unlike [actor_embedding()], which reports a single
#'   per-(actor, time) point estimate, this helper exposes the per-draw
#'   trajectory so you can plot ribbons.
#'
#' @param fit A dynamic / piecewise / lowrank / hmm fit with operator
#'   draws stored under `fit$A`, `fit$B`.
#' @param actors Optional integer vector of actor indices, or character
#'   vector of actor names. Default: all actors.
#' @param level Credible level (default 0.95).
#' @return A data frame with columns `actor`, `t`, `mean`, `lo`, `hi`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'            nscan = 200, burn = 100, verbose = FALSE)
#' ct <- coupling_trajectory(fit)
#' head(ct)
#' }
coupling_trajectory <- function(fit, actors = NULL, level = 0.95) {
	if (!inherits(fit, "dbn"))
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	A_list <- fit$A; B_list <- fit$B
	if (!is.list(A_list) || length(A_list) == 0L)
		cli::cli_abort(c(
			"{.fun coupling_trajectory} requires per-draw operators on the fit.",
			"i" = "Got {.val {fit$model}} with no posterior draws of {.code A}."
		))
	# resolve actor labels
	n_row <- dim(A_list[[1L]])[1L]
	if (length(dim(A_list[[1L]])) == 3L) Tt <- dim(A_list[[1L]])[3L] else Tt <- 1L
	actor_names <- dimnames(fit$Y)[[1L]] %||% paste0("a", seq_len(n_row))
	if (!is.null(actors)) {
		if (is.character(actors)) {
			idx <- match(actors, actor_names)
			if (anyNA(idx)) cli::cli_abort("Some {.arg actors} not found in {.code dimnames(fit$Y)[[1]]}.")
			actors <- idx
		}
		actors <- as.integer(actors)
	} else {
		actors <- seq_len(n_row)
	}
	# per-draw, per-time coupling for selected actors
	n_draws <- length(A_list)
	C_draws <- array(NA_real_, dim = c(n_draws, length(actors), Tt))
	for (m in seq_len(n_draws)) {
		Am <- A_list[[m]]
		Bm <- B_list[[m]]
		for (t in seq_len(Tt)) {
			At <- if (length(dim(Am)) == 3L) Am[, , min(t, dim(Am)[3L])] else Am
			Bt <- if (length(dim(Bm)) == 3L) Bm[, , min(t, dim(Bm)[3L])] else Bm
			Wt <- At %*% t(Bt)
			C_draws[m, , t] <- sqrt(rowSums(Wt[actors, , drop = FALSE]^2))
		}
	}
	alpha <- (1 - level) / 2
	out <- expand.grid(
		actor = actor_names[actors],
		t = seq_len(Tt),
		stringsAsFactors = FALSE
	)
	out$mean <- NA_real_
	out$lo <- NA_real_
	out$hi <- NA_real_
	for (i in seq_along(actors)) {
		for (t in seq_len(Tt)) {
			vals <- C_draws[, i, t]
			row_idx <- (t - 1L) * length(actors) + i
			out$mean[row_idx] <- mean(vals, na.rm = TRUE)
			out$lo[row_idx]   <- stats::quantile(vals, alpha, na.rm = TRUE)
			out$hi[row_idx]   <- stats::quantile(vals, 1 - alpha, na.rm = TRUE)
		}
	}
	out
}

####
#' Compact one-line representation of a dbn fit
#'
#' @description Returns a single-line description suitable for
#'   stakeholder-facing tables. Use with `knitr::kable(format(fit))`
#'   for inline reporting.
#'
#' @param x A `dbn` fit.
#' @param ... Unused.
#' @return A character scalar.
#' @author Tosin Salau and Shahryar Minhas
#' @export
format.dbn <- function(x, ...) {
	mdl <- x$model %||% "?"
	fam <- if (is.list(x$family)) x$family$name else x$family %||% "?"
	dims <- x$meta$dims %||% x$dims
	su <- x$meta$sampler_used %||% "?"
	sprintf("dbn(model=%s, family=%s, n=%d, T=%d, sampler=%s)",
	        mdl, fam,
	        dims$n_row %||% NA_integer_,
	        dims$Tt %||% NA_integer_,
	        su)
}

####
#' Posterior communities from a low-rank dbn fit
#'
#' @description Clusters actors by their latent-space positions
#'   `node_embedding(fit)` into `k` communities using k-means. The
#'   columns of `U` are sign-aligned across draws before clustering so
#'   the result is invariant to the per-draw sign indeterminacy.
#'
#'   Closer in spirit to a low-rank spectral embedding than to a
#'   stochastic block model: actors in the same community share
#'   similar influence-row patterns rather than identical edge
#'   probabilities. Use `dynsbm` or `blockmodels` for SBM-proper
#'   inference; this helper is for descriptive partitioning.
#'
#' @param fit A `model = "lowrank"` fit.
#' @param k Integer number of communities. Defaults to `ceiling(sqrt(n))`.
#' @param ... Forwarded to [stats::kmeans()].
#' @return A list with `membership` (named integer vector of community
#'   labels), `centers` (k x r matrix of cluster centroids), and `U`
#'   (the n x r embedding clustered).
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_lowrank_dbn(n = 10, p = 1, time = 8, r = 2, seed = 1)
#' fit <- dbn(sim$Z, model = "lowrank", family = "gaussian", r = 2,
#'            nscan = 200, burn = 100, verbose = FALSE)
#' pc <- posterior_communities(fit, k = 3)
#' pc$membership
#' }
posterior_communities <- function(fit, k = NULL, ...) {
	if (!inherits(fit, "dbn"))
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	if (!identical(fit$model, "lowrank"))
		cli::cli_abort(c(
			"{.fun posterior_communities} requires a low-rank fit.",
			"i" = "For dynamic / static fits, cluster on {.code dbn_operator(fit)$W_abs_mean} directly."
		))
	U <- node_embedding(fit, align = TRUE)
	n <- nrow(U)
	if (is.null(k)) k <- max(2L, as.integer(ceiling(sqrt(n))))
	if (k < 2L || k > n)
		cli::cli_abort("{.arg k} must satisfy {.code 2 <= k <= n_actors}.")
	km <- stats::kmeans(U, centers = k, ...)
	# preserve actor labels
	mem <- as.integer(km$cluster)
	if (!is.null(rownames(U))) names(mem) <- rownames(U)
	list(membership = mem, centers = km$centers, U = U)
}
