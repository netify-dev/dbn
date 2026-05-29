####
#' Build Shock Matrix for IRF Analysis
#'
#' @description Creates structured perturbation matrices for impulse response
#'   analysis. The shock matrix S is added to the latent network at a single
#'   time point, and [compute_irf()] tracks how the perturbation propagates.
#'
#' @details
#' **When to use each shock type:**
#' - `"unit_edge"`: Shock a single bilateral tie (e.g., "what if USA-Russia
#'   relations improve?"). Use when you care about a specific dyad.
#' - `"node_out"`: Shock all outgoing ties from one actor (e.g., "what if
#'   China engages all partners at once?"). Use when you care about one
#'   actor's overall engagement.
#' - `"node_in"`: Shock all incoming ties to one actor.
#' - `"density"`: Apply a uniform shock to all edges. Use to study how
#'   a system-wide shift propagates.
#'
#' For **symmetric (undirected) networks**, symmetrize the shock after
#' building it: `S = S + t(S)`.
#'
#' @param m Number of actors (senders), or a fitted `dbn` object. When a
#'   fitted object is passed, the network size and actor names are read off
#'   the fit so `i` and `j` may be supplied as character strings.
#' @param type Type of shock:
#'   \itemize{
#'     \item `"unit_edge"`: perturb a single directed edge (i -> j)
#'     \item `"node_out"`: perturb all outgoing edges from actor i
#'     \item `"node_in"`: perturb all incoming edges to actor i
#'     \item `"density"`: uniform perturbation to all edges
#'   }
#' @param i Source actor (sender). Integer index or character name; the
#'   latter is matched against `dimnames(fit$Y)[[1]]` when `m` is a fit.
#' @param j Target actor (receiver), used only for `"unit_edge"`. Integer
#'   index or character name (matched against `dimnames(fit$Y)[[2]]`).
#' @param magnitude Size of the perturbation. Scale relative to your data:
#'   e.g., if data ranges from -1 to 1, a magnitude of 0.5 is a large shock.
#' @param n_col Number of receiver nodes (default: same as `m`). Ignored
#'   when `m` is a fitted `dbn` object.
#' @return An `m x n_col` matrix of zeros with the specified entries set to
#'   `magnitude`. Pass this to [compute_irf()] as the `shock` argument.
#' @seealso \code{\link{compute_irf}}, \code{\link{plot.dbn_irf}}
#' @examples
#' # Integer indices (legacy interface)
#' S <- build_shock(m = 5, type = "unit_edge", i = 1, j = 2)
#'
#' # Named actors via a fitted object
#' \dontrun{
#'   Y <- simulate_dynamic_dbn(n = 5, time = 8)$Z
#'   dimnames(Y) <- list(c("USA","UK","CHN","RUS","DEU"),
#'                       c("USA","UK","CHN","RUS","DEU"), "tie", NULL)
#'   fit <- dbn(Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#'   S <- build_shock(fit, type = "unit_edge", i = "USA", j = "CHN")
#' }
#'
#' # Node-level shock: activate all outgoing edges from node 1
#' S_out <- build_shock(m = 5, type = "node_out", i = 1)
#'
#' # Density shock: small shock spread across the network
#' S_dens <- build_shock(m = 5, type = "density", magnitude = 0.1)
#'
#' # Bipartite shock: 4 senders, 6 receivers
#' S_bip <- build_shock(m = 4, n_col = 6, type = "unit_edge", i = 1, j = 3)
#' @author Tosin Salau and Shahryar Minhas
#' @export
build_shock <- function(m, type = c("unit_edge", "node_out", "node_in", "density"),
						i = 1, j = 2, magnitude = 1, n_col = m) {
	type <- match.arg(type)
	# accept a fitted dbn so i / j can be character actor names
	actor_names_row <- NULL
	actor_names_col <- NULL
	if (inherits(m, "dbn")) {
		fit <- m
		dims <- fit$meta$dims %||% fit$dims
		m <- dims$n_row %||% dim(fit$Y)[1L]
		n_col <- dims$n_col %||% dim(fit$Y)[2L]
		actor_names_row <- dimnames(fit$Y)[[1]]
		actor_names_col <- dimnames(fit$Y)[[2]]
	}
	# resolve character actor names against fit dimnames
	resolve_actor <- function(x, axis_names, axis_label) {
		if (is.character(x)) {
			if (is.null(axis_names)) {
				cli::cli_abort(c(
					"{.arg {axis_label}} is a character string but no actor names are available.",
					"i" = "Pass {.arg m = fit} (a fitted dbn) so actor names can be resolved, or use integer indices."
				))
			}
			idx <- match(x, axis_names)
			if (anyNA(idx)) {
				avail <- if (length(axis_names) > 6) paste0(paste(head(axis_names, 6), collapse = ", "), ", ...")
				         else paste(axis_names, collapse = ", ")
				cli::cli_abort(c(
					"{.arg {axis_label}} actor name(s) not found in fit dimnames.",
					"x" = "Unmatched: {.val {x[is.na(idx)]}}",
					"i" = "Available: {avail}"
				))
			}
			return(idx)
		}
		x
	}
	i <- resolve_actor(i, actor_names_row, "i")
	j <- resolve_actor(j, actor_names_col, "j")

	if (length(m) != 1L || !is.numeric(m) || !is.finite(m) || m < 1 || m != round(m)) {
		cli::cli_abort("{.arg m} (number of senders) must be a single positive integer or a fitted dbn.")
	}
	if (length(n_col) != 1L || !is.numeric(n_col) || !is.finite(n_col) || n_col < 1 || n_col != round(n_col)) {
		cli::cli_abort("{.arg n_col} (number of receivers) must be a single positive integer.")
	}
	if (length(magnitude) != 1L || !is.numeric(magnitude) || !is.finite(magnitude)) {
		cli::cli_abort("{.arg magnitude} must be a single finite number, got {.val {magnitude}}.")
	}
	n_row <- m
	S <- matrix(0, n_row, n_col)
	is_bipartite <- (n_row != n_col)

	switch(type,
		unit_edge = {
			if (i < 1 || i > n_row || j < 1 || j > n_col) {
				cli::cli_abort(c(
					"Node indices out of range.",
					"i" = "{.arg i} must be between 1 and {.val {n_row}}, got {.val {i}}.",
					"i" = "{.arg j} must be between 1 and {.val {n_col}}, got {.val {j}}."
				))
			}
			S[i, j] <- magnitude
		},
		node_out = {
			if (i < 1 || i > n_row) cli::cli_abort(c(
				"Sender index out of range.",
				"i" = "{.arg i} must be between 1 and {.val {n_row}}, got {.val {i}}."
			))
			S[i, ] <- magnitude
		},
		node_in = {
			if (i < 1 || i > n_col) cli::cli_abort(c(
				"Receiver index out of range.",
				"i" = "{.arg i} must be between 1 and {.val {n_col}}, got {.val {i}}."
			))
			S[, i] <- magnitude
		},
		density = {
			n_cells <- if (is_bipartite) n_row * n_col else n_row * (n_row - 1)
			S[, ] <- magnitude / n_cells
			if (!is_bipartite) diag(S) <- 0
		}
	)
	S
}
####

####
#' Network Statistic: Density
#'
#' @description Mean of off-diagonal elements.
#'
#' **Note:** \code{ggplot2} also exports a function called \code{stat_density}.
#' If both packages are loaded, use \code{dbn::stat_density} to refer to this
#' function unambiguously (e.g., when passing to \code{\link{compute_irf}}).
#'
#' @param X Network matrix
#' @return Scalar density value
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_in_degree}}, \code{\link{stat_out_degree}},
#'   \code{\link{stat_reciprocity}}, \code{\link{stat_transitivity}}
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' dbn::stat_density(X)
#' @author Tosin Salau and Shahryar Minhas
#' @export
stat_density <- function(X) {
	# for square (unipartite) networks the diagonal is excluded as self-ties;
	# for bipartite (rectangular) networks every cell is a valid sender-receiver
	# pair and the off-diagonal mask would incorrectly drop n_col valid ties.
	if (nrow(X) != ncol(X)) return(mean(X, na.rm = TRUE))
	mean(X[row(X) != col(X)], na.rm = TRUE)
}

#' Network Statistic: In-Degree
#'
#' @description Column sums of network matrix.  Includes the diagonal entry
#'   for each node.  In typical DBN usage the diagonal is `NA` or zero, so
#'   the result matches the conventional in-degree.  If your matrix has
#'   non-zero diagonal entries, subtract them manually.
#' @param X Network matrix
#' @return Vector of in-degrees
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_density}}, \code{\link{stat_out_degree}}
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' stat_in_degree(X)
#' @author Tosin Salau and Shahryar Minhas
#' @export
stat_in_degree <- function(X) {
	# require a 2D matrix; a 4D [n,n,p,T] input would otherwise pass through
	# colSums and silently return an n*p*T-vector.
	if (!is.matrix(X) && !(is.array(X) && length(dim(X)) == 2L)) {
		cli::cli_abort(c(
			"{.fun stat_in_degree} expects a 2D matrix (a single network snapshot).",
			"x" = "Got {.cls {class(X)[1]}} with dim {.val {paste(dim(X) %||% length(X), collapse = 'x')}}.",
			"i" = "If you have a 4D `[n_row, n_col, p, T]` array, apply per slice: {.code apply(Y, c(3, 4), stat_in_degree)}."
		))
	}
	colSums(X, na.rm = TRUE)
}

#' Network Statistic: Out-Degree
#'
#' @description Row sums of network matrix.  Includes the diagonal entry
#'   for each node.  In typical DBN usage the diagonal is `NA` or zero, so
#'   the result matches the conventional out-degree.  If your matrix has
#'   non-zero diagonal entries, subtract them manually.
#' @param X Network matrix
#' @return Vector of out-degrees
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_density}}, \code{\link{stat_in_degree}}
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' stat_out_degree(X)
#' @author Tosin Salau and Shahryar Minhas
#' @export
stat_out_degree <- function(X) {
	if (!is.matrix(X) && !(is.array(X) && length(dim(X)) == 2L)) {
		cli::cli_abort(c(
			"{.fun stat_out_degree} expects a 2D matrix (a single network snapshot).",
			"x" = "Got {.cls {class(X)[1]}} with dim {.val {paste(dim(X) %||% length(X), collapse = 'x')}}.",
			"i" = "If you have a 4D `[n_row, n_col, p, T]` array, apply per slice: {.code apply(Y, c(3, 4), stat_out_degree)}."
		))
	}
	rowSums(X, na.rm = TRUE)
}

#' Network Statistic: Reciprocity
#'
#' @description Edgewise matched-edge ratio
#'   \eqn{\sum_{i \ne j} X_{ij} X_{ji} / \sum_{i \ne j} X_{ij}^2}, which
#'   equals 1 for a fully symmetric network and 0 for a fully
#'   antisymmetric one. Matches \code{igraph::reciprocity(mode = "default")}
#'   and \code{sna::grecip(g, measure = "edgewise")}; differs from the
#'   dyad-census reciprocity (\code{sna::grecip(g, measure = "dyadic")}).
#'   Only defined for square (unipartite) networks. For bipartite networks,
#'   use \code{\link{stat_density}}, \code{\link{stat_in_degree}}, or
#'   \code{\link{stat_out_degree}} instead.
#' @param X Square network matrix
#' @return Scalar reciprocity value
#' @seealso \code{\link{stat_reciprocity_cor}} for the triangle-correlation
#'   form; \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_density}}, \code{\link{stat_transitivity}}
#' @references Garlaschelli, D., & Loffredo, M. I. (2004). Patterns of
#'   link reciprocity in directed networks. \emph{Physical Review
#'   Letters}, 93(26), 268701.
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' stat_reciprocity(X)
#' @author Tosin Salau and Shahryar Minhas
#' @export
stat_reciprocity <- function(X) {
	# edgewise reciprocity for a (possibly weighted) directed network:
	# sum(X * t(X)) / sum(X^2). returns 1 for fully symmetric networks,
	# 0 for fully antisymmetric. matches igraph::reciprocity default and
	# sna::grecip(measure="edgewise"); see stat_reciprocity_cor() for the
	# triangle-correlation form.
	if (nrow(X) != ncol(X)) {
		cli::cli_abort(c(
			"Reciprocity is undefined for bipartite (rectangular) networks.",
			"i" = "Use {.fn stat_density}, {.fn stat_in_degree}, or {.fn stat_out_degree} instead."
		))
	}
	# exclude diagonal (self-ties don't enter reciprocity).
	off_diag <- row(X) != col(X)
	X_off <- X
	X_off[!off_diag] <- 0
	X_off[is.na(X_off)] <- 0
	num <- sum(X_off * t(X_off))
	den <- sum(X_off^2)
	# undefined on an empty network
	if (den < .Machine$double.eps) return(NA_real_)
	# trivially 1 on networks that are exactly symmetric by construction
	if (isTRUE(all.equal(X_off, t(X_off), tolerance = sqrt(.Machine$double.eps),
	                     check.attributes = FALSE)))
		return(NA_real_)
	num / den
}

#' Network Statistic: Reciprocity (triangle correlation form)
#'
#' Correlation between the upper and lower triangles of a square matrix.
#' This is the *correlation* between forward and reverse ties, NOT the
#' standard SNA reciprocity (which is `sum(X * t(X)) / sum(X^2)`; see
#' \code{\link{stat_reciprocity}}). Provided for users who want the
#' previous behaviour explicitly.
#'
#' @param X Network matrix.
#' @return Scalar correlation, or 0 if either triangle has zero variance.
#' @seealso \code{\link{stat_reciprocity}} for the standard form.
#' @author Tosin Salau and Shahryar Minhas
#' @export
stat_reciprocity_cor <- function(X) {
	if (nrow(X) != ncol(X)) {
		cli::cli_abort(c(
			"Reciprocity is undefined for bipartite (rectangular) networks.",
			"i" = "Use {.fn stat_density}, {.fn stat_in_degree}, or {.fn stat_out_degree} instead."
		))
	}
	upper_idx <- upper.tri(X)
	upper_vals <- X[upper_idx]
	lower_vals <- t(X)[upper_idx]
	tri_eps <- .Machine$double.eps ^ 0.5
	if (length(upper_vals) > 0L && max(abs(upper_vals - lower_vals), na.rm = TRUE) < tri_eps) {
		cli::cli_warn(c(
			"{.fun stat_reciprocity_cor} on a symmetric matrix is trivially 1 (the two triangles are identical).",
			"i" = "If you fit with {.code symmetric = TRUE}, reciprocity carries no information."
		))
		return(1)
	}
	eps <- .Machine$double.eps ^ 0.5
	if (length(upper_vals) < 2L ||
	    sd(upper_vals, na.rm = TRUE) < eps ||
	    sd(lower_vals, na.rm = TRUE) < eps) return(0)
	cor(upper_vals, lower_vals, use = "complete.obs")
}

#' Network Statistic: Transitivity
#'
#' @description Clustering coefficient (fraction of closed triangles).
#'   Only defined for square (unipartite) networks. For bipartite networks,
#'   use \code{\link{stat_density}}, \code{\link{stat_in_degree}}, or
#'   \code{\link{stat_out_degree}} instead.
#'
#'   Implementation note: the input matrix is binarised at \code{threshold}
#'   and then symmetrised by \code{X_bin <- pmax(X_bin, t(X_bin))} before
#'   the clustering coefficient is computed. This is the standard
#'   undirected definition (3 x triangles / connected-triples). It will
#'   not match \code{sna::gtrans(..., measure = "weak"/"strong")} or
#'   \code{ergm::ttriple}, which use directed-triangle definitions on the
#'   un-symmetrised graph.
#'
#'   Transitivity is a property of a binary graph. Continuous-valued latent
#'   states (e.g. `Theta` from a Gaussian fit) must be thresholded into edges
#'   before transitivity is defined. Pass \code{threshold} (the cut for
#'   "this is an edge") or pre-binarize the input.
#' @param X Square network matrix.
#' @param threshold Numeric scalar: an edge is recorded when
#'   \code{abs(X) > threshold}. Defaults to \code{0} for a binary input
#'   (i.e., any non-zero is treated as an edge). For continuous inputs, pass
#'   an explicit threshold (e.g., the off-diagonal median, or a domain-derived
#'   value); the default of 0 would otherwise binarize every non-zero entry
#'   into an edge and return 1 trivially.
#' @return Scalar transitivity value in `[0, 1]`.
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_density}}, \code{\link{stat_reciprocity}}
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' stat_transitivity(X)
#' # continuous input: choose a threshold
#' set.seed(1)
#' Y <- matrix(rnorm(36), 6, 6); diag(Y) <- 0
#' stat_transitivity(Y, threshold = quantile(abs(Y), 0.7))
#' @author Tosin Salau and Shahryar Minhas
#' @export
stat_transitivity <- function(X, threshold = 0) {
	if (nrow(X) != ncol(X)) {
		cli::cli_abort(c(
			"Transitivity is undefined for bipartite (rectangular) networks.",
			"i" = "Use {.fn stat_density}, {.fn stat_in_degree}, or {.fn stat_out_degree} instead."
		))
	}
	# warn loudly when the input is continuous and threshold stays at 0:
	# any non-degenerate continuous matrix would binarise to the all-ones
	# graph and return 1. throttled to once per session to avoid emitting
	# T identical warnings from gof() over T time slices.
	if (threshold == 0 && !all(X %in% c(0, 1, -1, NA))) {
		cli::cli_warn(c(
			"{.fun stat_transitivity} was called with continuous input and {.code threshold = 0}.",
			"i" = "Every non-zero entry counts as an edge, so the binarized graph is fully connected and transitivity is trivially 1.",
			"i" = "Pass an explicit {.arg threshold} (e.g., {.code quantile(abs(X), 0.7)}) or pre-binarize the input."
		), .frequency = "once", .frequency_id = "dbn_stat_transitivity_continuous")
	}
	X_binary <- (abs(X) > threshold) * 1
	# treat as undirected: a directed pair with even one edge counts as a tie
	X_binary <- pmax(X_binary, t(X_binary))
	diag(X_binary) <- 0
	triangles <- sum(diag(X_binary %*% X_binary %*% X_binary)) / 6
	degrees <- rowSums(X_binary)
	triples <- sum(degrees * (degrees - 1)) / 2
	if (triples == 0) return(0)
	# standard global clustering coefficient: 3 * triangles / connected-triples
	3 * triangles / triples
}
####

####
#' Number of posterior draws stored in a static fit
#'
#' @description
#' Robustly determines the posterior draw count for a static `dbn` fit. The
#' static fit stores `B` as a list whose first element is an
#' `[n_col, n_col, draws]` array; older or alternative layouts store `B` as a
#' `[draws, n, n]` array or expose only `params`. Returns `NA_integer_` if no
#' usable component is found.
#'
#' @param fit A static `dbn` fit object.
#' @return Integer draw count, or `NA_integer_`.
#' @keywords internal
#' @noRd
static_draw_count <- function(fit) {
	if (!is.null(fit$B)) {
		if (is.list(fit$B) && length(fit$B) >= 1L && !is.null(dim(fit$B[[1]]))) {
			d <- dim(fit$B[[1]])
			if (length(d) >= 3L) return(as.integer(d[3]))
		}
		if (!is.null(dim(fit$B)) && length(dim(fit$B)) >= 1L) {
			return(as.integer(dim(fit$B)[1]))
		}
	}
	if (!is.null(fit$params) && !is.null(nrow(fit$params))) {
		return(as.integer(nrow(fit$params)))
	}
	NA_integer_
}
####

####
#' Compute IRF for a Single Posterior Draw
#'
#' @param fit dbn model fit object
#' @param draw_idx Index of posterior draw
#' @param shock Shock matrix
#' @param H Number of horizons
#' @param t0 Shock time (for dynamic models, 1-based)
#' @param stat_fun Network statistic function
#' @return Vector of IRF values at each horizon
#' @keywords internal
compute_irf_single <- function(fit, draw_idx, shock, H, t0 = 1,
                                stat_fun = stat_density, stat_n_elem = 1L) {
	dims <- fit$dims
	n_row <- as.integer(dims$n_row)
	n_col <- as.integer(dims$n_col)
	p <- as.integer(dims$p)

	# handle dynamic model with time-varying A, B
	if (fit$model == "dynamic") {
		if (is.null(fit$A) || is.null(fit$B)) {
			cli::cli_abort(c(
				"Dynamic model requires {.code A} and {.code B} matrices.",
				"i" = "Ensure the model was fit with {.code model = \"dynamic\"}."
			))
		}

		time_thin <- fit$time_thin %||% 1
		T_stored <- dim(fit$A[[draw_idx]])[3]
		t0_stored <- if (time_thin > 1) ceiling(t0 / time_thin) else t0

		# allow out-of-sample IRFs (t0 + H > T) by holding A_T, B_T constant
		# beyond T -- a conservative "operator stays at its last observed
		# value" extension. predict(fit, H = N) propagates the RW
		# innovation for a fully honest forecast.
		A_array <- fit$A[[draw_idx]]
		B_array <- fit$B[[draw_idx]]
		if (t0_stored + H > T_stored) {
			n_pad <- t0_stored + H - T_stored
			A_pad <- array(NA_real_, dim = c(dim(A_array)[1], dim(A_array)[2], n_pad))
			B_pad <- array(NA_real_, dim = c(dim(B_array)[1], dim(B_array)[2], n_pad))
			for (kk in seq_len(n_pad)) {
				A_pad[, , kk] <- A_array[, , T_stored]
				B_pad[, , kk] <- B_array[, , T_stored]
			}
			A_array <- abind::abind(A_array, A_pad, along = 3)
			B_array <- abind::abind(B_array, B_pad, along = 3)
		}
		Delta <- impulse_response_dynamic(A_array, B_array, shock, t0_stored - 1, H)
	####

	# handle static model: Theta_t = M + A (Theta_{t-1} - M) B' + eps
	# Delta_h = A^h S (B^h)' is the closed-form propagation under constant
	# A, B (sender / receiver factors from the static Tucker decomp)
	} else if (fit$model == "static") {
		if (is.null(fit$B)) cli::cli_abort(c(
			"Static model requires {.code B} matrix.",
			"i" = "Ensure the model was fit with {.code model = \"static\"}."
		))

		if (is.list(fit$B) && length(fit$B) >= 2L) {
			# fit$B[[1]] is the sender factor (analogue of dynamic A)
			# fit$B[[2]] is the receiver factor (analogue of dynamic B)
			A <- fit$B[[1]][, , draw_idx]
			B <- fit$B[[2]][, , draw_idx]
		} else if (is.list(fit$B)) {
			# legacy single-factor storage: treat as receiver and use A = I
			A <- diag(n_row)
			B <- fit$B[[1]][, , draw_idx]
		} else {
			A <- diag(n_row)
			B <- matrix(fit$B[draw_idx, , ], n_col, n_col)
		}
		Delta <- impulse_response_const(A, B, shock, H)
	####

	} else {
		cli::cli_abort(c(
			"IRF computation not yet implemented for model type {.val {fit$model}}.",
			"i" = "Supported model types are {.val static} and {.val dynamic}."
		))
	}

	# evaluate the network statistic at each horizon. stat_fun may return a
	# scalar (length 1) or a vector (per-element, e.g. per-actor); the caller
	# pre-probed the expected length as `stat_n_elem`.
	if (stat_n_elem == 1L) {
		irf_vals <- numeric(H + 1)
		for (h in 0:H) {
			baseline <- extract_baseline(fit, draw_idx, n_row, n_col, p)
			shocked_net <- baseline + Delta[, , h + 1]
			irf_vals[h + 1] <- stat_fun(shocked_net) - stat_fun(baseline)
		}
		irf_vals
	} else {
		irf_mat <- matrix(NA_real_, nrow = stat_n_elem, ncol = H + 1)
		for (h in 0:H) {
			baseline <- extract_baseline(fit, draw_idx, n_row, n_col, p)
			shocked_net <- baseline + Delta[, , h + 1]
			irf_mat[, h + 1] <- stat_fun(shocked_net) - stat_fun(baseline)
		}
		irf_mat
	}
	####
}
####

####
#' Extract Baseline Network for IRF
#'
#' @param fit dbn model fit object
#' @param draw_idx Posterior draw index
#' @param n_row Number of sender nodes
#' @param n_col Number of receiver nodes
#' @param p Number of relations
#' @return Baseline network matrix
#' @keywords internal
extract_baseline <- function(fit, draw_idx, n_row, n_col, p) {
	if (fit$model == "dynamic" && !is.null(fit$M)) {
		if (length(dim(fit$M)) == 4) {
			baseline <- matrix(0, n_row, n_col)
			for (rel in 1:p) baseline <- baseline + fit$M[, , rel, draw_idx]
			return(baseline / p)
		} else if (length(dim(fit$M)) == 3 && p == 1) {
			return(fit$M[, , draw_idx])
		}
	} else if (fit$model == "static" && !is.null(fit$M)) {
		if (is.matrix(fit$M)) return(fit$M)
		if (length(dim(fit$M)) == 3) return(fit$M[, , 1])
	}
	matrix(0, n_row, n_col)
}
####

####
#' Compute Network-Level Impulse Response Functions
#'
#' @description Computes IRFs for network-level statistics given a shock.
#'
#' **Important:** \code{ggplot2} exports its own \code{stat_density} function,
#' which masks \code{dbn::stat_density} when both packages are loaded.
#' If you use \code{ggplot2}, pass the statistic explicitly:
#' \code{stat_fun = dbn::stat_density}.
#'
#' @param fit A dbn model fit object
#' @param shock Shock matrix or shock type (see \code{\link{build_shock}})
#' @param H Number of horizons to compute. For a dynamic fit the response is
#'   propagated through the estimated operators within the observed series, so
#'   `t0 + H` must not exceed the number of time points `T`.
#' @param t0 Shock time for dynamic models (1-based). Must satisfy
#'   `t0 + H <= T`. The `horizon = 0` row of the result is the shock itself
#'   (the statistic of the unpropagated perturbation), not a response.
#' @param stat_fun Network statistic function (default: \code{\link{stat_density}}).
#'   Must accept a matrix and return a scalar. When \code{ggplot2} is loaded,
#'   use \code{dbn::stat_density} explicitly to avoid name masking.
#' @param n_draws Number of posterior draws to use
#' @param shock_pars List of parameters for \code{\link{build_shock}} if
#'   \code{shock} is a character string
#' @param seed Optional integer seed for reproducible per-draw sampling.
#' @param ... Additional arguments
#' @return Data frame with IRF results including posterior summaries
#' @seealso \code{\link{build_shock}}, \code{\link{plot.dbn_irf}},
#'   \code{\link{print.dbn_irf}}, \code{\link{stat_density}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' S <- build_shock(m = 6, type = "unit_edge", i = 1, j = 2)
#' irf <- compute_irf(fit, shock = S, H = 5, stat_fun = dbn::stat_density)
#' print(irf)
#' }
#' @importFrom stats complete.cases
#' @author Tosin Salau and Shahryar Minhas
#' @export
compute_irf <- function(fit, shock, H = 20, t0 = 1,
						stat_fun = stat_density, n_draws = NULL,
						shock_pars = list(), seed = NULL, ...) {
	# alias the common econometric synonyms for the horizon argument:
	# `H` is the canonical name (MTS / impulse-response convention) but
	# `horizon`, `horizons`, and `n.ahead` (vars / forecast convention)
	# all work too.
	extra <- list(...)
	for (alias in c("horizon", "horizons", "n.ahead")) {
		if (alias %in% names(extra)) {
			# only honor the alias when H was left at its default; if both
			# were specified explicitly, prefer the alias and warn.
			missing_H <- missing(H)
			if (!missing_H && !identical(H, 20)) {
				cli::cli_warn(c(
					"Both {.arg H} and {.arg {alias}} were specified; using {.arg {alias}} = {extra[[alias]]}.",
					"i" = "Pass only one."
				))
			}
			H <- extra[[alias]]
			extra[[alias]] <- NULL
			break
		}
	}
	if (length(extra) > 0L) {
		cli::cli_warn(c(
			"Ignoring {length(extra)} unrecognized argument{?s} to {.fun compute_irf}: {.arg {names(extra)}}.",
			"i" = "Allowed: {.arg fit}, {.arg shock}, {.arg H}, {.arg t0}, {.arg stat_fun}, {.arg n_draws}, {.arg shock_pars}, {.arg seed} (with {.arg horizon} / {.arg horizons} / {.arg n.ahead} accepted as aliases for {.arg H})."
		))
	}
	# honor `seed` for reproducible draw subsampling
	if (!is.null(seed)) {
		if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed)) {
			cli::cli_abort("{.arg seed} must be a single finite number.")
		}
		.restore <- .use_seed_locally(as.integer(seed))
		on.exit(.restore(), add = TRUE)
	}
	# capture the name of the statistic as the caller wrote it (e.g.
	# "stat_density"), so the result records the actual statistic used
	stat_label <- paste(deparse(substitute(stat_fun)), collapse = " ")
	if (length(stat_label) != 1L || nchar(stat_label) > 40L) {
		stat_label <- "custom_function"
	}
	if (!inherits(fit, "dbn")) cli::cli_abort(c(
		"{.arg fit} must be a {.cls dbn} object.",
		"i" = "Use {.code dbn()} to fit a model first."
	))
	# bootstrap-aware path: an ALS fit with `fit$bootstrap` carries R draws of
	# A, B, M in vectorised form. Reshape into per-draw arrays so the standard
	# per-draw IRF loop below sees R draws and the resulting CIs are genuine
	# bootstrap intervals (the static operator is tiled across time).
	su <- fit$meta$sampler_used
	is_als <- length(su) == 1L && !is.na(su) && su %in% c("als", "als_tv")
	if (is_als && !is.null(fit$bootstrap) && !isTRUE(fit$meta$bootstrap_expanded)) {
		fit <- .bootstrap_expand(fit)
	} else if (is_als || isFALSE(fit$meta$uncertainty_available)) {
		cli::cli_warn(c(
			"This is a point-estimate fit from {.fun dbn_als}.",
			"i" = "IRF credible-interval columns are degenerate (all equal to the mean).",
			"i" = "Refit with {.code bootstrap = N} or with {.fun dbn} for genuine intervals."
		))
	}
	# piecewise has block-constant operators that we tile into a per-time
	# cube and route through the dynamic IRF path. HMM has per-regime
	# operators (no regime path here) and lowrank uses a factored A; both
	# get a directive abort instead of the generic "not implemented".
	if (identical(fit$model, "piecewise")) {
		fit <- .piecewise_expand_to_dynamic(fit)
	}
	if (identical(fit$model, "hmm")) {
		cli::cli_abort(c(
			"IRF is not supported for {.val hmm} fits.",
			"i" = "HMM has regime-specific operators; an IRF would require conditioning on a regime path.",
			"i" = "If you have a posterior over regime assignments {.code fit$draws$S}, you can compute per-regime IRFs by extracting the per-regime operator and constructing a dummy piecewise fit; this is on the roadmap.",
			"i" = "For overall structural summaries on HMM, see {.fn compare_blocks}."
		))
	}
	if (identical(fit$model, "lowrank")) {
		cli::cli_abort(c(
			"IRF is not supported for {.val lowrank} fits.",
			"i" = "Low-rank parameterises `A = U diag(alpha) U^T`; the IRF reduces to a rank-r geometric series that requires a different computation than the full-rank path.",
			"i" = "On the roadmap; for stability inspection use {.fn dbn_operator} or the eigenvalue structure of {.code coef(fit)$A}."
		))
	}
	if (!fit$model %in% c("static", "dynamic")) {
		cli::cli_abort(c(
			"IRF currently only implemented for {.val static} and {.val dynamic} models.",
			"i" = "Got model type {.val {fit$model}}."
		))
	}
	# static fits have a time-invariant IRF (same operator at every
	# horizon -> pure geometric series of a fixed A). warn loudly so
	# users don't over-interpret it as a "dynamic" finding.
	if (identical(fit$model, "static") && !isTRUE(getOption("dbn.suppress_static_irf_warning", FALSE))) {
		cli::cli_warn(c(
			"Computing IRF on a {.val static} fit.",
			"i" = "Static fits apply the same operators at every horizon, so the IRF profile is purely a geometric series of fixed factor matrices (sender {.var A} and receiver {.var B}); it carries no time-varying information.",
			"i" = "If you want a time-varying IRF, refit with {.code model = \"dynamic\"}, {.val piecewise}, or {.val hmm}.",
			"i" = "Set {.code options(dbn.suppress_static_irf_warning = TRUE)} to silence."
		))
	}
	if (length(H) != 1L || !is.numeric(H) || !is.finite(H) || H < 1 || H != round(H)) {
		cli::cli_abort("{.arg H} (impulse-response horizon) must be a single positive integer.")
	}
	H <- as.integer(H)
	# stat_fun can be scalar- or vector-valued. vector-valued
	# (e.g. stat_in_degree, stat_out_degree, or a user-defined per-actor
	# stat) triggers the long-form output path: one row per
	# (horizon, element), which answers "who absorbs the shock from i?".
	stat_n_elem <- 1L
	stat_elem_names <- NULL
	if (!is.null(fit$dims$n_row)) {
		# probe stat_fun on a zero matrix to detect vector- vs scalar-valued
		# output. surface any error directly (except the ggplot2-masking
		# case, which has its own dedicated abort downstream), so users
		# see custom-stat-fn bugs here instead of via the per-draw aggregator.
		probe <- tryCatch(stat_fun(matrix(0, fit$dims$n_row, fit$dims$n_col)),
		                  error = function(e) {
			if (grepl("mapping.*aes|ggplot|stat_density", conditionMessage(e),
			          ignore.case = TRUE)) return(NULL)
			# the ggplot2-masking case is handled below; everything else
			# we surface here.
			cli::cli_warn(c(
				"Stat probe failed on a zero matrix; defaulting to scalar treatment.",
				"x" = "{conditionMessage(e)}",
				"i" = "Per-draw failures (if any) will be reported with the same error class at the end of the run."
			))
			NULL
		})
		if (!is.null(probe)) {
			stat_n_elem <- length(probe)
			if (!is.numeric(probe)) {
				cli::cli_abort(c(
					"{.arg stat_fun} must return a numeric (scalar or vector).",
					"x" = "Your {.arg stat_fun} returned a {.cls {class(probe)[1]}}."
				))
			}
			stat_elem_names <- names(probe)
		}
	}
	if (length(t0) != 1L || !is.numeric(t0) || !is.finite(t0) || t0 < 1 || t0 != round(t0)) {
		cli::cli_abort("{.arg t0} must be a single positive integer.")
	}
	t0 <- as.integer(t0)

	dims <- fit$dims
	n_row <- dims$n_row
	n_col <- dims$n_col

	# detect if ggplot2::stat_density is masking dbn::stat_density
	test_mat <- matrix(0, n_row, n_col)
	tryCatch(stat_fun(test_mat), error = function(e) {
		if (grepl("mapping.*aes|ggplot|stat_density", conditionMessage(e),
				  ignore.case = TRUE)) {
			cli::cli_abort(c(
				"{.fn ggplot2::stat_density} was passed instead of {.fn dbn::stat_density}.",
				"i" = "ggplot2 masks {.fn dbn::stat_density} when both packages are loaded.",
				"i" = "Use {.code stat_fun = dbn::stat_density} explicitly."
			))
		}
		cli::cli_abort(c(
			"{.arg stat_fun} must accept a matrix and return a scalar.",
			"i" = "Error when testing: {conditionMessage(e)}"
		))
	})

	# construct shock matrix from string shorthand if needed
	if (is.character(shock)) {
		shock_args <- c(list(m = n_row, n_col = n_col, type = shock), shock_pars)
		shock <- do.call(build_shock, shock_args)
	} else if (!is.matrix(shock) || nrow(shock) != n_row || ncol(shock) != n_col) {
		cli::cli_abort(c(
			"{.arg shock} must be an {.val {n_row}} x {.val {n_col}} matrix or a valid shock type.",
			"i" = "Valid shock types: {.val unit_edge}, {.val node_out}, {.val node_in}, {.val density}."
		))
	}
	####

	# determine how many posterior draws to use
	if (is.null(n_draws)) {
		if (fit$model == "dynamic") {
			if (!is.null(fit$A) && is.list(fit$A)) {
				n_draws <- length(fit$A)
			} else if (!is.null(fit$sigma2)) {
				n_draws <- length(fit$sigma2)
			} else {
				cli::cli_abort(c(
					"Cannot determine number of posterior draws for dynamic model.",
					"i" = "The fit object is missing both {.code A} and {.code sigma2} components."
				))
			}
		} else if (fit$model == "static") {
			n_draws <- static_draw_count(fit)
			if (is.na(n_draws)) {
				cli::cli_abort(c(
					"Cannot determine number of posterior draws for static model.",
					"i" = "The fit object is missing usable {.code B} and {.code params} components."
				))
			}
		}
	}

	available_draws <- if (fit$model == "dynamic" && !is.null(fit$A)) {
		length(fit$A)
	} else if (fit$model == "static") {
		static_draw_count(fit)
	} else {
		n_draws
	}
	if (is.na(available_draws)) available_draws <- n_draws
	n_draws <- suppressWarnings(min(n_draws, available_draws))
	if (!is.finite(n_draws) || n_draws < 1) {
		cli::cli_abort(c(
			"Could not determine a valid number of posterior draws for the IRF.",
			"i" = "Refit the model, or pass {.arg n_draws} explicitly."
		))
	}
	####

	# validate time bounds for dynamic model. t0 must be >= 1; t0 + H may
	# exceed T (out-of-sample IRF), in which case the per-draw path
	# extends the operator with constant padding.
	if (fit$model == "dynamic") {
		T_total <- dims$Tt
		if (is.null(T_total)) cli::cli_abort(c(
			"Cannot determine time dimension {.code T} from model.",
			"i" = "The fit object is missing {.code T}, {.code Tt}, and {.code TT} in {.code dims}."
		))
		if (t0 < 1) {
			cli::cli_abort(c(
				"{.arg t0} must be at least 1 for dynamic models.",
				"i" = "Got {.code t0 = {t0}}."
			))
		}
		if (t0 + H > T_total) {
			cli::cli_inform(c(
				"i" = "{.code t0 + H = {t0 + H}} extends past the last observed period ({.code T = {T_total}}); the operator is held constant from t = T forward."
			))
		}
	}
	####

	# loop over posterior draws and compute per-draw IRF
	# storage: 3D array [n_draws, stat_n_elem, H+1]. For scalar stat_fun the
	# middle dim is 1; for vector stat_fun (e.g. per-actor) it's > 1.
	irf_array <- array(NA_real_, dim = c(n_draws, stat_n_elem, H + 1))
	cli::cli_progress_bar("Computing IRFs", total = n_draws)

	# capture per-draw failures so the user sees the actual error class, not
	# just a draw count.
	draw_errors <- character(n_draws)
	for (s in 1:n_draws) {
		cli::cli_progress_update()
		tryCatch({
			vals <- compute_irf_single(fit, s, shock, H, t0, stat_fun,
			                            stat_n_elem = stat_n_elem)
			# vals is either length (H+1) (scalar stat) or matrix [stat_n_elem, H+1]
			if (stat_n_elem == 1L) {
				irf_array[s, 1, ] <- as.numeric(vals)
			} else {
				irf_array[s, , ] <- vals
			}
		}, error = function(e) {
			draw_errors[s] <<- conditionMessage(e)
			if (getOption("dbn.debug", FALSE)) {
				cat("Traceback for draw", s, ":\n")
				traceback()
			}
		})
	}
	cli::cli_progress_done()
	####

	# surface per-draw failures with the actual error classes.
	failed_mask <- nzchar(draw_errors)
	if (any(failed_mask)) {
		err_table <- sort(table(draw_errors[failed_mask]), decreasing = TRUE)
		n_shown <- min(3L, length(err_table))
		cli::cli_warn(c(
			"{sum(failed_mask)} of {n_draws} IRF draws failed.",
			"i" = "Top unique error{?s} (showing {n_shown}):",
			setNames(paste0("(", as.integer(err_table[1:n_shown]), "x) ",
				substr(names(err_table[1:n_shown]), 1, 80)),
				rep("x", n_shown))
		))
	}

	# drop draws where every (element, horizon) is NA (all-failed draws).
	draw_valid <- vapply(seq_len(n_draws), function(s) any(!is.na(irf_array[s, , ])),
	                    logical(1))
	if (sum(draw_valid) == 0L) {
		cli::cli_abort(c(
			"All IRF computations failed.",
			"i" = "Run {.fn debug_irf} for per-draw diagnostics (the function is now exported).",
			"i" = "Common causes: degenerate / singular operator at the target time, or an unusual {.arg shock} shape."
		))
	}
	irf_array <- irf_array[draw_valid, , , drop = FALSE]

	if (stat_n_elem == 1L) {
		# scalar stat: backward-compatible wide format (one row per horizon).
		mat <- matrix(irf_array, nrow = dim(irf_array)[1], ncol = H + 1)
		result <- data.frame(
			horizon = 0:H,
			mean   = colMeans(mat),
			median = apply(mat, 2, median),
			sd     = apply(mat, 2, sd),
			q025   = apply(mat, 2, quantile, 0.025),
			q975   = apply(mat, 2, quantile, 0.975),
			q10    = apply(mat, 2, quantile, 0.10),
			q90    = apply(mat, 2, quantile, 0.90)
		)
		attr(result, "irf_draws") <- mat
	} else {
		# vector stat: long format with per-element rows.
		elem_lbl <- stat_elem_names %||% paste0("elem_", seq_len(stat_n_elem))
		# irf_array now is [n_valid, stat_n_elem, H+1]
		# aggregate per (element, horizon)
		dims_v <- dim(irf_array)
		rows <- list()
		for (k in seq_len(stat_n_elem)) {
			# [n_valid x (H+1)]
			mat_k <- matrix(irf_array[, k, ], nrow = dims_v[1], ncol = dims_v[3])
			rows[[k]] <- data.frame(
				horizon = 0:H,
				element = elem_lbl[k],
				mean   = colMeans(mat_k),
				median = apply(mat_k, 2, median),
				sd     = apply(mat_k, 2, sd),
				q025   = apply(mat_k, 2, quantile, 0.025),
				q975   = apply(mat_k, 2, quantile, 0.975),
				q10    = apply(mat_k, 2, quantile, 0.10),
				q90    = apply(mat_k, 2, quantile, 0.90),
				stringsAsFactors = FALSE
			)
		}
		result <- do.call(rbind, rows)
		rownames(result) <- NULL
		attr(result, "irf_draws") <- irf_array
		attr(result, "stat_n_elem") <- stat_n_elem
	}

	attr(result, "shock") <- shock
	attr(result, "stat_fun") <- stat_label
	attr(result, "model") <- fit$model
	attr(result, "t0") <- t0
	class(result) <- c("dbn_irf", "data.frame")
	result
	####
}
####

####
#' Actor-level leverage (impulse-response propagation magnitude)
#'
#' @description
#' Leverage measures how widely a standardized shock involving an actor is
#' predicted to spread through the network over a finite horizon. For each
#' actor, `dbn_leverage()` applies a unit out-shock to that actor, propagates
#' it through the fitted lag operator with the package's impulse-response
#' machinery, and reports the average absolute network-wide effect over
#' horizons 1 to `H`.
#'
#' Leverage is distinct from coupling. Coupling (the row magnitude of the lag
#' operator \eqn{W_t}, see [dbn_coupling_rank_probs()]) asks whether the
#' previous network helps predict an actor's next-period involvement. Leverage
#' asks how far a shock involving the actor propagates once it occurs. An actor
#' can be high coupling and low leverage, or the reverse; the cases where the
#' two diverge are often the substantive payoff.
#'
#' @param fit A fitted dynamic `dbn` object. Static fits are rejected: with a
#'   constant operator a node shock propagates identically for every actor, so
#'   per-actor leverage is not identified.
#' @param H Propagation horizon (number of steps; default 3).
#' @param t0 The time index at which the shock is applied (default 1). Must
#'   satisfy `t0 + H <= T`.
#' @param magnitude Size of the unit out-shock (default 1).
#' @param actors Optional character labels for the actors.
#' @return A data frame of class `dbn_leverage`, one row per actor, with the
#'   posterior `mean`, `sd`, and 95% credible interval (`q025`, `q975`) of
#'   actor leverage. For a point-estimate ([dbn_als()]) fit the interval
#'   collapses to the point estimate.
#' @seealso [compute_irf()], [build_shock()], [dbn_coupling_rank_probs()]
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 8, time = 12, seed = 1)
#' fit <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'            nscan = 200, burn = 100, verbose = FALSE)
#' dbn_leverage(fit, H = 3)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_leverage <- function(fit, H = 3, t0 = 1, magnitude = 1, actors = NULL) {
	if (!inherits(fit, "dbn")) cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	# leverage requires a time-varying operator. Under the static model the
	# operator is constant and (in this package's IRF path) A = I, so a unit
	# node-out shock propagates identically for every actor and per-actor
	# leverage is not identified -- reject static fits rather than return a
	# vector of identical, meaningless values.
	if (fit$model != "dynamic") {
		cli::cli_abort(c(
			"{.fun dbn_leverage} requires a dynamic fit.",
			"x" = "Got model type {.val {fit$model}}.",
			"i" = "With a constant operator, a node shock propagates identically for every actor, so per-actor leverage is not identified. Use {.code model = \"dynamic\"}."
		))
	}
	if (!is.numeric(H) || length(H) != 1L || H < 1) {
		cli::cli_abort("{.arg H} must be a single integer of at least 1.")
	}
	H <- as.integer(H)
	dims <- fit$dims
	n_row <- as.integer(dims$n_row)
	n_col <- as.integer(dims$n_col)

	if (is.null(fit$A) || !is.list(fit$A)) {
		cli::cli_abort("Dynamic fit is missing time-varying {.code A} draws.")
	}
	n_draws <- length(fit$A)
	Tt <- dim(fit$A[[1]])[3]
	if (t0 < 1 || t0 + H > Tt) {
		cli::cli_abort(c(
			"{.code t0 + H} must not exceed the number of time points ({Tt}).",
			"i" = "Got {.code t0 = {t0}}, {.code H = {H}}."
		))
	}

	# per-draw, per-actor leverage: mean |Delta| over horizons 1..H of a unit
	# out-shock applied to the actor (slice 1 of Delta is the shock itself)
	lev <- matrix(NA_real_, n_draws, n_row)
	for (m in seq_len(n_draws)) {
		for (i in seq_len(n_row)) {
			shock <- build_shock(n_row, type = "node_out", i = i,
			                     magnitude = magnitude, n_col = n_col)
			Delta <- impulse_response_dynamic(fit$A[[m]], fit$B[[m]],
			                                  shock, t0 - 1, H)
			lev[m, i] <- mean(abs(Delta[, , seq.int(2L, H + 1L), drop = FALSE]))
		}
	}

	if (is.null(actors)) {
		actors <- paste0("actor_", seq_len(n_row))
	} else if (length(actors) != n_row) {
		cli::cli_abort(c(
			"{.arg actors} must have one label per actor.",
			"x" = "Got {length(actors)} label{?s} for {n_row} actor{?s}."
		))
	}
	# a single-draw (point-estimate) fit has no posterior spread: report
	# sd = 0 and a degenerate interval rather than NA
	sd_col <- if (n_draws < 2L) rep(0, n_row) else apply(lev, 2, stats::sd, na.rm = TRUE)
	result <- data.frame(
		actor = actors,
		mean  = colMeans(lev, na.rm = TRUE),
		sd    = sd_col,
		q025  = apply(lev, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
		q975  = apply(lev, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
		row.names = NULL, stringsAsFactors = FALSE
	)
	class(result) <- c("dbn_leverage", "data.frame")
	result
}
####

####
#' Print method for actor-level leverage
#'
#' @param x A `dbn_leverage` object from [dbn_leverage()].
#' @param digits Number of digits (default 3).
#' @param top Number of top-leverage actors to show (default 8).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @method print dbn_leverage
print.dbn_leverage <- function(x, digits = 3, top = 8, ...) {
	cli::cli_h2("DBN actor leverage")
	cli::cli_text("Average absolute network-wide response to a unit out-shock, by actor.")
	df <- as.data.frame(x)
	ord <- order(df$mean, decreasing = TRUE)
	show <- ord[seq_len(min(top, nrow(df)))]
	for (i in show) {
		cli::cli_text("  {df$actor[i]}: {round(df$mean[i], digits)} [{round(df$q025[i], digits)}, {round(df$q975[i], digits)}]")
	}
	if (nrow(df) > top) cli::cli_text("  {cli::symbol$ellipsis} {nrow(df) - top} more actor{?s}")
	invisible(x)
}
####

####
#' Plot IRF
#'
#' @description IRF plot with credible intervals
#' @param x A dbn_irf object from compute_irf
#' @param ci_level Credible interval level (0.90 or 0.95)
#' @param title Plot title
#' @param ... Ignored
#' @return A ggplot2 object
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{print.dbn_irf}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' S <- build_shock(m = 6, type = "unit_edge", i = 1, j = 2)
#' irf <- compute_irf(fit, shock = S, H = 5)
#' plot(irf)
#' }
#' @importFrom rlang .data
#' @author Tosin Salau and Shahryar Minhas
#' @export
plot.dbn_irf <- function(x, ci_level = 0.95, title = NULL, ...) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort(c(
		"Package {.pkg ggplot2} is required for this function.",
		"i" = "Install with {.code install.packages(\"ggplot2\")}."
	))

	if (ci_level == 0.95) {
		q_low <- "q025"; q_high <- "q975"
	} else if (ci_level == 0.90) {
		q_low <- "q10"; q_high <- "q90"
	} else {
		cli::cli_abort(c(
			"{.arg ci_level} must be {.val 0.90} or {.val 0.95}.",
			"i" = "Got {.val {ci_level}}."
		))
	}

	if (is.null(title)) {
		title <- paste0("Impulse Response Function: ", attr(x, "stat_fun"))
	}

	ggplot2::ggplot(x, ggplot2::aes(x = .data$horizon)) +
		ggplot2::geom_ribbon(
			ggplot2::aes(ymin = .data[[q_low]], ymax = .data[[q_high]]),
			alpha = 0.3, fill = "gray") +
		ggplot2::geom_line(ggplot2::aes(y = mean), linewidth = 1) +
		ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
		ggplot2::labs(title = title, x = "Horizon", y = "Response") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
}
####

####
#' Print IRF
#'
#' @param x A dbn_irf object
#' @param digits Number of digits
#' @param ... Ignored
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{plot.dbn_irf}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' S <- build_shock(m = 6, type = "unit_edge", i = 1, j = 2)
#' irf <- compute_irf(fit, shock = S, H = 5)
#' print(irf)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_irf <- function(x, digits = 3, ...) {
	cat("Network Impulse Response Function\n")
	cat("Model:", attr(x, "model"), "\n")
	cat("Statistic:", attr(x, "stat_fun"), "\n")
	if (attr(x, "model") == "dynamic") cat("Shock time:", attr(x, "t0"), "\n")
	cat("\nSummary:\n")
	# round only the numeric columns; vector-valued stat_fun output has a
	# character `element` column that round.data.frame cannot handle.
	df <- as.data.frame(x)
	num_cols <- vapply(df, is.numeric, logical(1))
	df[num_cols] <- lapply(df[num_cols], round, digits)
	print(df)
	invisible(x)
}
####

####
#' Debug IRF
#'
#' @description Diagnostic output for troubleshooting IRF computation.
#'   Internal helper retained for regression testing.
#' @param fit A dbn model fit object
#' @param draw_idx Posterior draw index to debug
#' @param shock_i Source node (default: 1)
#' @param shock_j Target node (default: 2)
#' @return Invisibly returns `NULL`; called for its printed diagnostic side
#'   effects.
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{plot.dbn_irf}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 15, time = 50, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' debug_irf(fit, draw_idx = 1, shock_i = 1, shock_j = 2)
#' }
#' @export
#' @keywords internal
debug_irf <- function(fit, draw_idx = 1, shock_i = 1, shock_j = 2) {
	if (!inherits(fit, "dbn")) {
		cli::cli_abort("{.arg fit} must be a {.cls dbn} object.")
	}
	# preserve the user-facing model label before any internal expansion
	model_label <- fit$model %||% "?"
	# static fits have no time-varying operator; debug_irf's C++ probe step
	# would receive a NULL A and crash the R session at the Rcpp layer.
	# refuse upfront with a friendly directive.
	if (identical(fit$model, "static")) {
		cli::cli_abort(c(
			"{.fn debug_irf} is not applicable to {.val static} fits.",
			"x" = "Static fits have no time-varying operator {.var A_t}; the IRF is a degenerate geometric series of a fixed {.var A}.",
			"i" = "For a static-fit IRF profile, call {.fn compute_irf} directly; it will warn and produce the constant-operator output.",
			"i" = "For time-varying diagnostics, refit with {.code model = \"dynamic\"} or {.code model = \"piecewise\"}."
		))
	}
	# auto-expand piecewise fits so the dynamic-shaped A / B cubes are
	# available; without this, fit$A is a list of per-block n_row x n_row
	# matrices and the C++ impulse_response_dynamic rejects the 2D input.
	if (identical(fit$model, "piecewise") &&
	    !isTRUE(fit$meta$piecewise_expanded)) {
		fit <- .piecewise_expand_to_dynamic(fit)
	}
	if (is.null(fit$A) || !is.list(fit$A) || length(fit$A) == 0L) {
		cli::cli_abort(c(
			"{.fn debug_irf} could not find time-varying operator draws on {.arg fit}.",
			"x" = "{.code fit$A} is {.val {if (is.null(fit$A)) NULL else class(fit$A)[1]}}.",
			"i" = "Make sure the fit was produced by {.code model = \"dynamic\"} or {.code model = \"piecewise\"}, and that {.code keep = ...} did not strip {.var A}."
		))
	}
	# validate draw_idx upfront so an out-of-bounds index errors with a
	# friendly directive rather than crashing mid-diagnostic.
	n_draws_fit <- length(fit$A)
	if (length(draw_idx) != 1L || !is.numeric(draw_idx) || !is.finite(draw_idx) ||
	    draw_idx < 1 || draw_idx > n_draws_fit || draw_idx != round(draw_idx)) {
		cli::cli_abort(c(
			"{.arg draw_idx} must be a single integer in 1..{n_draws_fit}.",
			"x" = "Got {.val {draw_idx}}; the fit has {.val {n_draws_fit}} stored draws."
		))
	}
	cat("=== IRF Debug Information ===\n")

	cat("\n1. Model Structure:\n")
	cat("   Model type:", model_label,
		if (!identical(model_label, fit$model %||% "?"))
			sprintf(" (expanded to %s for IRF)", fit$model %||% "?") else "",
		"\n", sep = "")
	cat("   n_row =", fit$dims$n_row, "\n")
	cat("   n_col =", fit$dims$n_col, "\n")
	cat("   p =", fit$dims$p, "\n")
	cat("   T =", fit$dims$Tt, "\n")

	cat("\n2. A and B matrices:\n")
	cat("   A is list:", is.list(fit$A), "\n")
	cat("   B is list:", is.list(fit$B), "\n")
	if (is.list(fit$A)) {
		cat("   Length of A:", length(fit$A), "\n")
		cat("   Dim of A[[1]]:", dim(fit$A[[1]]), "\n")
	}
	if (is.list(fit$B)) {
		cat("   Length of B:", length(fit$B), "\n")
		cat("   Dim of B[[1]]:", dim(fit$B[[1]]), "\n")
	}

	cat("\n3. M matrix:\n")
	cat("   M dimensions:", dim(fit$M), "\n")

	cat("\n4. Attempting to extract matrices for draw", draw_idx, ":\n")
	tryCatch({
		n_row <- fit$dims$n_row
		n_col <- fit$dims$n_col

		A_array <- fit$A[[draw_idx]]
		B_array <- fit$B[[draw_idx]]
		cat("   A_array dim:", dim(A_array), "\n")
		cat("   B_array dim:", dim(B_array), "\n")

		if (length(dim(fit$M)) == 4) {
			M_draw <- fit$M[, , 1, draw_idx]
			cat("   M for draw", draw_idx, "extracted, dim:", dim(M_draw), "\n")
		}

		S <- matrix(0, n_row, n_col)
		S[shock_i, shock_j] <- 1

		Tt_A <- dim(A_array)[3]
		H_dbg <- min(5, Tt_A - 1)
		t0_dbg <- max(0, Tt_A - H_dbg - 1)

		cat("\n5. Testing C++ impulse_response_dynamic:\n")
		Delta <- impulse_response_dynamic(A_array, B_array, S, t0_dbg, H_dbg)
		cat("   Delta dim:", dim(Delta), "\n")
		cat("   Delta[,,2] non-zero elements:", sum(Delta[, , 2] != 0), "\n")

		cat("\n6. Testing stat_density:\n")
		test_mat <- matrix(runif(n_row * n_col), n_row, n_col)
		cat("   stat_density on random matrix:", stat_density(test_mat), "\n")

		baseline <- M_draw
		shocked <- baseline + Delta[, , 1]
		cat("   Baseline density:", stat_density(baseline), "\n")
		cat("   Shocked density:", stat_density(shocked), "\n")
		cat("   Difference:", stat_density(shocked) - stat_density(baseline), "\n")
	}, error = function(e) {
		cat("   ERROR:", e$message, "\n")
	})

	cat("\n=== End Debug ===\n")
}
####
