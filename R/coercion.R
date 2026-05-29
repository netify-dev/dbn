####
#' Coerce common network formats to the dbn 4D array
#'
#' @description
#' `dbn()` and friends consume a 4D `[n_row, n_col, p, T]` array (with
#' diagonal `NA` on unipartite slices). Most R users start with an
#' `igraph` panel, a `statnet` `networkDynamic`, or a tidy edgelist.
#' These helpers do the coercion so callers don't have to write the
#' `array(NA_real_, dim = ...)` skeleton themselves.
#'
#' @details
#' All coercions assume a single relation (`p = 1`) by default. For
#' multi-relation data, call the helper once per relation and stack with
#' `abind::abind` along the `p` axis. Diagonal entries on unipartite
#' (square) outputs are set to `NA` to match the dbn convention.
#'
#' @name dbn_coercion
#' @param ... Method-specific arguments forwarded to the dispatched
#'   coercion method.
#' @return A 4D numeric array of shape `[n_row, n_col, 1, T]` with
#'   `dimnames` populated from the source object when available, suitable
#'   for direct use as the first argument to `dbn()`.
NULL

#' Coerce an object to a dbn 4D array (S3 generic)
#'
#' @param x Object to coerce. Dispatches on class; methods exist for
#'   `network` (statnet), `networkDynamic`, `igraph_panel` (list of
#'   igraphs), and a default list-of-networks. Edgelist users should call
#'   \code{\link{as_dbn_array_edgelist}} directly.
#' @param ... Method-specific arguments.
#' @return A 4D numeric array.
#' @author Tosin Salau and Shahryar Minhas
#' @export
as_dbn_array <- function(x, ...) UseMethod("as_dbn_array")

#' @rdname dbn_coercion
#' @param x A single `network` object (statnet). Returned as
#'   `[n, n, 1, 1]`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
as_dbn_array.network <- function(x, ...) {
	if (!requireNamespace("network", quietly = TRUE))
		cli::cli_abort(c("Package {.pkg network} is required.",
		                 "i" = "Install with {.code install.packages(\"network\")}"))
	A <- as.matrix(x)
	if (nrow(A) == ncol(A)) diag(A) <- NA_real_
	v <- network::network.vertex.names(x)
	array(A, dim = c(nrow(A), ncol(A), 1L, 1L),
		dimnames = list(v, v, "tie", "1"))
}

#' @rdname dbn_coercion
#' @param x A list of `network` objects sharing the same vertex set.
#' @author Tosin Salau and Shahryar Minhas
#' @export
as_dbn_array.list <- function(x, ...) {
	if (length(x) == 0L)
		cli::cli_abort("{.arg x} is empty.")
	if (all(vapply(x, inherits, logical(1L), what = "igraph"))) {
		# delegate to the igraph_panel method
		class(x) <- c("igraph_panel", class(x))
		return(as_dbn_array.igraph_panel(x, ...))
	}
	if (!all(vapply(x, inherits, logical(1L), what = "network")))
		cli::cli_abort(c(
			"{.arg x} must be a list of {.cls network} or {.cls igraph} objects.",
			"i" = "For tidy edgelists use {.fn as_dbn_array_edgelist}."
		))
	if (!requireNamespace("network", quietly = TRUE))
		cli::cli_abort("Package {.pkg network} is required.")
	v <- network::network.vertex.names(x[[1L]])
	n <- length(v)
	Tt <- length(x)
	Y <- array(NA_real_, dim = c(n, n, 1L, Tt),
		dimnames = list(v, v, "tie", as.character(seq_len(Tt))))
	for (t in seq_along(x)) {
		A <- as.matrix(x[[t]])
		Y[, , 1, t] <- A
	}
	for (t in seq_len(Tt)) diag(Y[, , 1, t]) <- NA_real_
	Y
}

#' @rdname dbn_coercion
#' @param x An `igraph` object panel, supplied as a `list` of `igraph`
#'   graphs (one per time point) sharing the same vertex set.
#' @param attr Edge attribute to extract (numeric weight). If `NULL` and
#'   the graph is unweighted, the binary adjacency is used.
#' @param time_labels Optional character vector of time labels (length
#'   `T`). Defaults to `seq_along(x)`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
as_dbn_array.igraph_panel <- function(x, attr = NULL, time_labels = NULL, ...) {
	if (!requireNamespace("igraph", quietly = TRUE))
		cli::cli_abort(c("Package {.pkg igraph} is required.", "i" = "Install with {.code install.packages(\"igraph\")}"))
	if (!is.list(x) || length(x) == 0L)
		cli::cli_abort("{.arg x} must be a non-empty list of igraph graphs.")
	if (!all(vapply(x, inherits, logical(1L), what = "igraph")))
		cli::cli_abort("Every element of {.arg x} must be an {.cls igraph} object.")
	# vertex sets must match across snapshots
	vsets <- lapply(x, function(g) igraph::V(g)$name %||% as.character(seq_len(igraph::vcount(g))))
	if (length(unique(lapply(vsets, sort))) != 1L)
		cli::cli_abort(c(
			"Vertex sets differ across snapshots in {.arg x}.",
			"i" = "Add isolates / drop missing actors so every snapshot shares the same labelled vertex set."
		))
	v <- vsets[[1L]]
	n <- length(v)
	Tt <- length(x)
	Y <- array(NA_real_, dim = c(n, n, 1L, Tt),
	           dimnames = list(v, v, "tie",
	                           time_labels %||% as.character(seq_len(Tt))))
	for (t in seq_len(Tt)) {
		g <- x[[t]]
		if (is.null(attr)) {
			adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
		} else {
			adj <- as.matrix(igraph::as_adjacency_matrix(g, attr = attr, sparse = FALSE))
		}
		# align rows and cols to canonical vertex order
		gv <- igraph::V(g)$name %||% as.character(seq_len(igraph::vcount(g)))
		adj <- adj[match(v, gv), match(v, gv), drop = FALSE]
		Y[, , 1, t] <- adj
	}
	for (t in seq_len(Tt)) diag(Y[, , 1, t]) <- NA_real_
	Y
}

#' @rdname dbn_coercion
#' @param x A `networkDynamic` object (from the `networkDynamic` package).
#' @param onsets Numeric vector of time points to snapshot at. Default
#'   uses `networkDynamic::get.change.times(x)`.
#' @author Tosin Salau and Shahryar Minhas
#' @export
as_dbn_array.networkDynamic <- function(x, onsets = NULL, ...) {
	if (!requireNamespace("networkDynamic", quietly = TRUE))
		cli::cli_abort(c("Package {.pkg networkDynamic} is required.",
		                 "i" = "Install with {.code install.packages(\"networkDynamic\")}"))
	if (!inherits(x, "networkDynamic"))
		cli::cli_abort("{.arg x} must be a {.cls networkDynamic} object.")
	if (is.null(onsets)) onsets <- networkDynamic::get.change.times(x)
	if (length(onsets) == 0L)
		cli::cli_abort("Could not derive time onsets from {.arg x}; pass {.arg onsets} explicitly.")
	# trim trailing onsets that correspond to the open-right endpoint of a
	# spell (a network.list of N discrete snapshots gives change.times of
	# 1..N+1; the last one has no active spell and as.matrix returns 0 rows)
	keep <- vapply(onsets, function(tt) {
		net_t <- tryCatch(networkDynamic::network.collapse(x, at = tt),
			error = function(e) NULL)
		!is.null(net_t) && network::network.size(net_t) > 0L
	}, logical(1L))
	onsets <- onsets[keep]
	if (length(onsets) == 0L)
		cli::cli_abort("No usable snapshots found; pass {.arg onsets} explicitly.")
	n <- network::network.size(x)
	v <- network::network.vertex.names(x)
	Tt <- length(onsets)
	Y <- array(NA_real_, dim = c(n, n, 1L, Tt),
	           dimnames = list(v, v, "tie", as.character(onsets)))
	for (t in seq_len(Tt)) {
		snap <- networkDynamic::network.collapse(x, at = onsets[t])
		Y[, , 1, t] <- as.matrix(snap)
	}
	for (t in seq_len(Tt)) diag(Y[, , 1, t]) <- NA_real_
	Y
}

#' @rdname dbn_coercion
#' @param x A tidy edgelist data frame with columns
#'   `from`, `to`, `time`, and optionally `weight`.
#' @param actors Optional character vector of canonical actor labels in
#'   the order they should appear in the array. Defaults to
#'   `sort(unique(c(x$from, x$to)))`.
#' @param time_levels Optional vector of time-axis levels. Defaults to
#'   `sort(unique(x$time))`.
#' @param weight Column name (string) holding the weight. If `NULL`,
#'   uses 1 for every observed edge (binary).
#' @param symmetric Logical; if `TRUE`, mirror across the diagonal.
#' @author Tosin Salau and Shahryar Minhas
#' @export
as_dbn_array_edgelist <- function(x, actors = NULL, time_levels = NULL,
                                  weight = NULL, symmetric = FALSE) {
	stopifnot(is.data.frame(x))
	required <- c("from", "to", "time")
	miss <- setdiff(required, names(x))
	if (length(miss) > 0L)
		cli::cli_abort(c("Required columns missing from {.arg x}.", "x" = "Missing: {.val {miss}}"))
	if (is.null(actors)) actors <- sort(unique(c(as.character(x$from), as.character(x$to))))
	if (is.null(time_levels)) time_levels <- sort(unique(x$time))
	n <- length(actors); Tt <- length(time_levels)
	Y <- array(NA_real_, dim = c(n, n, 1L, Tt),
	           dimnames = list(actors, actors, "tie", as.character(time_levels)))
	# missing edge = 0, diagonal = NA
	for (t in seq_along(time_levels)) Y[, , 1, t] <- 0
	w <- if (is.null(weight)) rep(1, nrow(x)) else x[[weight]]
	row_idx <- match(as.character(x$from), actors)
	col_idx <- match(as.character(x$to), actors)
	t_idx   <- match(x$time, time_levels)
	bad <- is.na(row_idx) | is.na(col_idx) | is.na(t_idx)
	if (any(bad))
		cli::cli_warn("Dropping {sum(bad)} row{?s} of {.arg x} whose from/to/time do not match the canonical levels.")
	for (k in which(!bad)) Y[row_idx[k], col_idx[k], 1, t_idx[k]] <- w[k]
	if (symmetric) {
		for (t in seq_len(Tt)) {
			S <- Y[, , 1, t]
			Y[, , 1, t] <- pmax(S, t(S), na.rm = TRUE)
		}
	}
	for (t in seq_len(Tt)) diag(Y[, , 1, t]) <- NA_real_
	Y
}

# peacesciencer_to_dbn_array was previously defined here. it has been
# removed: the package is data-source-agnostic, and users with dyad-year
# panels should construct their dbn array via netify or call
# as_dbn_array_edgelist() directly. that keeps the dbn surface focused
# on modelling rather than on any particular upstream data toolchain.
