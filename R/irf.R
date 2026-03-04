####
#' Build Shock Matrix for IRF Analysis
#'
#' @description Creates shock matrices for different types of network interventions
#' @param m Number of sender nodes
#' @param type Type of shock: "unit_edge", "node_out", "node_in", or "density"
#' @param i Source node index
#' @param j Target node index (for unit_edge)
#' @param magnitude Shock magnitude
#' @param n_col Number of receiver nodes (default: m)
#' @return n_row x n_col shock matrix
#' @seealso \code{\link{compute_irf}}, \code{\link{plot.dbn_irf}}
#' @examples
#' # Unit edge shock: activate edge from node 1 to node 2
#' S <- build_shock(m = 5, type = "unit_edge", i = 1, j = 2)
#' str(S)
#'
#' # Node-level shock: activate all outgoing edges from node 1
#' S_out <- build_shock(m = 5, type = "node_out", i = 1)
#'
#' # Density shock: small shock spread across the network
#' S_dens <- build_shock(m = 5, type = "density", magnitude = 0.1)
#'
#' # Bipartite shock: 4 senders, 6 receivers
#' S_bip <- build_shock(m = 4, n_col = 6, type = "unit_edge", i = 1, j = 3)
#' @export
build_shock <- function(m, type = c("unit_edge", "node_out", "node_in", "density"),
						i = 1, j = 2, magnitude = 1, n_col = m) {
	type <- match.arg(type)
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
#' @export
stat_density <- function(X) {
	mean(X[row(X) != col(X)])
}

#' Network Statistic: In-Degree
#'
#' @description Column sums of network matrix
#' @param X Network matrix
#' @return Vector of in-degrees
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_density}}, \code{\link{stat_out_degree}}
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' stat_in_degree(X)
#' @export
stat_in_degree <- function(X) {
	colSums(X)
}

#' Network Statistic: Out-Degree
#'
#' @description Row sums of network matrix
#' @param X Network matrix
#' @return Vector of out-degrees
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_density}}, \code{\link{stat_in_degree}}
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' stat_out_degree(X)
#' @export
stat_out_degree <- function(X) {
	rowSums(X)
}

#' Network Statistic: Reciprocity
#'
#' @description Correlation between X\[i,j\] and X\[j,i\]
#' @param X Network matrix
#' @return Scalar reciprocity value
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_density}}, \code{\link{stat_transitivity}}
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' stat_reciprocity(X)
#' @export
stat_reciprocity <- function(X) {
	upper_idx <- upper.tri(X)
	upper_vals <- X[upper_idx]
	lower_vals <- t(X)[upper_idx]
	if (length(unique(c(upper_vals, lower_vals))) == 1) return(0)
	cor(upper_vals, lower_vals)
}

#' Network Statistic: Transitivity
#'
#' @description Clustering coefficient
#' @param X Network matrix
#' @return Scalar transitivity value
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{stat_density}}, \code{\link{stat_reciprocity}}
#' @examples
#' X <- matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), 3, 3)
#' stat_transitivity(X)
#' @export
stat_transitivity <- function(X) {
	X_binary <- (X != 0) * 1
	diag(X_binary) <- 0
	triangles <- sum(diag(X_binary %*% X_binary %*% X_binary)) / 6
	degrees <- rowSums(X_binary)
	triples <- sum(degrees * (degrees - 1)) / 2
	if (triples == 0) return(0)
	triangles / triples
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
compute_irf_single <- function(fit, draw_idx, shock, H, t0 = 1, stat_fun = stat_density) {
	dims <- fit$dims
	n_row <- as.integer(dims$n_row)
	n_col <- as.integer(dims$n_col)
	p <- as.integer(dims$p)

	# dynamic model
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

		if (t0_stored + H > T_stored) {
			cli::cli_abort(c(
				"{.arg t0} + {.arg H} exceeds stored time points.",
				"i" = "With {.code time_thin = {time_thin}}, maximum {.code t0 + H} is {.val {T_stored * time_thin}}."
			))
		}

		A_array <- fit$A[[draw_idx]]
		B_array <- fit$B[[draw_idx]]
		Delta <- impulse_response_dynamic(A_array, B_array, shock, t0_stored - 1, H)
	####

	# static model
	} else if (fit$model == "static") {
		if (is.null(fit$B)) cli::cli_abort(c(
			"Static model requires {.code B} matrix.",
			"i" = "Ensure the model was fit with {.code model = \"static\"}."
		))

		if (is.list(fit$B)) {
			B <- fit$B[[1]][, , draw_idx]
		} else {
			B <- matrix(fit$B[draw_idx, , ], n_col, n_col)
		}
		A <- diag(n_row)
		Delta <- impulse_response_const(A, B, shock, H)
	####

	} else {
		cli::cli_abort(c(
			"IRF computation not yet implemented for model type {.val {fit$model}}.",
			"i" = "Supported model types are {.val static} and {.val dynamic}."
		))
	}

	# compute network statistic at each horizon
	irf_vals <- numeric(H + 1)
	for (h in 0:H) {
		baseline <- extract_baseline(fit, draw_idx, n_row, n_col, p)
		shocked_net <- baseline + Delta[, , h + 1]
		irf_vals[h + 1] <- stat_fun(shocked_net) - stat_fun(baseline)
	}
	irf_vals
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
#' @param H Number of horizons to compute
#' @param t0 Shock time for dynamic models
#' @param stat_fun Network statistic function (default: \code{\link{stat_density}}).
#'   Must accept a matrix and return a scalar. When \code{ggplot2} is loaded,
#'   use \code{dbn::stat_density} explicitly to avoid name masking.
#' @param n_draws Number of posterior draws to use
#' @param shock_pars List of parameters for \code{\link{build_shock}} if
#'   \code{shock} is a character string
#' @param ... Additional arguments
#' @return Data frame with IRF results including posterior summaries
#' @seealso \code{\link{build_shock}}, \code{\link{plot.dbn_irf}},
#'   \code{\link{print.dbn_irf}}, \code{\link{debug_irf}},
#'   \code{\link{stat_density}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' S <- build_shock(m = 6, type = "unit_edge", i = 1, j = 2)
#' irf <- compute_irf(fit, shock = S, H = 5, stat_fun = dbn::stat_density)
#' print(irf)
#' }
#' @importFrom stats complete.cases
#' @export
compute_irf <- function(fit, shock, H = 20, t0 = 1,
						stat_fun = stat_density, n_draws = NULL,
						shock_pars = list(), ...) {
	if (!inherits(fit, "dbn")) cli::cli_abort(c(
		"{.arg fit} must be a {.cls dbn} object.",
		"i" = "Use {.code dbn()} to fit a model first."
	))
	if (!fit$model %in% c("static", "dynamic")) {
		cli::cli_abort(c(
			"IRF currently only implemented for {.val static} and {.val dynamic} models.",
			"i" = "Got model type {.val {fit$model}}."
		))
	}

	dims <- fit$dims
	n_row <- dims$n_row
	n_col <- dims$n_col

	# guard against ggplot2::stat_density masking dbn::stat_density
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

	# build shock matrix if needed
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

	# determine number of draws
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
			if (!is.null(fit$B)) {
				n_draws <- dim(fit$B)[1]
			} else if (!is.null(fit$params)) {
				n_draws <- nrow(fit$params)
			} else {
				cli::cli_abort(c(
					"Cannot determine number of posterior draws for static model.",
					"i" = "The fit object is missing both {.code B} and {.code params} components."
				))
			}
		}
	}

	available_draws <- if (fit$model == "dynamic" && !is.null(fit$A)) {
		length(fit$A)
	} else if (fit$model == "static" && !is.null(fit$B)) {
		dim(fit$B)[1]
	} else {
		n_draws
	}
	n_draws <- min(n_draws, available_draws)
	####

	# check time bounds for dynamic
	if (fit$model == "dynamic") {
		T_total <- dims$Tt
		if (is.null(T_total)) cli::cli_abort(c(
			"Cannot determine time dimension {.code T} from model.",
			"i" = "The fit object is missing {.code T}, {.code Tt}, and {.code TT} in {.code dims}."
		))
		if (t0 < 1 || t0 + H > T_total) {
			cli::cli_abort(c(
				"{.code t0 + H} must not exceed {.code T} for dynamic models.",
				"i" = "Got {.code t0 = {t0}}, {.code H = {H}}, but {.code T = {T_total}}."
			))
		}
	}
	####

	# compute IRF for each posterior draw
	irf_array <- matrix(NA, n_draws, H + 1)
	cli::cli_progress_bar("Computing IRFs", total = n_draws)

	for (s in 1:n_draws) {
		cli::cli_progress_update()
		tryCatch({
			irf_array[s, ] <- compute_irf_single(fit, s, shock, H, t0, stat_fun)
		}, error = function(e) {
			warning("Error in draw ", s, ": ", e$message)
			if (getOption("dbn.debug", FALSE)) {
				cat("Traceback for draw", s, ":\n")
				traceback()
			}
		})
	}
	cli::cli_progress_done()
	####

	# posterior summaries
	valid_draws <- complete.cases(irf_array)
	if (sum(valid_draws) < n_draws) {
		warning(n_draws - sum(valid_draws), " draws failed and were removed")
		irf_array <- irf_array[valid_draws, , drop = FALSE]
	}
	if (nrow(irf_array) == 0) {
		cli::cli_abort(c(
			"All IRF computations failed.",
			"i" = "Run {.code debug_irf()} for diagnostics."
		))
	}

	result <- data.frame(
		horizon = 0:H,
		mean = colMeans(irf_array),
		median = apply(irf_array, 2, median),
		sd = apply(irf_array, 2, sd),
		q025 = apply(irf_array, 2, quantile, 0.025),
		q975 = apply(irf_array, 2, quantile, 0.975),
		q10 = apply(irf_array, 2, quantile, 0.10),
		q90 = apply(irf_array, 2, quantile, 0.90)
	)

	attr(result, "irf_draws") <- irf_array
	attr(result, "shock") <- shock
	attr(result, "stat_fun") <- if (is.function(stat_fun)) "custom_function" else deparse(substitute(stat_fun))
	attr(result, "model") <- fit$model
	attr(result, "t0") <- t0
	class(result) <- c("dbn_irf", "data.frame")
	result
	####
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
#' irf <- compute_irf(fit, S = S, H = 5)
#' plot(irf)
#' }
#' @importFrom rlang .data
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
		ggplot2::theme_minimal()
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
#' irf <- compute_irf(fit, S = S, H = 5)
#' print(irf)
#' }
#' @export
print.dbn_irf <- function(x, digits = 3, ...) {
	cat("Network Impulse Response Function\n")
	cat("Model:", attr(x, "model"), "\n")
	cat("Statistic:", attr(x, "stat_fun"), "\n")
	if (attr(x, "model") == "dynamic") cat("Shock time:", attr(x, "t0"), "\n")
	cat("\nSummary:\n")
	print(round(as.data.frame(x), digits))
	invisible(x)
}
####

####
#' Debug IRF
#'
#' @description Diagnostic output for troubleshooting IRF computation
#' @param fit A dbn model fit object
#' @param draw_idx Posterior draw index to debug
#' @param shock_i Source node
#' @param shock_j Target node
#' @seealso \code{\link{compute_irf}}, \code{\link{build_shock}},
#'   \code{\link{plot.dbn_irf}}
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 15, time = 50, seed = 1)
#' fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
#' debug_irf(fit, draw_idx = 1, shock_i = 1, shock_j = 2)
#' }
#' @export
debug_irf <- function(fit, draw_idx = 1, shock_i = 8, shock_j = 12) {
	cat("=== IRF Debug Information ===\n")

	cat("\n1. Model Structure:\n")
	cat("   Model type:", fit$model, "\n")
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

		cat("\n5. Testing C++ impulse_response_dynamic:\n")
		Delta <- impulse_response_dynamic(A_array, B_array, S, 46, 5)
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
