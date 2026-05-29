####
# exogenous covariate engine for dbn
#
# the bilinear dynamic model with covariates is
#   theta_t = L_t + A_t (theta_{t-1} - L_{t-1}) B_t^T + eps_t
# with linear predictor L_t = D_t %*% beta where D_t is the per-time
# design matrix assembled from dyad / row / col / actor components.
# beta has a conjugate gaussian full conditional so the gibbs update
# is closed form.
####

####
#' Construct a dbn covariate specification
#'
#' @description Builds the covariate object consumed by \code{dbn()} and
#'   \code{dbn_als()}. Accepts dyadic, sender-actor, receiver-actor, and
#'   shared-actor (symmetric) covariates as named lists of arrays. Stores
#'   the input arrays compactly (no \code{n_row x n_col x T} expansion of
#'   time-invariant covariates) and records the centering / scaling so
#'   downstream prediction can reapply it to new data.
#'
#' @details
#' Component semantics:
#' \itemize{
#'   \item \code{dyad}: dyadic covariates. Each element of the list is
#'     either an \code{n_row x n_col} matrix (time-invariant) or an
#'     \code{n_row x n_col x T} array (time-varying).
#'   \item \code{row}: sender-actor covariates. Each element is a length-
#'     \code{n_row} vector or \code{n_row x T} matrix. When entered into
#'     the design, the value is broadcast across receivers.
#'   \item \code{col}: receiver-actor covariates. Length \code{n_col} or
#'     \code{n_col x T}. Broadcast across senders.
#'   \item \code{actor}: shared actor covariates for symmetric square
#'     networks (\code{symmetric = TRUE}). Length \code{n_actor} or
#'     \code{n_actor x T}. The same coefficient applies to the sender and
#'     receiver positions, producing a symmetric contribution to the
#'     latent mean.
#' }
#'
#' Covariates are centered and scaled by default before the prior is
#' applied; the prior is then meaningfully N(0, 2.5^2) on the
#' standardized scale.
#'
#' @param dyad Named list of dyadic-covariate arrays.
#' @param row Named list of sender-actor covariate vectors / matrices.
#' @param col Named list of receiver-actor covariate vectors / matrices.
#' @param actor Named list of shared-actor covariates (symmetric fits).
#' @param standardize Logical (default TRUE). Center each covariate
#'   over the free dyad set and divide by its standard deviation.
#' @param coef_by One of \code{"shared"} (default; one coefficient per
#'   term shared across relations), \code{"relation"} (one per relation),
#'   or \code{"partial_pool"} (relation coefficients with hierarchical
#'   prior; deferred, currently aborts with a directive).
#' @return A \code{dbn_covariates} object.
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_covariates <- function(dyad = NULL,
                            row = NULL,
                            col = NULL,
                            actor = NULL,
                            standardize = TRUE,
                            coef_by = c("shared", "relation", "partial_pool")) {
	coef_by <- match.arg(coef_by)
	if (identical(coef_by, "partial_pool"))
		cli::cli_abort(c(
			"{.code coef_by = \"partial_pool\"} is not yet implemented.",
			"i" = "Use {.code \"shared\"} (default) or {.code \"relation\"}; partial pooling will land in a follow-up."
		))
	check_named_list <- function(x, label) {
		if (is.null(x)) return(NULL)
		if (!is.list(x))
			cli::cli_abort("{.arg {label}} must be a named list of arrays.")
		nm <- names(x)
		if (is.null(nm) || any(nm == "") || anyDuplicated(nm) > 0L)
			cli::cli_abort("{.arg {label}} must have unique non-empty names.")
		x
	}
	dyad <- check_named_list(dyad, "dyad")
	row <- check_named_list(row, "row")
	col <- check_named_list(col, "col")
	actor <- check_named_list(actor, "actor")
	if (is.null(dyad) && is.null(row) && is.null(col) && is.null(actor))
		cli::cli_abort("Pass at least one of {.arg dyad}, {.arg row}, {.arg col}, {.arg actor}.")
	if (!is.null(actor) && (!is.null(row) || !is.null(col)))
		cli::cli_abort("{.arg actor} is for symmetric square fits and must not be combined with {.arg row} or {.arg col}.")

	# wrap each entry in a small descriptor: keep the raw array and record
	# its detected shape so callers can broadcast time-invariant inputs
	# without duplicating the array
	wrap_entry <- function(value, name, component) {
		if (is.null(value))
			cli::cli_abort("Entry {.val {name}} in {.arg {component}} is NULL.")
		if (!is.numeric(value))
			cli::cli_abort("Entry {.val {name}} in {.arg {component}} must be numeric.")
		d <- dim(value)
		if (is.null(d)) d <- length(value)
		list(name = name, component = component, value = value, dim = d)
	}
	wrap_list <- function(xs, component) {
		if (is.null(xs)) return(list())
		lapply(seq_along(xs), function(i) wrap_entry(xs[[i]], names(xs)[i], component))
	}
	entries <- c(
		wrap_list(dyad, "dyad"),
		wrap_list(row, "row"),
		wrap_list(col, "col"),
		wrap_list(actor, "actor")
	)
	# detect duplicate term names across components
	all_names <- vapply(entries, function(e) e$name, character(1L))
	if (anyDuplicated(all_names) > 0L)
		cli::cli_abort(c(
			"Duplicate covariate name across components.",
			"i" = "Each covariate must have a unique name across {.arg dyad}, {.arg row}, {.arg col}, {.arg actor}."
		))
	# pre-extract standardization parameters as NA placeholders (filled when validated against data)
	terms <- data.frame(
		index = seq_along(entries),
		name = vapply(entries, function(e) e$name, character(1L)),
		component = vapply(entries, function(e) e$component, character(1L)),
		time_varying = NA,
		mean = NA_real_,
		sd = NA_real_,
		stringsAsFactors = FALSE
	)
	out <- list(
		version = 1L,
		entries = entries,
		terms = terms,
		standardize = isTRUE(standardize),
		coef_by = coef_by,
		validated = FALSE
	)
	class(out) <- c("dbn_covariates", "list")
	out
}

#' @author Tosin Salau and Shahryar Minhas
#' @export
print.dbn_covariates <- function(x, ...) {
	cli::cli_h2("dbn covariate specification (v{x$version})")
	cli::cli_text("{nrow(x$terms)} term{?s}; standardize = {.val {x$standardize}}; coef_by = {.val {x$coef_by}}")
	if (nrow(x$terms) > 0L) {
		cli::cli_text("Terms:")
		comps <- split(x$terms$name, x$terms$component)
		for (cmp in names(comps))
			cli::cli_li("{.field {cmp}}: {.val {comps[[cmp]]}}")
	}
	if (!isTRUE(x$validated))
		cli::cli_alert_info("Not yet validated against a data array; centering / scaling deferred until fit time.")
	invisible(x)
}

####
# validate a dbn_covariates object against a target data array shape
# fills in time_varying, mean, sd, standardized entries
.dbn_validate_covariates <- function(covars, n_row, n_col, p, Tt,
                                      symmetric = FALSE, mask = NULL) {
	stopifnot(inherits(covars, "dbn_covariates"))
	# mask is logical [n_row, n_col, p, Tt] or NULL (no masking beyond
	# diagonal for unipartite)
	for (i in seq_along(covars$entries)) {
		e <- covars$entries[[i]]
		comp <- e$component
		v <- e$value
		d <- e$dim
		# shape checks
		if (comp == "dyad") {
			if (length(d) == 2L) {
				if (!identical(as.integer(d), as.integer(c(n_row, n_col))))
					cli::cli_abort(c(
						"Dyadic covariate {.val {e$name}} has shape {paste(d, collapse = 'x')}; expected {n_row}x{n_col}{if (Tt > 1L) sprintf(' or %dx%dx%d', n_row, n_col, Tt) else ''}."
					))
				time_varying <- FALSE
			} else if (length(d) == 3L) {
				if (!identical(as.integer(d), as.integer(c(n_row, n_col, Tt))))
					cli::cli_abort("Dyadic covariate {.val {e$name}} has shape {paste(d, collapse = 'x')}; expected {n_row}x{n_col}x{Tt}.")
				time_varying <- TRUE
			} else {
				cli::cli_abort("Dyadic covariate {.val {e$name}} must be 2D ({n_row}x{n_col}) or 3D ({n_row}x{n_col}xT).")
			}
		} else if (comp == "row") {
			if (length(d) == 1L) {
				if (d[1] != n_row)
					cli::cli_abort("Row covariate {.val {e$name}} length {d[1]}; expected {n_row}.")
				time_varying <- FALSE
			} else if (length(d) == 2L) {
				if (!identical(as.integer(d), as.integer(c(n_row, Tt))))
					cli::cli_abort("Row covariate {.val {e$name}} shape {paste(d, collapse = 'x')}; expected {n_row} or {n_row}x{Tt}.")
				time_varying <- TRUE
			} else {
				cli::cli_abort("Row covariate {.val {e$name}} must be vector or 2D matrix.")
			}
		} else if (comp == "col") {
			if (length(d) == 1L) {
				if (d[1] != n_col)
					cli::cli_abort("Col covariate {.val {e$name}} length {d[1]}; expected {n_col}.")
				time_varying <- FALSE
			} else if (length(d) == 2L) {
				if (!identical(as.integer(d), as.integer(c(n_col, Tt))))
					cli::cli_abort("Col covariate {.val {e$name}} shape {paste(d, collapse = 'x')}; expected {n_col} or {n_col}x{Tt}.")
				time_varying <- TRUE
			} else {
				cli::cli_abort("Col covariate {.val {e$name}} must be vector or 2D matrix.")
			}
		} else if (comp == "actor") {
			if (!symmetric || n_row != n_col)
				cli::cli_abort("{.code component = \"actor\"} requires {.code symmetric = TRUE} with a square network.")
			n_a <- n_row
			if (length(d) == 1L) {
				if (d[1] != n_a)
					cli::cli_abort("Actor covariate {.val {e$name}} length {d[1]}; expected {n_a}.")
				time_varying <- FALSE
			} else if (length(d) == 2L) {
				if (!identical(as.integer(d), as.integer(c(n_a, Tt))))
					cli::cli_abort("Actor covariate {.val {e$name}} shape {paste(d, collapse = 'x')}; expected {n_a} or {n_a}x{Tt}.")
				time_varying <- TRUE
			} else {
				cli::cli_abort("Actor covariate {.val {e$name}} must be vector or 2D matrix.")
			}
		}
		covars$terms$time_varying[i] <- time_varying
		# symmetry check for dyadic covariates under symmetric model
		if (comp == "dyad" && isTRUE(symmetric)) {
			if (n_row != n_col)
				cli::cli_abort("symmetric = TRUE requires a square network; got {n_row}x{n_col}.")
			if (time_varying) {
				asym <- vapply(seq_len(Tt), function(t) {
					M <- v[, , t]
					max(abs(M - t(M)), na.rm = TRUE)
				}, numeric(1L))
				if (any(asym > sqrt(.Machine$double.eps)))
					cli::cli_abort(c(
						"Dyadic covariate {.val {e$name}} is not symmetric under {.code symmetric = TRUE}.",
						"x" = "Max asymmetry across t: {.val {signif(max(asym), 4)}}.",
						"i" = "Symmetrize via {.code X <- (X + t(X)) / 2} or pass an asymmetric model."
					))
			} else {
				asym <- max(abs(v - t(v)), na.rm = TRUE)
				if (asym > sqrt(.Machine$double.eps))
					cli::cli_abort(c(
						"Dyadic covariate {.val {e$name}} is not symmetric under {.code symmetric = TRUE}.",
						"x" = "Max asymmetry: {.val {signif(asym, 4)}}.",
						"i" = "Symmetrize via {.code X <- (X + t(X)) / 2} or pass an asymmetric model."
					))
			}
		}
	}
	# centering and scaling over the free dyad set
	# for now: free dyad set excludes diagonal on square networks
	# (mask handling for explicit user masks can extend later)
	if (isTRUE(covars$standardize)) {
		for (i in seq_along(covars$entries)) {
			e <- covars$entries[[i]]
			v <- e$value
			comp <- e$component
			tv <- covars$terms$time_varying[i]
			if (comp == "dyad") {
				if (tv) {
					vals_all <- as.numeric(v)
					# exclude diagonal for unipartite (n_row == n_col)
					if (n_row == n_col) {
						mask_diag <- rep(diag(n_row) == 1, Tt)
						vals_all <- vals_all[!mask_diag]
					}
				} else {
					vals_all <- as.numeric(v)
					if (n_row == n_col)
						vals_all <- vals_all[diag(n_row) == 0]
				}
			} else {
				vals_all <- as.numeric(v)
			}
			mu_k <- mean(vals_all, na.rm = TRUE)
			sd_k <- stats::sd(vals_all, na.rm = TRUE)
			if (!is.finite(sd_k) || sd_k < .Machine$double.eps) {
				cli::cli_abort(c(
					"Covariate {.val {e$name}} has zero variance over the free dyad set.",
					"x" = "Constant covariates are not identifiable; drop the term."
				))
			}
			covars$terms$mean[i] <- mu_k
			covars$terms$sd[i] <- sd_k
			# replace stored value with standardized version
			covars$entries[[i]]$value <- (v - mu_k) / sd_k
		}
	} else {
		# even without standardization, fill in mu and sd as 0 / 1 for the
		# accessor table; this lets predict() apply the same transformation
		covars$terms$mean <- 0
		covars$terms$sd <- 1
	}
	# rank check on the assembled per-time design at t = T (representative)
	# to detect aliased columns up front. assembled later in build_design_t.
	covars$n_row <- n_row
	covars$n_col <- n_col
	covars$p <- p
	covars$Tt <- Tt
	covars$symmetric <- symmetric
	covars$K <- nrow(covars$terms)
	covars$validated <- TRUE
	covars
}

####
# build the per-time design matrix at time t
# returns D_t with shape (n_row * n_col, K) in column-major vec order
.dbn_build_design_t <- function(covars, t) {
	stopifnot(inherits(covars, "dbn_covariates"), isTRUE(covars$validated))
	n_row <- covars$n_row
	n_col <- covars$n_col
	nm <- n_row * n_col
	K <- covars$K
	D <- matrix(0, nrow = nm, ncol = K)
	for (i in seq_along(covars$entries)) {
		e <- covars$entries[[i]]
		comp <- e$component
		tv <- covars$terms$time_varying[i]
		v <- e$value
		if (comp == "dyad") {
			Xt <- if (tv) v[, , t] else v
			D[, i] <- as.numeric(Xt)
		} else if (comp == "row") {
			rt <- if (tv) v[, t] else v
			Xt <- matrix(rt, nrow = n_row, ncol = n_col, byrow = FALSE)
			D[, i] <- as.numeric(Xt)
		} else if (comp == "col") {
			ct <- if (tv) v[, t] else v
			Xt <- matrix(ct, nrow = n_row, ncol = n_col, byrow = TRUE)
			D[, i] <- as.numeric(Xt)
		} else if (comp == "actor") {
			at <- if (tv) v[, t] else v
			Xt <- outer(at, rep(1, n_col)) + outer(rep(1, n_row), at)
			D[, i] <- as.numeric(Xt)
		}
	}
	colnames(D) <- covars$terms$name
	D
}

####
# build the linear predictor L_t at time t given coefficient vector beta
# returns an n_row x n_col matrix
.dbn_build_linear_predictor_t <- function(covars, beta, t) {
	D_t <- .dbn_build_design_t(covars, t)
	matrix(D_t %*% beta, nrow = covars$n_row, ncol = covars$n_col)
}

####
# build the transition design H_t = D_t - kron(B_t, A_t) D_{t-1}
# applies the kronecker product as matrix multiplication on reshape
# A_t is n_row x n_row, B_t is n_col x n_col
# D_{t-1} columns are vec(n_row x n_col matrices)
.dbn_build_transition_design <- function(covars, t, A_t, B_t, D_prev = NULL) {
	D_t <- .dbn_build_design_t(covars, t)
	if (is.null(D_prev))
		D_prev <- .dbn_build_design_t(covars, t - 1L)
	n_row <- covars$n_row
	n_col <- covars$n_col
	K <- ncol(D_t)
	H <- matrix(0, nrow = n_row * n_col, ncol = K)
	for (k in seq_len(K)) {
		Xk <- matrix(D_prev[, k], nrow = n_row, ncol = n_col)
		H[, k] <- D_t[, k] - as.numeric(A_t %*% Xk %*% t(B_t))
	}
	H
}

####
# gibbs update for beta given current state, operators, sigma^2
# returns a list with the new beta and its conditional mean/precision for tests
.dbn_update_beta <- function(covars, Theta, A_arr, B_arr, sigma2,
                              prior_mean = NULL, prior_sd = 2.5,
                              Tt = NULL) {
	stopifnot(inherits(covars, "dbn_covariates"), isTRUE(covars$validated))
	K <- covars$K
	if (is.null(prior_mean)) prior_mean <- rep(0, K)
	prior_prec <- 1 / prior_sd ^ 2
	# A_arr is n_row x n_row x Tt and B_arr is n_col x n_col x Tt
	# Theta is n_row x n_col x Tt
	dT <- dim(Theta)
	if (is.null(Tt)) Tt <- dT[3]
	# accumulate Λ = V_0^{-1} I + (1/σ²) Σ H_t' H_t, h = (1/σ²) Σ H_t' c_t
	HtH <- matrix(0, K, K)
	Htc <- numeric(K)
	get_A <- function(t) {
		if (length(dim(A_arr)) == 3L) A_arr[, , min(t, dim(A_arr)[3])] else A_arr
	}
	get_B <- function(t) {
		if (length(dim(B_arr)) == 3L) B_arr[, , min(t, dim(B_arr)[3])] else B_arr
	}
	D_prev <- .dbn_build_design_t(covars, 1L)
	for (t in 2:Tt) {
		D_t <- .dbn_build_design_t(covars, t)
		A_t <- get_A(t)
		B_t <- get_B(t)
		# transition design H_t
		H <- matrix(0, nrow = nrow(D_t), ncol = K)
		for (k in seq_len(K)) {
			Xk_prev <- matrix(D_prev[, k], covars$n_row, covars$n_col)
			H[, k] <- D_t[, k] - as.numeric(A_t %*% Xk_prev %*% t(B_t))
		}
		# transition residual c_t = vec(Theta_t - A_t Theta_{t-1} B_t')
		c_t <- as.numeric(Theta[, , t] - A_t %*% Theta[, , t - 1L] %*% t(B_t))
		HtH <- HtH + crossprod(H)
		Htc <- Htc + crossprod(H, c_t)
		D_prev <- D_t
	}
	Lambda <- diag(prior_prec, K) + HtH / sigma2
	rhs <- prior_prec * prior_mean + Htc / sigma2
	chol_L <- chol(Lambda)
	mu_post <- backsolve(chol_L, forwardsolve(t(chol_L), rhs))
	z <- rnorm(K)
	beta_draw <- mu_post + backsolve(chol_L, z)
	list(
		beta = as.numeric(beta_draw),
		mean = as.numeric(mu_post),
		cov = chol2inv(chol_L),
		HtH = HtH,
		Htc = Htc
	)
}

####
# als-style coordinate-descent update for beta (no draw, just argmax)
.dbn_als_update_beta <- function(covars, Theta, A_arr, B_arr,
                                  prior_mean = NULL, prior_sd = 2.5) {
	# same machinery as gibbs but return mu_post without adding noise
	res <- .dbn_update_beta(covars, Theta, A_arr, B_arr, sigma2 = 1,
		prior_mean = prior_mean, prior_sd = prior_sd, Tt = dim(Theta)[3])
	res$mean
}

####
# observation-residual gibbs update for beta
# uses the full conditional p(beta | R, Z, sigma2_obs) which arises from
# the observation eq Z_t = L_t + R_t + obs_noise. unlike the transition
# update this carries strong signal for time-invariant covariates (where
# the transition design H_t = D_t - A_t D_{t-1} B_t' collapses near zero
# when A_t, B_t are close to identity)
.dbn_update_beta_obs <- function(covars, Z, R, sigma2_obs,
                                  prior_mean = NULL, prior_sd = 2.5,
                                  mask = NULL) {
	stopifnot(inherits(covars, "dbn_covariates"), isTRUE(covars$validated))
	K <- covars$K
	if (is.null(prior_mean)) prior_mean <- rep(0, K)
	prior_prec <- 1 / prior_sd ^ 2
	n_row <- covars$n_row
	n_col <- covars$n_col
	# accept R, Z either as 3D (n_row x n_col x Tt) or as 4D with p=1
	if (length(dim(Z)) == 4L) Z <- Z[, , 1, , drop = TRUE]
	if (length(dim(R)) == 4L) R <- R[, , 1, , drop = TRUE]
	Tt <- dim(Z)[3]
	# off-diagonal cell mask for unipartite; default include all cells
	if (is.null(mask)) {
		mask <- if (n_row == n_col) as.logical(diag(n_row) == 0) else
			rep(TRUE, n_row * n_col)
	}
	mask_idx <- which(mask)
	DtD <- matrix(0, K, K)
	Dtr <- numeric(K)
	for (t in seq_len(Tt)) {
		D_t <- .dbn_build_design_t(covars, t)
		# residual: Z_t - R_t (the linear predictor target)
		res_t <- as.numeric(Z[, , t] - R[, , t])
		D_t_m <- D_t[mask_idx, , drop = FALSE]
		res_m <- res_t[mask_idx]
		DtD <- DtD + crossprod(D_t_m)
		Dtr <- Dtr + crossprod(D_t_m, res_m)
	}
	Lambda <- diag(prior_prec, K) + DtD / sigma2_obs
	rhs <- prior_prec * prior_mean + Dtr / sigma2_obs
	chol_L <- chol(Lambda)
	mu_post <- backsolve(chol_L, forwardsolve(t(chol_L), rhs))
	z <- stats::rnorm(K)
	beta_draw <- mu_post + backsolve(chol_L, z)
	list(
		beta = as.numeric(beta_draw),
		mean = as.numeric(mu_post),
		cov  = chol2inv(chol_L),
		DtD  = DtD,
		Dtr  = Dtr
	)
}

####
# time-varying beta FFBS
# state eq: beta_t = beta_{t-1} + N(0, diag(tau_beta2)), beta_1 ~ N(prior_mean, diag(prior_sd^2))
# obs eq:  (Z_t - R_t) | mask  = D_t |mask %*% beta_t + N(0, sigma2_obs I)
# returns beta_mat with shape K x Tt
.dbn_update_beta_tv <- function(covars, Z, R, sigma2_obs, tau_beta2,
                                  prior_mean = NULL, prior_sd = 2.5) {
	stopifnot(inherits(covars, "dbn_covariates"), isTRUE(covars$validated))
	K <- covars$K
	if (is.null(prior_mean)) prior_mean <- rep(0, K)
	n_row <- covars$n_row
	n_col <- covars$n_col
	if (length(dim(Z)) == 4L) Z <- Z[, , 1, , drop = TRUE]
	if (length(dim(R)) == 4L) R <- R[, , 1, , drop = TRUE]
	Tt <- dim(Z)[3]
	mask <- if (n_row == n_col) as.logical(diag(n_row) == 0) else
		rep(TRUE, n_row * n_col)
	mask_idx <- which(mask)
	V0_inv <- diag(1 / prior_sd^2, K)
	V0_inv_mu <- as.numeric(V0_inv %*% prior_mean)
	S_inv <- diag(1 / tau_beta2, K)

	# forward filter
	m_filt <- matrix(0, K, Tt)
	C_filt <- array(0, dim = c(K, K, Tt))
	C_pred_inv <- V0_inv
	C_pred_inv_mu <- V0_inv_mu

	for (t in seq_len(Tt)) {
		D_t <- .dbn_build_design_t(covars, t)
		D_m <- D_t[mask_idx, , drop = FALSE]
		resid_t <- as.numeric(Z[, , t] - R[, , t])
		r_m <- resid_t[mask_idx]
		DtD <- crossprod(D_m)
		Q_post <- C_pred_inv + DtD / sigma2_obs
		Dty <- as.numeric(crossprod(D_m, r_m))
		eta_post <- C_pred_inv_mu + Dty / sigma2_obs
		chol_Q <- chol(Q_post)
		m_t_t <- backsolve(chol_Q, forwardsolve(t(chol_Q), eta_post))
		C_t_t <- chol2inv(chol_Q)
		m_filt[, t] <- m_t_t
		C_filt[, , t] <- C_t_t
		if (t < Tt) {
			C_pred <- C_t_t + diag(tau_beta2, K)
			C_pred_inv <- solve(C_pred)
			C_pred_inv_mu <- as.numeric(C_pred_inv %*% m_t_t)
		}
	}

	# backward sample
	beta_mat <- matrix(0, K, Tt)
	chol_T <- chol(C_filt[, , Tt])
	beta_mat[, Tt] <- m_filt[, Tt] +
		as.numeric(crossprod(chol_T, stats::rnorm(K)))
	if (Tt >= 2L) {
		for (t in (Tt - 1L):1L) {
			C_tt_inv <- solve(C_filt[, , t])
			Q_back <- C_tt_inv + S_inv
			eta_back <- as.numeric(C_tt_inv %*% m_filt[, t]) +
				as.numeric(S_inv %*% beta_mat[, t + 1L])
			chol_Q <- chol(Q_back)
			m_back <- backsolve(chol_Q, forwardsolve(t(chol_Q), eta_back))
			C_back <- chol2inv(chol_Q)
			chol_C_back <- chol(C_back)
			beta_mat[, t] <- m_back +
				as.numeric(crossprod(chol_C_back, stats::rnorm(K)))
		}
	}
	rownames(beta_mat) <- covars$terms$name
	list(beta_mat = beta_mat)
}
