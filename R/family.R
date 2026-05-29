####
#' Family Objects for DBN Models
#'
#' @description Internal constructors for likelihood family objects
#' @name family
#' @keywords internal
NULL
####

####
#' Family Constructor
#'
#' @description Builds a dbn_family structure holding all family-specific callbacks
#' @param name Family name string
#' @param draw_latent Latent Z sampler
#' @param ffbs_wrapper FFBS dispatch wrapper
#' @param loglik Log-likelihood function
#' @param linkinv Inverse link function
#' @param rgen_obs Observation generator
#' @param init_pars Initial parameter list
#' @return A dbn_family object
#' @keywords internal
dbn_make_family <- function(name, draw_latent, ffbs_wrapper, loglik,
							linkinv, rgen_obs, init_pars = list()) {
	# pin closure environments to the package namespace and strip srcref so
	# the family list is identical across consecutive calls and survives a
	# saveRDS/readRDS round-trip under identical(). without this, each
	# family_*() call captures a fresh enclosing frame and the srcref
	# attribute carries a non-roundtrippable srcfile env, so identical()
	# always returns FALSE even on functionally-equivalent fits.
	stable_env <- topenv()
	.strip_src <- function(x) {
		# walk language objects and strip srcref-related attributes everywhere
		# they hide. R attaches these at parse time on nested { } blocks too,
		# and they carry a non-roundtrippable srcfile env that breaks
		# identical() after saveRDS/readRDS.
		if (is.call(x) || is.pairlist(x)) {
			attr(x, "srcref") <- NULL
			attr(x, "srcfile") <- NULL
			attr(x, "wholeSrcref") <- NULL
			for (i in seq_along(x)) {
				if (!is.null(x[[i]])) x[[i]] <- .strip_src(x[[i]])
			}
		}
		x
	}
	.stabilize <- function(f) {
		if (!is.function(f)) return(f)
		environment(f) <- stable_env
		attr(f, "srcref") <- NULL
		body(f) <- .strip_src(body(f))
		f
	}
	draw_latent  <- .stabilize(draw_latent)
	ffbs_wrapper <- .stabilize(ffbs_wrapper)
	loglik       <- .stabilize(loglik)
	linkinv      <- .stabilize(linkinv)
	rgen_obs     <- .stabilize(rgen_obs)
	structure(
		list(
			name = name,
			draw_latent = draw_latent,
			ffbs_wrapper = ffbs_wrapper,
			loglik = loglik,
			linkinv = linkinv,
			rgen_obs = rgen_obs,
			init_pars = init_pars
		),
		class = "dbn_family"
	)
}
####

####
#' Ordinal Family
#'
#' @description Rank likelihood family for ordinal outcomes.
#'   Observation variance fixed at 1 for identifiability.
#' @return A dbn_family object
#' @keywords internal
family_ordinal <- function() {
	dbn_make_family(
		name = "ordinal",

		draw_latent = function(pre, ...) {
			for (j in seq_len(pre$dims$p)) {
				EZ <- pre$Theta[, , j, ]
				for (t in 1:pre$dims$Tt) {
					EZ[, , t] <- EZ[, , t] + pre$M[, , j]
				}
				pre$Z[, , j, ] <- rz_fc(pre$R[, , j, ], pre$Z[, , j, ], EZ, pre$IR[[j]])
			}
			pre
		},

		ffbs_wrapper = function(Z, mu, Aarray, Barray, sigma2_proc, sigma2_obs, ...) {
			if (!missing(sigma2_obs) && sigma2_obs != 1) {
				warning("observation variance sigma2_obs is fixed at 1 for ordinal family (identifiability constraint)")
			}
			ffbs_theta(Z, mu, Aarray, Barray, sigma2_proc)
		},

		loglik = step_ll,

		linkinv = function(theta, misc) {
			if (is.null(misc$cuts)) {
				return(theta)
			}
			K <- length(misc$cuts) + 1
			cts <- seq_len(K)
			expected <- array(NA, dim = dim(theta))
			for (i in seq_len(length(theta))) {
				cum_probs <- c(0, pnorm(misc$cuts - theta[i]), 1)
				probs <- diff(cum_probs)
				expected[i] <- sum(cts * probs)
			}
			expected
		},

		rgen_obs = function(theta, misc) {
			theta_mean <- add_baseline_mean(theta, misc$M)
			if (is.null(misc$cuts)) {
				misc$cuts <- qnorm(seq(0.2, 0.8, by = 0.2))
			}

			u <- array(runif(length(theta_mean)), dim = dim(theta_mean))
			ranks <- array(NA, dim = dim(theta_mean))
			for (i in seq_len(length(theta_mean))) {
				cum_probs <- pnorm(misc$cuts - theta_mean[i])
				ranks[i] <- findInterval(u[i], c(0, cum_probs, 1))
			}

			n_cats <- length(misc$cuts) + 1
			pmin(pmax(ranks, 1L), n_cats)
		}
	)
}
####

####
#' Gaussian Family
#'
#' @description Identity-link family for continuous outcomes
#' @return A dbn_family object
#' @keywords internal
family_gaussian <- function() {
	dbn_make_family(
		name = "gaussian",

		draw_latent = function(pre, ...) pre,

		ffbs_wrapper = function(Z, mu, Aarray, Barray, sigma2_proc, sigma2_obs, ...) {
			ffbs_theta_struct_cpp(Z, mu, Aarray, Barray, sigma2_proc, sigma2_obs)
		},

		loglik = function(A, B, Theta_prev, Theta_curr, fam_pars) {
			sigma2_obs <- fam_pars$sigma2_obs
			resid <- Theta_curr - A %*% Theta_prev %*% t(B)
			-0.5 * sum(resid^2) / sigma2_obs -
				0.5 * length(resid) * log(2 * pi * sigma2_obs)
		},

		linkinv = function(theta, misc, ...) {
			theta
		},

		rgen_obs = function(theta, misc, sigma2_obs = 1) {
			theta_mean <- add_baseline_mean(theta, misc$M)
			if (!is.null(misc$sigma2_obs)) {
				sigma2_obs <- misc$sigma2_obs
			}
			theta_mean + array(rnorm(length(theta_mean), sd = sqrt(sigma2_obs)), dim(theta_mean))
		},

		init_pars = list(sigma2_obs = 1)
	)
}
####

####
#' Binary Family
#'
#' @description Probit-link family for binary outcomes.
#'   Observation variance fixed at 1 for identifiability.
#' @return A dbn_family object
#' @keywords internal
family_binary <- function() {
	dbn_make_family(
		name = "binary",

		draw_latent = function(pre, Aarray, Barray, ...) {
			if (!requireNamespace("truncnorm", quietly = TRUE)) {
				cli::cli_abort(c(
					"Package {.pkg truncnorm} is required for binary outcomes.",
					"i" = "Install with {.code install.packages(\"truncnorm\")}"
				))
			}
			for (t in 1:pre$dims$Tt) {
				for (rel in 1:pre$dims$p) {
					eta_t <- pre$Theta[, , rel, t] + pre$M[, , rel]
					Y_t <- pre$R[, , rel, t]
					Z_t <- pre$Z[, , rel, t]

					pos <- which(Y_t == 1)
					neg <- which(Y_t == 0)
					if (length(pos) > 0) {
						Z_t[pos] <- truncnorm::rtruncnorm(length(pos),
							a = 0, b = Inf, mean = eta_t[pos], sd = 1)
					}
					if (length(neg) > 0) {
						Z_t[neg] <- truncnorm::rtruncnorm(length(neg),
							a = -Inf, b = 0, mean = eta_t[neg], sd = 1)
					}
					pre$Z[, , rel, t] <- Z_t
				}
			}
			pre
		},

		ffbs_wrapper = function(Z, mu, Aarray, Barray, sigma2_proc, sigma2_obs = 1, ...) {
			if (!missing(sigma2_obs) && sigma2_obs != 1) {
				warning("observation variance sigma2_obs is fixed at 1 for binary family (probit identifiability)")
			}
			ffbs_theta(Z, mu, Aarray, Barray, sigma2_proc)
		},

		loglik = step_ll,

		linkinv = function(theta, misc) {
			pnorm(theta)
		},

		rgen_obs = function(theta, misc) {
			theta_mean <- add_baseline_mean(theta, misc$M)
			probs <- pnorm(theta_mean)
			array(rbinom(length(theta_mean), 1, probs), dim = dim(theta_mean))
		},

		init_pars = list()
	)
}
####

####
#' Add Baseline Mean to Theta
#'
#' @description Broadcasts M across time dimension of theta
#' @param theta Latent array
#' @param M Baseline mean (may be NULL)
#' @return theta + M with appropriate broadcasting
#' @keywords internal
add_baseline_mean <- function(theta, M) {
	if (is.null(M)) return(theta)

	dim_theta <- dim(theta)
	dim_M <- dim(M)

	if (length(dim_M) == 3 && length(dim_theta) == 4) {
		sweep(theta, 1:3, M, "+")
	} else if (length(dim_M) == length(dim_theta) && all(dim_M == dim_theta)) {
		theta + M
	} else {
		cli::cli_abort("Incompatible dimensions: theta is {paste(dim_theta, collapse = 'x')}, M is {paste(dim_M, collapse = 'x')}.")
	}
}
####
