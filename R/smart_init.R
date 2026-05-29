####
# smart initialization for the MCMC sampler
#
# motivation: the default starting values (A = B = I, sigma2 = 1, etc.)
# are mediocre when the data are at any scale other than ~unit variance.
# the chain walks slowly from identity toward the truth, which inflates
# burn-in. running a few ALS iterations gives initial values close to the
# posterior mode at near-zero cost; the user opts in via init = "smart"
# (default) and falls back to identity via init = "default".
####

####
#' Smart starting values for the MCMC sampler
#'
#' @description Runs a fast (~5 ALS iteration) point estimate on the data
#'   and uses the result to initialise the MCMC chain. Cuts burn-in by an
#'   order of magnitude on signal-rich data. Used internally by
#'   \code{\link{dbn}} when \code{init = "smart"} (the default).
#'
#' @param Y The 4D data array.
#' @param family Family name.
#' @param model One of \code{"static"}, \code{"dynamic"}, \code{"piecewise"},
#'   \code{"hmm"}, \code{"lowrank"}.
#' @param symmetric Logical; passed to the inner ALS run.
#' @param Tt Time dimension.
#' @param n_row,n_col,p Data dimensions.
#' @param verbose Logical.
#' @return A list with \code{A_init} (n_row x n_row x Tt), \code{B_init}
#'   (n_col x n_col x Tt), \code{M_init} (n_row x n_col x p),
#'   \code{sigma2_init}, \code{sigma2_obs_init}, \code{tauA2_init},
#'   \code{tauB2_init}, \code{g2_init}.
#' @keywords internal
#' @noRd
.dbn_smart_init <- function(Y, family, model, symmetric,
                              Tt, n_row, n_col, p,
                              verbose = FALSE) {
	# fast identity fallback if ALS would be heavier than the eventual
	# MCMC (small data) or if ALS isn't well-defined for the family
	if (Tt < 3L || min(n_row, n_col) < 2L) {
		return(.dbn_default_init(Tt, n_row, n_col, p))
	}

	# call dbn_als briefly. wrapped because dbn_als has heavy preprocessing;
	# any failure falls through to the default identity init so a bad warm
	# start never blocks an MCMC run
	als_fit <- tryCatch(
		suppressMessages(suppressWarnings(
			dbn_als(Y, family = family, symmetric = symmetric,
				bootstrap = 0, max_iter = 5L, verbose = FALSE)
		)),
		error = function(e) NULL
	)
	if (is.null(als_fit) || !is.list(als_fit$A) || !is.list(als_fit$B)) {
		return(.dbn_default_init(Tt, n_row, n_col, p))
	}

	# dbn_als returns A/B as lists of cubes (one cube per "draw" - here
	# just one point estimate); extract the first cube
	A_als <- als_fit$A[[1L]]
	B_als <- als_fit$B[[1L]]
	M_als <- als_fit$M
	# dbn_als M is [n_row, n_col, p, 1]; squeeze the last dim
	if (length(dim(M_als)) == 4L) {
		M_init <- M_als[, , , 1L, drop = TRUE]
		dim(M_init) <- c(n_row, n_col, p)
	} else {
		M_init <- M_als
	}

	# A_als / B_als are [n_row, n_row, Tt] cubes. for the static path the
	# cube has the same matrix at every t; for dynamic it varies per t
	# (TV-ALS). either way we tile / pass it through directly
	if (length(dim(A_als)) == 3L && dim(A_als)[3] == Tt) {
		A_init <- A_als
		B_init <- B_als
	} else if (length(dim(A_als)) == 3L) {
		# tile a static estimate across t
		A_init <- array(rep(c(A_als[, , 1L]), Tt), dim = c(n_row, n_row, Tt))
		B_init <- array(rep(c(B_als[, , 1L]), Tt), dim = c(n_col, n_col, Tt))
	} else {
		A_init <- array(rep(c(A_als), Tt), dim = c(n_row, n_row, Tt))
		B_init <- array(rep(c(B_als), Tt), dim = c(n_col, n_col, Tt))
	}

	# pin A_t = B_t = I at t = 1 (the identifiability anchor for the RW prior)
	A_init[, , 1L] <- diag(n_row)
	B_init[, , 1L] <- diag(n_col)

	# variance components: estimate from residuals
	# state-eq residual r_t = Theta_t - A_t Theta_{t-1} B_t' (approximated via Y)
	# obs-eq variance: rough cross-sectional variance of Y
	Y_finite <- as.numeric(Y[is.finite(Y)])
	y_var <- if (length(Y_finite) > 1L) stats::var(Y_finite) else 1
	# state innovation variance: half of total
	sigma2_init <- max(0.05, y_var / 2)
	sigma2_obs_init <- max(0.05, y_var / 2)
	# operator RW innovation variance: small relative to state noise
	tauA2_init <- max(0.001, sigma2_init / 100)
	tauB2_init <- tauA2_init
	# M variance: from the magnitude of M_init
	g2_init <- max(0.1, mean(M_init^2, na.rm = TRUE))

	if (isTRUE(verbose)) {
		cli::cli_alert_info(
			"Smart init: sigma2={sprintf('%.3f', sigma2_init)}, tauA2={sprintf('%.4f', tauA2_init)}, g2={sprintf('%.3f', g2_init)}"
		)
	}

	list(
		A_init           = A_init,
		B_init           = B_init,
		M_init           = M_init,
		sigma2_init      = sigma2_init,
		sigma2_obs_init  = sigma2_obs_init,
		tauA2_init       = tauA2_init,
		tauB2_init       = tauB2_init,
		g2_init          = g2_init,
		source           = "als_warm_start"
	)
}

####
# default identity init (legacy behavior). Used as a fallback when smart
# init fails or when the user requests `init = "default"`.
.dbn_default_init <- function(Tt, n_row, n_col, p) {
	A_init <- array(0, dim = c(n_row, n_row, Tt))
	B_init <- array(0, dim = c(n_col, n_col, Tt))
	for (t in seq_len(Tt)) {
		A_init[, , t] <- diag(n_row)
		B_init[, , t] <- diag(n_col)
	}
	list(
		A_init           = A_init,
		B_init           = B_init,
		M_init           = array(0, dim = c(n_row, n_col, p)),
		sigma2_init      = 1,
		sigma2_obs_init  = 1,
		tauA2_init       = 1,
		tauB2_init       = 1,
		g2_init          = 1,
		source           = "identity_default"
	)
}
