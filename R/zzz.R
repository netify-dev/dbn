####
# package initialization
####

#' @keywords internal
.onLoad <- function(libname, pkgname) {
	# compiled code is loaded via useDynLib(dbn, .registration = TRUE) in
	# the NAMESPACE; a broken shared library fails loudly at namespace load.

	# create the package env while the namespace is still mutable; later
	# assignments would hit a locked binding error.
	ns <- getNamespace(pkgname)
	if (!exists(".pkg_env", envir = ns, inherits = FALSE)) {
		assign(".pkg_env", new.env(parent = emptyenv()), envir = ns)
		assign("models", list(), envir = get(".pkg_env", envir = ns))
	}

	# default options
	op <- options()
	op.dbn <- list(
		dbn.use_cpp_stability = TRUE,
		dbn.use_ffbs_dlm_cpp = TRUE,
		dbn.use_ffbs_cpp = TRUE,
		dbn.use_cpp_ranklik = TRUE,
		dbn.use_batch_ffbs = FALSE,
		dbn.n_threads = 1L
	)
	toset <- !(names(op.dbn) %in% names(op))
	if (any(toset)) options(op.dbn[toset])

	# register the local gof generic + dbn method so dispatch works under
	# devtools::load_all and after install.
	registerS3method("gof", "dbn", gof.dbn, envir = asNamespace(pkgname))

	# register as_draws methods on the posterior generic when posterior is
	# available, so posterior::as_draws(fit) dispatches to our methods.
	if (requireNamespace("posterior", quietly = TRUE)) {
		registerS3method("as_draws", "dbn", as_draws.dbn,
		                 envir = asNamespace("posterior"))
		registerS3method("as_draws", "dbn_multichain", as_draws.dbn_multichain,
		                 envir = asNamespace("posterior"))
	}

	# register the forecast method on forecast::forecast when forecast is
	# available, so forecast(dbn_fit, h = ...) dispatches to our shim.
	if (requireNamespace("forecast", quietly = TRUE)) {
		registerS3method("forecast", "dbn", forecast.dbn,
		                 envir = asNamespace("forecast"))
	}

	# register autoplot methods on ggplot2::autoplot when ggplot2 is
	# available. Keeps ggplot2 in Suggests but lets `autoplot(fc)` dispatch
	# to our plot.dbn_forecast{,_ci} / plot.dbn_rank_probs wrappers.
	if (requireNamespace("ggplot2", quietly = TRUE)) {
		registerS3method("autoplot", "dbn_forecast", autoplot.dbn_forecast,
		                 envir = asNamespace("ggplot2"))
		registerS3method("autoplot", "dbn_forecast_ci", autoplot.dbn_forecast_ci,
		                 envir = asNamespace("ggplot2"))
		registerS3method("autoplot", "dbn_rank_probs", autoplot.dbn_rank_probs,
		                 envir = asNamespace("ggplot2"))
	}

	# register coda's as.mcmc when coda is available.
	if (requireNamespace("coda", quietly = TRUE)) {
		registerS3method("as.mcmc", "dbn", as.mcmc.dbn,
		                 envir = asNamespace("coda"))
	}

	# register posterior_predict and log_lik on rstantools
	if (requireNamespace("rstantools", quietly = TRUE)) {
		registerS3method("posterior_predict", "dbn", posterior_predict.dbn,
		                 envir = asNamespace("rstantools"))
		registerS3method("log_lik", "dbn", log_lik.dbn,
		                 envir = asNamespace("rstantools"))
	}

	# register tidybayes::tidy_draws methods
	if (requireNamespace("tidybayes", quietly = TRUE)) {
		registerS3method("tidy_draws", "dbn", tidy_draws.dbn,
		                 envir = asNamespace("tidybayes"))
		registerS3method("tidy_draws", "dbn_multichain", tidy_draws.dbn_multichain,
		                 envir = asNamespace("tidybayes"))
	}

	invisible()
}
