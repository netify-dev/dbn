####
#' Dynamic Bilinear Network Analysis
#'
#' @description
#' Main entry point for fitting Dynamic Bilinear Network (DBN) models. These
#' models estimate how past network interactions predict future interactions,
#' recovering time-varying influence structures from temporal relational data.
#'
#' For a fast point-estimate (seconds, no MCMC) call [dbn_als()] -- the
#' alternating-least-squares estimator -- and pair it with `bootstrap = N`
#' for entry-wise standard errors. `dbn()` itself runs MCMC and returns
#' full posterior draws; reach for it when you want credible intervals,
#' model comparison, or predictive checks.
#'
#' The core model is the deviation-form bilinear AR
#' \eqn{\Theta_t = M + A_t (\Theta_{t-1} - M) B_t' + \varepsilon_t},
#' where \eqn{A_t} captures sender influence, \eqn{B_t} captures receiver
#' influence, and \eqn{M} captures stable dyad-specific tendencies.
#'
#' @details
#' **Choosing a family:**
#' - `"gaussian"`: Use when your relational data are continuous measurements
#'   (e.g., trade volumes, similarity scores). This is the simplest family and
#'   converges fastest.
#' - `"ordinal"`: Use when your data are ordered categories or counts (e.g.,
#'   conflict severity 1-5, event counts) and you trust the ordering but not
#'   the exact values. Uses a rank likelihood.
#' - `"binary"`: Use for presence/absence data (0/1), e.g., whether a tie
#'   exists. Uses a probit link with data augmentation.
#'
#' **Choosing a model:**
#' - `"static"`: Simplest. Influence structure is fixed over time. Good
#'   starting point and for short time series.
#' - `"dynamic"`: Influence structure changes over time. Use when you expect
#'   shifting alliances, evolving trade patterns, etc.
#' - `"piecewise"`: Influence is constant within known regimes but differs
#'   across them. Use when you know when structural breaks occurred (e.g.,
#'   before/after a crisis).
#' - `"hmm"`: Like piecewise but discovers regimes from data. Use when breaks
#'   are unknown.
#' - `"lowrank"`: Like dynamic but with dimensionality reduction for large
#'   networks (50+ actors).
#'
#' **Fast path (point estimate, no MCMC):**
#' For exploratory work or when you need results in seconds rather than
#' minutes, use [dbn_als()]. ALS returns the same operator structure as
#' `dbn()` but as a single point estimate; pair it with `bootstrap = N`
#' for entry-wise standard errors. The MCMC `dbn()` entry point below
#' returns full posterior draws and is the right choice for inference,
#' credible intervals, model comparison, and predictive checks.
#'
#' **MCMC settings:**
#' The sampler draws `nscan` posterior samples after discarding the first
#' `burn` as warm-up. Setting `odens > 1` thins the output by saving every
#' k-th sample. For initial exploration, `nscan = 5000, burn = 2000, odens = 5`
#' is a reasonable starting point. For final inference use longer
#' chains (`nscan = 10000+`) and verify convergence with [check_convergence()].
#' For multi-chain Rhat / posterior-package interop, see [dbn_multichain()].
#'
#' @param data Numeric array of network data, or a file path to an `.RData`
#'   file that contains an object named `Y`.  The array should be
#'   3-dimensional `[actors, actors, time]` for a single relation type, or
#'   4-dimensional `[actors, actors, relations, time]` for multiple relation
#'   types.  Diagonal entries (self-ties) should be `NA` for unipartite
#'   networks.  For bipartite networks, pass a rectangular array where the
#'   first dimension (senders) differs from the second (receivers).
#'
#'   Missing values are supported: `NA` entries (including the `NA` diagonal)
#'   are treated as missing and imputed within the sampler, so partially
#'   observed networks can be fit directly. The data must, however, contain at
#'   least some finite observations; `NaN` and infinite values are anomalies
#'   (use `NA` to mark missing data) and trigger a warning.
#' @param family Character string specifying the outcome distribution.
#'   See **Details** for guidance on choosing:
#'   \itemize{
#'     \item `"ordinal"`: For ordinal/ranked data (positive integers)
#'     \item `"gaussian"`: For continuous data (any real numbers)
#'     \item `"binary"`: For binary data (0/1 or logical)
#'   }
#' @param model Character string specifying the model type.
#'   See **Details** for guidance on choosing:
#'   \itemize{
#'     \item `"static"`: Fixed sender/receiver effects across time
#'     \item `"dynamic"`: Time-varying sender/receiver effects
#'     \item `"lowrank"`: Low-rank factorization of sender effects (large networks)
#'     \item `"hmm"`: Regime-switching with data-driven regime discovery
#'     \item `"piecewise"`: Block-constant influence with known break points
#'   }
#' @param method Character string selecting the estimation method.
#'   \itemize{
#'     \item `"mcmc"` (default): Full Bayesian posterior via MCMC. Slow but
#'       gives credible intervals from the posterior directly.
#'     \item `"als"`: Fast alternating-least-squares point estimate via
#'       \code{\link{dbn_als}}. Routes to the dynamic-model ALS regardless of
#'       the requested `model`. MCMC-only arguments (`nscan`, `burn`, `odens`,
#'       `sampler`) are ignored on this path. Combine with `bootstrap = N` for
#'       resampling-based CIs in the same call.
#'   }
#' @param bootstrap Integer. Only used when `method = "als"`. If `> 0`, runs
#'   that many bootstrap replicates after the point-estimate fit and attaches
#'   CIs to the returned object (see \code{\link{dbn_als}}). Default `0`.
#' @param nscan Number of posterior samples to draw after burn-in (MCMC only)
#' @param burn Number of initial MCMC samples to discard (warm-up period)
#' @param odens Thinning interval: save every odens-th sample (reduces autocorrelation and memory)
#' @param verbose Logical or numeric. If TRUE, show progress. If numeric, print detailed info every n iterations (default: TRUE)
#' @param symmetric Logical. If TRUE, enforce B = A (symmetric / undirected network). Requires a square network (n_row == n_col). Supported on the static, dynamic, piecewise, and HMM model paths; rejected with an informative error on the low-rank model. Default: FALSE.
#' @param sampler Character string specifying the inference sampler (dynamic model only).
#'   \itemize{
#'     \item `"auto"` (default): Smart choice -- exact PCG for symmetric, FFBS for asymmetric.
#'     \item `"exact"`: Exact posterior via preconditioned conjugate gradient (PCG).
#'       Currently only implemented for symmetric fits (`symmetric = TRUE`); on an
#'       asymmetric fit a warning is issued and the sampler silently falls back
#'       to the approximate FFBS path.
#'     \item `"approx"`: Approximate posterior via forward-filtering/backward-sampling (FFBS).
#'   }
#' @param ... Additional model-specific parameters:
#'   \describe{
#'     \item{\code{r}}{Rank for lowrank model (default: 2)}
#'     \item{\code{R}}{Number of regimes for HMM model (default: 3)}
#'     \item{\code{blocks}}{Block specification for piecewise model: integer (number of equal blocks),
#'       numeric vector (block boundaries), named vector (labeled boundaries), or "auto" for automatic selection}
#'     \item{\code{ar1}}{Use AR(1) dynamics for dynamic model (default: FALSE)}
#'     \item{\code{update_rho}}{Update AR coefficient in dynamic model (default: FALSE)}
#'     \item{\code{seed}}{Random seed for reproducibility (default: 6886)}
#'     \item{\code{previous}}{Previous fit object to continue MCMC from}
#'     \item{\code{init}}{List of initial values for parameters}
#'     \item{\code{time_thin}}{Time thinning factor for dynamic/lowrank/HMM (default: auto for dynamic, 1 for others)}
#'     \item{\code{store_z}}{Store Z draws for dynamic model (default: auto based on memory)}
#'     \item{\code{store_theta}}{Store full Theta trajectory draws for piecewise model (default: TRUE).
#'       \strong{Critical for large networks:} Set to FALSE for networks with 100+ actors to avoid
#'       memory issues. Theta storage scales as O(n^2 * T * draws) -- a 200-actor network with 50 time
#'       points and 500 draws requires ~40 GB. With \code{store_theta = FALSE}, you retain posterior
#'       draws for A, B, M and variance parameters, \code{compare_blocks()} functionality, and
#'       convergence diagnostics. You lose full posterior uncertainty on individual Theta entries
#'       and \code{posterior_predict_dbn()} with uncertainty propagation.}
#'   }
#' @return A list of class `"dbn"` with model-specific contents. Common elements:
#'   \item{model}{Character string indicating which model was fit}
#'   \item{family}{Character string indicating the outcome family}
#'   \item{dims}{List of data dimensions (n_row, n_col, p, Tt)}
#'   \item{settings}{List of MCMC settings used}
#'   \item{Y}{Original data array}
#'   \item{M}{Posterior draws for the baseline mean M, an array with the draw
#'     index as the last dimension}
#'   \item{Theta}{Posterior draws for the latent network state, an array with
#'     the draw index as the last dimension}
#'
#'   Model-specific elements include:
#'   \item{A, B}{Posterior draws for the sender and receiver influence
#'     operators. For the dynamic model these are \emph{lists} of length
#'     `draws`, each element an `[n_row, n_row, T]` (respectively
#'     `[n_col, n_col, T]`) array. The same A and B are shared across all
#'     `p` relations: shapes are `[n_row, n_row, T]` whether `p = 1` or
#'     `p > 1`, and each relation evolves through the same operator pair
#'     (this is the bilinear sharing assumption -- relations differ in M
#'     and Theta, not in A/B). The posterior-mean operator at time `t` is,
#'     for example,
#'     `Reduce("+", lapply(fit$A, function(a) a[, , t])) / length(fit$A)`.
#'     The \strong{static} model has no time dimension and stores its operator
#'     draws differently: `fit$A` is `NULL`, and the draws are held together in
#'     `fit$B`, a three-element list of mode-wise operators -- `fit$B[[1]]` the
#'     `[n_row, n_row, draws]` sender operator (A), `fit$B[[2]]` the
#'     `[n_col, n_col, draws]` receiver operator (B), and `fit$B[[3]]` the
#'     `[p, p, draws]` relation-mode operator (`1 x 1` for single-relation data).}
#'   \item{sigma2, tau_A2, tau_B2, g2, sigma2_obs}{Posterior draws for the
#'     variance parameters. `sigma2` is the process (state-innovation)
#'     variance of the latent state; `sigma2_obs` is the observation-noise
#'     variance that separates `Y` from `Theta`; `tau_A2` and `tau_B2` are
#'     the random-walk innovation variances of the operators `A_t` and `B_t`
#'     (under the default RW prior on both symmetric and asymmetric paths),
#'     with `A_0 = B_0 = I` as the structural anchor. For the iid prior on
#'     the asymmetric path, set `options(dbn.prior_kind = "iid")`.
#'     `g2` is the variance of the baseline mean `M`.}
#'   \item{rhoA, rhoB}{AR(1) persistence parameters (dynamic model with ar1=TRUE)}
#'   \item{A_blocks}{List of regime-specific posterior mean A matrices (piecewise)}
#'   \item{time_kept}{Which time indices are stored (dynamic/lowrank/HMM)}
#'   \item{meta}{List of fit metadata, including `meta$runtime_sec` (wall-clock
#'     fit time in seconds) and the sampler/family/dimension records.}
#'
#'   Use [summary()], [plot()], [param_summary()], and [check_convergence()]
#'   to inspect results. See model-specific vignettes for full workflows.
#' @author Tosin Salau and Shahryar Minhas
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 10, seed = 6886)
#'
#' # static model with gaussian family
#' fit <- dbn(sim$Z, model = "static", family = "gaussian",
#'     nscan = 200, burn = 100, verbose = FALSE)
#'
#' # dynamic model
#' fit_dyn <- dbn(sim$Z, model = "dynamic", family = "gaussian",
#'     nscan = 200, burn = 100, verbose = FALSE)
#'
#' # piecewise model with 2 blocks
#' fit_pw <- dbn(sim$Y, model = "piecewise", blocks = 2,
#'     nscan = 200, burn = 100, verbose = FALSE)
#' }
####
dbn <- function(data,
				family = c("ordinal", "gaussian", "binary"),
				model = c("static", "dynamic", "lowrank", "hmm", "piecewise"),
				method = c("mcmc", "als"),
				nscan = 10000,
				burn = 1000,
				odens = 1,
				verbose = TRUE,
				symmetric = FALSE,
				sampler = "auto",
				bootstrap = 0,
				...) {
	# auto-detect gaussian when family was not specified and data looks continuous
	family_was_default <- missing(family)
	if (family_was_default && is.numeric(data)) {
		obs <- data[!is.na(data)]
		if (length(obs) > 0L) {
			looks_continuous <- !all(obs == round(obs)) || length(unique(obs)) > 7L
			looks_binary <- length(unique(obs)) <= 2L && all(obs %in% c(0, 1))
			if (looks_continuous && !looks_binary) {
				if (isTRUE(verbose > 0))
					cli::cli_inform(c(
						"i" = "{.arg family} not specified; defaulting to {.val gaussian} based on data values.",
						"i" = "Pass {.code family = \"ordinal\"} or {.code \"binary\"} explicitly to override."
					))
				family <- "gaussian"
			}
		}
	}
	family <- match.arg(family)
	# catch the common `model = "als"` mistake before match.arg fires so
	# the user gets a directive hint instead of a cryptic enum error.
	if (length(model) == 1L && is.character(model) && identical(model, "als")) {
		cli::cli_abort(c(
			"{.code model = \"als\"} is not a valid model type.",
			"i" = "ALS is an estimator, not a model. Either call {.fn dbn_als} directly,",
			"i" = "or use {.code dbn(..., method = \"als\")} on top of a {.code model = \"dynamic\"} / {.val static} / {.val piecewise} specification."
		))
	}
	# inform the user when `model` defaults: the static default is a common
	# beginner trap on dynamic-shaped data, so make the choice visible.
	# gated on `verbose` because production pipelines that set verbose=FALSE
	# explicitly want this kind of chatter suppressed; the user still gets
	# what the static default produces, which is reflected in fit$model.
	model_was_default <- missing(model)
	model <- match.arg(model)
	if (model_was_default && isTRUE(verbose > 0)) {
		# nudge harder toward dynamic when there are several time periods
		Tt_data <- if (is.array(data)) tail(dim(data), 1L) else NA_integer_
		if (!is.na(Tt_data) && Tt_data >= 8L) {
			cli::cli_inform(c(
				"!" = "{.arg model} not specified; defaulting to {.val static} but the panel has {Tt_data} time periods.",
				"i" = "Strongly consider {.code model = \"dynamic\"} for panels with this many periods; static is rarely the right choice."
			))
		} else {
			cli::cli_inform(c(
				"i" = "{.arg model} not specified; defaulting to {.val static}.",
				"i" = "For time-varying operators (most common for panels with multiple periods), pass {.code model = \"dynamic\"}."
			))
		}
	}
	method <- match.arg(method)
	sampler <- match.arg(sampler, choices = c("auto", "exact", "approx"))
	# validate verbose upfront so character / NA inputs error directively.
	if (!is.logical(verbose) && !is.numeric(verbose)) {
		cli::cli_abort(c(
			"{.arg verbose} must be {.code TRUE}/{.code FALSE} or a non-negative integer.",
			"x" = "Got {.cls {class(verbose)[1]}} {.val {verbose}}."
		))
	}
	if (is.numeric(verbose) && (length(verbose) != 1L || !is.finite(verbose) || verbose < 0)) {
		cli::cli_abort("{.arg verbose} must be a single non-negative number, got {.val {verbose}}.")
	}

	# fast point-estimate path: `method = "als"` routes to `dbn_als()`, which
	# also handles `bootstrap = N` for CIs. MCMC-only args (nscan/burn/odens/
	# sampler) are silently ignored on this path; warn if the user passed
	# values that suggest they expected MCMC.
	if (method == "als") {
		mc <- match.call()
		mcmc_only <- intersect(names(mc), c("nscan", "burn", "odens", "sampler"))
		if (length(mcmc_only) > 0L) {
			cli::cli_warn(c(
				"{.code method = \"als\"} ignores MCMC-only argument{?s} {.arg {mcmc_only}}.",
				"i" = "{.fun dbn_als} controls fitting via {.arg max_iter}, {.arg tol}, {.arg ridge}; pass those via {.code ...}."
			))
		}
		dots <- list(...)
		base <- list(data = data, family = family, symmetric = symmetric,
			bootstrap = bootstrap, verbose = verbose)
		return(do.call(dbn_als, c(base, dots)))
	}

	if (bootstrap != 0) {
		cli::cli_warn(c(
			"{.arg bootstrap} is only meaningful with {.code method = \"als\"}; ignoring.",
			"i" = "For MCMC uncertainty use the posterior draws directly (the chain itself supplies intervals)."
		))
	}

	# catch brms-style names (n_iter, thin, prior, chains) directively and
	# point at the dbn equivalent; otherwise they fall into `...` and the
	# user gets a fit that didn't honor what they asked for.
	dots_dbn <- list(...)
	if ("n_iter" %in% names(dots_dbn)) {
		cli::cli_abort(c(
			"{.fun dbn} uses {.arg nscan} (not {.arg n_iter}) for the MCMC iteration count.",
			"i" = "Replace {.code n_iter = N} with {.code nscan = N}."
		))
	}
	if ("thin" %in% names(dots_dbn)) {
		cli::cli_abort(c(
			"{.fun dbn} uses {.arg odens} (not {.arg thin}) for output thinning.",
			"i" = "Replace {.code thin = N} with {.code odens = N} (keep every N-th iteration)."
		))
	}
	if ("prior" %in% names(dots_dbn)) {
		cli::cli_abort(c(
			"{.fun dbn} does not yet accept a {.arg prior} argument (brms-style).",
			"i" = "Variance hyperparameters are controlled via {.arg a_tau}, {.arg b_tau}, {.arg kappa_Abar2} on the symmetric path.",
			"i" = "A general prior interface is on the roadmap."
		))
	}
	if ("chains" %in% names(dots_dbn)) {
		cli::cli_abort(c(
			"{.fun dbn} runs a single Gibbs chain.",
			"i" = "For multi-chain Rhat-style convergence checks, use {.fun dbn_multichain}.",
			"i" = "Example: {.code dbn_multichain(Y, chains = 4, seeds = 1:4, ...)}."
		))
	}
	# brms-style covariate alias names: redirect to the canonical {.arg covariates}
	# argument so a user typing `X = ...` or `x_covariates = ...` does not
	# silently get a fit that ignores them
	cov_alias <- intersect(c("X", "x_covariates", "Z_covariates"), names(dots_dbn))
	if (length(cov_alias) > 0L) {
		cli::cli_abort(c(
			"{.fun dbn} accepts exogenous covariates via the {.arg covariates} argument.",
			"i" = "Replace {.code {cov_alias[1]} = ...} with {.code covariates = dbn_covariates(...)}; see {.fn dbn_covariates}."
		))
	}
	# pull covariates and actor_effects out of `...`; either non-NULL
	# (covariates) or non-"none" (actor_effects) triggers the covariate-aware
	# sampler. currently supported on {.code model = "dynamic"} only
	covariates    <- dots_dbn$covariates
	actor_effects <- dots_dbn$actor_effects %||% "none"
	uses_covariate_path <- !is.null(covariates) ||
		(is.character(actor_effects) && length(actor_effects) == 1L &&
		 actor_effects != "none")
	if (uses_covariate_path) {
		if (!is.null(covariates) && !inherits(covariates, "dbn_covariates"))
			cli::cli_abort(c(
				"{.arg covariates} must be a {.cls dbn_covariates} object.",
				"i" = "Build one with {.fn dbn_covariates}, e.g. {.code dbn_covariates(dyad = list(distance = D))}."
			))
		if (!model %in% c("dynamic", "static", "piecewise", "hmm", "lowrank"))
			cli::cli_abort(c(
				"Covariate / actor-effects support requires a recognised {.arg model}.",
				"i" = "Got {.val {model}}."
			))
		if (isTRUE(symmetric))
			cli::cli_abort("Covariate / actor-effects support for {.code symmetric = TRUE} is on the roadmap.")
	}
	dots_dbn$covariates <- NULL

	# validate sampler vs model
	if (sampler != "auto" && model != "dynamic") {
		cli::cli_warn(c(
			"Sampler {.val {sampler}} is only supported for {.code model = \"dynamic\"}.",
			"i" = "Switching to {.code sampler = \"auto\"} for {.val {model}} model."
		))
		sampler <- "auto"
	}

	# load data from file path or use the array directly. A length-1
	# character is a path; a character array is malformed data, not a path.
	if (is.character(data) && length(data) == 1L) {
		cli::cli_inform("Loading data from: {.path {data}}")
		env <- new.env()
		load(data, envir = env)
		Y <- env$Y
		if (is.null(Y)) cli::cli_abort("Data file must contain object {.var Y}")
	} else {
		if (inherits(data, "formula")) {
			cli::cli_abort(c(
				"{.fun dbn} does not accept a formula interface ({.code dbn(Y ~ x)}).",
				"x" = "Got a {.cls formula} for {.arg data}.",
				"i" = "Pass the network data directly as a 3D ({.code [actors, actors, time]}) or 4D ({.code [actors, actors, relations, time]}) numeric array.",
				"i" = "DBN models the network's own dynamics; exogenous covariates are not currently supported. Include actor-level covariates by augmenting the network array if needed."
			))
		}
		if (is.character(data)) {
			cli::cli_abort(c(
				"{.arg data} must be a numeric array or a single file path.",
				"x" = "Got a character array of length {length(data)}.",
				"i" = "Convert factor/character data to numeric codes before fitting."
			))
		}
		if (inherits(data, c("igraph", "network"))) {
			cli::cli_abort(c(
				"{.fun dbn} does not accept {.cls {class(data)[1]}} objects directly.",
				"i" = "Coerce to an adjacency array first: for one time slice, {.code array(as_adjacency_matrix(g, sparse = FALSE), c(n, n, 1))}; for a panel, stack along the time axis."
			))
		}
		# catch the two natural panel layouts (list of network/igraph
		# objects, list of plain matrices) before they fall through to
		# the generic "must be numeric array" error
		if (is.list(data) && length(data) >= 1L &&
		    all(vapply(data, function(x) inherits(x, c("network", "igraph")), logical(1L)))) {
			cli::cli_abort(c(
				"{.fun dbn} does not accept lists of {.cls network} / {.cls igraph} objects directly.",
				"i" = "Stack them into a 3D adjacency array first:",
				"i" = "  {.code Y <- simplify2array(lapply(net_list, as.matrix))}",
				"i" = "If actors enter/exit across periods, align node labels first and pass {.code NA} for absent ties."
			))
		}
		if (is.list(data) && length(data) >= 2L &&
		    all(vapply(data, function(x) is.matrix(x) && is.numeric(x), logical(1L)))) {
			cli::cli_abort(c(
				"{.fun dbn} does not accept a list of plain matrices.",
				"i" = "Stack them into a 3D adjacency array first:",
				"i" = "  {.code Y <- simplify2array(mat_list)}",
				"i" = "This gives an {.code [n_row, n_col, time]} array which {.fun dbn} accepts directly."
			))
		}
		Y <- data
	}
	####

	# validate dimensions and reshape 3D to 4D if needed
	if (length(dim(Y)) == 3) {
		dim_orig <- dim(Y)
		# heuristic: a [time, actor, actor] array (a common panel-data
		# mistake) has its last two dimensions equal and unequal to the
		# first; the expected layout is [actor, actor, time]
		if (dim_orig[2] == dim_orig[3] && dim_orig[1] != dim_orig[2]) {
			cli::cli_warn(c(
				"Input array dimensions {dim_orig[1]} x {dim_orig[2]} x {dim_orig[3]} look transposed.",
				"i" = "Expected {.code [actor, actor, time]}; this looks like {.code [time, actor, actor]}.",
				"i" = "If the first dimension is time, use {.code aperm(Y, c(2, 3, 1))} first."
			))
		}
		# preserve dimnames across the 3D->4D auto-conversion so downstream
		# tidy() / coef() output keeps the actor labels
		dnames_orig <- dimnames(Y)
		Y <- array(Y, dim = c(dim_orig[1], dim_orig[2], 1, dim_orig[3]))
		if (!is.null(dnames_orig)) {
			dimnames(Y) <- list(dnames_orig[[1]], dnames_orig[[2]], NULL, dnames_orig[[3]])
		}
		# silence the conversion note under verbose = FALSE
		if (isTRUE(verbose > 0)) cli::cli_inform("Converting 3D array to 4D array with single relation")
	} else if (length(dim(Y)) != 4) {
		cli::cli_abort("Data must be a 3D array [actors x actors x time] or 4D array [actors x actors x relations x time]")
	}

	# warn (do not fail) if a square network has a non-NA diagonal
	warn_filled_diagonal(Y)

	# MCMC controls must be single whole numbers
	for (ctrl_nm in c("nscan", "burn", "odens")) {
		ctrl_v <- get(ctrl_nm)
		if (length(ctrl_v) != 1L || !is.numeric(ctrl_v) || !is.finite(ctrl_v) ||
		    ctrl_v != round(ctrl_v)) {
			cli::cli_abort(c(
				"{.arg {ctrl_nm}} must be a single whole number.",
				"x" = "Got {.val {ctrl_v}}."
			))
		}
		assign(ctrl_nm, as.integer(round(ctrl_v)))
	}
	if (nscan <= 0) cli::cli_abort("{.arg nscan} must be positive.")
	if (burn < 0) cli::cli_abort("{.arg burn} must be non-negative.")
	if (odens < 1) cli::cli_abort("{.arg odens} must be at least 1.")
	if (burn >= nscan) cli::cli_abort("{.arg burn} ({burn}) must be less than {.arg nscan} ({nscan}).")
	if (odens > nscan) cli::cli_abort("{.arg odens} ({odens}) too large: no iterations would be saved ({.arg nscan} = {nscan}).")
	####

	# enforce static model for single time point
	Tt <- dim(Y)[4]
	if (Tt < 2) {
		if (model != "static") {
			cli::cli_abort(c(
				"Model {.val {model}} requires at least 2 time points.",
				"i" = "Your data has {Tt} time point{?s}.",
				"i" = "Use {.code model = \"static\"} for cross-sectional data."
			))
		}
		cli::cli_inform("Single time point detected -- using static model.")
	}

	# piecewise model requires sufficient time points
	if (model == "piecewise" && Tt < 4) {
		cli::cli_abort(c(
			"Piecewise model requires at least 4 time points.",
			"i" = "Your data has {Tt} time point{?s}.",
			"i" = "Use {.code model = \"static\"} for short time series."
		))
	}
	####

	# print dimension summary
	if (verbose) {
		n_row <- dim(Y)[1]
		n_col <- dim(Y)[2]
		cli::cli_h3("Data dimensions")
		if (n_row != n_col) {
			cli::cli_bullets(c(
				" " = "Senders: {n_row}",
				" " = "Receivers: {n_col}",
				" " = "Relations: {dim(Y)[3]}",
				" " = "Time points: {Tt}",
				"i" = "Bipartite network detected"
			))
		} else {
			cli::cli_bullets(c(
				" " = "Nodes: {n_row}",
				" " = "Relations: {dim(Y)[3]}",
				" " = "Time points: {Tt}"
			))
		}
	}
	####

	# validate symmetric constraint
	if (symmetric) {
		n_r <- dim(Y)[1]
		n_c <- dim(Y)[2]
		if (n_r != n_c) {
			cli::cli_abort(c(
				"Symmetric networks require equal sender and receiver dimensions.",
				"i" = "Your data has {n_r} senders and {n_c} receivers.",
				"i" = "Symmetric networks are not compatible with bipartite data."
			))
		}
		if (model == "lowrank") {
			cli::cli_abort(c(
				"Symmetric networks are not yet supported for low-rank models.",
				"i" = "Use {.code model = \"dynamic\"}, {.val static}, {.val piecewise}, or {.val hmm} for symmetric fits.",
				"i" = "The constraint is not yet implemented for the low-rank Tucker factorisation."
			))
		}
		# warn unconditionally (correctness warnings ignore verbosity) when
		# symmetric = TRUE is requested on visibly directed data, since the
		# directed information collapses onto the upper triangle. routed via
		# the dispatcher so every model (dynamic / piecewise / hmm / ...)
		# inherits the same check.
		asym_mag <- .dbn_asymmetry_magnitude(Y, n_r, n_c)
		if (!is.na(asym_mag) && asym_mag > 1e-6) {
			cli::cli_warn(c(
				"{.code symmetric = TRUE} was requested but the input data shows directed asymmetry (max |Y_t - Y_t^T| = {sprintf('%.3g', asym_mag)}).",
				"i" = "The symmetric specification collapses directed information onto the upper triangle of {.var A_t}.",
				"i" = "If this is intentional (e.g. data are pre-symmetrized), ignore. Otherwise refit with {.code symmetric = FALSE}."
			))
		}
		if (verbose) {
			cli::cli_inform(c("i" = "Symmetric constraint active: B will be set equal to A."))
		}
	}
	####

	# suggest lowrank for large networks
	if (model == "dynamic" && max(dim(Y)[1], dim(Y)[2]) > 80 && verbose) {
		n_suggest <- max(dim(Y)[1], dim(Y)[2])
		r_suggest <- min(max(2L, as.integer(ceiling(log2(n_suggest)))), n_suggest - 1L)
		cli::cli_inform(c(
			"i" = "Large network ({n_suggest} nodes). Consider {.code model = \"lowrank\"} for better scalability.",
			"i" = "Suggested rank: {.code r = {r_suggest}} (log2(n), increase if fit is poor)."
		))
	}
	####

	# dispatch to model-specific sampler
	dots <- list(...)

	# warn on ... arguments the chosen model does not accept (e.g. an HMM's
	# `R` passed to a lowrank model), which would otherwise be swallowed
	# silently and produce a different model than the user intended
	model_fn <- switch(model,
		static = "dbn_static", dynamic = "dbn_dynamic",
		lowrank = "dbn_lowrank", hmm = "dbn_hmm",
		piecewise = "dbn_piecewise", NULL)
	# when uses_covariate_path is TRUE, the actual callee is a covariate
	# dispatcher (e.g. .dbn_with_covariates for dynamic), not dbn_dynamic
	# itself. add the dispatcher's formals to `known` so covariate-specific
	# args (prior_beta_scale, tau_beta2_init, prior_kind, etc.) are not
	# flagged as unknown and dropped before dispatch.
	cov_fn <- if (isTRUE(uses_covariate_path)) switch(model,
		dynamic = ".dbn_with_covariates",
		static = ".dbn_with_covariates_static",
		lowrank = ".dbn_with_covariates_lowrank",
		hmm = ".dbn_with_covariates_hmm",
		piecewise = ".dbn_with_covariates_piecewise",
		NULL) else NULL
	if (length(dots) > 0 && !is.null(model_fn) && exists(model_fn, mode = "function")) {
		known <- names(formals(get(model_fn)))
		if (!is.null(cov_fn) && exists(cov_fn, mode = "function")) {
			known <- unique(c(known, names(formals(get(cov_fn)))))
		}
		# `blocks` is consumed by dbn() itself for the piecewise model;
		# `covariates`, `actor_effects`, and `time_varying_beta` are
		# consumed for the dynamic covariate path
		unknown <- setdiff(names(dots),
			c(known, "blocks", "covariates", "actor_effects",
				"time_varying_beta"))
		unknown <- unknown[nzchar(unknown)]
		if (length(unknown) > 0) {
			cli::cli_warn(c(
				"Ignoring {length(unknown)} argument{?s} not used by the {.val {model}} model: {.arg {unknown}}.",
				"i" = "Check for a typo or an argument meant for a different model."
			))
			# drop the unknown args before forwarding. Without this, partial
			# matching in the callee can resolve an unknown arg to a real
			# formal -- e.g. `ar1` partial-matches both `ar1_alpha` and `ar1_B`
			# in `dbn_lowrank`, crashing with
			# "argument N matches multiple formal arguments".
			dots[unknown] <- NULL
		}
	}

	.fit_t0 <- proc.time()[["elapsed"]]
	# use the filtered `dots` rather than `...` so dropped-unknown arguments
	# do not partial-match a formal in the callee (see comment above).
	.base_args <- list(Y = Y, family = family, nscan = nscan, burn = burn,
		odens = odens, verbose = verbose, symmetric = symmetric)
	# drop covariates from dots: it is passed explicitly to .dbn_with_covariates
	# in the dynamic branch below and is not a formal of any other model fn
	dots$covariates <- NULL
	results <- switch(model,
		static = if (uses_covariate_path) {
			do.call(.dbn_with_covariates_static, c(
				list(Y = Y, covariates = covariates),
				.base_args[!names(.base_args) %in% "Y"], dots
			))
		} else {
			do.call(dbn_static, c(.base_args, dots))
		},
		dynamic = if (uses_covariate_path) {
			do.call(.dbn_with_covariates, c(
				list(Y = Y, covariates = covariates),
				.base_args[!names(.base_args) %in% "Y"],
				dots
			))
		} else {
			do.call(dbn_dynamic, c(.base_args, list(sampler = sampler), dots))
		},
		lowrank = if (uses_covariate_path) {
			do.call(.dbn_with_covariates_lowrank, c(
				list(Y = Y, covariates = covariates),
				.base_args[!names(.base_args) %in% "Y"], dots
			))
		} else {
			do.call(dbn_lowrank, c(.base_args, dots))
		},
		hmm = if (uses_covariate_path) {
			do.call(.dbn_with_covariates_hmm, c(
				list(Y = Y, covariates = covariates),
				.base_args[!names(.base_args) %in% "Y"], dots
			))
		} else {
			do.call(dbn_hmm, c(.base_args, dots))
		},
		piecewise = if (uses_covariate_path) {
			# resolve blocks (including blocks = "auto") *before* dispatch so
			# the covariate driver receives an integer / numeric specification
			blocks_cov <- dots$blocks
			if (is.null(blocks_cov)) {
				cli::cli_abort(c(
					"Piecewise model requires {.arg blocks} parameter.",
					"i" = "Use {.code blocks = 4} for 4 equal blocks,",
					"i" = "or {.code blocks = c(25, 50, 75)} for custom boundaries,",
					"i" = "or {.code blocks = \"auto\"} for automatic selection."
				))
			}
			if (is.character(blocks_cov) && length(blocks_cov) == 1L && blocks_cov == "auto") {
				if (isTRUE(verbose > 0)) cli::cli_h2("Automatic Block Selection")
				auto_K_result <- select_K_auto(Y, family = family,
											   K_min = dots$K_min %||% 1L,
											   K_max = dots$K_max %||% NULL,
											   verbose = (verbose > 0))
				# select_K_auto returns length-1 boundaries for K=1 which
				# parse_blocks misreads as K = Tt; use the integer K directly
				blocks_cov <- if (auto_K_result$selected_K == 1L) 1L else
					auto_K_result$selected_boundaries
				Tt_auto <- dim(Y)[4]
				n_blocks_auto <- if (length(blocks_cov) == 1L) blocks_cov else
					length(unique(blocks_cov[blocks_cov >= 1 & blocks_cov <= Tt_auto]))
				if (n_blocks_auto >= Tt_auto - 1L) {
					fallback_k <- max(2L, min(4L, Tt_auto %/% 4L))
					cli::cli_warn(c(
						"Automatic block selection returned a degenerate segmentation ({n_blocks_auto} blocks for {Tt_auto} time points).",
						"i" = "Falling back to {fallback_k} equal blocks; set {.arg blocks} explicitly for a considered segmentation."
					))
					blocks_cov <- fallback_k
				}
			}
			do.call(.dbn_with_covariates_piecewise, c(
				list(Y = Y, covariates = covariates, blocks = blocks_cov),
				.base_args[!names(.base_args) %in% "Y"],
				dots[!names(dots) %in% c("blocks", "K_min", "K_max")]
			))
		} else {
			blocks <- dots$blocks
			if (is.null(blocks)) {
				cli::cli_abort(c(
					"Piecewise model requires {.arg blocks} parameter.",
					"i" = "Use {.code blocks = 4} for 4 equal blocks,",
					"i" = "or {.code blocks = c(25, 50, 75)} for custom boundaries,",
					"i" = "or {.code blocks = \"auto\"} for automatic selection."
				))
			}

			# handle automatic block selection
			auto_K_result <- NULL
			if (is.character(blocks) && blocks == "auto") {
				if (verbose) cli::cli_h2("Automatic Block Selection")
				auto_K_result <- select_K_auto(Y, family = family,
											   K_min = dots$K_min %||% 1L,
											   K_max = dots$K_max %||% NULL,
											   verbose = (verbose > 0))
				# select_K_auto returns length-1 boundaries for K=1 which
				# parse_blocks misreads as K = Tt; use the integer K directly
				blocks <- if (auto_K_result$selected_K == 1L) 1L else
					auto_K_result$selected_boundaries
				# guard against a degenerate auto-segmentation (one block per
				# time point is maximally overfit and not a useful default)
				Tt_auto <- dim(Y)[4]
				n_blocks_auto <- if (length(blocks) == 1L) blocks else
					length(unique(blocks[blocks >= 1 & blocks <= Tt_auto]))
				if (n_blocks_auto >= Tt_auto - 1L) {
					fallback_k <- max(2L, min(4L, Tt_auto %/% 4L))
					cli::cli_warn(c(
						"Automatic block selection returned a degenerate segmentation ({n_blocks_auto} blocks for {Tt_auto} time points).",
						"i" = "Falling back to {fallback_k} equal blocks; set {.arg blocks} explicitly for a considered segmentation."
					))
					blocks <- fallback_k
				}
			}

			# filter out blocks parameter from dots
			piecewise_dots <- dots[!names(dots) %in% c("blocks", "K_min", "K_max")]

			fit <- do.call(dbn_piecewise, c(
				list(Y = Y, family = family, blocks = blocks,
					 nscan = nscan, burn = burn, odens = odens,
					 verbose = verbose, symmetric = symmetric),
				piecewise_dots
			))

			if (!is.null(auto_K_result)) {
				fit$auto_K <- auto_K_result
			}
			fit
		}
	)

	results$model <- model

	# record wall-clock fit time
	runtime_sec <- round(proc.time()[["elapsed"]] - .fit_t0, 2)
	if (is.null(results$meta)) results$meta <- list()
	results$meta$runtime_sec <- runtime_sec

	# preserve covariate-fit subclass; the default "dbn" class is the base
	prior_class <- class(results)
	class(results) <- if ("dbn_covariates_fit" %in% prior_class) {
		c("dbn_covariates_fit", "dbn")
	} else {
		"dbn"
	}

	# flag posterior variance chains pinned at the safe-IG ceiling (1e8) and
	# locally explosive posterior-mean operators. These are post-fit
	# DIAGNOSTICS; failing to compute them shouldn't kill the fit the user
	# just paid for. But silent failure also hides genuine bugs in the
	# diagnostics themselves, so surface the error class via cli_warn.
	tryCatch(.warn_variance_saturation(results), error = function(e) {
		# different models store scalars in different slots
		var_hint <- if (identical(results$model, "static") ||
		                identical(results$model, "piecewise")) {
			"{.code fit$params}"
		} else {
			"{.code fit$sigma2}, {.code fit$tau_A2} etc."
		}
		cli::cli_warn(c(
			"Could not run the variance-saturation diagnostic on this fit.",
			"x" = "{conditionMessage(e)}",
			"i" = paste0("The fit itself is fine; you can inspect ", var_hint, " directly.")
		))
	})
	tryCatch(.warn_operator_stability(results), error = function(e) {
		cli::cli_warn(c(
			"Could not run the operator-stability diagnostic on this fit.",
			"x" = "{conditionMessage(e)}",
			"i" = "Check stability manually with {.fun dbn_operator} on {.code fit}."
		))
	})

	# mirror the symmetric flag to a top-level slot so downstream code
	# can branch on isTRUE(fit$symmetric) without reaching into fit$dims.
	if (is.null(results$symmetric)) {
		results$symmetric <- isTRUE(results$dims$is_symmetric)
	}

	# attach actor / time dimnames from the original Y onto the operator
	# arrays and the latent state so accessors return labeled output
	results <- .attach_actor_names(results)

	return(results)
	####
}

####
#' Static DBN MCMC
#'
#' @description Fits a DBN with operators that do not vary over time
#'   (\eqn{A_t = A}, \eqn{B_t = B}). Good for small panels, short
#'   time series, and as a fast baseline. For the broader model
#'   structure (the deviation-form bilinear AR, priors, and the
#'   identifiability anchor \eqn{A_0 = B_0 = I}), see
#'   [dbn-package].
#'
#' @details The static model has the Tucker decomposition
#'   \eqn{\Theta_t = M + A (\Theta_{t-1} - M) B^T + \varepsilon_t} with
#'   constant \eqn{A} and \eqn{B}. When \code{symmetric = TRUE} the
#'   sampler enforces \eqn{B = A} and a symmetric \eqn{M} by averaging
#'   the corresponding Tucker factor with its transpose after each
#'   update. Use [dbn_dynamic()] when sender or receiver influence
#'   evolves over time.
#'
#' @param Y Data array (nodes x nodes x relations x time)
#' @param family Character string specifying the data family/distribution:
#'   - "ordinal": Ordinal data (ordered categories). Data should be positive integers.
#'   - "gaussian": Continuous data with Gaussian errors. Data can be any real numbers.
#'   - "binary": Binary (0/1) data. Data should be 0/1 or logical values.
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param seed Random seed for reproducibility
#' @param verbose Logical or numeric. If TRUE the progress bar is updated on every iteration (its display rate is throttled by the terminal); if numeric `n`, the bar is updated every n-th iteration; FALSE suppresses all progress output.
#' @param previous Previous dbn_static results to continue from (optional)
#' @param init List of initial values for parameters: B, s2, t2, g2, M, Z (optional)
#' @param symmetric Logical. If TRUE, enforce a symmetric (undirected) operator: averages each Tucker factor with its transpose at every update and forces a symmetric baseline mean M. Requires a unipartite square network (n_row == n_col). Default: FALSE.
#' @param a_sig,b_sig Shape and rate of the inverse-gamma prior on the
#'   process variance (default 1, 1).
#' @param a_g,b_g Shape and rate of the inverse-gamma prior on the
#'   baseline-mean variance (default 0.5, 0.5).
#' @return A list with class `dbn`. The static fit stores its factor draws
#'   under a single name, `B`, as a length-3 list -- one Tucker factor per
#'   mode -- rather than as separate `A` (senders) and `B` (receivers) matrices:
#'   \describe{
#'     \item{`B[[1]]`}{Sender factor draws, dimension `[n_row, n_row, draws]`.
#'       This is the static analogue of `A` in the dynamic model; static fits
#'       have no top-level `fit$A` field.}
#'     \item{`B[[2]]`}{Receiver factor draws, dimension `[n_col, n_col, draws]`.}
#'     \item{`B[[3]]`}{Relation factor draws, dimension `[p, p, draws]`.}
#'   }
#'   This layout reflects the static model's symmetric Tucker decomposition
#'   (all three modes are treated uniformly). For the per-mode access pattern
#'   that the dynamic model exposes (`fit$A` for senders, `fit$B` for
#'   receivers), use `fit$B[[1]]` and `fit$B[[2]]` respectively.
#' @seealso \code{\link{dbn}} for the main dispatcher, \code{\link{param_summary}} for posterior summaries
#' @examples
#' \donttest{
#' sim <- simulate_static_dbn(n = 8, time = 5, seed = 1)
#' fit <- dbn_static(sim$Y, nscan = 200, burn = 100, verbose = FALSE)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_static <- function(Y, family = c("ordinal", "gaussian", "binary"),
					   nscan = 10000, burn = 1000, odens = 1,
					   seed = NULL, verbose = TRUE,
					   previous = NULL, init = NULL,
					   symmetric = FALSE,
					   a_sig = 1, b_sig = 1,
					   a_g = 0.5, b_g = 0.5) {

	# set up family and preprocess data
	family <- match.arg(family)
	# symmetric static: require unipartite (n_row == n_col), enforce symmetric
	# tucker factor B[[1]] and symmetric baseline mean M by symmetrizing after
	# each update. This is the static analog of `dbn_dynamic(symmetric = TRUE)`.
	if (isTRUE(symmetric) && length(dim(Y)) >= 2L && dim(Y)[1] != dim(Y)[2]) {
		cli::cli_abort(c(
			"{.code symmetric = TRUE} requires a unipartite (square) network.",
			"x" = "Got data with shape {.val {paste(dim(Y), collapse = 'x')}}.",
			"i" = "Symmetric static is undefined for bipartite networks."
		))
	}
	# validate verbose: accept logical or non-negative integer; reject strings.
	if (!is.logical(verbose) && !is.numeric(verbose)) {
		cli::cli_abort(c(
			"{.arg verbose} must be {.code TRUE}/{.code FALSE} or a non-negative integer.",
			"x" = "Got {.cls {class(verbose)[1]}} {.val {verbose}}."
		))
	}
	if (is.numeric(verbose) && (length(verbose) != 1L || !is.finite(verbose) || verbose < 0)) {
		cli::cli_abort("{.arg verbose} must be a single non-negative number, got {.val {verbose}}.")
	}
	if (isTRUE(verbose)) verbose <- 100L
	if (isFALSE(verbose)) verbose <- 0L
	FAM <- switch(family,
		ordinal = family_ordinal(),
		gaussian = family_gaussian(),
		binary = family_binary()
	)

	.dbn_restore_seed <- .use_seed_locally(seed)
	on.exit(.dbn_restore_seed(), add = TRUE)

	pre <- shared_preprocess(Y, family = family)
	Z <- pre$Z
	R <- pre$R
	M <- pre$M
	dims <- pre$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	m <- n_row
	p <- dims$p
	n <- Tt <- dims$Tt
	is_bipartite <- dims$is_bipartite
	nc <- n_row * n_col

	if (n_row == 1 || n_col == 1) {
		cli::cli_abort("Network must have at least 2 nodes on each side.")
	}

	is_large_network <- (max(n_row, n_col) > 15) || (p > 1) || (Tt > 20) || (nc * p * Tt > 10000)

	if (family == "ordinal") {
		R_flat <- array(R, c(n_row, n_col, p * Tt))
		IR <- precompute_rank_structure(R_flat, n_row, n_col, p, Tt)
	} else {
		IR <- pre$IR
	}
	K <- 3
	d <- c(n_row, n_col, p)
	####

	# initialize parameters (warm-start from previous fit if available).
	# refuse a cross-model warm-start directively: continuing a dynamic /
	# piecewise / hmm / lowrank fit under model = "static" would hand the
	# C++ sampler operator slots it cannot consume.
	if (!is.null(previous)) {
		if (!inherits(previous, "dbn")) {
			cli::cli_abort(c(
				"{.arg previous} must be a {.cls dbn} fit.",
				"x" = "Got {.cls {class(previous)[1]}}."
			))
		}
		prev_model <- previous$model %||% NA_character_
		if (!is.na(prev_model) && !identical(prev_model, "static")) {
			cli::cli_abort(c(
				"Cross-model warm-start is not supported.",
				"x" = "{.arg previous} was fit with {.code model = {.val {prev_model}}}, but you are now requesting {.code model = \"static\"}.",
				"i" = "Refit the new model from scratch, or continue with the same model type."
			))
		}
		s2 <- tail(previous$draws$pars$s2, 1)
		t2 <- tail(previous$draws$pars$t2, 1)
		g2 <- tail(previous$draws$pars$g2, 1)
		M <- previous$draws$misc$M[[length(previous$draws$misc$M)]]
		B <- lapply(previous$draws$misc$B, function(b) {
			mat <- b[, , dim(b)[3], drop = FALSE]
			dim(mat) <- dim(b)[1:2]
			mat
		})
	} else {
		# validate prior hyperparameters
		for (nm in c("a_sig", "b_sig", "a_g", "b_g")) {
			v <- get(nm)
			if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v <= 0)
				cli::cli_abort("{.arg {nm}} must be a positive finite scalar.")
		}
		# smart init via ALS warm-start when init = "smart" (default if NULL)
		init_mode <- if (is.character(init) && length(init) == 1L) init else "smart"
		if (identical(init_mode, "smart") || is.null(init)) {
			init_pkg <- .dbn_smart_init(Y, family = FAM$name, model = "static",
				symmetric = symmetric, Tt = Tt, n_row = n_row, n_col = n_col,
				p = p, verbose = (is.numeric(verbose) && verbose > 1))
			# pick a single A0, B0 from the time-tiled ALS cube
			A0 <- init_pkg$A_init[, , 1L]
			B0 <- init_pkg$B_init[, , 1L]
			B <- list(A0, B0, diag(p))
			s2 <- init_pkg$sigma2_init
			t2 <- init_pkg$tauA2_init
			g2 <- init_pkg$g2_init
			M  <- init_pkg$M_init
		} else {
			s2 <- 1
			t2 <- 1
			g2 <- max(0.1, mean(M^2, na.rm = TRUE))
			B <- list(diag(n_row), diag(n_col), diag(p))
		}
	}
	####

	# allocate MCMC storage
	n_iter <- burn + nscan
	keep_idx <- seq(burn + 1, n_iter, by = odens)
	n_keep <- length(keep_idx)

	B_samples <- list()
	for (k in 1:K) B_samples[[k]] <- array(NA, c(d[k], d[k], n_keep))
	param_samples <- matrix(NA, n_keep, 3)
	colnames(param_samples) <- c("s2", "t2", "g2")
	Msave <- vector("list", n_keep)
	if (FAM$name %in% c("ordinal", "binary")) {
		Zsave <- vector("list", n_keep)
	}

	Z_flat <- array(Z, c(n_row, n_col, p * Tt))

	if (isTRUE(verbose > 0)) {
		st_notes <- character()
		if (.dbn_data_looks_undirected(Y, n_row, n_col) && !isTRUE(symmetric)) {
			st_notes <- c(st_notes,
				"Data looks undirected (symmetric adjacency) but {.code symmetric = FALSE}; {.code symmetric = TRUE} is better identified for undirected networks.")
		}
		.dbn_preflight("static", n_row, n_col, p, Tt, burn, nscan, odens,
					   n_keep, notes = st_notes)
	}

	t_start <- proc.time()[[3]]
	if (verbose) {
		cli::cli_alert_info("Running static DBN MCMC ({n_iter} iterations)")
		cli::cli_progress_bar("MCMC iterations", total = n_iter)
	}

	if (is_large_network) {
		Z_cube <- reshape_Z_to_cube_parallel(Z, n_row, n_col, p, Tt)
	} else {
		Z_cube <- reshape_Z_to_cube(Z, n_row, n_col, p, Tt)
	}

	if (!is.array(M) || length(dim(M)) != 3) {
		M <- array(as.numeric(M), dim = c(n_row, n_col, p))
	}
	storage.mode(Z_cube) <- "double"
	storage.mode(M) <- "double"

	if (is_large_network && family == "ordinal") {
		Z_flat_mat <- matrix(0, nrow = nc, ncol = p * Tt)
		EZ_flat_mat <- matrix(0, nrow = nc, ncol = p * Tt)
	}

	# precompute ordinal flags/arrays once (R doesn't change)
	if (family == "ordinal") {
		static_use_approx <- should_use_gaussian_approximation(R) ||
							(nc * p * Tt > 500)
		static_R_flat <- array(R, c(n_row, n_col, p * Tt))
		static_has_rz_cpp <- exists("rz_gaussian_approx_cpp", mode = "function")
	}
	####

	# timing diagnostics (activated by verbose >= 200, e.g., verbose = 200L)
	static_do_timing <- is.numeric(verbose) && verbose >= 200L
	if (static_do_timing) {
		stiming <- list(b_update = 0, z_update = 0, m_update = 0, s2_update = 0, storage = 0)
		stiming_iters <- 0L
	}

	# main MCMC loop
	for (iter in 1:n_iter) {

		# update B matrices
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (is_large_network) {
			B[[1]] <- update_B_static_tiled(Z_cube, M, s2, t2, n_row, n_col, p, Tt)
		} else {
			B[[1]] <- update_B_static(Z_cube, M, s2, t2, n_row, n_col, p, Tt)
		}
		# symmetric static: enforce the symmetric tucker decomposition.
		# B[[1]] (sender) is made symmetric by averaging with its transpose;
		# B[[2]] (receiver) is tied to B[[1]] so the dispatcher promise
		# "B = A" actually holds at the factor level, not just on B[[1]].
		if (isTRUE(symmetric)) {
			B[[1]] <- 0.5 * (B[[1]] + t(B[[1]]))
			if (length(B) >= 2L) B[[2]] <- B[[1]]
		}

		# update t2 (B / Tucker-factor precision). prior is IG(0.5, 0.5);
		# we leave this hardcoded because t2 controls the operator factor
		# variance separately from the observation-error variance s2
		sse <- compute_diagonal_sse(B, K)
		t2 <- safe_rinv_gamma((sum(d) + 1) / 2, (sse + 1) / 2)
		if (static_do_timing) stiming$b_update <- stiming$b_update + (proc.time()[[3]] - .t0)
		####

		# update latent Z (ordinal/binary families)
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (FAM$name == "binary") {
			for (j in 1:p) {
				for (t in 1:Tt) {
					eta <- M[, , j]
					Y_jt <- R[, , j, t]
					pos <- which(Y_jt == 1)
					neg <- which(Y_jt == 0)
					if (length(pos) > 0) {
						Z[, , j, t][pos] <- truncnorm::rtruncnorm(
							length(pos), a = 0, b = Inf, mean = eta[pos], sd = 1)
					}
					if (length(neg) > 0) {
						Z[, , j, t][neg] <- truncnorm::rtruncnorm(
							length(neg), a = -Inf, b = 0, mean = eta[neg], sd = 1)
					}
				}
			}
			if (is_large_network) {
				Z_cube <- reshape_Z_to_cube_parallel(Z, n_row, n_col, p, Tt)
			} else {
				Z_cube <- reshape_Z_to_cube(Z, n_row, n_col, p, Tt)
			}
		} else if (FAM$name == "ordinal") {
			Z_flat <- array(Z, c(n_row, n_col, p * Tt))
			EZ_cube <- broadcast_M_and_compute_EZ(M, s2, n_row, n_col, p, Tt)

			if (static_use_approx) {
				if (static_has_rz_cpp) {
					Z_flat <- rz_gaussian_approx_cpp(static_R_flat, Z_flat, EZ_cube)
				} else {
					Z_flat <- rz_gaussian_approx(static_R_flat, Z_flat, EZ_cube)
				}
			} else {
				Z_flat <- rz_fc_batch(static_R_flat, Z_flat, EZ_cube, IR, n_row, n_col, p, Tt)
			}

			Z <- array(Z_flat, c(n_row, n_col, p, Tt))

			if (is_large_network) {
				Z_cube <- reshape_Z_to_cube_parallel(Z, n_row, n_col, p, Tt)
			} else {
				Z_cube <- reshape_Z_to_cube(Z, n_row, n_col, p, Tt)
			}
		}
		if (static_do_timing) stiming$z_update <- stiming$z_update + (proc.time()[[3]] - .t0)
		####

		# update baseline mean M every iteration: prior code only updated every
		# 5th iteration for non-ordinal as a perf optimization, but that left
		# `confint(fit)$M` with lo == hi == mean (each kept sample is a copy of
		# the most-recent update). Every-iter update restores posterior diversity.
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (FAM$name != "ordinal") {
			Z_flat <- array(Z, c(n_row, n_col, p * Tt))
		}
		Z_flat_mat <- matrix(Z_flat, nrow = nc, ncol = p * Tt)
		if (is_large_network) {
			M <- compute_M_static_blocked(Z_flat_mat, n_row, n_col, p, Tt)
		} else {
			M <- compute_M_static(Z_flat_mat, n_row, n_col, p, Tt)
		}
		# symmetric static: enforce M[,,r] symmetric per relation.
		if (isTRUE(symmetric)) {
			for (r in seq_len(p)) M[, , r] <- 0.5 * (M[, , r] + t(M[, , r]))
		}
		M_sum_sq <- sum(M^2, na.rm = TRUE)
		g2 <- (2 * b_g + M_sum_sq) /
			(2 * rgamma(1, shape = (2 * a_g + nc * p) / 2, rate = 1))
		if (static_do_timing) stiming$m_update <- stiming$m_update + (proc.time()[[3]] - .t0)
		####

		# update observation variance s2 (gaussian only)
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (FAM$name == "gaussian") {
			if (is_large_network) {
				rss <- compute_rss_static_parallel(Z_cube, M, n_row, n_col, p, Tt)
			} else {
				rss <- compute_rss_static(Z, M, n_row, n_col, p, Tt)
			}
			s2 <- safe_rinv_gamma(a_sig + nc * p * Tt / 2, b_sig + rss / 2)
		} else {
			s2 <- 1
		}
		if (static_do_timing) stiming$s2_update <- stiming$s2_update + (proc.time()[[3]] - .t0)
		####

		# store thinned samples
		if (static_do_timing) .t0 <- proc.time()[[3]]
		if (iter %in% keep_idx) {
			idx <- which(keep_idx == iter)
			B_samples[[1]][, , idx] <- B[[1]]
			if (K > 1) {
				if (length(B) >= 2 && !is.null(B[[2]])) {
					B_samples[[2]][, , idx] <- B[[2]]
				} else {
					B_samples[[2]][, , idx] <- diag(d[2])
				}
			}
			if (K > 2) {
				if (length(B) >= 3 && !is.null(B[[3]])) {
					B_samples[[3]][, , idx] <- B[[3]]
				} else {
					id_mat <- diag(d[3])
					if (d[3] == 1) {
						id_mat <- matrix(1, 1, 1)
					}
					B_samples[[3]][, , idx] <- id_mat
				}
			}
			param_samples[idx, ] <- c(s2, t2, g2)
			Msave[[idx]] <- M
			if (FAM$name %in% c("ordinal", "binary")) {
				Zsave[[idx]] <- Z
			}
		}
		if (static_do_timing) stiming$storage <- stiming$storage + (proc.time()[[3]] - .t0)
		####

		if (static_do_timing) {
			stiming_iters <- stiming_iters + 1L
			if (iter == 10L) {
				total <- Reduce(`+`, stiming)
				if (total > 0) {
					cli::cli_alert_info("Static MCMC timing (first 10 iterations, avg per iter):")
					for (nm in names(stiming)) {
						pct <- round(100 * stiming[[nm]] / total, 1)
						ms <- round(1000 * stiming[[nm]] / stiming_iters, 1)
						cli::cli_alert("{nm}: {ms}ms ({pct}%)")
					}
					cli::cli_alert("Total: {round(1000 * total / stiming_iters, 1)}ms/iter")
				}
			}
		}

		if (verbose) {
			if (iter %% verbose == 0 || iter == n_iter) {
				cli::cli_progress_update(set = iter)
			}
			.dbn_heartbeat(iter, n_iter, t_start, verbose)
		}
	}
	####

	if (verbose) cli::cli_progress_done()
	if (verbose) {
		.dbn_closing("static", t_start,
					 conv = .dbn_mixing_note(param_samples[, "s2"], "s2"))
	}

	# assemble output list
	draws <- list(
		theta = NULL,
		z = if (FAM$name %in% c("ordinal", "binary")) Zsave else NULL,
		pars = data.frame(
			s2 = param_samples[, "s2"],
			t2 = param_samples[, "t2"],
			g2 = param_samples[, "g2"]
		),
		misc = list(
			M = Msave,
			B = B_samples
		)
	)

	dic <- NA
	pd <- NA
	deviance_mean <- NA
	if (FAM$name == "gaussian") {
		dev_samples <- numeric(n_keep)
		for (idx in seq_len(n_keep)) {
			M_idx <- Msave[[idx]]
			s2_idx <- param_samples[idx, "s2"]
			dev <- 0
			for (j in 1:p) {
				for (t in 1:Tt) {
					resid <- Y[, , j, t] - M_idx[, , j]
					dev <- dev + sum(resid^2, na.rm = TRUE)
				}
			}
			n_obs <- sum(!is.na(Y))
			dev_samples[idx] <- n_obs * log(2 * pi * s2_idx) + dev / s2_idx
		}
		deviance_mean <- mean(dev_samples)
		pd <- var(dev_samples) / 2
		dic <- deviance_mean + pd
	}

	out <- list(
		model = "static",
		family = family,
		Y = Y,
		R = R,
		dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, n = n,
					 Tt = Tt, is_bipartite = is_bipartite, is_symmetric = symmetric),
		settings = list(
			nscan = nscan,
			burn = burn,
			odens = odens,
			draws = n_keep
		),
		meta = list(
			family = family,
			dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
						is_bipartite = is_bipartite, is_symmetric = symmetric),
			draws = n_keep,
			settings = list(nscan = nscan, burn = burn, odens = odens),
			Omega = pre$Omega
		),
		params = param_samples,
		M = M,
		B = B_samples,
		draws = draws,
		diagnostics = list(
			deviance = deviance_mean,
			pd = pd,
			dic = dic
		),
		sigma2_obs = if (FAM$name == "gaussian") param_samples[, "s2"] else NULL,
		Z_final = if (FAM$name %in% c("ordinal", "binary")) Z else NULL,
		M_final = M
	)

	if (!is.null(previous)) {
		prev_total <- previous$total_iter %||% (previous$settings$burn + previous$settings$nscan)
		out$total_iter <- prev_total + nscan
		out$continued_from <- prev_total
	}

	class(out) <- "dbn"
	return(out)
	####
}
####

####
#' Dynamic DBN MCMC
#'
#' @description Fits a DBN with time-varying operators \eqn{A_t}, \eqn{B_t}
#'   under a Gaussian random-walk prior anchored at the identity
#'   (\eqn{A_0 = B_0 = I}, see [dbn-package] for the full prior /
#'   identifiability story). Use [dbn_static()] for time-invariant
#'   operators, or [dbn_piecewise()] when the operator is constant
#'   within known regimes.
#'
#' @details The deviation-form bilinear AR is
#'   \eqn{\Theta_t = M + A_t (\Theta_{t-1} - M) B_t^T + \varepsilon_t},
#'   with \eqn{A_t = A_{t-1} + \varepsilon^A_t} and
#'   \eqn{\varepsilon^A_t \sim N(0, \tau_A^2 I)} (and analogously for
#'   \eqn{B_t}). The \eqn{A_0 = B_0 = I} anchor resolves the
#'   \eqn{(A_t, B_t) \to (c A_t, c^{-1} B_t)} scale ambiguity; under
#'   \code{symmetric = TRUE} the sampler enforces \eqn{B_t = A_t} and a
#'   symmetric \eqn{A_t} at every iteration. See [actor_embedding()] and
#'   [compute_irf()] for downstream summaries.
#'
#' @param Y Data array (nodes x nodes x relations x time)
#' @param family Character string specifying the data family/distribution:
#'   - "ordinal": Ordinal data (ordered categories). Data should be positive integers.
#'   - "gaussian": Continuous data with Gaussian errors. Data can be any real numbers.
#'   - "binary": Binary (0/1) data. Data should be 0/1 or logical values.
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param ar1 Use AR(1) dynamics instead of random walk (default: FALSE)
#' @param update_rho Update AR coefficients (default: FALSE)
#' @param seed Random seed
#' @param verbose Logical or numeric. If TRUE the progress bar is updated on every iteration (its display rate is throttled by the terminal); if numeric `n`, the bar is updated every n-th iteration; FALSE suppresses all progress output.
#' @param time_thin Save every nth time point to reduce memory (default: NULL = auto)
#' @param store_z Whether to store Z draws (default: NULL = auto based on memory)
#' @param previous Previous dbn_dynamic results to continue from (optional)
#' @param init List of initial values: A, B, sigma2, tau_A2, tau_B2, g2, rho_A, rho_B, Theta, M, Z (optional)
#' @param symmetric Logical. If TRUE, enforce `B_t = A_t` and use the
#'   symmetric-specification sampler (matrix-free PCG on the symmetric
#'   zero-diagonal latent state, per-entry M-H on the symmetric operator).
#'   Requires unipartite square data (`n_row == n_col`, `p == 1`). Default:
#'   FALSE.
#' @param lambda_diag Diagonal-penalty weight for the symmetric specification.
#'   `0` (default) disables the diagonal penalty; positive values activate
#'   the per-entry M-H sampler with prior `A_t = Abar + Delta_t` and
#'   stationary AR(1) on `Delta_t`. Only used when `symmetric = TRUE`.
#' @param phi_diag Persistence parameter on `Delta_t` when `lambda_diag > 0`
#'   (default: 0.999, near random-walk).
#' @param kappa_Abar2 Per-entry prior variance on `Abar` when
#'   `lambda_diag > 0` (default: 0.5).
#' @param tau_A_fixed If non-NULL, hold `tau_A^2` fixed at this value and
#'   skip the conjugate IG update. Honored on both the symmetric and the
#'   asymmetric specifications.
#' @param tau_B_fixed If non-NULL, hold `tau_B^2` fixed at this value and
#'   skip the conjugate IG update. Asymmetric specification only; ignored
#'   when `symmetric = TRUE` (where `tau_B^2 = tau_A^2` by construction).
#' @param a_tau,b_tau Shape and rate of the inverse-gamma prior on the
#'   random-walk innovation variances `tau_A^2` (and `tau_B^2` in the
#'   asymmetric specification), default 0.5 each, a weakly informative prior.
#'   Exposed for prior-sensitivity analysis; applies to both the symmetric and
#'   asymmetric specifications.
#' @param a_sig,b_sig Shape and rate of the inverse-gamma prior on the
#'   process variance `sigma^2` (default 2, 2). Exposed for prior-sensitivity
#'   analysis.
#' @param a_sig_obs,b_sig_obs Shape and rate of the inverse-gamma prior on
#'   the observation variance `sigma^2_obs` (Gaussian family only). `NULL`
#'   (default) reuses the process-variance prior `a_sig`, `b_sig`.
#' @param a_g,b_g Shape and rate of the inverse-gamma prior on the baseline
#'   mean variance `g^2` (default 2, 2).
#' @param sampler Sampler choice. `"auto"` (default) picks the best
#'   sampler for the model; `"approx"` forces the FFBS path; `"exact"`
#'   forces the exact PCG path where available.
#' @param kappa_A2 If non-NULL, decouple the row-1 prior variance from
#'   `tau_A^2` by setting each row's `A_1` prior variance to this value.
#'   NULL (default) ties the prior to `tau_A^2`.
#' @param rho_max Optional cap on the per-time spectral radius of `A_t`
#'   via row-trajectory rejection sampling in the suffstat path. NULL
#'   (default) disables.
#' @param rho_max_rejects Per-row resampling attempts before leaving the
#'   row unchanged for this sweep (default: 50). Only used when
#'   `rho_max > 0`.
#' @param shrink_rho Contractivity prior on the operator, on by default. The
#'   per-time `A_t`, `B_t` are under-identified (one transition informs each
#'   slice), and the random-walk prior does not bound their spectral radius,
#'   so impulse responses and forecasts can drift into the explosive region.
#'   `shrink_rho` adds a Gaussian ridge prior that pulls each operator slice
#'   toward zero, keeping the spectral radius contractive and multi-step
#'   propagation stable. It defaults to `0.9`; smaller values shrink harder,
#'   and `shrink_rho = NULL` disables it entirely (recovering the raw
#'   random-walk prior). The realized spectral radius sits at or (for short
#'   panels) well below `shrink_rho`, since the random-walk coupling shares
#'   the same precision. This regularizes the magnitude of the operator; to
#'   recover its structure, pool transitions with `model = "piecewise"` or
#'   `model = "lowrank"`.
#' @param keep Optional character vector of draw-indexed components to retain
#'   in the returned object, controlling its size on disk and in memory.
#'   Recognized: `"Theta"` (latent-state draws -- by far the largest
#'   component), `"A"`, `"B"`, `"M"`, `"Z"`. `NULL` (default) keeps all.
#'   Coupling, leverage, and rank-probability analysis need only `A`, `B`,
#'   and `M`; passing `keep = c("A", "B", "M")` drops the latent-state draws
#'   and substantially shrinks the saved object. See [estimate_memory()].
#' @return List containing MCMC results
#' @seealso \code{\link{dbn}} for the main dispatcher, \code{\link{param_summary}} for posterior summaries
#' @examples
#' \donttest{
#' sim <- simulate_dynamic_dbn(n = 6, time = 5, seed = 1)
#' fit <- dbn_dynamic(sim$Y, nscan = 200, burn = 100, verbose = FALSE)
#' }
#' @author Tosin Salau and Shahryar Minhas
#' @export
dbn_dynamic <- function(Y,
						family = c("ordinal", "gaussian", "binary"),
						nscan = 10000,
						burn = 1000,
						odens = 1,
						ar1 = FALSE,
						update_rho = FALSE,
						seed = NULL,
						verbose = TRUE,
						time_thin = NULL,
						store_z = NULL,
						previous = NULL,
						init = NULL,
						symmetric = FALSE,
						lambda_diag = 0,
						phi_diag = 0.999,
						kappa_Abar2 = 0.5,
						tau_A_fixed = NULL,
						tau_B_fixed = NULL,
						a_tau = 0.5,
						b_tau = 0.5,
						a_sig = 2,
						b_sig = 2,
						a_sig_obs = NULL,
						b_sig_obs = NULL,
						a_g = 2,
						b_g = 2,
						kappa_A2 = NULL,
						rho_max = NULL,
						rho_max_rejects = 50L,
						shrink_rho = 0.9,
						keep = NULL,
						sampler = "auto") {

	# set up family and preprocess data
	.dbn_restore_seed <- .use_seed_locally(seed)
	on.exit(.dbn_restore_seed(), add = TRUE)
	family <- match.arg(family)
	# validate verbose: accept logical or non-negative integer; reject strings.
	if (!is.logical(verbose) && !is.numeric(verbose)) {
		cli::cli_abort(c(
			"{.arg verbose} must be {.code TRUE}/{.code FALSE} or a non-negative integer.",
			"x" = "Got {.cls {class(verbose)[1]}} {.val {verbose}}."
		))
	}
	if (is.numeric(verbose) && (length(verbose) != 1L || !is.finite(verbose) || verbose < 0)) {
		cli::cli_abort("{.arg verbose} must be a single non-negative number, got {.val {verbose}}.")
	}
	if (isTRUE(verbose)) verbose <- 100L
	if (isFALSE(verbose)) verbose <- 0L
	# validate ar1 / update_rho coupling. ar1 must be a single logical;
	# numeric values are rejected so a typo like 0.8 cannot silently coerce
	# to TRUE. ar1 = TRUE with update_rho = FALSE pins rho at zero, which is
	# the random-walk default; auto-promote update_rho when ar1 = TRUE so the
	# caller's request is honoured (explicit update_rho = FALSE opts out).
	if (!is.logical(ar1) || length(ar1) != 1L || is.na(ar1)) {
		cli::cli_abort(c(
			"{.arg ar1} must be a single logical (TRUE or FALSE).",
			"x" = "Got {.val {ar1}}.",
			"i" = "For a fixed persistence value, fit with {.code ar1 = TRUE} and inspect {.code fit$rhoA}."
		))
	}
	if (ar1 && !update_rho) {
		ar1_call <- match.call()
		# only promote when the user did not pass `update_rho` explicitly
		if (!"update_rho" %in% names(ar1_call)) {
			update_rho <- TRUE
			if (verbose > 0) cli::cli_inform(c(
				"i" = "{.code ar1 = TRUE} requires {.code update_rho = TRUE} to actually estimate persistence -- promoting automatically.",
				"i" = "Pass {.code update_rho = FALSE} explicitly to pin rho at 0 (equivalent to {.code ar1 = FALSE})."
			))
		} else {
			cli::cli_warn(c(
				"{.code ar1 = TRUE} with {.code update_rho = FALSE} pins rhoA = rhoB = 0; this is identical to the random-walk default ({.code ar1 = FALSE}).",
				"i" = "Drop {.code update_rho = FALSE} (or set it to {.val TRUE}) to actually estimate AR(1) persistence."
			))
		}
	}
	# capture the user-facing `keep` argument under a distinct name: the
	# sampler reuses the bare name `keep` internally for the vector of
	# retained MCMC iteration indices.
	keep_components <- NULL
	if (!is.null(keep)) {
		keep_components <- as.character(keep)
		valid_keep <- c("Theta", "A", "B", "M", "Z", "draws")
		bad <- setdiff(keep_components, valid_keep)
		if (length(bad) > 0) {
			cli::cli_abort(c(
				"Unrecognized {.arg keep} value{?s}: {.val {bad}}.",
				"i" = "Valid components are {.val {valid_keep}}."
			))
		}
	}
	FAM <- switch(family,
		ordinal = family_ordinal(),
		gaussian = family_gaussian(),
		binary = family_binary()
	)

	# validate and resolve sampler parameter
	sampler <- match.arg(sampler, choices = c("auto", "exact", "approx"))
	if (sampler == "auto") {
		sampler_resolved <- if (symmetric) "exact" else "approx"
	} else {
		sampler_resolved <- sampler
		if (sampler == "approx" && symmetric && verbose > 0) {
			cli::cli_inform("i" = "Sampler {.val \"approx\"} (symmetric FFBS) not yet implemented; using exact PCG.")
			sampler_resolved <- "exact"
		}
		# exact PCG is implemented only for the symmetric specification;
		# fall back rather than crash in the asymmetric exact path
		if (sampler == "exact" && !symmetric) {
			cli::cli_warn(c(
				"Exact PCG is not yet implemented for asymmetric (directed) models.",
				"i" = "Falling back to the approximate FFBS sampler.",
				"i" = "Use {.code symmetric = TRUE} for exact PCG, or {.code sampler = \"auto\"}."
			))
			sampler_resolved <- "approx"
		}
	}

	if (dim(Y)[4] < 2) {
		cli::cli_abort("Dynamic model requires at least 2 time points. Use {.code model = \"static\"} for cross-sectional data.")
	}
	# reject single-actor "networks": the operator A_t collapses to a
	# 1x1 scalar with no edge structure to estimate
	if (dim(Y)[1] < 2 || dim(Y)[2] < 2) {
		cli::cli_abort(c(
			"Dynamic model requires at least 2 actors on each side.",
			"x" = "Got data with dim {.val {paste(dim(Y), collapse = 'x')}}.",
			"i" = "A single-actor network has no edge structure; the bilinear operator is a 1x1 scalar."
		))
	}

	# the IG(a_tau, b_tau) prior on the random-walk innovation variance now
	# applies to both the symmetric and asymmetric specifications
	if (length(a_tau) != 1 || !is.numeric(a_tau) || !is.finite(a_tau) ||
		a_tau <= 0 || length(b_tau) != 1 || !is.numeric(b_tau) ||
		!is.finite(b_tau) || b_tau <= 0) {
		cli::cli_abort("{.arg a_tau} and {.arg b_tau} must be positive finite scalars.")
	}

	# small binary networks can be numerically unstable
	if (family == "binary" && min(dim(Y)[1], dim(Y)[2]) < 15) {
		cli::cli_warn(c(
			"Dynamic binary models with small networks (n < 15) may encounter numerical singularities.",
			"i" = "Consider using {.code model = \"static\"} or a larger network.",
			"i" = "The model will attempt to run, but may produce unreliable results."
		))
	}

	pre <- shared_preprocess(Y, family = family)
	Z <- pre$Z
	R <- pre$R
	IR <- pre$IR
	M <- pre$M
	dims <- pre$dims
	n_row <- dims$n_row
	n_col <- dims$n_col
	m <- n_row
	p <- dims$p
	Tt <- dims$Tt
	is_bipartite <- dims$is_bipartite
	nc <- n_row * n_col
	d <- nc
	####

	# symmetric-specification validation and knob bookkeeping. when
	# symmetric = TRUE the package activates the unipartite-square sampler
	# (matrix-free PCG on the symmetric, zero-diagonal latent state and
	# per-entry M-H on the symmetric operator).
	# tau_A_fixed / tau_B_fixed apply to both the symmetric and asymmetric
	# paths; validate up front and honor them at the tau update site below.
	if (!is.null(tau_A_fixed)) {
		if (length(tau_A_fixed) != 1 || !is.numeric(tau_A_fixed) ||
			!is.finite(tau_A_fixed) || tau_A_fixed <= 0) {
			cli::cli_abort("{.arg tau_A_fixed} must be a single positive finite value.")
		}
	}
	if (!is.null(tau_B_fixed)) {
		if (length(tau_B_fixed) != 1 || !is.numeric(tau_B_fixed) ||
			!is.finite(tau_B_fixed) || tau_B_fixed <= 0) {
			cli::cli_abort("{.arg tau_B_fixed} must be a single positive finite value.")
		}
		if (isTRUE(symmetric)) {
			cli::cli_warn(c(
				"{.arg tau_B_fixed} is ignored when {.code symmetric = TRUE}.",
				"i" = "Under the symmetric specification {.code tau_B^2 = tau_A^2}; pass {.arg tau_A_fixed} instead."
			))
		}
	}
	if (symmetric) {
		if (n_row != n_col || p != 1) {
			cli::cli_abort("symmetric = TRUE requires unipartite square data (n_row == n_col, p == 1).")
		}
		# directed-asymmetry warning has been moved up into dbn() so it
		# fires for every model (piecewise / hmm / static), not just dynamic.
		# the dispatcher always sees Y first.
		if (length(lambda_diag) != 1 || !is.numeric(lambda_diag) ||
			!is.finite(lambda_diag) || lambda_diag < 0) {
			cli::cli_abort("{.arg lambda_diag} must be a non-negative finite scalar.")
		}
		if (lambda_diag > 0) {
			if (length(phi_diag) != 1 || !is.numeric(phi_diag) ||
				!is.finite(phi_diag) || abs(phi_diag) >= 1) {
				cli::cli_abort("{.arg phi_diag} must be a finite scalar in (-1, 1).")
			}
			if (length(kappa_Abar2) != 1 || !is.numeric(kappa_Abar2) ||
				!is.finite(kappa_Abar2) || kappa_Abar2 <= 0) {
				cli::cli_abort("{.arg kappa_Abar2} must be a positive finite scalar.")
			}
		}
		if (!is.null(kappa_A2)) {
			if (length(kappa_A2) != 1 || !is.numeric(kappa_A2) ||
				!is.finite(kappa_A2) || kappa_A2 <= 0) {
				cli::cli_abort("{.arg kappa_A2} must be a positive finite scalar (or NULL).")
			}
		}
		if (!is.null(rho_max)) {
			if (length(rho_max) != 1 || !is.numeric(rho_max) ||
				!is.finite(rho_max) || rho_max <= 0) {
				cli::cli_abort("{.arg rho_max} must be a positive finite scalar (or NULL).")
			}
		}
	}
	# shrink_rho governs the directed (asymmetric) operator only; validate it
	# on every path. under symmetric = TRUE it is a no-op (that path uses
	# kappa_A2).
	if (!is.null(shrink_rho)) {
		if (length(shrink_rho) != 1 || !is.numeric(shrink_rho) ||
			!is.finite(shrink_rho) || shrink_rho <= 0 || shrink_rho > 1) {
			cli::cli_abort("{.arg shrink_rho} must be a scalar in (0, 1] (or NULL).")
		}
	}
	####

	# configure time-thinning and estimate memory usage
	if (is.null(time_thin)) {
		time_thin <- max(1L, Tt %/% 20L)
		if (time_thin > 1 && isTRUE(verbose > 0)) {
			cli::cli_inform(c(
				"i" = "Auto time-thinning: storing every {time_thin}th time point ({length(seq(1, Tt, by = time_thin))} of {Tt}).",
				"i" = "Override with {.code time_thin = 1} to store all."
			))
		}
	}

	if (is.null(store_z)) {
		n_keep_est <- floor(nscan / odens)
		n_time_est <- length(seq(1, Tt, by = time_thin))
		z_gb <- nc * p * n_time_est * n_keep_est * 8 / 1024^3
		store_z <- (z_gb < 2)
		if (!store_z && family %in% c("ordinal", "binary") && isTRUE(verbose > 0)) {
			cli::cli_inform(c(
				"i" = "Z array would be {round(z_gb, 1)} GB. Skipping Z storage to save memory.",
				"i" = "Override with {.code store_z = TRUE} to force storage."
			))
		}
	}

	preflight_est <- c(object_gb = NA_real_, peak_gb = NA_real_)
	if (isTRUE(verbose > 0)) {
		preflight_est <- estimate_memory(n_row, n_col, p, Tt, nscan, burn,
										 odens, time_thin, family,
										 store_z = isTRUE(store_z),
										 keep = keep_components, quiet = TRUE)
		n_keep_pf <- length(seq.int(burn + 1, burn + nscan, by = odens))
		pf_notes <- character()
		if (.dbn_data_looks_undirected(Y, n_row, n_col) && !isTRUE(symmetric)) {
			pf_notes <- c(pf_notes,
				"Data looks undirected (symmetric adjacency) but {.code symmetric = FALSE}. The directed model is supported and its scale drift is auto-corrected, but {.code symmetric = TRUE} is faster and better identified for undirected networks.")
		}
		if (preflight_est[["object_gb"]] > 4 && is.null(keep_components)) {
			pf_notes <- c(pf_notes, paste0(
				"Large output (~", round(preflight_est[["object_gb"]], 1),
				" GB). Pass {.code keep = c(\"A\", \"B\", \"M\")} to drop the Theta draws, or raise {.code time_thin} / {.code odens}."))
		}
		.dbn_preflight("dynamic", n_row, n_col, p, Tt, burn, nscan, odens,
					   n_keep_pf, preflight_est[["object_gb"]],
					   preflight_est[["peak_gb"]], pf_notes)
	}

	is_large_network <- (max(n_row, n_col) > 15) || (p > 1) || (Tt > 20) || (nc * p * Tt > 10000)

	if (family == "ordinal" && !is.null(IR)) {
		IR_time_indices <- precompute_time_indices(IR, n_row, n_col, p, Tt)
	}
	####

	# flatten 4D arrays to matrices for vectorized operations
	Z_4d <- matrix(Z, nrow = nc, ncol = p * Tt)
	R_4d <- matrix(R, nrow = nc, ncol = p * Tt)
	####

	# initialize parameters (warm-start from previous fit if available).
	# refuse a cross-model warm-start directively.
	if (!is.null(previous)) {
		if (!inherits(previous, "dbn")) {
			cli::cli_abort(c(
				"{.arg previous} must be a {.cls dbn} fit.",
				"x" = "Got {.cls {class(previous)[1]}}."
			))
		}
		prev_model <- previous$model %||% NA_character_
		if (!is.na(prev_model) && !identical(prev_model, "dynamic")) {
			cli::cli_abort(c(
				"Cross-model warm-start is not supported.",
				"x" = "{.arg previous} was fit with {.code model = {.val {prev_model}}}, but you are now requesting {.code model = \"dynamic\"}.",
				"i" = "Refit the new model from scratch, or continue with the same model type."
			))
		}
		if (is.null(previous$A) || is.null(previous$B) || is.null(previous$sigma2)) {
			cli::cli_abort("{.arg previous} must be results from {.fun dbn_dynamic}.")
		}
		last_idx <- length(previous$sigma2)

		if (!is.null(previous$Theta)) {
			n_prev_iter <- dim(previous$Theta)[5]
			last_theta_idx <- n_prev_iter
			prev_time_thin <- previous$time_thin %||% 1
			n_time_stored <- dim(previous$Theta)[4]

			if (n_time_stored < Tt) {
				time_indices <- seq(1, Tt, by = prev_time_thin)
				Theta_all <- array(0, dim = c(n_row, n_col, p, Tt))
				for (i in 1:min(length(time_indices), n_time_stored)) {
					if (time_indices[i] <= Tt) {
						slice <- previous$Theta[,,,i,last_theta_idx, drop = FALSE]
						dim(slice) <- c(n_row, n_col, p)
						Theta_all[,,,time_indices[i]] <- slice
					}
				}
				for (t in 1:Tt) {
					if (all(Theta_all[,,,t] == 0)) {
						available_times <- time_indices[time_indices <= Tt]
						nearest_idx <- which.min(abs(available_times - t))
						nearest_stored <- nearest_idx
						if (nearest_stored <= n_time_stored) {
							slice <- previous$Theta[,,,nearest_stored,last_theta_idx, drop = FALSE]
							dim(slice) <- c(n_row, n_col, p)
							Theta_all[,,,t] <- slice
						}
					}
				}
			} else {
				Theta_all <- previous$Theta[,,,1:Tt,last_theta_idx, drop = FALSE]
				dim(Theta_all) <- c(n_row, n_col, p, Tt)
			}
		} else {
			Theta_all <- pre$Theta
		}

		sigma2 <- previous$sigma2[last_idx]
		sigma2_obs <- if (!is.null(previous$sigma2_obs)) previous$sigma2_obs[last_idx] else 1
		tauA2 <- previous$tau_A2[last_idx]
		tauB2 <- previous$tau_B2[last_idx]
		g2 <- previous$g2[last_idx]
		rhoA <- if (!is.null(previous$rhoA)) previous$rhoA[last_idx] else 0
		rhoB <- if (!is.null(previous$rhoB)) previous$rhoB[last_idx] else 0

		last_A <- previous$A[[last_idx]]
		last_B <- previous$B[[last_idx]]

		if (dim(last_A)[3] < Tt) {
			prev_time_thin <- previous$time_thin %||% 1
			time_indices <- seq(1, Tt, by = prev_time_thin)
			Aarray <- array(0, dim = c(n_row, n_row, Tt))
			Barray <- array(0, dim = c(n_col, n_col, Tt))
			for (i in 1:length(time_indices)) {
				if (time_indices[i] <= Tt && i <= dim(last_A)[3]) {
					Aarray[,,time_indices[i]] <- last_A[,,i]
					Barray[,,time_indices[i]] <- last_B[,,i]
				}
			}
			for (t in 1:Tt) {
				if (all(Aarray[,,t] == 0)) {
					available_times <- time_indices[time_indices <= Tt]
					nearest_idx <- which.min(abs(available_times - t))
					nearest_time <- available_times[nearest_idx]
					nearest_stored <- which(time_indices == nearest_time)
					Aarray[,,t] <- last_A[,,nearest_stored]
					Barray[,,t] <- last_B[,,nearest_stored]
				}
			}
		} else {
			Aarray <- last_A
			Barray <- last_B
		}
	} else {
		Theta_all <- pre$Theta
		# smart init: a fast ALS warm-start gives initial A_t / B_t / M /
		# variance components close to the posterior mode. defaults to
		# enabled; falls back to identity init if ALS fails or `init = "default"`
		init_mode <- if (is.character(init) && length(init) == 1L) init else "smart"
		if (identical(init_mode, "smart") || is.null(init)) {
			init_pkg <- .dbn_smart_init(Y, family = FAM$name, model = "dynamic",
				symmetric = symmetric, Tt = Tt, n_row = n_row, n_col = n_col,
				p = p, verbose = (is.numeric(verbose) && verbose > 1))
			Aarray <- init_pkg$A_init
			Barray <- init_pkg$B_init
			sigma2     <- init_pkg$sigma2_init
			sigma2_obs <- FAM$init_pars$sigma2_obs %||% init_pkg$sigma2_obs_init
			tauA2      <- init_pkg$tauA2_init
			tauB2      <- init_pkg$tauB2_init
			g2         <- init_pkg$g2_init
		} else {
			Aarray <- array(0, dim = c(n_row, n_row, Tt))
			Barray <- array(0, dim = c(n_col, n_col, Tt))
			for (t in 1:Tt) {
				Aarray[, , t] <- diag(n_row)
				Barray[, , t] <- diag(n_col)
			}
			sigma2 <- 1
			sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
			tauA2 <- tauB2 <- 1
			g2 <- 1
		}
		rhoA <- rhoB <- 0
	}

	if (symmetric && !is.null(tau_A_fixed)) tauA2 <- tau_A_fixed

	# symmetric-path state: per-entry M-H proposal scale (adapted in R) and
	# the operator mean Abar (lambda_diag > 0 only). Warm-start from
	# `previous` if available; otherwise initialize Abar to zero and the
	# proposal sd to a modest starting scale.
	if (symmetric) {
		if (!is.null(previous) && !is.null(previous$mh_proposal_sd) &&
			is.finite(previous$mh_proposal_sd) && previous$mh_proposal_sd > 0) {
			mh_proposal_sd <- previous$mh_proposal_sd
		} else {
			mh_proposal_sd <- 0.05
		}
		if (lambda_diag > 0) {
			if (!is.null(previous) && !is.null(previous$Abar) &&
				length(previous$Abar) > 0) {
				Abar_state <- previous$Abar[[length(previous$Abar)]]
				if (!all(dim(Abar_state) == c(n_row, n_row))) {
					Abar_state <- matrix(0, n_row, n_row)
				}
			} else {
				Abar_state <- matrix(0, n_row, n_row)
			}
		} else {
			Abar_state <- NULL
		}
		# enforce symmetric initial A and B = A
		for (t in seq_len(Tt)) {
			At <- 0.5 * (Aarray[, , t] + t(Aarray[, , t]))
			Aarray[, , t] <- At
			Barray[, , t] <- At
		}
		# 0-based upper-triangle pairs for the PCG sampler
		pairs_mat <- matrix(0L, nrow = n_row * (n_row - 1) / 2, ncol = 2)
		d_idx <- 0L
		for (i_p in 0L:(n_row - 2L)) {
			for (j_p in (i_p + 1L):(n_row - 1L)) {
				d_idx <- d_idx + 1L
				pairs_mat[d_idx, ] <- c(i_p, j_p)
			}
		}
		# acceptance history for adaptive proposal scaling during burn
		mh_accept_history <- numeric(0)
	}

	Theta_4d <- matrix(Theta_all, nrow = nc, ncol = p * Tt)
	####

	# allocate MCMC storage
	n_iter <- burn + nscan
	keep <- seq(burn + 1, n_iter, by = odens)
	n_keep <- length(keep)
	time_keep <- seq(1, Tt, by = time_thin)
	n_time_keep <- length(time_keep)

	Theta_store <- array(NA, dim = c(n_row, n_col, p, n_time_keep, n_keep))
	if (FAM$name %in% c("ordinal", "binary") && store_z) {
		Z_store <- array(NA, dim = c(n_row, n_col, p, n_time_keep, n_keep))
	}

	A_store <- B_store <- vector("list", n_keep)
	M_store <- array(NA, dim = c(n_row, n_col, p, n_keep))
	Abar_store <- if (symmetric && lambda_diag > 0) vector("list", n_keep) else NULL

	sigma2_store <- numeric(n_keep)
	sigma2_obs_store <- if (FAM$name == "gaussian") numeric(n_keep) else NULL
	tau_A2_store <- tau_B2_store <- g2_store <- numeric(n_keep)
	rhoA_store <- rhoB_store <- if (ar1) numeric(n_keep) else NULL

	keep_id <- 0
	# user-supplied IG hyperparams; default sig_obs to sig if not specified
	for (nm in c("a_sig", "b_sig", "a_g", "b_g")) {
		v <- get(nm)
		if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v <= 0)
			cli::cli_abort("{.arg {nm}} must be a positive finite scalar.")
	}
	if (is.null(a_sig_obs)) a_sig_obs <- a_sig
	if (is.null(b_sig_obs)) b_sig_obs <- b_sig
	for (nm in c("a_sig_obs", "b_sig_obs")) {
		v <- get(nm)
		if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v <= 0)
			cli::cli_abort("{.arg {nm}} must be a positive finite scalar.")
	}
	eye_nr <- diag(n_row)
	eye_nc <- diag(n_col)

	if (is_large_network) {
		# IG(a_tau, b_tau) prior; default a_tau = 0.5 gives the 2*0.5 = 1 term
		shape_tauA <- (2 * a_tau + n_row * n_row * (Tt - 1)) / 2
		shape_tauB <- (2 * a_tau + n_col * n_col * (Tt - 1)) / 2
		shape_sigma_proc <- (a_sig + nc * (Tt - 1) * p) / 2.0
		shape_sigma_obs <- (1.0 + nc * Tt * p) / 2.0
	}

	use_approx <- FALSE
	if (FAM$name == "ordinal") {
		use_approx <- should_use_gaussian_approximation(R) ||
					 (nc * p * Tt > 5000)
		if (use_approx) {
			EZ_cube <- array(0, c(n_row, n_col, p * Tt))
			Z_cube  <- array(0, c(n_row, n_col, p * Tt))
			R_cube  <- array(R_4d, c(n_row, n_col, p * Tt))
			M_expanded <- rep(as.vector(M), times = Tt)
		}
	}

	# precompute function availability flags (avoid per-iteration exists() calls)
	has_rz_gaussian_approx_cpp <- exists("rz_gaussian_approx_cpp", mode = "function")
	has_IR_time_indices <- exists("IR_time_indices")
	has_shape_sigma_proc <- exists("shape_sigma_proc")

	t_start <- proc.time()[[3]]
	if (verbose) {
		cli::cli_alert_info("Running dynamic DBN MCMC ({n_iter} iterations)")
		cli::cli_progress_bar("MCMC iterations", total = n_iter)
	}
	####

	# timing diagnostics (activated by verbose >= 200, e.g., verbose = 200L)
	do_timing <- is.numeric(verbose) && verbose >= 200L
	if (do_timing) {
		timing_accum <- list(z_update = 0, mu_update = 0, ffbs = 0, ab_update = 0,
							 tau_update = 0, var_update = 0, rho_update = 0, storage = 0)
		timing_iters <- 0L
	}

	# contractivity ridge precision on A_t, B_t. per-entry prior variance
	# shrink_rho^2 / n gives a prior operator spectral radius near shrink_rho
	# (circular law); pass the reciprocal to c++. null shrink_rho = off.
	kappaA_inv <- if (is.null(shrink_rho)) 0 else n_row / (shrink_rho^2)
	kappaB_inv <- if (is.null(shrink_rho)) 0 else n_col / (shrink_rho^2)

	# main MCMC loop
	for (g in 1:n_iter) {

		# update latent Z (ordinal/binary families)
		if (do_timing) .t0 <- proc.time()[[3]]
		if (FAM$name == "ordinal") {
			if (use_approx) {
				EZ_cube[] <- Theta_4d + M_expanded
				Z_cube[]  <- Z_4d
				R_cube[]  <- R_4d

				if (has_rz_gaussian_approx_cpp) {
					Z_cube <- rz_gaussian_approx_cpp(R_cube, Z_cube, EZ_cube)
				} else {
					Z_cube <- rz_gaussian_approx(R_cube, Z_cube, EZ_cube)
				}
				Z_4d <- matrix(Z_cube, nrow = nc, ncol = p * Tt)
			} else {
				if (has_IR_time_indices) {
					Z_4d <- batch_update_Z_ordinal_fast(R_4d, Z_4d, Theta_4d, M, IR, IR_time_indices, n_row, n_col, p, Tt)
				} else {
					Z_4d <- batch_update_Z_ordinal(R_4d, Z_4d, Theta_4d, M, IR, n_row, n_col, p, Tt)
				}
			}
			if (g %in% keep) {
				Z <- array(Z_4d, dim = c(n_row, n_col, p, Tt))
			}
		} else if (FAM$name == "binary") {
			Z <- update_Z_optimized(R, Z, Theta_all, M, IR = NULL, family = "binary")
			Z_4d <- matrix(Z, nrow = nc, ncol = p * Tt)
		}
		if (do_timing) timing_accum$z_update <- timing_accum$z_update + (proc.time()[[3]] - .t0)
		####

		# update baseline mean M and g2
		if (do_timing) .t0 <- proc.time()[[3]]
		mu_result <- update_mu_dynamic(Z_4d, Theta_4d, M, g2, a_g, b_g, n_row, n_col, p, Tt)
		M <- mu_result$M
		g2 <- mu_result$g2
		# under symmetric specification, project M onto the symmetric
		# subspace; no-op when the data are pre-symmetrized
		if (symmetric && p == 1L) {
			Mtmp <- M[, , 1]
			M[, , 1] <- 0.5 * (Mtmp + t(Mtmp))
		}
		if (use_approx) M_expanded <- rep(as.vector(M), times = Tt)
		if (do_timing) timing_accum$mu_update <- timing_accum$mu_update + (proc.time()[[3]] - .t0)
		####

		if (symmetric) {
			# symmetric path: exact PCG or approximate FFBS on the symmetric, zero-diag
			# latent state; per-entry M-H on the symmetric operator
			# (random-walk or AR(1)-around-Abar with the diagonal penalty);
			# triangle-only RSS for sigma^2 and sigma^2_obs; conjugate IG
			# for tau_A^2; Gibbs draw of Abar when lambda_diag > 0.

			# 1. Theta: exact PCG or approximate FFBS given current A
			if (do_timing) .t0 <- proc.time()[[3]]
			Y_3d <- array(Z_4d[, 1:Tt, drop = FALSE], dim = c(n_row, n_col, Tt))
			M_2d <- matrix(M[, , 1], n_row, n_col)
			if (sampler_resolved == "exact") {
				# exact PCG sampler (default for symmetric)
				Theta_3d <- theta_pcg_sampler_cpp(Y_3d, M_2d, Aarray, sigma2, sigma2_obs, pairs_mat)
			} else {
				# approximate FFBS sampler (symmetric FFBS with B_arr = A_arr)
				Barray_tmp <- Aarray  # Set B = A for symmetric case
				Theta_cube_ffbs <- batch_ffbs_all_relations(Z_4d, M, Aarray, Barray_tmp, sigma2, n_row, n_col, p, Tt)
				# theta_cube_ffbs is n_row*n_col x p*Tt, extract first relation (p=1)
				Theta_3d <- array(Theta_cube_ffbs, dim = c(n_row, n_col, Tt))
			}
			Theta_all <- array(0, dim = c(n_row, n_col, p, Tt))
			Theta_all[, , 1, ] <- Theta_3d
			Theta_4d <- matrix(Theta_all, nrow = nc, ncol = p * Tt)
			Theta_cube <- Theta_all
			if (do_timing) {
				if (sampler_resolved == "exact") {
					timing_accum$pcg <- timing_accum$pcg + (proc.time()[[3]] - .t0)
				} else {
					timing_accum$ffbs <- timing_accum$ffbs + (proc.time()[[3]] - .t0)
				}
			}

			# 2. A update. lambda_diag = 0 -> sufficient-statistic row-FFBS
			# (matrix-free O(n^4 T) Gibbs scan). lambda_diag > 0 -> per-entry
			# M-H scan with diagonal-penalty likelihood; adapt mh_proposal_sd
			# in R from the running acceptance rate over the burn window.
			if (do_timing) .t0 <- proc.time()[[3]]
			if (lambda_diag > 0) {
				sigma2_diag <- sigma2 / lambda_diag
				mh_result <- update_A_ar1_mh_diag_cpp(
					Theta_3d, Aarray, Abar_state, phi_diag,
					sigma2, sigma2_diag, tauA2,
					mh_proposal_sd,
					n_sweeps = 1L, burn = 0L, odens = 1L, adapt = FALSE
				)
				Aarray <- mh_result$A_final
				sweep_ar <- as.numeric(mh_result$accept_rates)[1]
				mh_accept_history <- c(mh_accept_history, sweep_ar)
				if (g <= burn && length(mh_accept_history) >= 25 && (g %% 25) == 0) {
					recent_ar <- mean(tail(mh_accept_history, 25))
					if (recent_ar < 0.2) mh_proposal_sd <- mh_proposal_sd * 0.85
					else if (recent_ar > 0.5) mh_proposal_sd <- mh_proposal_sd * 1.15
				}
			} else {
				tau2_init_per_row_arg <- if (is.null(kappa_A2)) numeric(0) else rep(kappa_A2, n_row)
				rho_max_arg <- if (is.null(rho_max)) 0 else rho_max
				suff_result <- update_A_suffstat_cpp(
					Theta_3d, Aarray, sigma2, tauA2,
					tau2_init_per_row = tau2_init_per_row_arg,
					rho_max = rho_max_arg,
					max_rejects = rho_max_rejects
				)
				Aarray <- suff_result$Aarray
			}
			Barray <- Aarray
			if (do_timing) timing_accum$ab_update <- timing_accum$ab_update + (proc.time()[[3]] - .t0)

			# 3. tau_A^2 via conjugate IG over unique upper-triangle
			# innovations (random walk when lambda_diag = 0; AR(1) around
			# abar when lambda_diag > 0). skip when tau_A_fixed is set.
			if (do_timing) .t0 <- proc.time()[[3]]
			if (is.null(tau_A_fixed)) {
				innov_ss_sym <- 0
				tri_mask <- upper.tri(matrix(0, n_row, n_row), diag = TRUE)
				if (Tt >= 2) {
					for (t in 2:Tt) {
						if (lambda_diag > 0) {
							innov_mat <- Aarray[, , t] -
								(1 - phi_diag) * Abar_state -
								phi_diag * Aarray[, , t - 1]
						} else {
							innov_mat <- Aarray[, , t] - Aarray[, , t - 1]
						}
						innov_ss_sym <- innov_ss_sym + sum(innov_mat[tri_mask]^2)
					}
				}
				n_unique_per_t <- n_row * (n_row + 1L) / 2L
				# IG(a_tau, b_tau) prior; written so the default a_tau = b_tau =
				# 0.5 reproduces the original (1 + .)/2 expressions exactly.
				shape_post_tauA <- (2 * a_tau + n_unique_per_t * (Tt - 1)) / 2
				tauA2 <- safe_rinv_gamma(shape_post_tauA, (2 * b_tau + innov_ss_sym) / 2)
			}
			tauB2 <- tauA2
			if (do_timing) timing_accum$tau_update <- timing_accum$tau_update + (proc.time()[[3]] - .t0)

			# 4. Abar Gibbs draw (lambda_diag > 0 only)
			if (lambda_diag > 0) {
				Abar_state <- sample_Abar_internal(Aarray, phi_diag, tauA2, kappa_Abar2)
			}

			# 5. sigma^2 and sigma^2_obs (triangle-only RSS via
			# `exclude_diagonal = TRUE`)
			if (do_timing) .t0 <- proc.time()[[3]]
			var_result <- update_variances_dynamic(
				Theta_4d, Z_4d, M, Aarray, Barray,
				a_sig, b_sig, n_row, n_col, p, Tt,
				is_gaussian = (FAM$name == "gaussian"),
				exclude_diagonal = TRUE
			)
			sigma2 <- var_result$sigma2
			if (FAM$name == "gaussian") {
				sigma2_obs <- var_result$sigma2_obs
			}
			if (do_timing) timing_accum$var_update <- timing_accum$var_update + (proc.time()[[3]] - .t0)
		} else {
			# asymmetric (general directed) path: exact PCG or approximate FFBS for Theta
			# and row-wise updates for A and B.

			# theta sampling: exact PCG or approximate FFBS
			if (do_timing) .t0 <- proc.time()[[3]]
			if (sampler_resolved == "exact") {
				# exact asymmetric PCG is on the roadmap but not yet implemented;
				# the symmetric path has the exact PCG sampler today.
				cli::cli_abort(c(
					"{.code sampler = \"exact\"} is not yet available for asymmetric dynamic fits.",
					"i" = "Use {.code sampler = \"approx\"} (default for asymmetric), or set {.code symmetric = TRUE} to get the exact PCG sampler."
				))
			} else {
				# approximate FFBS sampler
				if (is_large_network && max(n_row, n_col) > 100) {
					Theta_cube <- batch_ffbs_all_relations_blocked(Z_4d, M, Aarray, Barray, sigma2, n_row, n_col, p, Tt)
				} else {
					Theta_cube <- batch_ffbs_all_relations(Z_4d, M, Aarray, Barray, sigma2, n_row, n_col, p, Tt)
				}
				if (g %in% keep) {
					Theta_all <- array(Theta_cube, dim = c(n_row, n_col, p, Tt))
				}
				if (do_timing) timing_accum$ffbs <- timing_accum$ffbs + (proc.time()[[3]] - .t0)
			}
			Theta_4d <- matrix(Theta_cube, nrow = nc, ncol = p * Tt)

			# prior_kind = "rw" is the random-walk prior anchored at
			# A_0 = B_0 = I (default). pass options(dbn.prior_kind = "iid")
			# for the iid alternative. validate so a typo errors loudly
			# instead of silently falling through.
			prior_kind <- getOption("dbn.prior_kind", "rw")
			if (!identical(prior_kind, "rw") && !identical(prior_kind, "iid")) {
				cli::cli_abort(c(
					"{.code getOption(\"dbn.prior_kind\")} must be {.val rw} or {.val iid}, got {.val {prior_kind}}.",
					"i" = "Reset with {.code options(dbn.prior_kind = \"rw\")} (default) or {.val iid}."
				))
			}
			if (do_timing) .t0 <- proc.time()[[3]]
			if (is_large_network && max(n_row, n_col) > 100) {
				AB_result <- update_AB_batch_large(
					Theta_4d, Aarray, Barray,
					sigma2, tauA2, tauB2,
					ar1, rhoA, rhoB,
					n_row, n_col, p, Tt,
					prior_kind = prior_kind,
					kappaA_inv = kappaA_inv, kappaB_inv = kappaB_inv
				)
			} else {
				AB_result <- update_AB_batch_extended(
					Theta_4d, Aarray, Barray,
					sigma2, tauA2, tauB2,
					ar1, rhoA, rhoB,
					n_row, n_col, p, Tt,
					prior_kind = prior_kind,
					kappaA_inv = kappaA_inv, kappaB_inv = kappaB_inv
				)
			}
			Aarray <- AB_result$Aarray
			Barray <- AB_result$Barray
			if (do_timing) timing_accum$ab_update <- timing_accum$ab_update + (proc.time()[[3]] - .t0)

			# scale renormalization: the asymmetric model has a (A_t, B_t) ->
			# (c A_t, c^-1 B_t) scale ambiguity. under the rw prior (default),
			# the A_0 = B_0 = I anchor pins this scale (the prior penalty on
			# A_1 - I is not invariant under rescaling), so the rebalancer is
			# skipped. under the iid prior the anchor carries no probabilistic
			# weight, so the rebalancer runs to prevent overflow.
			if (!identical(prior_kind, "rw")) {
				if (all(is.finite(Aarray)) && all(is.finite(Barray))) {
					rmsA <- sqrt(mean(Aarray^2))
					rmsB <- sqrt(mean(Barray^2))
					if (is.finite(rmsA) && is.finite(rmsB) && rmsA > 0 && rmsB > 0) {
						c_bal  <- sqrt(rmsB / rmsA)
						Aarray <- Aarray * c_bal
						Barray <- Barray / c_bal
					}
				}
			}

			# divergence guard: the asymmetric specification has a scale
			# ambiguity ((A_t, B_t) and (c A_t, c^-1 B_t) give the same
			# fit), so the MCMC can drift along the scale orbit until the
			# influence matrices overflow. Catch it here and fail with an
			# informative error rather than letting non-finite values
			# reach the next iteration's compiled samplers.
			if (!all(is.finite(Aarray)) || !all(is.finite(Barray))) {
				cli::cli_abort(c(
					"Dynamic sampler diverged at iteration {g} of {n_iter}: the influence matrices contain non-finite values.",
					"i" = "The asymmetric specification ({.code symmetric = FALSE}) has a scale ambiguity: {.code (A_t, B_t)} and {.code (c A_t, c^-1 B_t)} produce the same fit, so the sampler can drift until {.var A_t} or {.var B_t} overflows.",
					"*" = "For undirected or symmetric networks, refit with {.code symmetric = TRUE} (recommended).",
					"*" = "Otherwise, shorten the chain, rescale the input data toward unit variance, or report this issue."
				))
			}

			# update innovation variances tauA2, tauB2
			if (do_timing) .t0 <- proc.time()[[3]]
			# IG(a_tau, b_tau) prior on tauA2/tauB2; the default a_tau = b_tau
			# = 0.5 reproduces the original (1 + .)/2 expressions exactly
			if (ar1) {
				innovA_ss <- compute_ar1_innovation_ss_cpp(Aarray, rhoA, n_row, Tt)
				innovB_ss <- compute_ar1_innovation_ss_cpp(Barray, rhoB, n_col, Tt)
				if (is_large_network) {
					tauA2 <- safe_rinv_gamma(shape_tauA, (2 * b_tau + innovA_ss)/2)
					tauB2 <- safe_rinv_gamma(shape_tauB, (2 * b_tau + innovB_ss)/2)
				} else {
					tauA2 <- safe_rinv_gamma((2 * a_tau + n_row * n_row * (Tt-1))/2, (2 * b_tau + innovA_ss)/2)
					tauB2 <- safe_rinv_gamma((2 * a_tau + n_col * n_col * (Tt-1))/2, (2 * b_tau + innovB_ss)/2)
				}
			} else {
				# tau update sufficient statistic depends on the prior:
				#   "rw"  : sum_{t=1..T-1} ||A_t - A_{t-1}||_F^2 (rw innovations).
				#           reuse compute_ar1_innovation_ss_cpp with rho = 1.0;
				#           the formula then becomes ||A_t - A_{t-1}||^2.
				#   "iid" : sum_{t=1..T-1} ||A_t||_F^2 (iid sufficient statistic).
				if (identical(prior_kind, "rw")) {
					A_sum <- compute_ar1_innovation_ss_cpp(Aarray, 1.0, n_row, Tt)
					B_sum <- compute_ar1_innovation_ss_cpp(Barray, 1.0, n_col, Tt)
				} else {
					A_sum <- compute_deviation_sum(Aarray, n_row, Tt)
					B_sum <- compute_deviation_sum(Barray, n_col, Tt)
				}
				if (is_large_network) {
					tauA2 <- safe_rinv_gamma(shape_tauA, (2 * b_tau + A_sum)/2)
					tauB2 <- safe_rinv_gamma(shape_tauB, (2 * b_tau + B_sum)/2)
				} else {
					tauA2 <- safe_rinv_gamma((2 * a_tau + n_row * n_row * (Tt-1))/2, (2 * b_tau + A_sum)/2)
					tauB2 <- safe_rinv_gamma((2 * a_tau + n_col * n_col * (Tt-1))/2, (2 * b_tau + B_sum)/2)
				}
			}
			# honor tau_A_fixed / tau_B_fixed on the asymmetric and AR(1)
			# branches. the symmetric path applies tau_A_fixed before the
			# loop and forces tauB2 = tauA2 inside its own tau branch.
			if (!isTRUE(symmetric)) {
				if (!is.null(tau_A_fixed)) tauA2 <- tau_A_fixed
				if (!is.null(tau_B_fixed)) tauB2 <- tau_B_fixed
			}
			if (do_timing) timing_accum$tau_update <- timing_accum$tau_update + (proc.time()[[3]] - .t0)

			# update process and observation variances
			if (do_timing) .t0 <- proc.time()[[3]]
			if (is_large_network && max(n_row, n_col) > 100) {
				proc_rss <- compute_process_variance_blocked(Theta_4d, Aarray, Barray, n_row, n_col, p, Tt)
				if (has_shape_sigma_proc) {
					sigma2 <- (b_sig + proc_rss / 2.0) / rgamma(1, shape = shape_sigma_proc, rate = 1)
				} else {
					sigma2 <- (b_sig + proc_rss / 2.0) / rgamma(1, shape = (a_sig + nc * (Tt - 1) * p) / 2.0, rate = 1)
				}

				if (FAM$name == "gaussian") {
					obs_rss <- compute_gaussian_obs_residuals_dynamic_cpp(Z_4d, Theta_4d, M, n_row, n_col, p, Tt)
					sigma2_obs <- (b_sig_obs + obs_rss / 2) /
						rgamma(1, shape = (a_sig_obs + nc * Tt * p) / 2, rate = 1)
				}
			} else {
				var_result <- update_variances_dynamic(
					Theta_4d, Z_4d, M, Aarray, Barray,
					a_sig, b_sig, n_row, n_col, p, Tt,
					is_gaussian = (FAM$name == "gaussian")
				)
				sigma2 <- var_result$sigma2
				if (FAM$name == "gaussian") {
					sigma2_obs <- var_result$sigma2_obs
				}
			}
			if (do_timing) timing_accum$var_update <- timing_accum$var_update + (proc.time()[[3]] - .t0)
		}
		####

		# update AR(1) coefficients rhoA, rhoB
		if (do_timing) .t0 <- proc.time()[[3]]
		if (ar1 && update_rho) {
			rhoA_result <- compute_rho_update_cpp(Aarray, n_row, Tt)
			rho_mean <- rhoA_result$num / (rhoA_result$denom + 1e-10)
			rho_var  <- tauA2 / (rhoA_result$denom + 1e-10)
			rhoA <- truncnorm::rtruncnorm(1, a = -0.99, b = 0.99, mean = rho_mean, sd = sqrt(rho_var))

			rhoB_result <- compute_rho_update_cpp(Barray, n_col, Tt)
			rho_mean <- rhoB_result$num / (rhoB_result$denom + 1e-10)
			rho_var  <- tauB2 / (rhoB_result$denom + 1e-10)
			rhoB <- truncnorm::rtruncnorm(1, a = -0.99, b = 0.99, mean = rho_mean, sd = sqrt(rho_var))
			if (symmetric) rhoB <- rhoA
		}
		if (do_timing) timing_accum$rho_update <- timing_accum$rho_update + (proc.time()[[3]] - .t0)
		####

		# store thinned samples
		if (do_timing) .t0 <- proc.time()[[3]]
		if (g %in% keep) {
			keep_id <- keep_id + 1

			if (!exists("Theta_all") || !is.array(Theta_all) || length(dim(Theta_all)) != 4) {
				Theta_all <- array(Theta_cube, dim = c(n_row, n_col, p, Tt))
			}
			Theta_store[,,,,keep_id] <- Theta_all[,,,time_keep, drop = FALSE]

			if (FAM$name %in% c("ordinal", "binary") && store_z) {
				if (!is.array(Z) || length(dim(Z)) != 4) {
					Z <- array(Z_4d, dim = c(n_row, n_col, p, Tt))
				}
				Z_store[,,,,keep_id] <- Z[,,,time_keep, drop = FALSE]
			}

			if (symmetric) {
				# project each retained A_t onto the symmetric subspace so
				# fit$A draws are exactly transpose-symmetric. internal
				# sampler state is unaffected.
				A_thin <- Aarray[, , time_keep, drop = FALSE]
				for (tt in seq_len(dim(A_thin)[3])) {
					Atmp <- A_thin[, , tt]
					A_thin[, , tt] <- 0.5 * (Atmp + t(Atmp))
				}
				A_store[[keep_id]] <- A_thin
				B_store[[keep_id]] <- A_thin
				if (lambda_diag > 0) Abar_store[[keep_id]] <- Abar_state
			} else {
				A_store[[keep_id]] <- Aarray[, , time_keep, drop = FALSE]
				B_store[[keep_id]] <- Barray[, , time_keep, drop = FALSE]
			}

			M_store[,,,keep_id] <- M
			sigma2_store[keep_id] <- sigma2
			if (FAM$name == "gaussian") {
				sigma2_obs_store[keep_id] <- sigma2_obs
			}

			tau_A2_store[keep_id] <- tauA2
			tau_B2_store[keep_id] <- tauB2
			g2_store[keep_id] <- g2

			if (ar1) {
				rhoA_store[keep_id] <- rhoA
				rhoB_store[keep_id] <- rhoB
			}
		}
		if (do_timing) timing_accum$storage <- timing_accum$storage + (proc.time()[[3]] - .t0)
		####

		if (do_timing) {
			timing_iters <- timing_iters + 1L
			if (g == 10L) {
				total <- Reduce(`+`, timing_accum)
				if (total > 0) {
					cli::cli_alert_info("MCMC timing (first 10 iterations, avg per iter):")
					for (nm in names(timing_accum)) {
						pct <- round(100 * timing_accum[[nm]] / total, 1)
						ms <- round(1000 * timing_accum[[nm]] / timing_iters, 1)
						cli::cli_alert("{nm}: {ms}ms ({pct}%)")
					}
					cli::cli_alert("Total: {round(1000 * total / timing_iters, 1)}ms/iter")
				}
			}
		}

		if (verbose) {
			if (g %% verbose == 0 || g == n_iter) {
				cli::cli_progress_update(set = g)
			}
			.dbn_heartbeat(g, n_iter, t_start, verbose)
		}
	}
	####

	if (verbose) cli::cli_progress_done()
	if (verbose) {
		.dbn_closing("dynamic", t_start,
					 object_gb = preflight_est[["object_gb"]],
					 conv = .dbn_mixing_note(sigma2_store, "sigma2"))
	}

	# convert arrays to lists for unified draws structure
	Theta_list <- lapply(seq_len(n_keep), function(i) Theta_store[, , , , i, drop = FALSE])
	for (i in seq_along(Theta_list)) {
		dim(Theta_list[[i]]) <- dim(Theta_store)[1:4]
	}
	M_list <- lapply(seq_len(n_keep), function(i) M_store[, , , i, drop = FALSE])
	for (i in seq_along(M_list)) {
		dim(M_list[[i]]) <- dim(M_store)[1:3]
	}

	# build unified draws structure (matching static/piecewise)
	draws <- list(
		theta = Theta_list,
		z = NULL,
		pars = data.frame(
			sigma2 = sigma2_store,
			tau_A2 = tau_A2_store,
			tau_B2 = tau_B2_store,
			g2 = g2_store
		),
		misc = list(
			A = A_store,
			B = B_store,
			M = M_list
		)
	)

	# add Z to draws if stored
	if (FAM$name %in% c("ordinal", "binary") && store_z) {
		Z_list <- lapply(seq_len(n_keep), function(i) Z_store[, , , , i, drop = FALSE])
		for (i in seq_along(Z_list)) {
			dim(Z_list[[i]]) <- dim(Z_store)[1:4]
		}
		draws$z <- Z_list
	}

	# add AR1 parameters if used
	if (ar1) {
		draws$pars$rhoA <- rhoA_store
		draws$pars$rhoB <- rhoB_store
	}

	# add sigma2_obs for gaussian
	if (FAM$name == "gaussian") {
		draws$pars$sigma2_obs <- sigma2_obs_store
	}

	# params data.frame for API consistency
	params <- draws$pars

	# assemble output list with unified structure
	out <- list(
		model = "dynamic",
		family = family,
		Y = Y,
		dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
					is_bipartite = is_bipartite, is_symmetric = symmetric),
		settings = list(
			nscan = nscan,
			burn = burn,
			odens = odens,
			time_thin = time_thin,
			store_z = store_z,
			shrink_rho = shrink_rho,
			draws = n_keep
		),
		meta = list(
			family = family,
			dims = list(m = n_row, n_row = n_row, n_col = n_col, p = p, Tt = Tt,
						is_bipartite = is_bipartite, is_symmetric = symmetric),
			draws = n_keep,
			settings = list(nscan = nscan, burn = burn, odens = odens, time_thin = time_thin),
			sampler_requested = sampler,
			sampler_used = if (sampler_resolved == "exact") "pcg" else "ffbs",
			uncertainty_available = TRUE,
			Omega = pre$Omega,
			approximation_note = if (sampler_resolved == "ffbs" && !symmetric)
			                     "Asymmetric FFBS approximates row covariance propagation."
			                     else NULL
		),
		params = params,
		draws = draws,
		time_kept = time_keep,

		# top-level accessors for convenience
		Theta = Theta_store,
		A = A_store,
		B = B_store,
		M = M_store,
		sigma2 = sigma2_store,
		tau_A2 = tau_A2_store,
		tau_B2 = tau_B2_store,
		g2 = g2_store,
		n_iter = n_iter,
		burn = burn,
		thin = odens,
		time_thin = time_thin
	)

	# gaussian-specific outputs
	if (FAM$name == "gaussian") {
		out$sigma2_obs <- sigma2_obs_store
	}
	if (FAM$name %in% c("ordinal", "binary") && store_z) {
		out$Z <- Z_store
		out$Z_final <- Z_4d  # match static/piecewise naming
	}

	if (ar1) {
		out$rhoA <- rhoA_store
		out$rhoB <- rhoB_store
		out$ar1 <- TRUE
	}

	if (symmetric) {
		out$mh_proposal_sd <- mh_proposal_sd
		out$lambda_diag <- lambda_diag
		out$tau_A_fixed <- tau_A_fixed
		out$kappa_A2 <- kappa_A2
		out$rho_max <- rho_max
		out$rho_max_rejects <- rho_max_rejects
		if (lambda_diag > 0) {
			out$Abar <- Abar_store
			out$phi_diag <- phi_diag
			out$kappa_Abar2 <- kappa_Abar2
		}
	}

	if (!is.null(previous)) {
		prev_total <- previous$n_iter %||% (previous$burn + length(previous$sigma2))
		out$total_iter <- prev_total + nscan
		out$continued_from <- prev_total
	}

	# drop draw-indexed components the caller does not want retained, to
	# control the size of the saved object. The latent-state draws `Theta`
	# are by far the largest component and are not needed for coupling /
	# leverage / rank-probability analysis, which use only A, B, and M.
	# `keep_components` is NULL (keep everything) by default.
	if (!is.null(keep_components)) {
		if (!("Theta" %in% keep_components)) {
			out$Theta <- NULL
			out$draws$theta <- NULL
		}
		if (!("Z" %in% keep_components)) {
			out$Z <- NULL
			out$draws$z <- NULL
		}
		if (!("M" %in% keep_components)) {
			out$M <- NULL
			if (!is.null(out$draws$misc)) out$draws$misc$M <- NULL
		}
		if (!("A" %in% keep_components)) {
			out$A <- NULL
			if (!is.null(out$draws$misc)) out$draws$misc$A <- NULL
		}
		if (!("B" %in% keep_components)) {
			out$B <- NULL
			if (!is.null(out$draws$misc)) out$draws$misc$B <- NULL
		}
		# also drop the bulky $draws bundle unless explicitly requested.
		# this is the slot that holds pars + misc + theta + z; it is the
		# largest storage component after Theta itself
		if (!("draws" %in% keep_components)) {
			out$draws <- NULL
		}
	}

	class(out) <- "dbn"
	out
	####
}
####
