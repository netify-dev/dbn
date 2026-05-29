####
# package metadata and imports
####

#' @keywords internal
#'
#'
#' @section Known limitations:
#' \describe{
#'   \item{Bipartite HMM/lowrank}{The HMM and low-rank models currently
#'     support only unipartite (square) networks. Bipartite data will
#'     produce an informative error. Use \code{model = "dynamic"} for
#'     bipartite networks.}
#'   \item{Dynamic binary with small networks}{The dynamic model with
#'     \code{family = "binary"} may encounter numerical singularities when
#'     the network has fewer than ~15 nodes. A warning is issued at model
#'     entry. Consider \code{model = "static"} or a larger network.}
#'   \item{HMM label switching}{Regime numbering (1, 2, ..., R) in the
#'     HMM model is arbitrary and may differ across MCMC chains.
#'     Compare regimes by their estimated A/B matrices, not by label.}
#'   \item{Lowrank Stiefel identifiability}{The orthonormal factor matrix
#'     U in the low-rank model is identified only up to orthogonal rotation.
#'     Factor loadings \eqn{\alpha_t} and U should be interpreted together.}
#'   \item{Baseline M identifiability}{The model
#'     \eqn{\Theta_t = M + A_t (\Theta_{t-1} - M) B_t^\top + \varepsilon}
#'     defines \eqn{M} as the long-run baseline of \eqn{\Theta_t}, but the
#'     current M update conditions on \eqn{Z - \Theta} (which carries no
#'     information about \eqn{M} given \eqn{\Theta}), so \eqn{M} shrinks
#'     toward the prior mean (0) and \eqn{\Theta} absorbs the baseline.
#'     Interpret \eqn{\Theta_t} as the un-centred latent state and ignore
#'     \eqn{M}.}
#'   \item{Asymmetric operator scale}{In the asymmetric (directed)
#'     specification \eqn{(A_t, B_t)} and \eqn{(c A_t, c^{-1} B_t)} give the
#'     same likelihood. Under the default RW prior, the \eqn{A_0 = B_0 = I}
#'     anchor pins the scale identification: the prior penalty on
#'     \eqn{A_1 - I} (and \eqn{B_1 - I}) is not invariant under rescaling.
#'     Both \eqn{A_t}, \eqn{B_t} and their product \eqn{W_t = A_t B_t^\top}
#'     are interpretable; \eqn{\tau_A^2} and \eqn{\tau_B^2} are honest
#'     innovation variances. Under the alternative iid prior
#'     (\code{options(dbn.prior_kind = "iid")}), only \eqn{W_t} is identified.}
#'   \item{Operator stability}{Dynamic fits apply a contractivity prior by
#'     default (\code{shrink_rho = 0.9}) that keeps the operator's spectral
#'     radius bounded, so impulse responses and forecasts stay stable. With
#'     \code{shrink_rho = NULL} the raw random-walk prior does not constrain
#'     contractivity and individual periods or draws can be locally explosive.
#'     Use \code{\link{dbn_operator}()} to check the per-period gain.}
#'   \item{Random-walk prior anchored at identity}{The asymmetric non-AR1 path
#'     implements a Gaussian random walk \eqn{A_t = A_{t-1} + \varepsilon^A_t},
#'     \eqn{\varepsilon^A_t \sim N(0, \tau_A^2 I)}, with \eqn{A_0 = I} pinned
#'     as the structural anchor (and analogously for \eqn{B_t}). The anchor
#'     at \eqn{t = 0} resolves the
#'     \eqn{(A_t, B_t) \to (c A_t, c^{-1} B_t)} scale ambiguity that breaks
#'     unanchored RW priors: the prior on \eqn{A_1} is \eqn{N(I, \tau_A^2)},
#'     which is NOT invariant under rescaling, so the chain does not drift
#'     in scale. \eqn{\tau_A^2} and \eqn{\tau_B^2} are honest innovation
#'     variances. \strong{Chain length}: the temporal coupling of the RW
#'     prior makes mixing slower than the iid prior (~20x lower ESS
#'     per A entry, ~3x lower for variance hyperparameters on small data).
#'     For credible intervals to be trustworthy, use at least
#'     \code{nscan = 2000} when \code{Tt} is small; the default 10000 is
#'     comfortable. Multi-chain runs (manual: \code{lapply(seeds, dbn, ...)}
#'     with different \code{seed=} values) recommended for Rhat checks.
#'     The default random-walk prior anchors at \eqn{A_0 = B_0 = I};
#'     pass \code{options(dbn.prior_kind = "iid")} for the iid prior
#'     \eqn{A_t \sim N(0, \tau_A^2 I)}. For explicit AR(1) persistence,
#'     use \code{ar1 = TRUE}.}
#' }
#'
#' @section Priors:
#' The variance parameters carry inverse-gamma priors. The random-walk
#' innovation variances \eqn{\tau_A^2} and \eqn{\tau_B^2} use
#' \eqn{\mathrm{IG}(a_\tau, b_\tau)} with a weakly informative default
#' \eqn{a_\tau = b_\tau = 0.5}, exposed through the \code{a_tau} and
#' \code{b_tau} arguments of \code{\link{dbn_dynamic}()} for prior-sensitivity
#' analysis. The process variance \eqn{\sigma^2}, the observation variance
#' \eqn{\sigma^2_{\mathrm{obs}}}, and the baseline-mean variance \eqn{g^2}
#' carry inverse-gamma priors with fixed weakly informative hyperparameters.
#' The baseline mean \eqn{M} has a Gaussian prior with variance \eqn{g^2}. The
#' operators \eqn{A_t}, \eqn{B_t} follow a Gaussian random walk
#' \eqn{A_t = A_{t-1} + \varepsilon^A_t}, \eqn{\varepsilon^A_t \sim N(0, \tau_A^2 I)}
#' with \eqn{A_0 = B_0 = I} as structural anchors (resolves the
#' \eqn{(A, B) \to (c A, c^{-1} B)} scale ambiguity). On the AR(1) path
#' (\code{ar1 = TRUE}), \eqn{A_t = \rho_A A_{t-1} + \varepsilon^A_t}. Variance
#' posteriors are strongly data-driven and only mildly sensitive to the
#' hyperparameters. An iid alternative is available via
#' \code{options(dbn.prior_kind = "iid")}; see "Known limitations" for
#' the identifiability trade-off.
#'
#' @section Glossary:
#' Terms that recur in the documentation and console output:
#' \describe{
#'   \item{Bilinear autoregression}{The model regresses the network on its own
#'     past through \emph{two} matrix operators, one acting on senders (rows)
#'     and one on receivers (columns): \eqn{\Theta_t = M + A_t (\Theta_{t-1} - M) B_t^\top}.
#'     "Bilinear" means linear in each operator separately.}
#'   \item{MCMC / posterior draws}{The sampler returns many plausible parameter
#'     values (draws) rather than one point estimate; together they form the
#'     posterior distribution. Summaries (mean, 95 percent interval) are computed
#'     across draws.}
#'   \item{Credible interval}{The Bayesian analogue of a confidence interval: a
#'     range containing the parameter with stated posterior probability (e.g.
#'     95 percent). Reported as \code{[lower, upper]} in \code{summary()}.}
#'   \item{Burn-in (\code{burn})}{Early MCMC iterations discarded before the
#'     chain settles; \code{nscan} is the number of post-burn-in iterations.}
#'   \item{Thinning (\code{odens})}{Store every \code{odens}-th draw to save
#'     memory. All iterations are still computed; raising \code{odens} alone
#'     does not improve effective sample size.}
#'   \item{Effective sample size (ESS)}{The number of \emph{independent} draws
#'     the autocorrelated chain is worth. ESS below ~200 signals poor mixing.}
#'   \item{Geweke diagnostic}{A z-score comparing the first and last portions of
#'     a chain; \eqn{|z| > 2} suggests the chain has not converged.}
#'   \item{FFBS}{Forward-filtering backward-sampling -- a Kalman-recursion-based
#'     sampler that draws the latent state of a linear-Gaussian state space
#'     model one time-sweep at a time.}
#'   \item{PCG}{Preconditioned conjugate gradient -- an iterative linear solver
#'     used to sample the latent state without forming large dense matrices.}
#'   \item{Operator gain / spectral radius}{The largest absolute eigenvalue of
#'     the lag operator. Gain below 1 means perturbations decay (stable);
#'     above 1 means they can grow ("locally explosive").}
#'   \item{Label switching}{In the HMM model, regime numbers (1, 2, ...) are
#'     arbitrary and may permute across MCMC draws; compare regimes by their
#'     estimated operators, not by label.}
#'   \item{ALS}{Alternating least squares -- the fast point-estimate algorithm
#'     behind \code{\link{dbn_als}()}. It returns no posterior uncertainty.}
#' }
#'
#' @importFrom grDevices adjustcolor colorRampPalette heat.colors rainbow
#' @importFrom graphics abline arrows barplot image legend lines matplot par plot.new polygon text
#' @importFrom methods as
#' @importFrom stats acf aggregate AIC BIC complete.cases cor density dnorm ecdf kmeans mad median na.pass optim pnorm predict qnorm quantile rbinom rgamma rnorm runif sd time var
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp
#' @useDynLib dbn, .registration = TRUE
"_PACKAGE"

####
# package options
####

#' DBN Package Options
#'
#' @description options controlling R vs C++ backend selection.
#'
#' @section Performance Options:
#' \describe{
#'   \item{dbn.use_cpp_update_ab}{Logical (default: TRUE). Use C++ for A/B updates in HMM.}
#'   \item{dbn.use_cpp_build_f}{Logical (default: TRUE). Use C++ for design matrix F in lowrank.}
#'   \item{dbn.use_batch_ffbs}{Logical (default: TRUE). Use batch FFBS updates.}
#'   \item{dbn.use_cpp_stability}{Logical (default: TRUE). Use C++ for stability functions.}
#'   \item{dbn.use_ffbs_dlm_cpp}{Logical (default: TRUE). Use C++ FFBS for DLMs.}
#'   \item{dbn.use_ffbs_cpp}{Logical (default: TRUE). Use C++ time-varying FFBS.}
#'   \item{dbn.use_cpp_ranklik}{Logical (default: TRUE). Use C++ rank likelihood sampling.}
#' }
#'
#' @section Setting Options:
#' \preformatted{
#' options(dbn.use_cpp_ranklik = FALSE)
#' }
#'
#' @name dbn-options
NULL

####
# global variable declarations
####

if (getRversion() >= "2.15.1") {
	utils::globalVariables(c(
		".", "aarray", "actor", "backward_sample_fast", "barray",
		"collect_dynamic", "collect_hmm", "collect_lowrank", "collect_static",
		"compute_bilinear_residuals_fast", "ecdf", "forward_hmm_fast", "freq", "from",
		"group", "hi", "i", "init_dynamic", "init_hmm", "init_static", "iter",
		"iteration", "j", "loading", "lo", "lower", "med", "posterior_mean",
		"prob", "quant", "receiver", "regime", "rhoA", "rhoB", "running_mean",
		"sel", "sender", "set", "time", "to", "type", "update_AB_dynamic",
		"update_AB_static", "update_Theta_dynamic", "update_Theta_hmm",
		"update_Theta_lowrank", "update_Theta_static", "update_Z_hmm",
		"update_Z_lowrank", "update_Z_static", "update_factor_lowrank",
		"update_hyper_dynamic", "update_hyper_hmm", "update_hyper_lowrank",
		"update_hyper_static", "update_regime_hmm", "update_state_hmm",
		"upper", "val", "value", "z", "k", "linewidth"
	))
}
