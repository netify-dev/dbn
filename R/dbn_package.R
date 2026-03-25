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
#' }
#'
#' @importFrom grDevices adjustcolor colorRampPalette heat.colors rainbow
#' @importFrom graphics abline arrows barplot image legend lines matplot par plot.new polygon text
#' @importFrom methods as
#' @importFrom stats acf aggregate complete.cases cor density dnorm ecdf kmeans median na.pass pnorm predict qnorm quantile rbinom rgamma rnorm runif sd time var
#' @importFrom utils tail
#' @importFrom Rcpp evalCpp
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
