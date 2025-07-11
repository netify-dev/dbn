#' @keywords internal
#' @importFrom grDevices adjustcolor heat.colors rainbow
#' @importFrom graphics abline barplot image legend lines matplot par polygon
#' @importFrom methods as
#' @importFrom stats aggregate cor density dnorm ecdf kmeans median na.pass pnorm predict qnorm quantile rbinom rgamma rnorm runif sd time var
#' @importFrom utils tail
#' @importFrom Rcpp evalCpp
"_PACKAGE"

#' DBN Package Options
#'
#' @description The dbn package uses several options to control behavior,
#'   particularly for toggling between R and C++ implementations.
#'
#' @section Performance Options:
#' \describe{
#'   \item{dbn.use_cpp_update_ab}{Logical (default: TRUE). Use C++ implementation
#'     for updating A and B matrices in HMM models.}
#'   \item{dbn.use_cpp_build_f}{Logical (default: TRUE). Use C++ implementation
#'     for building design matrix F in low-rank models.}
#'   \item{dbn.use_batch_ffbs}{Logical (default: TRUE). Use batch FFBS updates
#'     when available for better performance.}
#'   \item{dbn.use_cpp_stability}{Logical (default: TRUE). Use C++ implementations
#'     for matrix stability functions (spectral radius, positive definite checks).}
#'   \item{dbn.use_ffbs_dlm_cpp}{Logical (default: TRUE). Use C++ implementation
#'     of Forward-Filter Backward-Sample for dynamic linear models.}
#'   \item{dbn.use_ffbs_cpp}{Logical (default: TRUE). Use C++ implementation
#'     of time-varying FFBS algorithm.}
#'   \item{dbn.use_cpp_ranklik}{Logical (default: TRUE). Use C++ implementation
#'     for rank likelihood sampling (significant speedup on large networks).}
#' }
#'
#' @section Setting Options:
#' Options can be set using \code{options()}:
#' \preformatted{
#' # Disable C++ rank likelihood (e.g., for debugging)
#' options(dbn.use_cpp_ranklik = FALSE)
#'
#' # Disable all C++ implementations
#' options(
#'   dbn.use_cpp_update_ab = FALSE,
#'   dbn.use_cpp_build_f = FALSE,
#'   dbn.use_batch_ffbs = FALSE,
#'   dbn.use_cpp_stability = FALSE,
#'   dbn.use_ffbs_dlm_cpp = FALSE,
#'   dbn.use_ffbs_cpp = FALSE,
#'   dbn.use_cpp_ranklik = FALSE
#' )
#' }
#'
#' @name dbn-options
NULL

# Global variables from ggplot2 and other packages
if(getRversion() >= "2.15.1") {
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
    "upper", "val", "value", "z"
  ))
}
