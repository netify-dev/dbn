#include <RcppArmadillo.h>
#include "dbn_stability.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// forward declaration
arma::cube ffbs_theta_struct_cpp(const arma::cube& Z,
                                 const arma::mat& mu,
                                 const arma::cube& A_array,
                                 const arma::cube& B_array,
                                 double sigma2_proc,
                                 double sigma2_obs);

// 5-argument version for backward compatibility (ordinal/binary)
//' Structured FFBS for Theta
//'
//' @description Memory-efficient FFBS that avoids creating n^2 x n^2 matrices.
//'   Supports rectangular (bipartite) Theta via row-covariance approximation.
//' @param Z Observations (n_row x n_col x T array as cube)
//' @param mu Baseline mean (n_row x n_col matrix)
//' @param A_array Time-varying A matrices (n_row x n_row x T)
//' @param B_array Time-varying B matrices (n_col x n_col x T)
//' @param sigma2 Innovation variance (scalar)
//' @return Sampled Theta array (n_row x n_col x T)
// [[Rcpp::export(name = "ffbs_theta_struct_5arg_cpp")]]
arma::cube ffbs_theta_struct_5arg_cpp(const arma::cube& Z,
                                      const arma::mat& mu,
                                      const arma::cube& A_array,
                                      const arma::cube& B_array,
                                      double sigma2) {
  // call 6-argument version with sigma2_obs = 1.0
  return ffbs_theta_struct_cpp(Z, mu, A_array, B_array, sigma2, 1.0);
}

// 6-argument version for gaussian family
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube ffbs_theta_struct_cpp(const arma::cube& Z,
                                 const arma::mat& mu,
                                 const arma::cube& A_array,
                                 const arma::cube& B_array,
                                 double sigma2_proc,
                                 double sigma2_obs) {

  int n_row = Z.n_rows;
  int n_col = Z.n_cols;
  int T = Z.n_slices;

  // make sure variances aren't too small
  sigma2_proc = std::max(sigma2_proc, 1e-6);
  sigma2_obs = std::max(sigma2_obs, 1e-6);

  // set up output array (n_row x n_col x T)
  arma::cube Theta(n_row, n_col, T);

  // Pre-compute centered observations for better cache usage
  arma::cube Y(n_row, n_col, T);
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int t = 0; t < T; ++t) {
    Y.slice(t) = Z.slice(t) - mu;
  }

  // storage for forward pass
  // m_filt: filtered means (n_row x n_col x T)
  // C_filt: row-covariance approximation (n_row x n_row x T)
  //   treats each column of Theta independently with shared row covariance
  arma::cube m_filt(n_row, n_col, T);
  arma::cube C_filt(n_row, n_row, T);

  // initial conditions
  m_filt.slice(0).zeros();
  C_filt.slice(0) = sigma2_proc * arma::eye(n_row, n_row);  // prior variance

  // forward pass
  for (int t = 0; t < T; ++t) {
    if (t > 0) {
      // prediction step: m_t|t-1 = a_t * m_{t-1} * b_t'
      // A(n_row x n_row) * m_filt(n_row x n_col) * B'(n_col x n_col) -> n_row x n_col
      m_filt.slice(t) = A_array.slice(t) * m_filt.slice(t-1) * B_array.slice(t).t();
      // for scalar variance: c_t|t-1 = c_{t-1} + sigma2_proc * I
      C_filt.slice(t) = C_filt.slice(t-1) + sigma2_proc * arma::eye(n_row, n_row);
    }

    // forward update using full matrix form
    arma::mat R_pred = C_filt.slice(t) + sigma2_obs * arma::eye(n_row, n_row);

    // add small regularization to ensure numerical stability
    double reg_eps = std::max(1e-8, 1e-8 * arma::norm(R_pred, "fro"));
    R_pred.diag() += reg_eps;

    // Kalman gain computation using Cholesky
    // K is n_row x n_row, applied to each column of (Y - m_filt) independently
    arma::mat K;
    arma::mat L_R;
    bool chol_success = arma::chol(L_R, force_sym(R_pred), "lower");

    if (chol_success) {
        arma::mat temp = arma::solve(arma::trimatl(L_R), C_filt.slice(t));
        K = arma::solve(arma::trimatu(L_R.t()), temp).t();
    } else {
        // Fallback with regularization
        arma::mat R_sym = 0.5 * (R_pred + R_pred.t());
        arma::mat R_inv;
        bool inv_ok = arma::inv_sympd(R_inv, R_sym);
        if (!inv_ok) {
            double reg = 1e-6 * arma::norm(R_sym, "fro") + 1e-8;
            R_sym.diag() += reg;
            inv_ok = arma::inv_sympd(R_inv, R_sym);
            if (!inv_ok) {
                R_inv = arma::inv(R_sym);
            }
        }
        K = C_filt.slice(t) * R_inv;
    }

    // K(n_row x n_row) * (Y - m_filt)(n_row x n_col) -> n_row x n_col
    m_filt.slice(t) += K * (Y.slice(t) - m_filt.slice(t));

    // Joseph form update for numerical stability
    arma::mat I_minus_K = arma::eye(n_row, n_row) - K;
    C_filt.slice(t) = I_minus_K * C_filt.slice(t) * I_minus_K.t() + K * sigma2_obs * K.t();
  }

  // backward sampling pass
  // sample from posterior at final time point
  arma::mat C_T_reg = C_filt.slice(T-1) + 1e-8 * arma::eye(n_row, n_row);
  // force exact symmetry before cholesky
  C_T_reg = 0.5 * (C_T_reg + C_T_reg.t());
  arma::mat L_T;
  bool chol_T_ok = arma::chol(L_T, C_T_reg, "lower");
  arma::mat noise_T;
  if (chol_T_ok) {
    // L_T(n_row x n_row) * randn(n_row x n_col) -> n_row x n_col noise
    noise_T = L_T * arma::randn(n_row, n_col);
  } else {
    // Fallback: diagonal noise
    arma::vec diag_sd = arma::sqrt(arma::abs(C_T_reg.diag()));
    noise_T = arma::diagmat(diag_sd) * arma::randn(n_row, n_col);
  }
  Theta.slice(T-1) = mu + m_filt.slice(T-1) + noise_T;

  // backward recursion
  for (int t = T-2; t >= 0; --t) {
    // prediction variance at t+1
    arma::mat R_pred = C_filt.slice(t) + sigma2_proc * arma::eye(n_row, n_row);

    // backward smoother update using full matrix form
    arma::mat R_sym = 0.5 * (R_pred + R_pred.t());
    arma::mat R_inv;
    bool inv_ok = arma::inv_sympd(R_inv, R_sym);
    if (!inv_ok) {
      double reg = 1e-6 * arma::norm(R_sym, "fro") + 1e-8;
      R_sym.diag() += reg;
      inv_ok = arma::inv_sympd(R_inv, R_sym);
      if (!inv_ok) {
        R_inv = arma::inv(R_sym);
      }
    }
    arma::mat J = C_filt.slice(t) * R_inv;    // smoother gain

    arma::mat m_smooth = m_filt.slice(t) +
                         J * (Theta.slice(t + 1) - mu -
                              A_array.slice(t + 1) * m_filt.slice(t) *
                              B_array.slice(t + 1).t());

    arma::mat C_smooth = C_filt.slice(t) -
                         J * (R_pred - C_filt.slice(t + 1)) * J.t();

    // sample smoothed state
    arma::mat C_smooth_reg = C_smooth + 1e-8 * arma::eye(n_row, n_row);
    // force exact symmetry before cholesky
    C_smooth_reg = 0.5 * (C_smooth_reg + C_smooth_reg.t());
    arma::mat L_smooth;
    bool chol_smooth_ok = arma::chol(L_smooth, C_smooth_reg, "lower");
    arma::mat noise_t;
    if (chol_smooth_ok) {
      noise_t = L_smooth * arma::randn(n_row, n_col);
    } else {
      arma::vec diag_sd = arma::sqrt(arma::abs(C_smooth_reg.diag()));
      noise_t = arma::diagmat(diag_sd) * arma::randn(n_row, n_col);
    }
    Theta.slice(t) = mu + m_smooth + noise_t;
  }

  return Theta;
}
