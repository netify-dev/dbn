#include <RcppArmadillo.h>

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
//' @description Memory-efficient FFBS that avoids creating n² × n² matrices
//' @param Z Observations (m × m × T array as cube)
//' @param mu Baseline mean (m × m matrix)
//' @param A_array Time-varying A matrices (m × m × T)
//' @param B_array Time-varying B matrices (m × m × T)
//' @param sigma2 Innovation variance (scalar)
//' @return Sampled Theta array (m × m × T)
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
  
  int m = Z.n_rows;
  int T = Z.n_slices;
  
  // make sure variances aren't too small
  sigma2_proc = std::max(sigma2_proc, 1e-6);
  sigma2_obs = std::max(sigma2_obs, 1e-6);
  
  // set up output array
  arma::cube Theta(m, m, T);
  
  // Pre-compute centered observations for better cache usage
  arma::cube Y(m, m, T);
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int t = 0; t < T; ++t) {
    Y.slice(t) = Z.slice(t) - mu;
  }
  
  // storage for forward pass
  arma::cube m_filt(m, m, T);  // filtered means
  arma::cube C_filt(m, m, T);  // filtered variances (diagonal stored as matrix)
  
  // initial conditions
  m_filt.slice(0).zeros();
  C_filt.slice(0) = sigma2_proc * arma::eye(m, m);  // prior variance
  
  // forward pass
  for (int t = 0; t < T; ++t) {
    if (t > 0) {
      // prediction step: m_t|t-1 = a_t * m_{t-1} * b_t'
      m_filt.slice(t) = A_array.slice(t) * m_filt.slice(t-1) * B_array.slice(t).t();
      // for scalar variance: c_t|t-1 = c_{t-1} + sigma2_proc * i
      C_filt.slice(t) = C_filt.slice(t-1) + sigma2_proc * arma::eye(m, m);
    }
    
    // forward update using full matrix form
    arma::mat R_pred = C_filt.slice(t) + sigma2_obs * arma::eye(m, m);  // observation variance
    
    // add small regularization to ensure numerical stability
    double reg_eps = std::max(1e-8, 1e-8 * arma::norm(R_pred, "fro"));
    R_pred.diag() += reg_eps;
    
    // Kalman gain computation using Cholesky
    arma::mat K;
    arma::mat L_R;
    bool chol_success = arma::chol(L_R, R_pred, "lower");
    
    if (chol_success) {
        arma::mat temp = arma::solve(arma::trimatl(L_R), C_filt.slice(t));
        K = arma::solve(arma::trimatu(L_R.t()), temp).t();
    } else {
        // Fallback to standard computation
        K = C_filt.slice(t) * arma::inv_sympd(R_pred);
    }
    
    m_filt.slice(t) += K * (Y.slice(t) - m_filt.slice(t));
    
    // form update for numerical stability
    arma::mat I_minus_K = arma::eye(m, m) - K;
    C_filt.slice(t) = I_minus_K * C_filt.slice(t) * I_minus_K.t() + K * sigma2_obs * K.t();
  }
  
  // backward sampling pass
  // sample from posterior at time t
  arma::mat C_T_reg = C_filt.slice(T-1) + 1e-8 * arma::eye(m, m);  // regularize
  // force exact symmetry before cholesky
  C_T_reg = 0.5 * (C_T_reg + C_T_reg.t());
  arma::mat L_T = arma::chol(C_T_reg, "lower");
  arma::mat noise_T = L_T * arma::randn(m, m);
  Theta.slice(T-1) = mu + m_filt.slice(T-1) + noise_T;
  
  // backward recursion
  for (int t = T-2; t >= 0; --t) {
    // prediction variance at t+1
    arma::mat R_pred = C_filt.slice(t) + sigma2_proc * arma::eye(m, m);
    
    // backward smoother update using full matrix form
    arma::mat J = C_filt.slice(t) * arma::inv_sympd(R_pred);    // smoother gain
    
    arma::mat m_smooth = m_filt.slice(t) +
                         J * (Theta.slice(t + 1) - mu -
                              A_array.slice(t + 1) * m_filt.slice(t) *
                              B_array.slice(t + 1).t());
    
    arma::mat C_smooth = C_filt.slice(t) -
                         J * (R_pred - C_filt.slice(t + 1)) * J.t();
    
    // sample smoothed state
    arma::mat C_smooth_reg = C_smooth + 1e-8 * arma::eye(m, m);  // regularize
    // force exact symmetry before cholesky
    C_smooth_reg = 0.5 * (C_smooth_reg + C_smooth_reg.t());
    arma::mat L_smooth = arma::chol(C_smooth_reg, "lower");
    arma::mat noise_t = L_smooth * arma::randn(m, m);
    Theta.slice(t) = mu + m_smooth + noise_t;
  }
  
  return Theta;
}