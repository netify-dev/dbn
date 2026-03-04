#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' Fast FFBS for Bilinear Model
//'
//' @description Efficient FFBS that avoids Kronecker products.
//'   Supports rectangular (bipartite) Theta via scalar variance approximation.
//' @param Z Observations (n_row x n_col x T)
//' @param mu Baseline mean (n_row x n_col)
//' @param A_array Time-varying A matrices (n_row x n_row x T)
//' @param B_array Time-varying B matrices (n_col x n_col x T)
//' @param sigma2_proc Process variance
//' @param sigma2_obs Observation variance (default 1.0)
//' @return Sampled Theta array (n_row x n_col x T)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube ffbs_bilinear(const arma::cube& Z,
                               const arma::mat& mu,
                               const arma::cube& A_array,
                               const arma::cube& B_array,
                               double sigma2_proc,
                               double sigma2_obs = 1.0) {

  int n_row = Z.n_rows;
  int n_col = Z.n_cols;
  int T = Z.n_slices;

  // floor variances to avoid numerical blowup
  sigma2_proc = std::max(sigma2_proc, 1e-8);
  sigma2_obs = std::max(sigma2_obs, 1e-8);

  // output array
  arma::cube Theta(n_row, n_col, T);

  // forward pass storage
  arma::cube m_fwd(n_row, n_col, T);    // filtered means
  arma::cube v_fwd(n_row, n_col, T);    // filtered variances (scalar approx)

  // first time point (centered: filter estimates Theta - mu)
  double v1 = sigma2_proc / (sigma2_proc + sigma2_obs);
  m_fwd.slice(0) = v1 * (Z.slice(0) - mu);
  v_fwd.slice(0).fill(v1);

  // forward filtering
  for (int t = 1; t < T; t++) {
    // prediction: m_pred = A_t * m_{t-1} * B_t'
    arma::mat m_pred = A_array.slice(t) * m_fwd.slice(t-1) * B_array.slice(t).t();

    // scalar prediction variance
    double v_pred = sigma2_proc + v_fwd(0, 0, t-1);

    // kalman update
    double kalman_gain = v_pred / (v_pred + sigma2_obs);
    m_fwd.slice(t) = m_pred + kalman_gain * (Z.slice(t) - mu - m_pred);
    v_fwd.slice(t).fill(kalman_gain * sigma2_obs);
  }

  // backward sampling
  // sample terminal state
  arma::mat noise_T(n_row, n_col, arma::fill::randn);
  Theta.slice(T-1) = m_fwd.slice(T-1) + std::sqrt(v_fwd(0, 0, T-1)) * noise_T;

  // backward recursion
  for (int t = T-2; t >= 0; t--) {
    // backward smoother gain
    double v_curr = v_fwd(0, 0, t);
    double v_pred_next = sigma2_proc + v_curr;
    double backward_gain = v_curr / v_pred_next;

    // one-step-ahead prediction from filtered mean
    arma::mat theta_pred_next = A_array.slice(t+1) * m_fwd.slice(t) * B_array.slice(t+1).t();

    // smoothed mean
    arma::mat m_back = m_fwd.slice(t) +
      backward_gain * (Theta.slice(t+1) - theta_pred_next);

    // smoothed variance
    double v_back = v_curr * (1.0 - backward_gain);

    // sample from smoothed distribution
    arma::mat noise_t(n_row, n_col, arma::fill::randn);
    Theta.slice(t) = m_back + std::sqrt(v_back) * noise_t;
  }

  // add back baseline mean
  for (int t = 0; t < T; t++) {
    Theta.slice(t) += mu;
  }

  return Theta;
}
