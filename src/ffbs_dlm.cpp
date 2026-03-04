#include <RcppArmadillo.h>
#include "dbn_stability.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "thread_control.h"

// regularize a matrix for numerical stability
arma::mat regularize_matrix(const arma::mat& M, 
                            double eps_scale = 5e-8, 
                            double min_eps = 1e-12, 
                            bool force_pd = false) {
  double eps = std::max(min_eps, eps_scale * arma::mean(M.diag()));
  arma::mat M_reg = M + eps * arma::eye(M.n_rows, M.n_cols);
  
  if (force_pd) {
    // symmetrize before eigendecomposition
    M_reg = 0.5 * (M_reg + M_reg.t());
    
    // eigendecomposition to check for non-PD
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, M_reg);
    
    if (eigval.min() <= 0) {
      // clamp eigenvalues to enforce positive definiteness
      double max_eig = eigval.max();
      if (max_eig <= min_eps) {
        // all eigenvalues tiny, fall back to scaled identity
        M_reg = arma::eye(M_reg.n_rows, M_reg.n_cols) * min_eps;
      } else {
        eigval = arma::clamp(eigval, min_eps, max_eig);
        M_reg = eigvec * arma::diagmat(eigval) * eigvec.t();
        // symmetrize after reconstruction
        M_reg = 0.5 * (M_reg + M_reg.t());
      }
    }
  }
  
  return M_reg;
}

//' Fast FFBS for Dynamic Linear Model in C++
//' 
//' @description C++ implementation of ffbs_dlm with OpenMP parallelization
//' @param y List of observation vectors (length T)
//' @param Flist List of design matrices (T x (r x p))
//' @param V Observation variance matrix (p x p)
//' @param W State innovation variance matrix (r x r)
//' @param m0 Prior mean vector (r x 1)
//' @param C0 Prior covariance matrix (r x r)
//' @param ar1 Logical: use AR(1) dynamics instead of random walk
//' @param rho AR(1) coefficient (used if ar1=TRUE)
//' @return Matrix of sampled state vectors (r x T)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat ffbs_dlm_cpp(const Rcpp::List& y,
                       const Rcpp::List& Flist,
                       const arma::mat& V,
                       const arma::mat& W,
                       const arma::vec& m0,
                       const arma::mat& C0,
                       bool ar1 = false,
                       double rho = 0.0) {
  
  int Tt = y.size();
  int r = m0.n_elem;
  
  // forward pass storage
  std::vector<arma::vec> a(Tt), m(Tt);
  std::vector<arma::mat> R(Tt), C(Tt);
  
  // initial conditions
  a[0] = m0;
  R[0] = C0;
  
  // transition matrix: identity or rho*I for AR(1)
  arma::mat G = ar1 ? arma::mat(arma::eye(r, r) * rho) : arma::mat(arma::eye(r, r));
  
  // forward filtering
  for (int t = 0; t < Tt; t++) {
    // observation and design matrix at time t
    arma::vec y_t = Rcpp::as<arma::vec>(y[t]);
    arma::mat F_t = Rcpp::as<arma::mat>(Flist[t]);
    
    // prediction step (skip at t=0)
    if (t > 0) {
      a[t] = G * m[t-1];
      R[t] = G * C[t-1] * G.t() + W;
    }
    
    // kalman update step
    arma::mat Qt = F_t * R[t] * F_t.t() + V;
    
    // regularize for numerical stability
    Qt = regularize_matrix(Qt);
    
    // kalman gain via safe inverse
    arma::mat Qt_inv;
    bool qt_ok = arma::inv_sympd(Qt_inv, Qt);
    if (!qt_ok) {
      double reg = 1e-6 * arma::norm(Qt, "fro") + 1e-8;
      Qt.diag() += reg;
      qt_ok = arma::inv_sympd(Qt_inv, Qt);
      if (!qt_ok) { Qt_inv = arma::inv(Qt); }
    }
    arma::mat Kt = R[t] * F_t.t() * Qt_inv;
    
    // filtered mean and covariance
    arma::vec et = y_t - F_t * a[t];
    m[t] = a[t] + Kt * et;
    C[t] = R[t] - Kt * F_t * R[t];
    
    // enforce symmetry on C[t]
    C[t] = 0.5 * (C[t] + C[t].t());
  }
  
  // backward sampling
  arma::mat theta(r, Tt);
  
  // sample terminal state
  arma::mat Ctt_reg = regularize_matrix(C[Tt-1], 5e-8, 1e-12, true);
  try {
    theta.col(Tt-1) = arma::mvnrnd(m[Tt-1], Ctt_reg);
  } catch (...) {
    theta.col(Tt-1) = m[Tt-1] + arma::sqrt(arma::abs(Ctt_reg.diag())) % arma::randn(r);
  }

  // backward recursion
  for (int t = Tt-2; t >= 0; t--) {
    // regularize R[t+1] to ensure PD
    arma::mat Rt1_reg = regularize_matrix(R[t+1], 5e-8, 1e-12, true);

    // backward smoother gain
    arma::mat Rt1_inv;
    bool inv_ok = arma::inv_sympd(Rt1_inv, Rt1_reg);
    if (!inv_ok) {
      double reg = 1e-6 * arma::norm(Rt1_reg, "fro") + 1e-8;
      Rt1_reg.diag() += reg;
      inv_ok = arma::inv_sympd(Rt1_inv, Rt1_reg);
      if (!inv_ok) { Rt1_inv = arma::inv(Rt1_reg); }
    }
    arma::mat J = C[t] * G.t() * Rt1_inv;

    // smoothed mean and covariance
    arma::vec mt = m[t] + J * (theta.col(t+1) - a[t+1]);
    arma::mat Ct = C[t] - J * (R[t+1] - C[t+1]) * J.t();

    // regularize to ensure PD
    Ct = regularize_matrix(Ct, 5e-8, 1e-12, true);

    // sample
    try {
      theta.col(t) = arma::mvnrnd(mt, Ct);
    } catch (...) {
      theta.col(t) = mt + arma::sqrt(arma::abs(Ct.diag())) % arma::randn(r);
    }
  }
  
  return theta;
}

//' Batch FFBS for multiple independent DLM problems
//' 
//' @description Parallelized FFBS for updating all rows of A or columns of B simultaneously
//' @param Y_batch Stacked observations (batch_size x T)
//' @param F_batch Design matrices for each problem (batch_size x T x state_dim x obs_dim)
//' @param V_batch Observation variances (can be scalar or vector of length batch_size)
//' @param W State innovation variance (common across batch)
//' @param m0 Prior mean (common)
//' @param C0 Prior covariance (common)
//' @param ar1 Use AR(1) dynamics
//' @param rho AR(1) coefficient
//' @return Matrix of sampled states (batch_size x state_dim x T)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube ffbs_dlm_batch_cpp(const arma::mat& Y_batch,
                              const arma::cube& F_batch,
                              const arma::vec& V_batch,
                              const arma::mat& W,
                              const arma::vec& m0,
                              const arma::mat& C0,
                              bool ar1 = false,
                              double rho = 0.0) {
  set_dbn_threads();
  
  int batch_size = Y_batch.n_rows;
  int Tt = Y_batch.n_cols;
  int r = m0.n_elem;
  
  // F_batch layout: batch_size x T x obs_dim x state_dim
  int obs_dim = F_batch.n_rows / (batch_size * Tt);
  int state_dim = r;
  
  // output: batch_size x r x Tt
  arma::cube theta_batch(batch_size, r, Tt);
  
  // sequential: arma::mvnrnd and Rcpp::List aren't thread-safe
  for (int b = 0; b < batch_size; b++) {
    // build observation and design-matrix lists for this element
    Rcpp::List y_list(Tt);
    Rcpp::List F_list(Tt);
    
    for (int t = 0; t < Tt; t++) {
      y_list[t] = Rcpp::wrap(arma::vec(1, arma::fill::value(Y_batch(b, t))));
      // extract design matrix from flattened F_batch
      arma::mat F_slice(obs_dim, state_dim);
      for (int i = 0; i < obs_dim; i++) {
        for (int j = 0; j < state_dim; j++) {
          F_slice(i, j) = F_batch(b, t, i * state_dim + j);
        }
      }
      F_list[t] = Rcpp::wrap(F_slice);
    }
    
    // scalar observation variance for this element
    arma::mat V_b = arma::eye(1, 1) * V_batch(b);
    
    // run FFBS
    arma::mat theta_b = ffbs_dlm_cpp(y_list, F_list, V_b, W, m0, C0, ar1, rho);
    
    // pack into output cube
    for (int t = 0; t < Tt; t++) {
      theta_batch.subcube(b, 0, t, b, r-1, t) = theta_b.col(t).t();
    }
  }
  
  return theta_batch;
}