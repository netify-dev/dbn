#include <RcppArmadillo.h>
#include "dbn_stability.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// safe helpers: delegate to the shared definitions in dbn_stability.h
static inline arma::mat safe_inv_sympd_ab(const arma::mat& M) {
    return dbn_safe_inv_sympd(M);
}
static inline arma::vec safe_mvnrnd_ab(const arma::vec& mu, const arma::mat& Sigma) {
    return dbn_safe_mvnrnd(mu, Sigma);
}

// update A via FFBS (forward-filter backward-sample)
// A is n_row x n_row (sender dynamics)
// for row k: y = row k of Theta_t, F = B_t * Theta_{t-1}', a = A(k,:)'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_A_batch(const arma::cube& Theta_all_1, const arma::cube& Theta_all_2,
                    const arma::cube& Theta_all_3, const arma::cube& Theta_all_4,
                    const arma::cube& Aarray, const arma::cube& Barray,
                    double sigma2, double tauA2, bool ar1 = false, double rhoA = 0.0,
                    int p = 1) {
  int n_row = Theta_all_1.n_rows;
  int n_col = Theta_all_1.n_cols;
  int Tt = Theta_all_1.n_slices;

  // output A array
  arma::cube A_new = arma::zeros(n_row, n_row, Tt);

  // observations per time step (n_col per relation)
  int V_dim = (p == 1) ? n_col : n_col * p;
  arma::mat V = arma::eye(V_dim, V_dim) * sigma2;
  arma::mat W = arma::eye(n_row, n_row) * tauA2;
  arma::vec m0 = arma::zeros(n_row);
  arma::mat C0 = arma::eye(n_row, n_row) * tauA2;

  // preallocate y and F storage
  arma::vec y_full(V_dim * Tt);
  arma::mat F_full(V_dim * Tt, n_row);

  // preallocate forward filter storage
  arma::cube m_fwd_store(n_row, Tt, 1);
  arma::cube C_fwd_store(n_row, n_row, Tt);

  // update each row of A
  for(int k = 0; k < n_row; k++) {
    // fill y and F for this row
    y_full.zeros();
    F_full.zeros();

    for(int t = 1; t < Tt; t++) {
      int start_idx = t * V_dim;

      if(p == 1) {
        // single relation case
        // y: row k of Theta_t (length n_col)
        y_full.subvec(start_idx, start_idx + n_col - 1) = Theta_all_1.slice(t).row(k).t();
        // F: (Theta_{t-1} * B_t')' = B_t * Theta_{t-1}' (n_col x n_row)
        F_full.submat(start_idx, 0, start_idx + n_col - 1, n_row - 1) =
          Barray.slice(t) * Theta_all_1.slice(t-1).t();
      } else {
        // multirelation case, stack them up
        for(int rel = 0; rel < p; rel++) {
          int rel_start = start_idx + rel * n_col;

          // select appropriate Theta cube based on relation
          const arma::cube* Theta_ptr = &Theta_all_1;
          if(rel == 1) Theta_ptr = &Theta_all_2;
          else if(rel == 2) Theta_ptr = &Theta_all_3;
          else if(rel == 3) Theta_ptr = &Theta_all_4;

          y_full.subvec(rel_start, rel_start + n_col - 1) =
            (*Theta_ptr).slice(t).row(k).t();
          F_full.submat(rel_start, 0, rel_start + n_col - 1, n_row - 1) =
            Barray.slice(t) * (*Theta_ptr).slice(t-1).t();
        }
      }
    }

    // forward filter
    arma::vec m_fwd = m0;
    arma::mat C_fwd = C0;
    m_fwd_store.zeros();
    C_fwd_store.zeros();

    m_fwd_store.slice(0).col(0) = m_fwd;
    C_fwd_store.slice(0) = C_fwd;

    for(int t = 1; t < Tt; t++) {
      // predict
      if(ar1) {
        m_fwd = rhoA * m_fwd;
      }
      C_fwd = (ar1 ? rhoA * rhoA : 1.0) * C_fwd + W;

      // update with observations
      int start_idx = t * V_dim;
      arma::vec y_t = y_full.subvec(start_idx, start_idx + V_dim - 1);
      arma::mat F_t = F_full.submat(start_idx, 0, start_idx + V_dim - 1, n_row - 1);

      arma::mat S = F_t * C_fwd * F_t.t() + V;
      arma::mat K = C_fwd * F_t.t() * safe_inv_sympd_ab(S);

      m_fwd = m_fwd + K * (y_t - F_t * m_fwd);
      C_fwd = C_fwd - K * F_t * C_fwd;

      m_fwd_store.slice(0).col(t) = m_fwd;
      C_fwd_store.slice(t) = C_fwd;
    }

    // backward sample
    // enforce symmetry before mvnrnd
    C_fwd = 0.5 * (C_fwd + C_fwd.t());
    arma::vec a_sample = safe_mvnrnd_ab(m_fwd, C_fwd);
    A_new.slice(Tt-1).row(k) = a_sample.t();

    for(int t = Tt-2; t >= 0; t--) {
      arma::vec m_t = m_fwd_store.slice(0).col(t);
      arma::mat C_t = C_fwd_store.slice(t);

      // backward recursion
      arma::mat C_pred = (ar1 ? rhoA * rhoA : 1.0) * C_t + W;
      arma::mat B_t = C_t * (ar1 ? rhoA : 1.0) * safe_inv_sympd_ab(C_pred);

      arma::vec h_t = m_t + B_t * (a_sample - (ar1 ? rhoA : 1.0) * m_t);
      arma::mat H_t = C_t - B_t * C_pred * B_t.t();
      // enforce symmetry before mvnrnd
      H_t = 0.5 * (H_t + H_t.t());

      a_sample = safe_mvnrnd_ab(h_t, H_t);
      A_new.slice(t).row(k) = a_sample.t();
    }
  }

  return List::create(Named("A") = A_new);
}

// update B via FFBS (forward-filter backward-sample)
// B is n_col x n_col (receiver dynamics)
// for column k: y = col k of Theta_t, F = A_t * Theta_{t-1}, b = B(:,k)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_B_batch(const arma::cube& Theta_all_1, const arma::cube& Theta_all_2,
                    const arma::cube& Theta_all_3, const arma::cube& Theta_all_4,
                    const arma::cube& Aarray, const arma::cube& Barray,
                    double sigma2, double tauB2, bool ar1 = false, double rhoB = 0.0,
                    int p = 1) {
  int n_row = Theta_all_1.n_rows;
  int n_col = Theta_all_1.n_cols;
  int Tt = Theta_all_1.n_slices;

  // output B array
  arma::cube B_new = arma::zeros(n_col, n_col, Tt);

  // observations per time step (n_row per relation)
  int V_dim = (p == 1) ? n_row : n_row * p;
  arma::mat V = arma::eye(V_dim, V_dim) * sigma2;
  arma::mat W = arma::eye(n_col, n_col) * tauB2;
  arma::vec m0 = arma::zeros(n_col);
  arma::mat C0 = arma::eye(n_col, n_col) * tauB2;

  arma::vec y_full(V_dim * Tt);
  arma::mat F_full(V_dim * Tt, n_col);

  // preallocate forward filter storage
  arma::cube m_fwd_store_B(n_col, Tt, 1);
  arma::cube C_fwd_store_B(n_col, n_col, Tt);

  // update each column of B
  for(int k = 0; k < n_col; k++) {
    y_full.zeros();
    F_full.zeros();

    for(int t = 1; t < Tt; t++) {
      int start_idx = t * V_dim;

      if(p == 1) {
        // y: column k of Theta_t (length n_row)
        y_full.subvec(start_idx, start_idx + n_row - 1) = Theta_all_1.slice(t).col(k);
        // F: A_t * Theta_{t-1} (n_row x n_col)
        // design: A_t * Theta_{t-1} (n_row x n_col)
        F_full.submat(start_idx, 0, start_idx + n_row - 1, n_col - 1) =
          (Aarray.slice(t) * Theta_all_1.slice(t-1));
      } else {
        for(int rel = 0; rel < p; rel++) {
          int rel_start = start_idx + rel * n_row;

          const arma::cube* Theta_ptr = &Theta_all_1;
          if(rel == 1) Theta_ptr = &Theta_all_2;
          else if(rel == 2) Theta_ptr = &Theta_all_3;
          else if(rel == 3) Theta_ptr = &Theta_all_4;

          y_full.subvec(rel_start, rel_start + n_row - 1) =
            (*Theta_ptr).slice(t).col(k);
          F_full.submat(rel_start, 0, rel_start + n_row - 1, n_col - 1) =
            Aarray.slice(t) * (*Theta_ptr).slice(t-1);
        }
      }
    }

    // forward filter
    arma::vec m_fwd = m0;
    arma::mat C_fwd = C0;
    m_fwd_store_B.zeros();
    C_fwd_store_B.zeros();

    m_fwd_store_B.slice(0).col(0) = m_fwd;
    C_fwd_store_B.slice(0) = C_fwd;

    for(int t = 1; t < Tt; t++) {
      if(ar1) {
        m_fwd = rhoB * m_fwd;
      }
      C_fwd = (ar1 ? rhoB * rhoB : 1.0) * C_fwd + W;

      int start_idx = t * V_dim;
      arma::vec y_t = y_full.subvec(start_idx, start_idx + V_dim - 1);
      arma::mat F_t = F_full.submat(start_idx, 0, start_idx + V_dim - 1, n_col - 1);

      arma::mat S = F_t * C_fwd * F_t.t() + V;
      arma::mat K = C_fwd * F_t.t() * safe_inv_sympd_ab(S);

      m_fwd = m_fwd + K * (y_t - F_t * m_fwd);
      C_fwd = C_fwd - K * F_t * C_fwd;

      m_fwd_store_B.slice(0).col(t) = m_fwd;
      C_fwd_store_B.slice(t) = C_fwd;
    }

    // enforce symmetry before mvnrnd
    C_fwd = 0.5 * (C_fwd + C_fwd.t());
    arma::vec b_sample = safe_mvnrnd_ab(m_fwd, C_fwd);
    B_new.slice(Tt-1).col(k) = b_sample;

    for(int t = Tt-2; t >= 0; t--) {
      arma::vec m_t = m_fwd_store_B.slice(0).col(t);
      arma::mat C_t = C_fwd_store_B.slice(t);

      arma::mat C_pred = (ar1 ? rhoB * rhoB : 1.0) * C_t + W;
      arma::mat B_t = C_t * (ar1 ? rhoB : 1.0) * safe_inv_sympd_ab(C_pred);

      arma::vec h_t = m_t + B_t * (b_sample - (ar1 ? rhoB : 1.0) * m_t);
      arma::mat H_t = C_t - B_t * C_pred * B_t.t();
      // enforce symmetry before mvnrnd
      H_t = 0.5 * (H_t + H_t.t());

      b_sample = safe_mvnrnd_ab(h_t, H_t);
      B_new.slice(t).col(k) = b_sample;
    }
  }

  return List::create(Named("B") = B_new);
}
