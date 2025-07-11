#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_A_batch(const arma::cube& Theta_all_1, const arma::cube& Theta_all_2,
                    const arma::cube& Theta_all_3, const arma::cube& Theta_all_4,
                    const arma::cube& Aarray, const arma::cube& Barray,
                    double sigma2, double tauA2, bool ar1 = false, double rhoA = 0.0,
                    int p = 1) {
  // theta_all_1,2,3,4 are the 4d array sliced by relation (if p > 1)
  // for p=1, only theta_all_1 is used
  // dimensions: m x m x tt
  
  int m = Theta_all_1.n_rows;
  int Tt = Theta_all_1.n_slices;
  
  // output a array
  arma::cube A_new = arma::zeros(m, m, Tt);
  
  // preallocate workspace for all rows
  int V_dim = (p == 1) ? m : m * p;
  arma::mat V = arma::eye(V_dim, V_dim) * sigma2;
  arma::mat W = arma::eye(m, m) * tauA2;
  arma::vec m0 = arma::zeros(m);
  arma::mat C0 = arma::eye(m, m) * tauA2;
  
  // preallocate y and f storage for maximum size
  arma::vec y_full(V_dim * Tt);
  arma::mat F_full(V_dim * Tt, m);
  
  // update each row of A
  for(int k = 0; k < m; k++) {
    // Fill y and F for this row
    y_full.zeros();
    F_full.zeros();
    
    for(int t = 1; t < Tt; t++) {
      int start_idx = t * V_dim;
      
      if(p == 1) {
        // single relation case
        y_full.subvec(start_idx, start_idx + m - 1) = Theta_all_1.slice(t).row(k).t();
        F_full.submat(start_idx, 0, start_idx + m - 1, m - 1) = 
          Barray.slice(t) * Theta_all_1.slice(t-1).t();
      } else {
        // multirelation case, stack them up
        for(int rel = 0; rel < p; rel++) {
          int rel_start = start_idx + rel * m;
          
          // select appropriate Theta cube based on relation
          const arma::cube* Theta_ptr = &Theta_all_1;
          if(rel == 1) Theta_ptr = &Theta_all_2;
          else if(rel == 2) Theta_ptr = &Theta_all_3;
          else if(rel == 3) Theta_ptr = &Theta_all_4;
          
          y_full.subvec(rel_start, rel_start + m - 1) = 
            (*Theta_ptr).slice(t).row(k).t();
          F_full.submat(rel_start, 0, rel_start + m - 1, m - 1) = 
            Barray.slice(t) * (*Theta_ptr).slice(t-1).t();
        }
      }
    }
    
    // fwd filter
    arma::vec m_fwd = m0;
    arma::mat C_fwd = C0;
    arma::cube m_fwd_store(m, Tt, 1);
    arma::cube C_fwd_store(m, m, Tt);
    
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
      arma::mat F_t = F_full.submat(start_idx, 0, start_idx + V_dim - 1, m - 1);
      
      arma::mat S = F_t * C_fwd * F_t.t() + V;
      arma::mat K = C_fwd * F_t.t() * inv_sympd(S);
      
      m_fwd = m_fwd + K * (y_t - F_t * m_fwd);
      C_fwd = C_fwd - K * F_t * C_fwd;
      
      // 
      m_fwd_store.slice(0).col(t) = m_fwd;
      C_fwd_store.slice(t) = C_fwd;
    }
    
    // backward samp
    // force exact symmetry before mvnrnd
    C_fwd = 0.5 * (C_fwd + C_fwd.t());
    arma::vec a_sample = mvnrnd(m_fwd, C_fwd);
    A_new.slice(Tt-1).row(k) = a_sample.t();
    
    for(int t = Tt-2; t >= 0; t--) {
      arma::vec m_t = m_fwd_store.slice(0).col(t);
      arma::mat C_t = C_fwd_store.slice(t);
      
      // backward recursion
      arma::mat C_pred = (ar1 ? rhoA * rhoA : 1.0) * C_t + W;
      arma::mat B_t = C_t * (ar1 ? rhoA : 1.0) * inv_sympd(C_pred);
      
      arma::vec h_t = m_t + B_t * (a_sample - (ar1 ? rhoA : 1.0) * m_t);
      arma::mat H_t = C_t - B_t * C_pred * B_t.t();
      // force exact symmetry before mvnrnd
      H_t = 0.5 * (H_t + H_t.t());
      
      a_sample = mvnrnd(h_t, H_t);
      A_new.slice(t).row(k) = a_sample.t();
    }
  }
  
  return List::create(Named("A") = A_new);
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_B_batch(const arma::cube& Theta_all_1, const arma::cube& Theta_all_2,
                    const arma::cube& Theta_all_3, const arma::cube& Theta_all_4,
                    const arma::cube& Aarray, const arma::cube& Barray,
                    double sigma2, double tauB2, bool ar1 = false, double rhoB = 0.0,
                    int p = 1) {
  // similar structure but for columns of B
  int m = Theta_all_1.n_rows;
  int Tt = Theta_all_1.n_slices;
  
  arma::cube B_new = arma::zeros(m, m, Tt);
  
  int V_dim = (p == 1) ? m : m * p;
  arma::mat V = arma::eye(V_dim, V_dim) * sigma2;
  arma::mat W = arma::eye(m, m) * tauB2;
  arma::vec m0 = arma::zeros(m);
  arma::mat C0 = arma::eye(m, m) * tauB2;
  
  arma::vec y_full(V_dim * Tt);
  arma::mat F_full(V_dim * Tt, m);
  
  // update each column of B
  for(int k = 0; k < m; k++) {
    y_full.zeros();
    F_full.zeros();
    
    for(int t = 1; t < Tt; t++) {
      int start_idx = t * V_dim;
      
      if(p == 1) {
        y_full.subvec(start_idx, start_idx + m - 1) = Theta_all_1.slice(t).col(k);
        F_full.submat(start_idx, 0, start_idx + m - 1, m - 1) = 
          (Aarray.slice(t) * Theta_all_1.slice(t-1)).t();
      } else {
        for(int rel = 0; rel < p; rel++) {
          int rel_start = start_idx + rel * m;
          
          const arma::cube* Theta_ptr = &Theta_all_1;
          if(rel == 1) Theta_ptr = &Theta_all_2;
          else if(rel == 2) Theta_ptr = &Theta_all_3;
          else if(rel == 3) Theta_ptr = &Theta_all_4;
          
          y_full.subvec(rel_start, rel_start + m - 1) = 
            (*Theta_ptr).slice(t).col(k);
          F_full.submat(rel_start, 0, rel_start + m - 1, m - 1) = 
            Aarray.slice(t) * (*Theta_ptr).slice(t-1);
        }
      }
    }
    
    // forward filter and backward sample (same structure as A update)
    arma::vec m_fwd = m0;
    arma::mat C_fwd = C0;
    arma::cube m_fwd_store(m, Tt, 1);
    arma::cube C_fwd_store(m, m, Tt);
    
    m_fwd_store.slice(0).col(0) = m_fwd;
    C_fwd_store.slice(0) = C_fwd;
    
    for(int t = 1; t < Tt; t++) {
      if(ar1) {
        m_fwd = rhoB * m_fwd;
      }
      C_fwd = (ar1 ? rhoB * rhoB : 1.0) * C_fwd + W;
      
      int start_idx = t * V_dim;
      arma::vec y_t = y_full.subvec(start_idx, start_idx + V_dim - 1);
      arma::mat F_t = F_full.submat(start_idx, 0, start_idx + V_dim - 1, m - 1);
      
      arma::mat S = F_t * C_fwd * F_t.t() + V;
      arma::mat K = C_fwd * F_t.t() * inv_sympd(S);
      
      m_fwd = m_fwd + K * (y_t - F_t * m_fwd);
      C_fwd = C_fwd - K * F_t * C_fwd;
      
      m_fwd_store.slice(0).col(t) = m_fwd;
      C_fwd_store.slice(t) = C_fwd;
    }
    
    // force exact symmetry before mvnrnd
    C_fwd = 0.5 * (C_fwd + C_fwd.t());
    arma::vec b_sample = mvnrnd(m_fwd, C_fwd);
    B_new.slice(Tt-1).col(k) = b_sample;
    
    for(int t = Tt-2; t >= 0; t--) {
      arma::vec m_t = m_fwd_store.slice(0).col(t);
      arma::mat C_t = C_fwd_store.slice(t);
      
      arma::mat C_pred = (ar1 ? rhoB * rhoB : 1.0) * C_t + W;
      arma::mat B_t = C_t * (ar1 ? rhoB : 1.0) * inv_sympd(C_pred);
      
      arma::vec h_t = m_t + B_t * (b_sample - (ar1 ? rhoB : 1.0) * m_t);
      arma::mat H_t = C_t - B_t * C_pred * B_t.t();
      // force exact symmetry before mvnrnd
      H_t = 0.5 * (H_t + H_t.t());
      
      b_sample = mvnrnd(h_t, H_t);
      B_new.slice(t).col(k) = b_sample;
    }
  }
  
  return List::create(Named("B") = B_new);
}