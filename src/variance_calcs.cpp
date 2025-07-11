#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_residual_sum_squares(const List& Z_field, 
                                    const arma::cube& M,
                                    int p) {
  // z_field contains p cubes, each of dimension m x m x tt
  // m is m x m x p
  
  double rss = 0.0;
  
  for(int j = 0; j < p; j++) {
    arma::cube Z_j = as<arma::cube>(Z_field[j]);
    int Tt = Z_j.n_slices;
    
    for(int t = 0; t < Tt; t++) {
      arma::mat diff = Z_j.slice(t) - M.slice(j);
      // only sum over non-na values
      arma::uvec finite_idx = find_finite(diff);
      rss += accu(diff.elem(finite_idx) % diff.elem(finite_idx));
    }
  }
  
  return rss;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_innovation_variance(const arma::cube& X_curr,
                                   const arma::cube& X_prev,
                                   bool ar1 = false,
                                   double rho = 0.0) {
  // compute sum of squared innovations for AR(1) or random walk
  // X_curr and X_prev are m x m x (Tt-1)
  
  arma::cube innovations;
  
  if(ar1) {
    innovations = X_curr - rho * X_prev;
  } else {
    innovations = X_curr - X_prev;
  }
  
  return accu(innovations % innovations);
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_M_update(const List& Z_field,
                            double g2, int m, int p) {
  // compute posterior mean of M given Z and g2
  
  arma::cube M_sum = arma::zeros(m, m, p);
  arma::cube counts = arma::zeros(m, m, p);
  
  for(int j = 0; j < p; j++) {
    arma::cube Z_j = as<arma::cube>(Z_field[j]);
    int Tt = Z_j.n_slices;
    
    for(int t = 0; t < Tt; t++) {
      arma::mat Z_slice = Z_j.slice(t);
      arma::uvec finite_idx = find_finite(Z_slice);
      
      for(arma::uword idx : finite_idx) {
        arma::uword i = idx % m;
        arma::uword k = idx / m;
        M_sum(i, k, j) += Z_slice(i, k);
        counts(i, k, j) += 1.0;
      }
    }
  }
  
  // compute posterior mean with prior centered at 0
  arma::cube M_new(m, m, p);
  double prior_prec = 1.0 / g2;
  
  for(int j = 0; j < p; j++) {
    for(int i = 0; i < m; i++) {
      for(int k = 0; k < m; k++) {
        double n = counts(i, k, j);
        if(n > 0) {
          double post_prec = n + prior_prec;
          double post_mean = M_sum(i, k, j) / post_prec;
          double post_var = 1.0 / post_prec;
          M_new(i, k, j) = post_mean + sqrt(post_var) * randn();
        } else {
          // No observations, sample from prior
          M_new(i, k, j) = sqrt(g2) * randn();
        }
      }
    }
  }
  
  return M_new;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List compute_AB_innovations(const arma::cube& Aarray,
                            const arma::cube& Barray,
                            bool ar1 = false,
                            double rhoA = 0.0,
                            double rhoB = 0.0) {
  // compute innovations for both A and B arrays
  int Tt = Aarray.n_slices;
  
  // extract relevant time slices
  arma::cube A_curr = Aarray.slices(1, Tt-1);
  arma::cube A_prev = Aarray.slices(0, Tt-2);
  arma::cube B_curr = Barray.slices(1, Tt-1);
  arma::cube B_prev = Barray.slices(0, Tt-2);
  
  double innovA_ss = compute_innovation_variance(A_curr, A_prev, ar1, rhoA);
  double innovB_ss = compute_innovation_variance(B_curr, B_prev, ar1, rhoB);
  
  return List::create(
    Named("innovA_ss") = innovA_ss,
    Named("innovB_ss") = innovB_ss
  );
}