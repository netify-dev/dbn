#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "thread_control.h"

using namespace Rcpp;
using namespace arma;

//' Fast parallel Gaussian approximation for ordinal data
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube rz_gaussian_approx_cpp(const arma::cube& R, const arma::cube& Z, 
                                  const arma::cube& EZ, double sigma = 1.0) {
    set_dbn_threads(); // set threads from R options
    
    arma::cube Z_new = Z; // copy
    int n_slices = R.n_slices;
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int s = 0; s < n_slices; s++) {
        arma::mat R_slice = R.slice(s);
        arma::mat Z_slice = Z.slice(s);
        arma::mat EZ_slice = EZ.slice(s);
        
        // find non-NA entries
        arma::uvec valid_idx = find_finite(R_slice);
        
        if (valid_idx.n_elem > 0) {
            // extract valid values
            arma::vec R_valid = R_slice.elem(valid_idx);
            arma::vec EZ_valid = EZ_slice.elem(valid_idx);
            
            // normalize ranks to [0,1]
            double R_min = R_valid.min();
            double R_max = R_valid.max();
            
            arma::vec Z_updates(valid_idx.n_elem);
            
            if (R_max > R_min) {
                arma::vec R_norm = (R_valid - R_min) / (R_max - R_min);
                
                // vectorized sampling with slight bias based on normalized rank
                arma::vec means = EZ_valid + 0.1 * (R_norm - 0.5);
                Z_updates = means + sigma * arma::randn(valid_idx.n_elem);
            } else {
                // all same rank - just sample from prior
                Z_updates = EZ_valid + sigma * arma::randn(valid_idx.n_elem);
            }
            
            // update Z_new
            Z_slice.elem(valid_idx) = Z_updates;
            Z_new.slice(s) = Z_slice;
        }
    }
    
    return Z_new;
}