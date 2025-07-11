#include <RcppArmadillo.h>
#include <omp.h>
#include "dbn_thread_utils.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double compute_gaussian_obs_residuals_cpp(const arma::cube& Z, 
                                         const arma::cube& Theta, 
                                         const arma::cube& M) {
    // set number of threads from r options
    set_dbn_threads();
    
    const int m = Z.n_rows;
    const int m2 = Z.n_cols;
    const int n_slices = Z.n_slices;
    
    double resid_obs = 0.0;
    
    #pragma omp parallel reduction(+:resid_obs)
    {
        #pragma omp for
        for (int k = 0; k < n_slices; k++) {
            mat diff = Z.slice(k) - Theta.slice(k) - M.slice(k % M.n_slices);
            resid_obs += accu(diff % diff);
        }
    }
    
    return resid_obs;
}

// [[Rcpp::export]]
double compute_gaussian_obs_residuals_4d_cpp(const arma::cube& Z_flat,
                                             const arma::cube& Theta_flat,
                                             const arma::cube& M_flat,
                                             int m, int p, int Tt) {
    // set number of threads from r options
    set_dbn_threads();
    
    double resid_obs = 0.0;
    
    #pragma omp parallel reduction(+:resid_obs)
    {
        #pragma omp for
        for (int j = 0; j < p; j++) {
            int offset = j * Tt;
            for (int t = 0; t < Tt; t++) {
                mat diff = Z_flat.slice(offset + t) - 
                          Theta_flat.slice(offset + t) - 
                          M_flat.slice(j);
                resid_obs += accu(diff % diff);
            }
        }
    }
    
    return resid_obs;
}

// vectorized version for better performance with large arrays
// [[Rcpp::export]]
double compute_gaussian_obs_residuals_batch_cpp(const arma::cube& Z,
                                               const arma::cube& Theta,
                                               const arma::mat& M,
                                               int m, int p, int Tt) {
    // set number of threads from r options
    set_dbn_threads();
    
    const int n_total = p * Tt;
    double resid_obs = 0.0;
    
    // parallel over time-relation combinations
    #pragma omp parallel reduction(+:resid_obs)
    {
        #pragma omp for
        for (int idx = 0; idx < n_total; idx++) {
            int j = idx / Tt;  // relation index
            int t = idx % Tt;  // time index
            int slice_idx = j * Tt + t;
            
            // compute residual for this slice
            mat diff = Z.slice(slice_idx) - Theta.slice(slice_idx);
            
            // subtract appropriate M slice
            if (M.n_cols == p) {
                // M is vectorized (m*m x p)
                mat M_j(m, m);
                for (int i = 0; i < m; i++) {
                    for (int k = 0; k < m; k++) {
                        M_j(i, k) = M(i * m + k, j);
                    }
                }
                diff -= M_j;
            }
            
            resid_obs += accu(diff % diff);
        }
    }
    
    return resid_obs;
}