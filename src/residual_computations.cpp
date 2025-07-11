#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// forward declaration for function defined in rank_likelihood_fast.cpp
arma::vec rz_fc_cpp(const arma::vec& R, const arma::vec& Z, const arma::vec& EZ, const List& iranks);

// fast bilinear product A * Theta * B^T avoiding intermediate matrices
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat bilinear_product_fast(const arma::mat& A, const arma::mat& Theta, const arma::mat& B) {
    int m = A.n_rows;
    arma::mat result(m, m);
    
    // (A * Theta) * B^T
    // but done element-wise to avoid large intermediate matrix
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            double sum = 0.0;
            for(int k = 0; k < m; k++) {
                for(int l = 0; l < m; l++) {
                    sum += A(i, k) * Theta(k, l) * B(j, l);
                }
            }
            result(i, j) = sum;
        }
    }
    
    return result;
}

// compute bilinear residuals for all relations and time points in one pass
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_bilinear_residuals(const arma::cube& Theta_flat,
                                 const arma::cube& Aarray,
                                 const arma::cube& Barray,
                                 int m, int p, int Tt) {
    double total_rss = 0.0;
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:total_rss)
    #endif
    for(int rel = 0; rel < p; rel++) {
        for(int t = 1; t < Tt; t++) {
            int idx_curr = rel * Tt + t;
            int idx_prev = rel * Tt + t - 1;
            
            arma::mat Theta_curr = Theta_flat.slice(idx_curr);
            arma::mat Theta_prev = Theta_flat.slice(idx_prev);
            arma::mat A_t = Aarray.slice(t);
            arma::mat B_t = Barray.slice(t);
            
            // compute residual: Theta_t - A_t * Theta_{t-1} * B_t^T
            arma::mat pred = A_t * Theta_prev * B_t.t();
            arma::mat resid = Theta_curr - pred;
            
            total_rss += accu(resid % resid);
        }
    }
    
    return total_rss;
}

// Bilinear residuals 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_bilinear_residuals_fast(const arma::cube& Theta,
                                      const arma::cube& Aarray,
                                      const arma::cube& Barray,
                                      int m, int p, int Tt) {
    double total_rss = 0.0;
    
    // process in blocks for better cache utilization
    const int block_size = std::min(10, Tt-1);
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:total_rss)
    #endif
    for(int rel = 0; rel < p; rel++) {
        for(int t_start = 1; t_start < Tt; t_start += block_size) {
            int t_end = std::min(t_start + block_size, Tt);
            
            // process block of time steps
            for(int t = t_start; t < t_end; t++) {
                // for 4D array, we need to flatten the index
                int idx_curr = rel * Tt + t;
                int idx_prev = rel * Tt + t - 1;
                
                arma::mat Theta_curr = Theta.slice(idx_curr);
                arma::mat Theta_prev = Theta.slice(idx_prev);
                
                const arma::mat& A_t = Aarray.slice(t);
                const arma::mat& B_t = Barray.slice(t);
                
                // bilinear product using temp variable
                arma::mat temp = A_t * Theta_prev;
                arma::mat pred = temp * B_t.t();
                
                // compute residual
                pred = Theta_curr - pred;
                total_rss += accu(pred % pred);
            }
        }
    }
    
    return total_rss;
}

// compute observation residuals for gaussian family
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_observation_residuals(const arma::cube& Z_flat,
                                   const arma::cube& Theta_flat,
                                   const arma::cube& M,
                                   int m, int p, int Tt) {
    double total_rss = 0.0;
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:total_rss)
    #endif
    for(int rel = 0; rel < p; rel++) {
        arma::mat M_rel = M.slice(rel);
        
        for(int t = 0; t < Tt; t++) {
            int idx = rel * Tt + t;
            
            arma::mat Z_rt = Z_flat.slice(idx);
            arma::mat Theta_rt = Theta_flat.slice(idx);
            
            // residual: Z - (Theta + M)
            arma::mat resid = Z_rt - (Theta_rt + M_rel);
            
            total_rss += accu(resid % resid);
        }
    }
    
    return total_rss;
}

// vectorized z-score computation for preprocessing
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_zscores_batch(const arma::cube& Y, const arma::vec& means, 
                                 const arma::vec& sds, int m, int p, int Tt) {
    arma::cube Z(m, m, p * Tt);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int idx = 0; idx < p * Tt; idx++) {
        int rel = idx / Tt;
        arma::mat Y_slice = Y.slice(idx);
        
        // vectorized z-score computation
        Z.slice(idx) = (Y_slice - means(rel)) / sds(rel);
        
        // preserve diagonal as NA
        for(int i = 0; i < m; i++) {
            Z(i, i, idx) = datum::nan;
        }
    }
    
    return Z;
}

// batch update Z for all relations (vectorized version)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube update_Z_batch(const arma::cube& R_flat,
                         const arma::cube& Theta_flat, 
                         const arma::cube& M,
                         const List& IR_list_flat,
                         int m, int p, int Tt) {
    int n_total = p * Tt;
    arma::cube Z_new(m, m, n_total);
    
    // process all relation-time pairs in parallel
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int idx = 0; idx < n_total; idx++) {
        int rel = idx / Tt;
        
        // get data for this slice
        arma::mat R_slice = R_flat.slice(idx);
        arma::mat Theta_slice = Theta_flat.slice(idx);
        arma::mat M_rel = M.slice(rel);
        
        // compute EZ = Theta + M
        arma::mat EZ = Theta_slice + M_rel;
        
        // vectorize and update using rank likelihood
        arma::vec R_vec = vectorise(R_slice);
        arma::vec EZ_vec = vectorise(EZ);
        
        // get rank indices for this slice
        List IR_idx = IR_list_flat[idx];
        
        // sample Z using rank likelihood
        arma::vec Z_vec = rz_fc_cpp(R_vec, EZ_vec, EZ_vec, IR_idx);
        
        // reshape and store
        Z_new.slice(idx) = reshape(Z_vec, m, m);
    }
    
    return Z_new;
}


// vectorized regime array construction for hmm
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List build_regime_arrays_vectorized(const IntegerVector& S,
                                   const List& A_list,
                                   const List& B_list,
                                   int m, int Tt) {
    arma::cube Aarray(m, m, Tt);
    arma::cube Barray(m, m, Tt);
    
    // vectorized assignment using advanced indexing
    for(int t = 0; t < Tt; t++) {
        int regime = S[t] - 1; // convert to 0-based
        Aarray.slice(t) = as<arma::mat>(A_list[regime]);
        Barray.slice(t) = as<arma::mat>(B_list[regime]);
    }
    
    return List::create(
        Named("Aarray") = Aarray,
        Named("Barray") = Barray
    );
}

// compute all outer products at once for low-rank model
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_outer_products_batch(const arma::mat& U) {
    int m = U.n_rows;
    int r = U.n_cols;
    arma::cube outer_prods(m, m, r);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int k = 0; k < r; k++) {
        outer_prods.slice(k) = U.col(k) * U.col(k).t();
    }
    
    return outer_prods;
}

// fast computation of all A matrices for low-rank model using preomputed outer products
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_A_lowrank_batch(const arma::cube& outer_prods,
                                  const arma::mat& alpha,
                                  int Tt) {
    int m = outer_prods.n_rows;
    int r = outer_prods.n_slices;
    arma::cube Aarray(m, m, Tt);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int t = 0; t < Tt; t++) {
        arma::mat A_t(m, m, fill::zeros);
        arma::vec alpha_t = alpha.col(t);
        
        // use precomputed outer products
        for(int k = 0; k < r; k++) {
            A_t += alpha_t(k) * outer_prods.slice(k);
        }
        
        Aarray.slice(t) = A_t;
    }
    
    return Aarray;
}