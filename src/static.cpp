#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "thread_control.h"

using namespace Rcpp;
using namespace arma;

//' Reshape 4D array to 3D for C++ processing
//' @description Efficiently reshape Z from m x m x p x Tt to m x m x (p*Tt)
//' @param Z_4d Input 4D array as R array
//' @param m Number of nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return 3D cube for C++ processing
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube reshape_Z_to_cube(const NumericVector& Z_4d, int m, int p, int Tt) {
    // create output cube 
    arma::cube Z_cube(m, m, p * Tt);
    
    // direct memory access with better cache locality
    const double* Z_ptr = Z_4d.begin();
    
    for (int slice = 0; slice < p * Tt; slice++) {
        int j = slice / Tt;
        int t = slice % Tt;
        
        double* cube_slice = Z_cube.slice_memptr(slice);
        const double* z_offset = Z_ptr + (t * m * m * p + j * m * m);
        
        // copy entire slice at once for better memory bandwidth
        std::memcpy(cube_slice, z_offset, m * m * sizeof(double));
    }
    
    return Z_cube;
}

//' Compute SSE for diagonal elements efficiently
//' @description Compute sum of squared errors for diagonal elements of B matrices
//' @param B_list List of B matrices
//' @param K Number of matrices
//' @return SSE value
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_diagonal_sse(const List& B_list, int K) {
    double sse = 0.0;
    
    for (int k = 0; k < K; k++) {
        arma::mat B = as<arma::mat>(B_list[k]);
        arma::vec diag_B = B.diag();
        arma::vec diff = diag_B - 1.0;
        sse += arma::dot(diff, diff);
    }
    
    return sse;
}


//' Compute sum of squared deviations from identity for A or B arrays
//' @description Efficiently compute sum of squared deviations from the identity matrix 
//'   for each matrix slice of A or B from time t = 2 to Tt
//' @param ABarray 3D array of A or B matrices (m x m x Tt)
//' @param m Number of nodes
//' @param Tt Number of time points
//' @return Sum of squared deviations
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_deviation_sum(const arma::cube& ABarray, int m, int Tt) {
    double sum_sq = 0.0;
    arma::mat I = arma::eye(m, m);
    
    // start from t=2 (index 1) as per R code
    for (int t = 1; t < Tt; t++) {
        arma::mat diff = ABarray.slice(t) - I;
        sum_sq += arma::accu(diff % diff);  // element-wise square and sum
    }
    
    return sum_sq;
}

//' Compute mean M for static model
//' @description Efficiently compute mean of Z across time for each relation
//' @param Z_flat Flattened Z array (m*m x p*Tt)
//' @param m Number of nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return M array (m x m x p)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_M_static(const arma::mat& Z_flat, int m, int p, int Tt) {
    set_dbn_threads(); // set threads from R options
    arma::cube M(m, m, p);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int j = 0; j < p; j++) {
        
        arma::vec ones_vec(Tt, arma::fill::ones);
        arma::mat Z_j_cols(m * m, Tt);
        
        // extract columns for this relation
        for (int t = 0; t < Tt; t++) {
            Z_j_cols.col(t) = Z_flat.col(j * Tt + t);
        }
        
        // compute sum using mat vector mult
        arma::vec sum_vec = Z_j_cols * ones_vec;
        M.slice(j) = arma::reshape(sum_vec / Tt, m, m);
    }
    
    return M;
}

//' Compute residual sum of squares for static model
//' @description Compute sum((Z - M)^2) across all times and relations
//' @param Z_4d Original Z array (m x m x p x Tt) as flattened vector
//' @param M Mean array (m x m x p)
//' @param m Number of nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return Sum of squared residuals
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_rss_static(const NumericVector& Z_4d, const arma::cube& M,
                         int m, int p, int Tt) {
    double rss = 0.0;
    
    for (int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);
        
        for (int t = 0; t < Tt; t++) {
            // extract Z[,,j,t]
            for (int i = 0; i < m; i++) {
                for (int k = 0; k < m; k++) {
                    int idx_4d = i + k * m + j * m * m + t * m * m * p;
                    double resid = Z_4d[idx_4d] - M_j(i, k);
                    rss += resid * resid;
                }
            }
        }
    }
    
    return rss;
}

//' Reshape for large networks with parallelization
//' @description Parallel version of reshape_Z_to_cube for large networks
//' @param Z_4d Input 4D array as R array
//' @param m Number of nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return 3D cube for C++ processing
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube reshape_Z_to_cube_parallel(const NumericVector& Z_4d, int m, int p, int Tt) {
    set_dbn_threads(); // Set threads from R options
    arma::cube Z_cube(m, m, p * Tt);
    
    // direct memory access with better cache locality
    const double* Z_ptr = Z_4d.begin();
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 4)
    #endif
    for (int slice = 0; slice < p * Tt; slice++) {
        int j = slice / Tt;
        int t = slice % Tt;
        
        double* cube_slice = Z_cube.slice_memptr(slice);
        const double* z_offset = Z_ptr + (t * m * m * p + j * m * m);
        
        // copy entire slice at once for better memory bandwidth
        std::memcpy(cube_slice, z_offset, m * m * sizeof(double));
    }
    
    return Z_cube;
}

//' M computation with blocked summation for numerical stability
//' @description Compute mean with blocked summation for large networks
//' @param Z_flat Flattened Z array (m*m x p*Tt)
//' @param m Number of nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return M array (m x m x p)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_M_static_blocked(const arma::mat& Z_flat, int m, int p, int Tt) {
    set_dbn_threads(); // set threads from R options
    arma::cube M(m, m, p, arma::fill::zeros);
    const int block_size = 64; // Cache-friendly block size
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int j = 0; j < p; j++) {
        arma::mat Z_j_sum(m, m, arma::fill::zeros);
        
        // blocked summation for better cache usage
        for (int t_block = 0; t_block < Tt; t_block += block_size) {
            int t_end = std::min(t_block + block_size, Tt);
            
            for (int t = t_block; t < t_end; t++) {
                int col_idx = j * Tt + t;
                arma::vec z_vec = Z_flat.col(col_idx);
                Z_j_sum += arma::reshape(z_vec, m, m);
            }
        }
        
        M.slice(j) = Z_j_sum / Tt;
    }
    
    return M;
}

//' RSS computation with parallel reduction
//' @description Parallel computation of residual sum of squares
//' @param Z_cube Z array as cube (m x m x p*Tt)
//' @param M Mean array (m x m x p)
//' @param m Number of nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return Sum of squared residuals
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_rss_static_parallel(const arma::cube& Z_cube, const arma::cube& M,
                                  int m, int p, int Tt) {
    double rss = 0.0;
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:rss)
    #endif
    for (int j = 0; j < p; j++) {
        double local_rss = 0.0;
        arma::mat M_j = M.slice(j);
        
        for (int t = 0; t < Tt; t++) {
            int idx = j * Tt + t;
            arma::mat diff = Z_cube.slice(idx) - M_j;
            local_rss += arma::accu(diff % diff);
        }
        rss += local_rss;
    }
    
    return rss;
}

//' Cache-efficient B update for static model
//' @description B update using tiled matrix operations
//' @param Z_cube Z array as cube
//' @param M Mean array
//' @param s2 Observation variance
//' @param t2 Prior variance
//' @param m Number of nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return Updated B matrix
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat update_B_static_tiled(const arma::cube& Z_cube, const arma::cube& M,
                                double s2, double t2, int m, int p, int Tt) {
    set_dbn_threads(); // set threads from R options
    
    // preallocate all workspace to avoid allocation in hot path
    arma::mat XtX(m, m, arma::fill::zeros);
    arma::mat XtY(m, m, arma::fill::zeros);
    
    arma::mat M_reshaped(m * p, m);
    for (int j = 0; j < p; j++) {
        M_reshaped.rows(j * m, (j + 1) * m - 1) = M.slice(j);
    }
    
    // compute XtX = sum_j (Tt * M_j' * M_j) 
    arma::mat temp = M_reshaped.t() * M_reshaped;
    XtX = Tt * temp;
    
    // compute XtY using batched operations 
    #ifdef _OPENMP
    #pragma omp parallel
    {
        arma::mat local_XtY(m, m, arma::fill::zeros);
        
        #pragma omp for schedule(static)
        for (int j = 0; j < p; j++) {
            arma::mat M_j = M.slice(j);
            arma::mat sum_Z(m, m, arma::fill::zeros);
            
            // Sum Z matrices for this relation
            for (int t = 0; t < Tt; t++) {
                sum_Z += Z_cube.slice(j * Tt + t);
            }
            
            // Single BLAS call per relation
            local_XtY += sum_Z * M_j.t();
        }
        
        #pragma omp critical
        XtY += local_XtY;
    }
    #else
    for (int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);
        arma::mat sum_Z(m, m, arma::fill::zeros);
        
        for (int t = 0; t < Tt; t++) {
            sum_Z += Z_cube.slice(j * Tt + t);
        }
        
        XtY += sum_Z * M_j.t();
    }
    #endif
    
    // posterior calculations with numerical stability
    double lambda = 1.0 / t2;
    double gamma = 1.0 / s2;
    
    // add prior precision to diagonal
    XtX.diag() += (lambda / gamma) * arma::ones(m);
    
    // efficient posterior covariance using Cholesky
    arma::mat L_post;
    bool chol_success = arma::chol(L_post, gamma * XtX, "lower");
    
    if (!chol_success) {
        // fallback to SVD for ill-conditioned matrices
        arma::mat U, V;
        arma::vec s;
        arma::svd_econ(U, s, V, XtX);
        
        // regularize small eigenvalues
        double eps = 1e-8 * s.max();
        s = arma::max(s, eps * arma::ones(s.n_elem));
        
        arma::mat post_cov = V * arma::diagmat(1.0 / (gamma * s + lambda)) * V.t();
        arma::mat post_mean = post_cov * (gamma * XtY + lambda * arma::eye(m, m));
        
        // generate sample
        arma::mat noise = arma::randn(m, m);
        return post_mean + V * arma::diagmat(arma::sqrt(1.0 / (gamma * s + lambda))) * V.t() * noise;
    }
    
    // solve for posterior mean efficiently
    arma::mat post_mean_part = arma::solve(arma::trimatl(L_post), gamma * XtY + lambda * arma::eye(m, m));
    arma::mat post_mean = arma::solve(arma::trimatu(L_post.t()), post_mean_part);
    
    // sample from posterior
    arma::mat Z_sample = arma::randn(m, m);
    arma::mat B_new = post_mean + arma::solve(arma::trimatu(L_post.t()), Z_sample) / sqrt(gamma);
    
    return B_new;
}

