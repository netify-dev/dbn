#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// Fast tensor product for static model B updates
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_XB_tensor(const arma::cube& X, 
                                  const List& B,
                                  int m, int p, int n) {
    // compute X tensor B using batched matrix operations
    // X is m x m x p x n, B is list of 3 matrices
    
    arma::mat B1 = as<arma::mat>(B[0]);
    arma::mat B2 = as<arma::mat>(B[1]);
    arma::mat B3 = as<arma::mat>(B[2]);
    
    arma::cube XB(m * m * p, n, 1);
    
    // use OpenMP for parallel computation if available
    #ifdef _OPENMP
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    #endif
    for(int t = 0; t < n; t++) {
        for(int r = 0; r < p; r++) {
            arma::mat X_slice(m, m);
            // Extract X[,,r,t]
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    X_slice(i, j) = X(i + j*m + r*m*m + t*m*m*p);
                }
            }
            
            // compute tensor product 
            arma::mat result = B1 * X_slice * B2.t();
            result *= B3(r, r); // Only diagonal of B3 matters for relation r
            
            // store result
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    XB(i + j*m + r*m*m, t, 0) = result(i, j);
                }
            }
        }
    }
    
    // reshape to match original dimensions
    XB.reshape(m, m, p*n);
    return XB;
}

// Static model B update - returns single matrix
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat update_B_static(const arma::cube& Z, const arma::cube& M,
                              double s2, double t2, int m, int p, int n) {
    // for static model, B is a single m x m matrix
    // prior: B ~ N(I, t2*I)
    // likelihood: Z[,,j,t] ~ N(B * M[,,j], s2*I)
    
    // preallocate all workspace to avoid allocation in hot path
    arma::mat XtX(m, m, arma::fill::zeros);
    arma::mat XtY(m, m, arma::fill::zeros);
    
    if (m >= 50) {
        arma::mat M_reshaped(m * p, m);
        for (int j = 0; j < p; j++) {
            M_reshaped.rows(j * m, (j + 1) * m - 1) = M.slice(j);
        }
        
        // compute XtX = sum_j (n * M_j' * M_j) 
        arma::mat temp = M_reshaped.t() * M_reshaped;
        XtX = n * temp;
        
        // compute XtY using batched operations
        #ifdef _OPENMP
        #pragma omp parallel
        {
            arma::mat local_XtY(m, m, arma::fill::zeros);
            
            #pragma omp for schedule(static)
            for (int j = 0; j < p; j++) {
                arma::mat M_j = M.slice(j);
                arma::mat sum_Z(m, m, arma::fill::zeros);
                
                // sum Z matrices for this relation
                for (int t = 0; t < n; t++) {
                    sum_Z += Z.slice(j * n + t);
                }
                
                local_XtY += sum_Z * M_j.t();
            }
            
            #pragma omp critical
            XtY += local_XtY;
        }
        #else
        for (int j = 0; j < p; j++) {
            arma::mat M_j = M.slice(j);
            arma::mat sum_Z(m, m, arma::fill::zeros);
            
            for (int t = 0; t < n; t++) {
                sum_Z += Z.slice(j * n + t);
            }
            
            XtY += sum_Z * M_j.t();
        }
        #endif
        
    } else {
        // small matrix path
        for (int j = 0; j < p; j++) {
            arma::mat M_j = M.slice(j);
            XtX += n * (M_j.t() * M_j);
            
            for (int t = 0; t < n; t++) {
                XtY += Z.slice(j * n + t) * M_j.t();
            }
        }
    }
    
    // posterior calculations with numerical stability
    double lambda = 1.0 / t2;
    double gamma = 1.0 / s2;
    
    // add prior precision to diagonal
    XtX.diag() += (lambda / gamma) * arma::ones(m);
    
    // posterior covariance using Cholesky
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
        
        // gen samp
        arma::mat noise = arma::randn(m, m);
        return post_mean + V * arma::diagmat(arma::sqrt(1.0 / (gamma * s + lambda))) * V.t() * noise;
    }
    
    // solve for posterior mean
    arma::mat post_mean_part = arma::solve(arma::trimatl(L_post), gamma * XtY + lambda * arma::eye(m, m));
    arma::mat post_mean = arma::solve(arma::trimatu(L_post.t()), post_mean_part);
    
    // sample from posterior
    arma::mat Z_sample = arma::randn(m, m);
    arma::mat B_new = post_mean + arma::solve(arma::trimatu(L_post.t()), Z_sample) / sqrt(gamma);
    
    return B_new;
}

// Fast broadcast and residual computation for static model
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube broadcast_M_and_compute_EZ(const arma::cube& M, double s2,
                                      int m, int p, int Tt) {
    // M is m x m x p, broadcast to m x m x p x Tt and add noise
    arma::cube EZ(m, m, p * Tt);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);
        for(int t = 0; t < Tt; t++) {
            int idx = j * Tt + t;
            EZ.slice(idx) = M_j + randn(m, m) * sqrt(s2);
        }
    }
    
    return EZ;
}

// Vectorized s2 update for static model
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_s2_update(const List& Z_field,
                             const arma::cube& M,
                             int m, int p, int Tt,
                             double a_prior, double b_prior) {
    // compute residual sum of squares 
    double rss = 0.0;
    int n_obs = 0;
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:rss,n_obs)
    #endif
    for(int j = 0; j < p; j++) {
        arma::cube Z_j = as<arma::cube>(Z_field[j]);
        arma::mat M_j = M.slice(j);
        
        for(int t = 0; t < Tt; t++) {
            arma::mat diff = Z_j.slice(t) - M_j;
            arma::uvec finite_idx = find_finite(diff);
            
            if(finite_idx.n_elem > 0) {
                arma::vec diff_vec = diff.elem(finite_idx);
                rss += dot(diff_vec, diff_vec);
                n_obs += finite_idx.n_elem;
            }
        }
    }
    
    // sample from inverse gamma
    double shape = (n_obs + a_prior) / 2.0;
    double scale = (rss + b_prior) / 2.0;
    
    return scale / randg(distr_param(shape, 1.0));
}

// Batch update for multiple variance parameters
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec update_variances_batch(const arma::vec& sum_squares,
                                const arma::vec& counts,
                                const arma::vec& a_priors,
                                const arma::vec& b_priors) {
    int n = sum_squares.n_elem;
    arma::vec variances(n);
    
    for(int i = 0; i < n; i++) {
        double shape = (counts(i) + a_priors(i)) / 2.0;
        double scale = (sum_squares(i) + b_priors(i)) / 2.0;
        variances(i) = scale / randg(distr_param(shape, 1.0));
    }
    
    return variances;
}

// Fast matrix stabilization for numerical stability
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat stabilize_matrix(const arma::mat& M, double min_eig = 1e-6) {
    arma::vec eigval;
    arma::mat eigvec;
    
    eig_sym(eigval, eigvec, M);
    
    // ensure all eigenvalues are at least min_eig
    eigval = clamp(eigval, min_eig, datum::inf);
    
    return eigvec * diagmat(eigval) * eigvec.t();
}

// Z update for dynamic model
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube update_Z_dynamic(const arma::cube& R, const arma::cube& Z_current,
                                const arma::cube& Theta, const arma::cube& M,
                                const List& IR, int m, int p, int Tt) {
    arma::cube Z_new = Z_current;
    
    // compute EZ = Theta + M
    arma::cube EZ(m, m, p * Tt);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);
        for(int t = 0; t < Tt; t++) {
            int idx = j * Tt + t;
            // Extract Theta[,,j,t]
            arma::mat Theta_jt(m, m);
            for(int i = 0; i < m; i++) {
                for(int k = 0; k < m; k++) {
                    Theta_jt(i, k) = Theta(i + k*m + j*m*m + t*m*m*p);
                }
            }
            EZ.slice(idx) = Theta_jt + M_j;
        }
    }
    
    return Z_new;
}

// Pre-compute and cache matrix products
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List precompute_products(const arma::cube& Aarray, const arma::cube& Barray,
                        int m, int Tt) {
    List products(Tt);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int t = 0; t < Tt; t++) {
        arma::mat AtBt = Aarray.slice(t).t() * Barray.slice(t);
        products[t] = AtBt;
    }
    
    return products;
}