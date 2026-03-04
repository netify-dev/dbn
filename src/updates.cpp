#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// tensor product for static model B updates
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_XB_tensor(const arma::cube& X, 
                                  const List& B,
                                  int m, int p, int n) {
    // X tensor B via batched matrix operations
    // X is m x m x p x n, B is a list of 3 matrices
    
    arma::mat B1 = as<arma::mat>(B[0]);
    arma::mat B2 = as<arma::mat>(B[1]);
    arma::mat B3 = as<arma::mat>(B[2]);
    
    arma::cube XB(m * m * p, n, 1);
    
    // parallelize over time if OpenMP available
    #ifdef _OPENMP
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    #endif
    for(int t = 0; t < n; t++) {
        for(int r = 0; r < p; r++) {
            arma::mat X_slice(m, m);
            // extract X[,,r,t]
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    X_slice(i, j) = X(i + j*m + r*m*m + t*m*m*p);
                }
            }
            
            // tensor product B1 * X * B2' scaled by B3
            arma::mat result = B1 * X_slice * B2.t();
            result *= B3(r, r); // only diagonal of B3 matters for relation r
            
            // pack result back into flat storage
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    XB(i + j*m + r*m*m, t, 0) = result(i, j);
                }
            }
        }
    }
    
    // reshape to original dimensions
    XB.reshape(m, m, p*n);
    return XB;
}

// static model B[[1]] (sender effects) update, returns n_row x n_row matrix
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat update_B_static(const arma::cube& Z, const arma::cube& M,
                              double s2, double t2, int n_row, int n_col, int p, int n) {
    // sender effects: n_row x n_row
    // model: Z[,,j,t] = B * M[,,j] + E
    // column-wise regression: z_c = B * m_c
    // XtX = sum M_j * M_j' (n_row x n_row), XtY = sum Z_jt * M_j' (n_row x n_row)

    arma::mat XtX(n_row, n_row, arma::fill::zeros);
    arma::mat XtY(n_row, n_row, arma::fill::zeros);

    // accumulate sufficient statistics over relations
    for (int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);
        // (n_row x n_col) * (n_col x n_row) = n_row x n_row
        XtX += n * (M_j * M_j.t());

        arma::mat sum_Z(n_row, n_col, arma::fill::zeros);
        for (int t = 0; t < n; t++) {
            sum_Z += Z.slice(j * n + t);
        }
        // (n_row x n_col) * (n_col x n_row) = n_row x n_row
        XtY += sum_Z * M_j.t();
    }

    // posterior precision and mean
    double lambda = 1.0 / t2;
    double gamma = 1.0 / s2;

    XtX.diag() += (lambda / gamma) * arma::ones(n_row);

    arma::mat gXtX = gamma * XtX;
    gXtX = 0.5 * (gXtX + gXtX.t());
    arma::mat L_post;
    bool chol_success = arma::chol(L_post, gXtX, "lower");

    // try cholesky-based solve if factorization succeeded
    bool solve_ok = false;
    if (chol_success) {
        arma::mat post_mean_part, post_mean, noise_part;
        solve_ok = arma::solve(post_mean_part, arma::trimatl(L_post),
                               gamma * XtY + lambda * arma::eye(n_row, n_row));
        if (solve_ok) {
            solve_ok = arma::solve(post_mean, arma::trimatu(L_post.t()), post_mean_part);
        }
        if (solve_ok) {
            arma::mat Z_sample = arma::randn(n_row, n_row);
            solve_ok = arma::solve(noise_part, arma::trimatu(L_post.t()), Z_sample);
        }
        if (solve_ok) {
            return post_mean + noise_part / sqrt(gamma);
        }
    }

    // SVD fallback when cholesky or triangular solve fails
    arma::mat U, V;
    arma::vec s;
    bool svd_ok = arma::svd_econ(U, s, V, XtX);

    if (!svd_ok || s.n_elem == 0) {
        return arma::eye(n_row, n_row) + 0.01 * arma::randn(n_row, n_row);
    }

    double eps = 1e-8 * s.max();
    s = arma::max(s, eps * arma::ones(s.n_elem));

    arma::mat post_cov = V * arma::diagmat(1.0 / (gamma * s + lambda)) * V.t();
    arma::mat post_mean_svd = post_cov * (gamma * XtY + lambda * arma::eye(n_row, n_row));

    arma::mat noise = arma::randn(n_row, n_row);
    return post_mean_svd + V * arma::diagmat(arma::sqrt(1.0 / (gamma * s + lambda))) * V.t() * noise;
}

// broadcast M across time and add noise for static model
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube broadcast_M_and_compute_EZ(const arma::cube& M, double s2,
                                      int n_row, int n_col, int p, int Tt) {
    // M is n_row x n_col x p; replicate across time with gaussian noise
    arma::cube EZ(n_row, n_col, p * Tt);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);
        for(int t = 0; t < Tt; t++) {
            int idx = j * Tt + t;
            EZ.slice(idx) = M_j + randn(n_row, n_col) * sqrt(s2);
        }
    }

    return EZ;
}

// observation variance update for static model
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_s2_update(const List& Z_field,
                             const arma::cube& M,
                             int m, int p, int Tt,
                             double a_prior, double b_prior) {
    // RSS across all relations and time points
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
    
    // sample from inverse-gamma posterior
    double shape = (n_obs + a_prior) / 2.0;
    double scale = (rss + b_prior) / 2.0;
    
    return scale / randg(distr_param(shape, 1.0));
}

// batch update for multiple variance parameters
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

// stabilize a symmetric matrix by clamping eigenvalues
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat stabilize_matrix(const arma::mat& M, double min_eig = 1e-6) {
    arma::vec eigval;
    arma::mat eigvec;
    
    eig_sym(eigval, eigvec, M);
    
    // floor eigenvalues at min_eig
    eigval = clamp(eigval, min_eig, datum::inf);
    
    return eigvec * diagmat(eigval) * eigvec.t();
}

// z update for dynamic model
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube update_Z_dynamic(const arma::cube& R, const arma::cube& Z_current,
                                const arma::cube& Theta, const arma::cube& M,
                                const List& IR, int n_row, int n_col, int p, int Tt) {
    int nc = n_row * n_col;
    arma::cube Z_new = Z_current;

    // expected z = Theta + M
    arma::cube EZ(n_row, n_col, p * Tt);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);
        for(int t = 0; t < Tt; t++) {
            int idx = j * Tt + t;
            // extract Theta[,,j,t]
            arma::mat Theta_jt(n_row, n_col);
            for(int i = 0; i < n_row; i++) {
                for(int k = 0; k < n_col; k++) {
                    Theta_jt(i, k) = Theta(i + k*n_row + j*nc + t*nc*p);
                }
            }
            EZ.slice(idx) = Theta_jt + M_j;
        }
    }

    return Z_new;
}

// precompute and cache A'*B products across time
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