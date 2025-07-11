#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// declarations for functions defined in other files
arma::mat ffbs_dlm_cpp(const List& y, const List& Flist, const arma::mat& V, 
                       const arma::mat& W, const arma::vec& m0, const arma::mat& C0, 
                       bool ar1, double rho);

arma::cube ffbs_bilinear(const arma::cube& Z, const arma::mat& mu, 
                             const arma::cube& A_array, const arma::cube& B_array,
                             double sigma2_proc, double sigma2_obs);

// A_t = U * diag(alpha_t) * U'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat compute_A_lowrank(const arma::mat& U, const arma::vec& alpha_t) {
    int m = U.n_rows;
    int r = U.n_cols;
    
    if(alpha_t.n_elem == 1) {
        // scalar: A = alpha * U * U'
        return alpha_t(0) * U * U.t();
    } else {
        // vector: A = U * diag(alpha) * U'
        arma::mat A(m, m);
        
        // computation using outer products
        A.zeros();
        for(int k = 0; k < r; k++) {
            A += alpha_t(k) * U.col(k) * U.col(k).t();
        }
        
        return A;
    }
}

// Batch computation of all A matrices
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_all_A_lowrank(const arma::mat& U, const arma::mat& alpha,
                                 int Tt) {
    int m = U.n_rows;
    int r = U.n_cols;
    arma::cube Aarray(m, m, Tt);
    
    // for large m, use blocked computation for better cache usage
    if (m > 100) {
        const int block_size = 32;  
        
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for(int t = 0; t < Tt; t++) {
            arma::vec alpha_t = alpha.col(t);
            arma::mat A_t(m, m, fill::zeros);
            
            // blocked computation for better cache usage
            for(int i_block = 0; i_block < m; i_block += block_size) {
                int i_end = std::min(i_block + block_size, m);
                
                for(int j_block = 0; j_block < m; j_block += block_size) {
                    int j_end = std::min(j_block + block_size, m);
                    
                    // compute block contribution
                    arma::mat block_contrib(i_end - i_block, j_end - j_block, fill::zeros);
                    
                    for(int k = 0; k < r; k++) {
                        arma::vec u_k_i = U.submat(i_block, k, i_end-1, k);
                        arma::vec u_k_j = U.submat(j_block, k, j_end-1, k);
                        block_contrib += alpha_t(k) * (u_k_i * u_k_j.t());
                    }
                    
                    A_t.submat(i_block, j_block, i_end-1, j_end-1) = block_contrib;
                }
            }
            
            Aarray.slice(t) = A_t;
        }
    } else {
        // implementation for small mats
        // pre-compute all outer products once
        arma::cube U_outer(m, m, r);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int k = 0; k < r; k++) {
            U_outer.slice(k) = U.col(k) * U.col(k).t();
        }
        
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int t = 0; t < Tt; t++) {
            arma::vec alpha_t = alpha.col(t);
            arma::mat A_t(m, m, fill::zeros);
            
            // linear combination
            for(int k = 0; k < r; k++) {
                A_t += alpha_t(k) * U_outer.slice(k);
            }
            
            Aarray.slice(t) = A_t;
        }
    }
    
    return Aarray;
}

// Build design matrix F for alpha update (single time/relation)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat build_F_alpha_vectorized(const arma::mat& U, 
                                   const arma::mat& Theta_prev,
                                   const arma::mat& B_t,
                                   int p, int Tt) {
    int m = U.n_rows;
    int r = U.n_cols;
    int m2 = m * m;
    
    // for a single time point and relation
    arma::mat F_block(m2, r);
    arma::mat BTheta = B_t * Theta_prev.t();
    
    // 
    for(int k = 0; k < r; k++) {
        arma::vec u_k = U.col(k);
        arma::mat Sk = u_k * u_k.t();
        arma::mat prod = Sk * BTheta;
        F_block.col(k) = vectorise(prod);
    }
    
    return F_block;
}

// Build observation vector y for alpha update
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec build_y_alpha_vectorized(const arma::cube& Theta_flat,
                                   const arma::cube& M,
                                   int m, int p, int Tt) {
    int m2 = m * m;
    arma::vec y_all(m2 * p * (Tt - 1));
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int t = 1; t < Tt; t++) {
        for(int rel = 0; rel < p; rel++) {
            int idx = rel * Tt + t;
            arma::mat Theta_t = Theta_flat.slice(idx);
            arma::mat M_rel = M.slice(rel);
            
            // vectorize Theta - M
            arma::vec y_block = vectorise(Theta_t - M_rel);
            
            // store in correct position
            int start = ((t - 1) * p + rel) * m2;
            y_all.subvec(start, start + m2 - 1) = y_block;
        }
    }
    
    return y_all;
}

// Log-likelihood for U update
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double loglik_U(const arma::mat& U, const arma::mat& alpha,
                    const arma::cube& Theta_avg, const arma::cube& Barray,
                    double sigma2) {
    int Tt = alpha.n_cols;
    double total_ss = 0.0;
    
    // precompute U'U
    arma::mat UtU = U.t() * U;
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:total_ss)
    #endif
    for(int t = 1; t < Tt; t++) {
        arma::vec alpha_t = alpha.col(t);
        arma::mat B_t = Barray.slice(t);
        arma::mat Theta_curr = Theta_avg.slice(t);
        arma::mat Theta_prev = Theta_avg.slice(t-1);
        
        // compute predicted value
        // A_t * Theta_{t-1} * B_t' = U * diag(alpha_t) * U' * Theta_{t-1} * B_t'
        arma::mat temp = U.t() * Theta_prev * B_t.t();
        arma::mat pred = U * (arma::diagmat(alpha_t) * temp);
        
        // compute residual sum of squares
        arma::mat resid = Theta_curr - pred;
        total_ss += arma::accu(resid % resid);
    }
    
    return -0.5 * total_ss / sigma2;
}

// Cayley transform for Stiefel manifold update 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat cayley_transform(const arma::mat& U, const arma::mat& W, double epsilon) {
    int m = U.n_rows;
    int r = U.n_cols;
    
    // for small epsilon, use Woodbury formula again
    if (std::abs(epsilon) < 0.1 && m > 50) {
        double eps_half = 0.5 * epsilon;
        double eps2_quarter = 0.25 * epsilon * epsilon;
        
        // Compute W^2
        arma::mat W2 = W * W;
        
        // 
        arma::mat inner = arma::eye(m, m) - eps2_quarter * W2;
        
        // solve using Cholesky if possible
        arma::mat inner_inv;
        arma::mat L;
        bool chol_success = arma::chol(L, inner, "lower");
        
        if (chol_success) {
            // cholesky decomp
            arma::mat temp = arma::solve(arma::trimatl(L), arma::eye(m, m));
            inner_inv = arma::solve(arma::trimatu(L.t()), temp);
        } else {
            // fallback to standard inverse
            inner_inv = arma::inv_sympd(inner);
        }
        
        // 
        arma::mat cayley_mat = arma::eye(m, m) + eps_half * W * inner_inv;
        
        // apply to U using matrix multiplication
        return cayley_mat * U;
    } else {
        // standard implementation for larger epsilon or smaller matrices
        arma::mat I_m = eye(m, m);
        arma::mat A = I_m - 0.5 * epsilon * W;
        arma::mat B = I_m + 0.5 * epsilon * W;
        
        // Solve AX = B for X, then multiply by U
        arma::mat X;
        bool success = solve(X, A, B);
        
        if(!success) {
            // Fall back to simple update
            return U + epsilon * W * U;
        }
        
        return X * U;
    }
}

// Generate skew-symmetric proposal for U update
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat generate_skew_proposal(int m, double norm_cap) {
    arma::mat W = randn(m, m);
    W = W - W.t(); // make skew-symmetric
    
    // 
    double W_norm = norm(W, "fro");
    if(W_norm > norm_cap) {
        W *= (norm_cap / W_norm);
    }
    
    return W;
}

// Batch B update using parallel FFBS
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube update_B_lowrank_batch(const arma::cube& Theta_flat,
                                 const arma::mat& U,
                                 const arma::mat& alpha,
                                 double sigma2_proc,
                                 double tau_B2,
                                 bool ar1, double rho,
                                 int m, int p, int Tt) {
    arma::cube Barray(m, m, Tt);
    
    // build all a mats
    arma::cube Aarray = compute_all_A_lowrank(U, alpha, Tt);
    
    // Update each column of B in parallel
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int k = 0; k < m; k++) {
        // build y and F for this column
        arma::vec y_all;
        arma::mat F_all;
        
        if(p == 1) {
            // single relation case
            y_all = arma::vec((Tt - 1) * m);
            F_all = arma::mat((Tt - 1) * m, m);
            
            for(int t = 1; t < Tt; t++) {
                arma::mat Theta_t = Theta_flat.slice(t);
                arma::mat Theta_prev = Theta_flat.slice(t - 1);
                
                y_all.subvec((t-1)*m, t*m - 1) = Theta_t.col(k);
                F_all.rows((t-1)*m, t*m - 1) = Aarray.slice(t) * Theta_prev;
            }
        } else {
            // multiple relations
            y_all = arma::vec((Tt - 1) * m * p);
            F_all = arma::mat((Tt - 1) * m * p, m);
            
            for(int t = 1; t < Tt; t++) {
                for(int rel = 0; rel < p; rel++) {
                    int idx = rel * Tt + t;
                    int idx_prev = rel * Tt + t - 1;
                    
                    arma::mat Theta_t = Theta_flat.slice(idx);
                    arma::mat Theta_prev = Theta_flat.slice(idx_prev);
                    
                    int start = ((t-1)*p + rel) * m;
                    y_all.subvec(start, start + m - 1) = Theta_t.col(k);
                    F_all.rows(start, start + m - 1) = Aarray.slice(t) * Theta_prev;
                }
            }
        }
        
        // add dummy for t=1
        int V_dim = (p == 1) ? m : m * p;
        arma::vec y_full = join_vert(arma::zeros(V_dim), y_all);
        arma::mat F_full = join_vert(arma::zeros(V_dim, m), F_all);
        
        // convert to list format for ffbs_dlm_cpp
        List y_list(Tt);
        List F_list(Tt);
        
        for(int t = 0; t < Tt; t++) {
            y_list[t] = y_full.subvec(t * V_dim, (t + 1) * V_dim - 1);
            F_list[t] = F_full.rows(t * V_dim, (t + 1) * V_dim - 1);
        }
        
        // 
        arma::mat b_path = ffbs_dlm_cpp(y_list, F_list,
                                       eye(V_dim, V_dim) * sigma2_proc,
                                       eye(m, m) * tau_B2,
                                       zeros(m), eye(m, m) * tau_B2,
                                       ar1, rho);
        
        // store result - b_path is m x Tt, we need to store in column k for all t
        for(int t = 0; t < Tt; t++) {
            Barray.slice(t).col(k) = b_path.col(t);
        }
    }
    
    return Barray;
}

// Compute all residuals for sigma2 update
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_lowrank_residuals(const arma::cube& Theta_flat,
                                const arma::cube& M,
                                const arma::mat& U,
                                const arma::mat& alpha,
                                const arma::cube& Barray,
                                int m, int p, int Tt) {
    double rss = 0.0;
    
    // precompute all A mats
    arma::cube Aarray = compute_all_A_lowrank(U, alpha, Tt);
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:rss)
    #endif
    for(int rel = 0; rel < p; rel++) {
        for(int t = 1; t < Tt; t++) {
            int idx = rel * Tt + t;
            int idx_prev = rel * Tt + t - 1;
            
            arma::mat Theta_t = Theta_flat.slice(idx);
            arma::mat Theta_prev = Theta_flat.slice(idx_prev);
            arma::mat M_rel = M.slice(rel);
            
            // Prediction
            arma::mat pred = M_rel + Aarray.slice(t) * (Theta_prev - M_rel) * Barray.slice(t).t();
            
            // Residual
            arma::mat resid = Theta_t - pred;
            rss += accu(resid % resid);
        }
    }
    
    return rss;
}

// batch residual computation using pre-computed arrays
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_sigma2_lowrank_batch(const arma::cube& Theta,
                                   const arma::cube& Aarray,
                                   const arma::cube& Barray,
                                   int m, int p, int Tt) {
    double total_ss = 0.0;
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:total_ss)
    #endif
    for(int t = 1; t < Tt; t++) {
        arma::mat A_t = Aarray.slice(t);
        arma::mat B_t = Barray.slice(t);
        arma::mat B_t_trans = B_t.t();
        
        for(int rel = 0; rel < p; rel++) {
            // residual computation without temporary matrices
            arma::mat pred = A_t * Theta.slice(rel * Tt + t - 1) * B_t_trans;
            arma::mat diff = Theta.slice(rel * Tt + t) - pred;
            total_ss += arma::dot(diff, diff);
        }
    }
    
    return total_ss;
}

// sigma2 computation
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_sigma2_simd(const arma::cube& Theta,
                          const arma::mat& U,
                          const arma::mat& alpha,
                          const arma::cube& Barray,
                          int m, int p, int Tt, int r) {
    double total_ss = 0.0;
    
    // pre-compute U*alpha for all time points
    arma::mat Ualpha(m, r * Tt);
    for(int t = 0; t < Tt; t++) {
        for(int k = 0; k < r; k++) {
            Ualpha.col(t * r + k) = alpha(k, t) * U.col(k);
        }
    }
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:total_ss)
    #endif
    for(int t = 1; t < Tt; t++) {
        arma::mat B_t = Barray.slice(t);
        arma::mat B_t_trans = B_t.t();
        
        // compute A_t using pre-computed values
        arma::mat A_t(m, m, arma::fill::zeros);
        for(int k = 0; k < r; k++) {
            A_t += Ualpha.col(t * r + k) * U.col(k).t();
        }
        
        for(int rel = 0; rel < p; rel++) {
            int curr_idx = rel * Tt + t;
            int prev_idx = rel * Tt + t - 1;
            
            // 
            arma::mat pred = A_t * Theta.slice(prev_idx) * B_t_trans;
            arma::mat diff = Theta.slice(curr_idx) - pred;
            
            // use vectorized dot product
            total_ss += arma::dot(diff, diff);
        }
    }
    
    return total_ss;
}

// batch ffbs for all theta relations simultaneously
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube ffbs_theta_all_relations(const arma::cube& Z_all,
                                   const arma::cube& M_all,
                                   const arma::cube& Aarray,
                                   const arma::cube& Barray,
                                   double sigma2_proc,
                                   double sigma2_obs,
                                   int m, int p, int Tt) {
    
    arma::cube Theta_all(m, m, p * Tt);
    
    // process each relation using fast bilinear ffbs
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int rel = 0; rel < p; rel++) {
        // extract z and mu for this relation
        arma::cube Z_rel(m, m, Tt);
        arma::mat mu_rel = M_all.slice(rel);
        
        for(int t = 0; t < Tt; t++) {
            Z_rel.slice(t) = Z_all.slice(rel * Tt + t);
        }
        
        // call bilinear ffbs
        arma::cube Theta_rel = ffbs_bilinear(Z_rel, mu_rel, Aarray, Barray, 
                                                  sigma2_proc, sigma2_obs);
        
        // store back
        for(int t = 0; t < Tt; t++) {
            Theta_all.slice(rel * Tt + t) = Theta_rel.slice(t);
        }
    }
    
    return Theta_all;
}

// theta update using blocked operations and cache-friendly access
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube ffbs_theta_blocked(const arma::cube& Z_all,
                             const arma::cube& M_all,
                             const arma::cube& Aarray,
                             const arma::cube& Barray,
                             double sigma2_proc,
                             double sigma2_obs,
                             int m, int p, int Tt) {
    
    arma::cube Theta_all(m, m, p * Tt);
    
    // process multiple relations in blocks for better cache usage
    const int block_size = std::min(4, p);  // process 4 relations at a time
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(int block_start = 0; block_start < p; block_start += block_size) {
        int block_end = std::min(block_start + block_size, p);
        
        for(int rel = block_start; rel < block_end; rel++) {
            // extract z and mu for this relation
            arma::cube Z_rel(m, m, Tt);
            arma::mat mu_rel = M_all.slice(rel);
            
            // vectorized copy
            for(int t = 0; t < Tt; t++) {
                Z_rel.slice(t) = Z_all.slice(rel * Tt + t);
            }
            
            // call bilinear ffbs
            arma::cube Theta_rel = ffbs_bilinear(Z_rel, mu_rel, Aarray, Barray, 
                                                      sigma2_proc, sigma2_obs);
            
            // store back
            for(int t = 0; t < Tt; t++) {
                Theta_all.slice(rel * Tt + t) = Theta_rel.slice(t);
            }
        }
    }
    
    return Theta_all;
}

// alpha update with reduced memory allocations
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat update_alpha_optimized(const arma::cube& Theta,
                                const arma::mat& U,
                                const arma::cube& Barray,
                                double sigma2_proc,
                                double tau_alpha2,
                                bool ar1_alpha,
                                double rho_alpha,
                                int m, int p, int Tt, int r) {
    
    // preallocate all working memory
    arma::mat alpha(r, Tt);
    arma::mat m_filt(r, Tt);
    arma::cube P_filt(r, r, Tt);
    
    // precompute U outer products once
    arma::cube U_outer(m, m, r);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int k = 0; k < r; k++) {
        U_outer.slice(k) = U.col(k) * U.col(k).t();
    }
    
    // initialize
    m_filt.col(0).zeros();
    P_filt.slice(0) = tau_alpha2 * arma::eye(r, r);
    
    // preallocate working arrays for the entire filter
    int obs_per_time = p * m * m;
    arma::vec y_t(obs_per_time);
    arma::mat H_t(obs_per_time, r);
    
    // forward filter
    for(int t = 1; t < Tt; t++) {
        arma::mat B_t = Barray.slice(t);
        
        // build observations and design matrix in parallel
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int rel = 0; rel < p; rel++) {
            int base_idx = rel * m * m;
            int prev_idx = rel * Tt + t - 1;
            int curr_idx = rel * Tt + t;
            
            arma::mat Theta_prev = Theta.slice(prev_idx);
            arma::mat Theta_curr = Theta.slice(curr_idx);
            arma::mat BTheta = B_t.t() * Theta_prev;
            
            // fill observations
            y_t.subvec(base_idx, base_idx + m*m - 1) = vectorise(Theta_curr);
            
            // fill design matrix
            for(int k = 0; k < r; k++) {
                arma::vec h_k = vectorise(U_outer.slice(k) * BTheta);
                H_t.col(k).subvec(base_idx, base_idx + m*m - 1) = h_k;
            }
        }
        
        // kalman update with numerical stability
        arma::vec m_pred;
        arma::mat P_pred;
        if (ar1_alpha) {
            m_pred = rho_alpha * m_filt.col(t-1);
            P_pred = rho_alpha * rho_alpha * P_filt.slice(t-1) + tau_alpha2 * arma::eye(r, r);
        } else {
            m_pred = m_filt.col(t-1);
            P_pred = P_filt.slice(t-1) + tau_alpha2 * arma::eye(r, r);
        }
        
        // use Cholesky decomposition for numerical stability
        arma::mat HtPH = H_t * P_pred * H_t.t();
        HtPH.diag() += sigma2_proc;  // add observation noise
        arma::mat S_chol = arma::chol(HtPH, "lower");
        
        // solve for Kalman gain using forward/backward substitution
        arma::mat PHt = P_pred * H_t.t();
        arma::mat K = arma::solve(arma::trimatl(S_chol), PHt.t()).t();
        K = arma::solve(arma::trimatu(S_chol.t()), K.t()).t();
        
        m_filt.col(t) = m_pred + K * (y_t - H_t * m_pred);
        P_filt.slice(t) = P_pred - K * H_t * P_pred;
        
        // ensure symmetry
        P_filt.slice(t) = 0.5 * (P_filt.slice(t) + P_filt.slice(t).t());
    }
    
    // backward sampling with numerical stability
    arma::mat P_final = P_filt.slice(Tt-1);
    P_final = 0.5 * (P_final + P_final.t());
    arma::mat P_chol = arma::chol(P_final + 1e-6 * arma::eye(r, r), "lower");
    alpha.col(Tt-1) = m_filt.col(Tt-1) + P_chol * arma::randn(r);
    
    for(int t = Tt-2; t >= 0; t--) {
        arma::mat P_pred;
        if (ar1_alpha) {
            P_pred = rho_alpha * rho_alpha * P_filt.slice(t) + tau_alpha2 * arma::eye(r, r);
        } else {
            P_pred = P_filt.slice(t) + tau_alpha2 * arma::eye(r, r);
        }
        
        // use Cholesky for stable inversion
        arma::mat P_pred_chol = arma::chol(P_pred + 1e-6 * arma::eye(r, r), "lower");
        arma::mat G = arma::solve(arma::trimatu(P_pred_chol.t()), 
                                  arma::solve(arma::trimatl(P_pred_chol), P_filt.slice(t).t())).t();
        
        arma::vec m_pred_t = ar1_alpha ? arma::vec(rho_alpha * m_filt.col(t)) : arma::vec(m_filt.col(t));
        arma::vec m_smooth = m_filt.col(t) + G * (alpha.col(t+1) - m_pred_t);
        arma::mat P_smooth = P_filt.slice(t) - G * (P_pred - P_filt.slice(t+1)) * G.t();
        
        // ensure symmetry and positive definiteness
        P_smooth = 0.5 * (P_smooth + P_smooth.t());
        arma::mat P_smooth_chol = arma::chol(P_smooth + 1e-6 * arma::eye(r, r), "lower");
        alpha.col(t) = m_smooth + P_smooth_chol * arma::randn(r);
    }
    
    return alpha;
}

// B update using blocked operations
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube update_B_blocked(const arma::cube& Theta,
                           const arma::cube& Aarray,
                           double sigma2_proc,
                           double tau_B2,
                           bool ar1_B,
                           double rho_B,
                           int m, int p, int Tt) {
    
    arma::cube Barray(m, m, Tt);
    
    // process columns in blocks for better cache usage
    const int block_size = std::min(8, m);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(int block_start = 0; block_start < m; block_start += block_size) {
        int block_end = std::min(block_start + block_size, m);
        
        // preallocate working memory for this block
        int block_cols = block_end - block_start;
        arma::mat m_filt(block_cols, Tt);
        arma::cube P_filt(block_cols, block_cols, Tt);
        
        // initialize
        m_filt.col(0).zeros();
        P_filt.slice(0) = tau_B2 * arma::eye(block_cols, block_cols);
        
        // forward filter for this block of columns
        for(int t = 1; t < Tt; t++) {
            arma::mat A_t = Aarray.slice(t);
            
            // build observations and design for this block
            arma::mat Y_block(p * m, block_cols);
            arma::mat F_block(p * m, block_cols);
            
            for(int rel = 0; rel < p; rel++) {
                int base_row = rel * m;
                int curr_idx = rel * Tt + t;
                int prev_idx = rel * Tt + t - 1;
                
                arma::mat Theta_curr = Theta.slice(curr_idx);
                arma::mat Theta_prev = Theta.slice(prev_idx);
                arma::mat ATheta = A_t * Theta_prev;
                
                Y_block.rows(base_row, base_row + m - 1) = 
                    Theta_curr.cols(block_start, block_end - 1);
                F_block.rows(base_row, base_row + m - 1) = 
                    ATheta.cols(block_start, block_end - 1);
            }
            
            // kalman update for block
            arma::vec y_vec = vectorise(Y_block);
            arma::mat F_mat = F_block;
            
            arma::vec m_pred;
            arma::mat P_pred;
            if (ar1_B) {
                m_pred = rho_B * m_filt.col(t-1);
                P_pred = rho_B * rho_B * P_filt.slice(t-1) + tau_B2 * arma::eye(block_cols, block_cols);
            } else {
                m_pred = m_filt.col(t-1);
                P_pred = P_filt.slice(t-1) + tau_B2 * arma::eye(block_cols, block_cols);
            }
            
            // update using Cholesky
            arma::mat S = F_mat.t() * F_mat + (sigma2_proc / P_pred(0,0)) * arma::eye(block_cols, block_cols);
            arma::mat S_chol = arma::chol(S, "lower");
            arma::mat K = arma::solve(arma::trimatu(S_chol.t()), 
                                     arma::solve(arma::trimatl(S_chol), F_mat.t())).t();
            
            m_filt.col(t) = m_pred + K.t() * (y_vec - F_mat * m_pred);
            P_filt.slice(t) = P_pred - K.t() * F_mat * P_pred;
        }
        
        // backward sample for block
        arma::mat P_final = P_filt.slice(Tt-1);
        P_final = 0.5 * (P_final + P_final.t());
        arma::mat P_chol = arma::chol(P_final + 1e-6 * arma::eye(block_cols, block_cols), "lower");
        arma::vec b_final = m_filt.col(Tt-1) + P_chol * arma::randn(block_cols);
        
        for(int k = 0; k < block_cols; k++) {
            Barray.slice(Tt-1)(block_start + k, block_start + k) = b_final(k);
        }
        
        // backward recursion
        for(int t = Tt-2; t >= 0; t--) {
            arma::mat P_pred;
            if (ar1_B) {
                P_pred = rho_B * rho_B * P_filt.slice(t) + tau_B2 * arma::eye(block_cols, block_cols);
            } else {
                P_pred = P_filt.slice(t) + tau_B2 * arma::eye(block_cols, block_cols);
            }
            
            arma::mat G = P_filt.slice(t) * arma::inv_sympd(P_pred);
            arma::vec m_pred_t = ar1_B ? arma::vec(rho_B * m_filt.col(t)) : arma::vec(m_filt.col(t));
            arma::vec b_next(block_cols);
            for(int k = 0; k < block_cols; k++) {
                b_next(k) = Barray.slice(t+1)(block_start + k, block_start + k);
            }
            
            arma::vec m_smooth = m_filt.col(t) + G * (b_next - m_pred_t);
            arma::mat P_smooth = P_filt.slice(t) - G * (P_pred - P_filt.slice(t+1)) * G.t();
            P_smooth = 0.5 * (P_smooth + P_smooth.t());
            
            arma::mat P_smooth_chol = arma::chol(P_smooth + 1e-6 * arma::eye(block_cols, block_cols), "lower");
            arma::vec b_t = m_smooth + P_smooth_chol * arma::randn(block_cols);
            
            for(int k = 0; k < block_cols; k++) {
                Barray.slice(t)(block_start + k, block_start + k) = b_t(k);
            }
        }
    }
    
    // fill off-diag elems
    for(int t = 0; t < Tt; t++) {
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < m; j++) {
                if(i != j) Barray.slice(t)(i, j) = 0;
            }
        }
    }
    
    return Barray;
}

// ordinal z update - handles global indices correctly
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
void update_Z_ordinal_global(arma::cube& Z_all,
                               const arma::cube& Theta_all,
                               const arma::cube& M_all,
                               const List& IR_list,
                               int m, int p, int Tt) {
    
    int m2 = m * m;  // elements per matrix
    
    // process each relation
    for(int rel = 0; rel < p; rel++) {
        List IR_rel = as<List>(IR_list[rel]);
        arma::mat M_rel = M_all.slice(rel);
        
        // create combined z and theta vectors for all time points
        arma::vec z_all_vec(m2 * Tt);
        arma::vec ez_all_vec(m2 * Tt);
        
        // fill combined vectors
        for(int t = 0; t < Tt; t++) {
            int idx = rel * Tt + t;
            arma::mat Z_t = Z_all.slice(idx);
            arma::mat Theta_t = Theta_all.slice(idx);
            arma::mat EZ = Theta_t + M_rel;
            
            // copy to combined vector at correct offset
            z_all_vec.subvec(t * m2, (t + 1) * m2 - 1) = vectorise(Z_t);
            ez_all_vec.subvec(t * m2, (t + 1) * m2 - 1) = vectorise(EZ);
        }
        
        // get rank names
        CharacterVector rank_names = IR_rel.names();
        
        // process each rank level
        for(int k = 0; k < IR_rel.size(); k++) {
            std::string rank_str = as<std::string>(rank_names[k]);
            if(rank_str == "NA") continue;
            
            // get indices for this rank (global indices across all time points)
            arma::vec idx_raw = as<arma::vec>(IR_rel[k]);
            if(idx_raw.n_elem == 0) continue;
            
            // convert to 0-based indices
            arma::uvec curr_idx = conv_to<arma::uvec>::from(idx_raw - 1);
            
            // determine bounds based on rank value
            int rank_val = std::atoi(rank_str.c_str());
            double lb = -arma::datum::inf;
            double ub = arma::datum::inf;
            
            // find bounds from adjacent ranks
            if(rank_val > 1) {
                // find max of previous rank
                std::string prev_rank_str = std::to_string(rank_val - 1);
                for(int j = 0; j < IR_rel.size(); j++) {
                    if(as<std::string>(rank_names[j]) == prev_rank_str) {
                        arma::vec prev_idx_raw = as<arma::vec>(IR_rel[j]);
                        if(prev_idx_raw.n_elem > 0) {
                            arma::uvec prev_idx = conv_to<arma::uvec>::from(prev_idx_raw - 1);
                            if(prev_idx.n_elem > 0) {
                                lb = z_all_vec.elem(prev_idx).max();
                            }
                        }
                        break;
                    }
                }
            }
            
            // check for next rank
            std::string next_rank_str = std::to_string(rank_val + 1);
            for(int j = 0; j < IR_rel.size(); j++) {
                if(as<std::string>(rank_names[j]) == next_rank_str) {
                    arma::vec next_idx_raw = as<arma::vec>(IR_rel[j]);
                    if(next_idx_raw.n_elem > 0) {
                        arma::uvec next_idx = conv_to<arma::uvec>::from(next_idx_raw - 1);
                        if(next_idx.n_elem > 0) {
                            ub = z_all_vec.elem(next_idx).min();
                        }
                    }
                    break;
                }
            }
            
            // sample from truncated normal for this rank
            arma::vec mu_vec = ez_all_vec.elem(curr_idx);
            for(unsigned int j = 0; j < curr_idx.n_elem; j++) {
                double mu = mu_vec(j);
                
                // truncated normal sampling
                double p_a = R::pnorm(lb, mu, 1.0, true, false);
                double p_b = R::pnorm(ub, mu, 1.0, true, false);
                
                if(p_b > p_a) {
                    double u = R::runif(p_a, p_b);
                    z_all_vec(curr_idx(j)) = R::qnorm(u, mu, 1.0, true, false);
                } else {
                    z_all_vec(curr_idx(j)) = mu;  // fallback
                }
            }
        }
        
        // copy back to individual time slices
        for(int t = 0; t < Tt; t++) {
            int idx = rel * Tt + t;
            arma::vec z_t_vec = z_all_vec.subvec(t * m2, (t + 1) * m2 - 1);
            Z_all.slice(idx) = reshape(z_t_vec, m, m);
        }
    }
}