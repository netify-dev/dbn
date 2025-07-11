#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// forward declarations
arma::mat ffbs_dlm_cpp(const List& y, const List& Flist, const arma::mat& V, 
                       const arma::mat& W, const arma::vec& m0, const arma::mat& C0, 
                       bool ar1, double rho);

// parallel b matrix update with memory access
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube update_B_parallel(const arma::cube& Theta,
                                   const arma::cube& Aarray,
                                   double sigma2_proc,
                                   double tau_B2,
                                   bool ar1_B,
                                   double rho_B,
                                   int m, int p, int Tt) {
    
    arma::cube Barray(m, m, Tt);
    
    // process all columns in parallel
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int k = 0; k < m; k++) {
        // build observations and design matrices for column k
        arma::mat Y_stacked((Tt - 1) * m * p, 1);
        arma::mat F_stacked((Tt - 1) * m * p, m);
        
        int row_idx = 0;
        for(int t = 1; t < Tt; t++) {
            arma::mat A_t = Aarray.slice(t);
            
            for(int rel = 0; rel < p; rel++) {
                // compute slice index
                int curr_idx = rel * Tt + t;
                int prev_idx = rel * Tt + t - 1;
                
                // check bounds
                if (curr_idx >= Theta.n_slices || prev_idx >= Theta.n_slices) {
                    Rcpp::stop("Index out of bounds in update_B_parallel: curr=%d, prev=%d, n_slices=%d", 
                               curr_idx, prev_idx, Theta.n_slices);
                }
                
                // observations: k-th column of theta at time t
                arma::vec y_tk = Theta.slice(curr_idx).col(k);
                Y_stacked.rows(row_idx, row_idx + m - 1) = y_tk;
                
                // design: A_t * Theta_{t-1}
                arma::mat Theta_prev = Theta.slice(prev_idx);
                F_stacked.rows(row_idx, row_idx + m - 1) = A_t * Theta_prev;
                
                row_idx += m;
            }
        }
        
        // run kalman filter and smoother
        int state_dim = m;
        int obs_dim = m * p;
        
        // forward filter
        arma::mat m_filt(state_dim, Tt);
        arma::cube P_filt(state_dim, state_dim, Tt);
        
        // initialize
        m_filt.col(0).zeros();
        P_filt.slice(0) = tau_B2 * arma::eye(state_dim, state_dim);
        
        // forward pass
        arma::vec m_t = m_filt.col(0);
        arma::mat P_t = P_filt.slice(0);
        
        for(int t = 1; t < Tt; t++) {
            // extract current observations
            arma::vec y_t = Y_stacked.rows((t-1)*obs_dim, t*obs_dim - 1);
            arma::mat F_t = F_stacked.rows((t-1)*obs_dim, t*obs_dim - 1);
            
            // prediction
            arma::vec m_pred;
            arma::mat P_pred;
            if (ar1_B) {
                m_pred = rho_B * m_t;
                P_pred = rho_B * rho_B * P_t + tau_B2 * arma::eye(state_dim, state_dim);
            } else {
                m_pred = m_t;
                P_pred = P_t + tau_B2 * arma::eye(state_dim, state_dim);
            }
            
            // update with numerical stability
            arma::mat S = F_t * P_pred * F_t.t() + sigma2_proc * arma::eye(obs_dim, obs_dim);
            
            // kalman gain cmp using cholesky
            arma::mat K;
            arma::mat L_S;
            bool chol_success = arma::chol(L_S, S, "lower");
            
            if (chol_success) {
                arma::mat temp = arma::solve(arma::trimatl(L_S), F_t * P_pred);
                K = arma::solve(arma::trimatu(L_S.t()), temp).t();
            } else {
                // Fallback to standard computation
                K = P_pred * F_t.t() * arma::inv(S);
            }
            
            m_t = m_pred + K * (y_t - F_t * m_pred);
            
            // joseph form update again so things dont blow up hopefully
            arma::mat I_minus_KF = arma::eye(state_dim, state_dim) - K * F_t;
            P_t = I_minus_KF * P_pred * I_minus_KF.t() + K * sigma2_proc * K.t();
            
            // force exact symmetry to avoid warnings
            P_t = 0.5 * (P_t + P_t.t());
            P_t = P_t + P_t.t();
            P_t = 0.5 * P_t;
            
            // store
            m_filt.col(t) = m_t;
            P_filt.slice(t) = P_t;
        }
        
        // backward sample
        arma::mat P_final = P_filt.slice(Tt-1);
        // ensure positive definiteness before sampling
        P_final = P_final + P_final.t();
        P_final = 0.5 * P_final;
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, P_final);
        eigval = arma::clamp(eigval, 1e-6, arma::datum::inf);
        P_final = eigvec * arma::diagmat(eigval) * eigvec.t();
        
        arma::vec b_T = mvnrnd(m_filt.col(Tt-1), P_final);
        Barray.slice(Tt-1).col(k) = b_T;
        
        for(int t = Tt-2; t >= 0; t--) {
            // smoother gain
            arma::mat P_pred;
            if (ar1_B) {
                P_pred = rho_B * rho_B * P_filt.slice(t) + tau_B2 * arma::eye(state_dim, state_dim);
            } else {
                P_pred = P_filt.slice(t) + tau_B2 * arma::eye(state_dim, state_dim);
            }
            arma::mat G = P_filt.slice(t) * arma::inv(P_pred);  // use regular inv
            
            // smooth
            arma::vec m_pred_t;
            if (ar1_B) {
                m_pred_t = rho_B * m_filt.col(t);
            } else {
                m_pred_t = m_filt.col(t);
            }
            arma::vec m_smooth = m_filt.col(t) + G * (Barray.slice(t+1).col(k) - m_pred_t);
            arma::mat P_smooth = P_filt.slice(t) - G * (P_pred - P_filt.slice(t+1)) * G.t();
            
            // force exact symmetry before sampling
            P_smooth = P_smooth + P_smooth.t();
            P_smooth = 0.5 * P_smooth;
            
            // ensure positive definiteness
            arma::vec eigval2;
            arma::mat eigvec2;
            arma::eig_sym(eigval2, eigvec2, P_smooth);
            eigval2 = arma::clamp(eigval2, 1e-6, arma::datum::inf);
            P_smooth = eigvec2 * arma::diagmat(eigval2) * eigvec2.t();
            
            // sample
            Barray.slice(t).col(k) = mvnrnd(m_smooth, P_smooth);
        }
    }
    
    return Barray;
}

// batch alpha update with block operations
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat update_alpha_batch(const arma::cube& Theta,
                             const arma::mat& U,
                             const arma::cube& Barray,
                             double sigma2_proc,
                             double tau_alpha2,
                             bool ar1_alpha,
                             double rho_alpha,
                             int m, int p, int Tt, int r) {
    
    // pre-compute u outer products
    arma::cube U_outer(m, m, r);
    for(int k = 0; k < r; k++) {
        U_outer.slice(k) = U.col(k) * U.col(k).t();
    }
    
    // build giant design matrix efficiently
    int total_rows = (Tt - 1) * p * m * m;
    arma::mat H_all(total_rows, r);
    arma::vec y_all(total_rows);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int t = 1; t < Tt; t++) {
        arma::mat B_t = Barray.slice(t);
        
        for(int rel = 0; rel < p; rel++) {
            // compute slice indices
            int prev_idx = rel * Tt + t - 1;
            int curr_idx = rel * Tt + t;
            
            // check bounds
            if (curr_idx >= Theta.n_slices || prev_idx >= Theta.n_slices) {
                Rcpp::stop("Index out of bounds in update_alpha_batch: curr=%d, prev=%d, n_slices=%d", 
                           curr_idx, prev_idx, Theta.n_slices);
            }
            
            arma::mat Theta_prev = Theta.slice(prev_idx);
            arma::mat Theta_curr = Theta.slice(curr_idx);
            
            // compute B * Theta' efficiently
            arma::mat BTheta = B_t.t() * Theta_prev;
            
            // fill design matrix
            int base_row = ((t-1) * p + rel) * m * m;
            
            for(int k = 0; k < r; k++) {
                arma::vec h_k = vectorise(U_outer.slice(k) * BTheta);
                H_all.col(k).rows(base_row, base_row + m*m - 1) = h_k;
            }
            
            // fill observations
            y_all.rows(base_row, base_row + m*m - 1) = vectorise(Theta_curr);
        }
    }
    
    // run kalman filter
    arma::mat alpha(r, Tt);
    
    // forward filter
    arma::mat m_filt(r, Tt);
    arma::cube P_filt(r, r, Tt);
    
    m_filt.col(0).zeros();
    P_filt.slice(0) = tau_alpha2 * arma::eye(r, r);
    
    for(int t = 1; t < Tt; t++) {
        // extract block
        int start_idx = (t-1) * p * m * m;
        int end_idx = t * p * m * m - 1;
        
        arma::mat H_t = H_all.rows(start_idx, end_idx);
        arma::vec y_t = y_all.rows(start_idx, end_idx);
        
        // kalman update
        arma::vec m_pred;
        arma::mat P_pred;
        if (ar1_alpha) {
            m_pred = rho_alpha * m_filt.col(t-1);
            P_pred = rho_alpha * rho_alpha * P_filt.slice(t-1) + tau_alpha2 * arma::eye(r, r);
        } else {
            m_pred = m_filt.col(t-1);
            P_pred = P_filt.slice(t-1) + tau_alpha2 * arma::eye(r, r);
        }
        
        // For large observation dimensions, use Woodbury formula
        arma::mat K;
        if (H_t.n_rows > 1000 && r < 50) {
            // Woodbury formula for (HPH' + σ²I)^{-1}
            double inv_sigma2 = 1.0 / sigma2_proc;
            arma::mat inner = arma::eye(r, r) + inv_sigma2 * H_t.t() * H_t * P_pred;
            arma::mat inner_inv = arma::inv_sympd(inner);
            K = P_pred * H_t.t() * (inv_sigma2 * arma::eye(H_t.n_rows, H_t.n_rows) - 
                                   inv_sigma2 * inv_sigma2 * H_t * P_pred * inner_inv * H_t.t());
        } else {
            // Standard approach for smaller problems
            arma::mat S = H_t * P_pred * H_t.t() + sigma2_proc * arma::eye(H_t.n_rows, H_t.n_rows);
            
            // Efficient computation using Cholesky
            arma::mat L_S;
            bool chol_success = arma::chol(L_S, S, "lower");
            
            if (chol_success) {
                arma::mat temp = arma::solve(arma::trimatl(L_S), H_t * P_pred);
                K = arma::solve(arma::trimatu(L_S.t()), temp).t();
            } else {
                K = P_pred * H_t.t() * arma::inv(S);
            }
        }
        
        m_filt.col(t) = m_pred + K * (y_t - H_t * m_pred);
        
        // Joseph form update for numerical stability
        arma::mat I_minus_KH = arma::eye(r, r) - K * H_t;
        P_filt.slice(t) = I_minus_KH * P_pred * I_minus_KH.t() + K * sigma2_proc * K.t();
    }
    
    // backward sample - force exact symmetry
    arma::mat P_final = P_filt.slice(Tt-1);
    P_final = P_final + P_final.t();
    P_final = 0.5 * P_final;
    
    // ensure positive definiteness
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, P_final);
    eigval = arma::clamp(eigval, 1e-6, arma::datum::inf);
    P_final = eigvec * arma::diagmat(eigval) * eigvec.t();
    
    alpha.col(Tt-1) = mvnrnd(m_filt.col(Tt-1), P_final);
    
    for(int t = Tt-2; t >= 0; t--) {
        arma::mat P_pred;
        if (ar1_alpha) {
            P_pred = rho_alpha * rho_alpha * P_filt.slice(t) + tau_alpha2 * arma::eye(r, r);
        } else {
            P_pred = P_filt.slice(t) + tau_alpha2 * arma::eye(r, r);
        }
        arma::mat G = P_filt.slice(t) * arma::inv(P_pred);  // use regular inv
        
        arma::vec m_pred_t;
        if (ar1_alpha) {
            m_pred_t = rho_alpha * m_filt.col(t);
        } else {
            m_pred_t = m_filt.col(t);
        }
        arma::vec m_smooth = m_filt.col(t) + G * (alpha.col(t+1) - m_pred_t);
        arma::mat P_smooth = P_filt.slice(t) - G * (P_pred - P_filt.slice(t+1)) * G.t();
        
        // force exact symmetry before sampling
        P_smooth = P_smooth + P_smooth.t();
        P_smooth = 0.5 * P_smooth;
        
        // ensure positive definiteness
        arma::vec eigval3;
        arma::mat eigvec3;
        arma::eig_sym(eigval3, eigvec3, P_smooth);
        eigval3 = arma::clamp(eigval3, 1e-6, arma::datum::inf);
        P_smooth = eigvec3 * arma::diagmat(eigval3) * eigvec3.t();
        
        alpha.col(t) = mvnrnd(m_smooth, P_smooth);
    }
    
    return alpha;
}

// vectorized ordinal z update
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
void update_Z_ordinal_vectorized(arma::cube& Z_all,
                           const arma::cube& Theta_all,
                           const arma::cube& M_all,
                           const List& IR_list,
                           int m, int p, int Tt) {
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int rel = 0; rel < p; rel++) {
        List IR_rel = IR_list[rel];
        arma::mat M_rel = M_all.slice(rel);
        
        for(int t = 0; t < Tt; t++) {
            int idx = rel * Tt + t;
            arma::mat Z_t = Z_all.slice(idx);
            arma::mat Theta_t = Theta_all.slice(idx);
            arma::mat EZ = Theta_t + M_rel;
            
            // vectorized rank constraint sampling
            arma::vec z_vec = vectorise(Z_t);
            arma::vec ez_vec = vectorise(EZ);
            
            // process each rank level
            CharacterVector rank_names = IR_rel.names();
            for(int k = 0; k < IR_rel.size(); k++) {
                std::string rank_str = as<std::string>(rank_names[k]);
                if(rank_str == "NA") continue;
                
                arma::vec idx_k_raw = as<arma::vec>(IR_rel[k]);
                arma::uvec idx_k = conv_to<arma::uvec>::from(idx_k_raw - 1);
                if(idx_k.n_elem == 0) continue;
                
                // determine bounds
                int rank_val = std::atoi(rank_str.c_str());
                double lb = -datum::inf;
                double ub = datum::inf;
                
                // find lower bound from previous rank
                if(rank_val > 1) {
                    std::string prev_rank = std::to_string(rank_val - 1);
                    for(int j = 0; j < IR_rel.size(); j++) {
                        if(as<std::string>(rank_names[j]) == prev_rank) {
                            arma::vec idx_prev_raw = as<arma::vec>(IR_rel[j]);
                            arma::uvec idx_prev = conv_to<arma::uvec>::from(idx_prev_raw - 1);
                            if(idx_prev.n_elem > 0) {
                                lb = z_vec.elem(idx_prev).max();
                            }
                            break;
                        }
                    }
                }
                
                // find upper bound from next rank
                if(rank_val < IR_rel.size() - 1) {
                    std::string next_rank = std::to_string(rank_val + 1);
                    for(int j = 0; j < IR_rel.size(); j++) {
                        if(as<std::string>(rank_names[j]) == next_rank) {
                            arma::vec idx_next_raw = as<arma::vec>(IR_rel[j]);
                            arma::uvec idx_next = conv_to<arma::uvec>::from(idx_next_raw - 1);
                            if(idx_next.n_elem > 0) {
                                ub = z_vec(idx_next).min();
                            }
                            break;
                        }
                    }
                }
                
                // sample from truncated normal
                for(unsigned int i = 0; i < idx_k.n_elem; i++) {
                    double mu = ez_vec(idx_k(i));
                    double p_lo = R::pnorm(lb, mu, 1.0, true, false);
                    double p_hi = R::pnorm(ub, mu, 1.0, true, false);
                    
                    if(p_lo < p_hi) {
                        double u = R::runif(p_lo, p_hi);
                        z_vec(idx_k(i)) = R::qnorm(u, mu, 1.0, true, false);
                    } else {
                        z_vec(idx_k(i)) = mu;
                    }
                }
            }
            
            // reshape back
            Z_all.slice(idx) = reshape(z_vec, m, m);
        }
    }
}