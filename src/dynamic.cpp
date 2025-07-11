#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "thread_control.h"

using namespace Rcpp;
using namespace arma;

// Forward declarations for existing functions
arma::vec rz_fc_cpp(const arma::vec& R, const arma::vec& Z, const arma::vec& EZ, const List& iranks);
arma::cube ffbs_theta_struct_5arg_cpp(const arma::cube& Z, const arma::mat& mu, 
                                     const arma::cube& A_array, const arma::cube& B_array, 
                                     double sigma2);

// Batch update z for ordinal data
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat batch_update_Z_ordinal(const arma::mat& R_4d,
                                 const arma::mat& Z_4d, 
                                 const arma::mat& Theta_4d,
                                 const arma::cube& M,
                                 const List& IR,
                                 int m, int p, int Tt) {
    // create output array
    arma::mat Z_new(m * m, p * Tt);
    
    // process all relations
    for(int j = 0; j < p; j++) {
        // get rank indices for this relation
        List IR_j = IR[j];
        
        // process all time points for this relation
        for(int t = 0; t < Tt; t++) {
            // extract columns
            int col_idx = j * Tt + t;
            arma::vec R_vec = R_4d.col(col_idx);
            arma::vec Z_vec = Z_4d.col(col_idx);
            arma::vec Theta_vec = Theta_4d.col(col_idx);
            
            // compute EZ = Theta + M for this relation
            arma::mat M_j = M.slice(j);
            arma::vec M_vec = vectorise(M_j);
            arma::vec EZ_vec = Theta_vec + M_vec;
            
            // create adjusted IR for this specific time point
            List IR_t;
            CharacterVector names = IR_j.names();
            for(int k = 0; k < IR_j.size(); k++) {
                std::string rank_name = Rcpp::as<std::string>(names[k]);
                arma::vec idx_vec = Rcpp::as<arma::vec>(IR_j[k]);
                
                // filter indices for current time point
                arma::uvec local_indices;
                for(unsigned int i = 0; i < idx_vec.n_elem; i++) {
                    int global_idx = idx_vec(i) - 1; // convert to 0-based
                    int time_idx = global_idx / (m * m);
                    if(time_idx == t) {
                        int spatial_idx = global_idx % (m * m);
                        local_indices.resize(local_indices.n_elem + 1);
                        local_indices(local_indices.n_elem - 1) = spatial_idx + 1; // convert back to 1-based
                    }
                }
                
                if(local_indices.n_elem > 0) {
                    IR_t[rank_name] = local_indices;
                } else {
                    IR_t[rank_name] = arma::uvec(); // empty vector
                }
            }
            
            // update Z using rank likelihood
            arma::vec Z_updated = rz_fc_cpp(R_vec, Z_vec, EZ_vec, IR_t);
            
            // store result
            Z_new.col(col_idx) = Z_updated;
        }
    }
    
    return Z_new;
}

// Z update with preallocated memory and better cache usage
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat batch_update_Z_ordinal_fast(const arma::mat& R_4d,
                                       const arma::mat& Z_4d, 
                                       const arma::mat& Theta_4d,
                                       const arma::cube& M,
                                       const List& IR,
                                       const List& IR_time_indices,  // precomputed time indices
                                       int m, int p, int Tt) {
    set_dbn_threads(); // set threads from R options
    
    // allocate output once
    arma::mat Z_new(m * m, p * Tt);
    
    // preallocate workspace for each thread
    #ifdef _OPENMP
    int n_threads = omp_get_max_threads();
    std::vector<arma::vec> workspace_R(n_threads, arma::vec(m * m));
    std::vector<arma::vec> workspace_Z(n_threads, arma::vec(m * m));
    std::vector<arma::vec> workspace_EZ(n_threads, arma::vec(m * m));
    #endif
    
    // parallel over relations with better load balancing
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for(int j = 0; j < p; j++) {
        #ifdef _OPENMP
        int tid = omp_get_thread_num();
        arma::vec& local_R = workspace_R[tid];
        arma::vec& local_Z = workspace_Z[tid];
        arma::vec& local_EZ = workspace_EZ[tid];
        #else
        arma::vec local_R(m * m);
        arma::vec local_Z(m * m);
        arma::vec local_EZ(m * m);
        #endif
        
        // get pre-computed time indices for this relation
        List IR_j_time = IR_time_indices[j];
        arma::mat M_j = M.slice(j);
        arma::vec M_vec = vectorise(M_j);
        
        // process time points with better memory access pattern
        for(int t = 0; t < Tt; t++) {
            int col_idx = j * Tt + t;
            
            // use pre-computed indices for this time point
            List IR_t = IR_j_time[t];
            
            // vectorized operations
            arma::vec EZ_vec = Theta_4d.col(col_idx) + M_vec;
            
            // update z
            Z_new.col(col_idx) = rz_fc_cpp(R_4d.col(col_idx), Z_4d.col(col_idx), EZ_vec, IR_t);
        }
    }
    
    return Z_new;
}

// Precompute time-specific rank indices to avoid repeated computation
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List precompute_time_indices(const List& IR, int m, int p, int Tt) {
    List IR_time_indices(p);
    
    for(int j = 0; j < p; j++) {
        List IR_j = IR[j];
        List IR_j_time(Tt);
        CharacterVector names = IR_j.names();
        
        // precompute indices for each time point
        for(int t = 0; t < Tt; t++) {
            List IR_t;
            
            for(int k = 0; k < IR_j.size(); k++) {
                std::string rank_name = Rcpp::as<std::string>(names[k]);
                arma::vec idx_vec = Rcpp::as<arma::vec>(IR_j[k]);
                
                // vectorized filtering
                arma::uvec time_mask = (arma::floor((idx_vec - 1) / (m * m)) == t);
                arma::vec filtered_global = idx_vec.elem(find(time_mask));
                
                if(filtered_global.n_elem > 0) {
                    // convert to local indices
                    arma::uvec local_indices(filtered_global.n_elem);
                    for(unsigned int i = 0; i < filtered_global.n_elem; i++) {
                        int global_idx = filtered_global(i) - 1;
                        local_indices(i) = (global_idx % (m * m)) + 1;
                    }
                    IR_t[rank_name] = local_indices;
                }
            }
            
            IR_j_time[t] = IR_t;
        }
        
        IR_time_indices[j] = IR_j_time;
    }
    
    return IR_time_indices;
}

// mu update with blocked computation for large m
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_mu_dynamic(const arma::mat& Z_4d,
                              const arma::mat& Theta_4d,
                              double g2,
                              double a_g, double b_g,
                              int m, int p, int Tt) {
    double mu_var = 1.0 / (Tt + 1.0 / g2);
    arma::cube M(m, m, p);
    
    // process relations in parallel
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int j = 0; j < p; j++) {
        // blocked summation for numerical stability with large Tt
        arma::mat sum_diff(m, m, fill::zeros);
        
        // sum over time 
        for(int t = 0; t < Tt; t++) {
            int col_idx = j * Tt + t;
            sum_diff += reshape(Z_4d.col(col_idx) - Theta_4d.col(col_idx), m, m);
        }
        
        // compute posterior mean
        arma::mat mu_hat = mu_var * sum_diff;
        
        // sample from posterior
        arma::mat noise(m, m);
        noise.randn();
        M.slice(j) = mu_hat + sqrt(mu_var) * noise;
    }
    
    // update g2 with stable computation
    double M_sum_sq = 0.0;
    for(int j = 0; j < p; j++) {
        M_sum_sq += accu(square(M.slice(j)));
    }
    
    double shape = (a_g + m * m * p) / 2.0;
    double rate = (b_g + M_sum_sq) / 2.0;
    double g2_new = rate / R::rgamma(shape, 1.0);
    
    return List::create(
        Named("M") = M,
        Named("g2") = g2_new
    );
}

// ffbs with pre-allocated workspace
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube batch_ffbs_all_relations(const arma::mat& Z_4d,
                           const arma::cube& M,
                           const arma::cube& Aarray,
                           const arma::cube& Barray,
                           double sigma2,
                           int m, int p, int Tt) {
    arma::cube Theta_4d_out(m, m, p * Tt);
    
    // parallel over relations with thread-local workspace
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int j = 0; j < p; j++) {
        // extract data for this relation
        arma::cube Z_j(m, m, Tt);
        
        // vectorized extraction
        for(int t = 0; t < Tt; t++) {
            Z_j.slice(t) = reshape(Z_4d.col(j * Tt + t), m, m);
        }
        
        // run ffbs
        arma::cube Theta_j = ffbs_theta_struct_5arg_cpp(Z_j, M.slice(j), Aarray, Barray, sigma2);
        
        // vectorized storage
        for(int t = 0; t < Tt; t++) {
            Theta_4d_out.slice(j * Tt + t) = Theta_j.slice(t);
        }
    }
    
    return Theta_4d_out;
}

// AB update for large networks
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_AB_batch_extended(const arma::mat& Theta_4d,
                    const arma::cube& Aarray_old,
                    const arma::cube& Barray_old,
                    double sigma2, double tauA2, double tauB2,
                    bool ar1, double rhoA, double rhoB,
                    int m, int p, int Tt) {
    arma::cube Aarray(m, m, Tt);
    arma::cube Barray(m, m, Tt);
    
    // precompute frequently used matrices
    arma::mat eye_m = eye(m, m);
    double inv_sigma2 = 1.0 / sigma2;
    double inv_tauA2 = 1.0 / tauA2;
    double inv_tauB2 = 1.0 / tauB2;
    
    // update a with blocked computation for large p
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for(int i = 0; i < m; i++) {
        for(int t = 1; t < Tt; t++) {
            // allocate workspace once per thread
            arma::mat F_it(p * m, m);
            arma::vec y_it(p * m);
            
            // build design matrix with better memory access
            for(int j = 0; j < p; j++) {
                int base_idx = j * m;
                int col_idx_prev = j * Tt + t - 1;
                int col_idx_curr = j * Tt + t;
                
                // extract and reshape in one operation
                arma::mat Theta_prev = reshape(Theta_4d.col(col_idx_prev), m, m);
                arma::mat Theta_curr = reshape(Theta_4d.col(col_idx_curr), m, m);
                
                F_it.rows(base_idx, base_idx + m - 1) = (Theta_prev * Barray_old.slice(t).t()).t();
                y_it.subvec(base_idx, base_idx + m - 1) = Theta_curr.row(i).t();
            }
            
            // compute posterior with pre-computed constants
            arma::mat V_inv = inv_sigma2 * (F_it.t() * F_it) + inv_tauA2 * eye_m;
            arma::mat V = inv_sympd(V_inv);
            // force exact symmetry before mvnrnd
            V = 0.5 * (V + V.t());
            arma::vec m_post = V * (inv_sigma2 * (F_it.t() * y_it));
            
            if(ar1 && t > 1) {
                m_post += (rhoA * inv_tauA2) * (V * Aarray.slice(t - 1).row(i).t());
            }
            
            // sample new row
            arma::vec a_new = mvnrnd(m_post, V);
            Aarray.slice(t).row(i) = a_new.t();
        }
    }
    
    // set a_1 = i
    Aarray.slice(0) = eye_m;
    
    // update b 
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for(int k = 0; k < m; k++) {
        for(int t = 1; t < Tt; t++) {
            arma::mat F_kt(p * m, m);
            arma::vec y_kt(p * m);
            
            for(int j = 0; j < p; j++) {
                int base_idx = j * m;
                int col_idx_prev = j * Tt + t - 1;
                int col_idx_curr = j * Tt + t;
                
                arma::mat Theta_prev = reshape(Theta_4d.col(col_idx_prev), m, m);
                arma::mat Theta_curr = reshape(Theta_4d.col(col_idx_curr), m, m);
                
                F_kt.rows(base_idx, base_idx + m - 1) = Aarray.slice(t) * Theta_prev;
                y_kt.subvec(base_idx, base_idx + m - 1) = Theta_curr.col(k);
            }
            
            arma::mat V_inv = inv_sigma2 * (F_kt.t() * F_kt) + inv_tauB2 * eye_m;
            arma::mat V = inv_sympd(V_inv);
            // force exact symmetry before mvnrnd
            V = 0.5 * (V + V.t());
            arma::vec m_post = V * (inv_sigma2 * (F_kt.t() * y_kt));
            
            if(ar1 && t > 1) {
                m_post += (rhoB * inv_tauB2) * (V * Barray.slice(t - 1).col(k));
            }
            
            arma::vec b_new = mvnrnd(m_post, V);
            Barray.slice(t).col(k) = b_new;
        }
    }
    
    // set b_1 = i
    Barray.slice(0) = eye_m;
    
    return List::create(
        Named("Aarray") = Aarray,
        Named("Barray") = Barray
    );
}

// Combined variance update with minimal memory allocation
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_variances_dynamic(const arma::mat& Theta_4d,
                           const arma::mat& Z_4d,
                           const arma::cube& M,
                           const arma::cube& Aarray,
                           const arma::cube& Barray,
                           double a_sig, double b_sig,
                           int m, int p, int Tt,
                           bool is_gaussian = false) {
    // compute process variance residuals with blocking
    double proc_rss = 0.0;
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:proc_rss)
    #endif
    for(int j = 0; j < p; j++) {
        double local_rss = 0.0;
        
        for(int t = 1; t < Tt; t++) {
            int idx_curr = j * Tt + t;
            int idx_prev = j * Tt + t - 1;
            
            // use in-place operations to minimize memory allocation
            arma::mat Theta_curr = reshape(Theta_4d.col(idx_curr), m, m);
            arma::mat Theta_prev = reshape(Theta_4d.col(idx_prev), m, m);
            
            // compute residual = Theta_curr - A_t * Theta_prev * B_t'
            arma::mat pred = Aarray.slice(t) * Theta_prev * Barray.slice(t).t();
            local_rss += accu(square(Theta_curr - pred));
        }
        
        proc_rss += local_rss;
    }
    
    double sigma2 = (b_sig + proc_rss / 2.0) / R::rgamma((a_sig + m * m * (Tt - 1) * p) / 2.0, 1.0);
    
    double sigma2_obs = 1.0;
    if(is_gaussian) {
        double obs_rss = 0.0;
        
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:obs_rss)
        #endif
        for(int j = 0; j < p; j++) {
            double local_rss = 0.0;
            arma::mat M_j = M.slice(j);
            
            for(int t = 0; t < Tt; t++) {
                int idx = j * Tt + t;
                
                arma::mat Z_jt = reshape(Z_4d.col(idx), m, m);
                arma::mat Theta_jt = reshape(Theta_4d.col(idx), m, m);
                
                local_rss += accu(square(Z_jt - (Theta_jt + M_j)));
            }
            
            obs_rss += local_rss;
        }
        
        sigma2_obs = (1.0 + obs_rss / 2.0) / R::rgamma((1.0 + m * m * Tt * p) / 2.0, 1.0);
    }
    
    return List::create(
        Named("sigma2") = sigma2,
        Named("sigma2_obs") = sigma2_obs
    );
}

// FFBS for very large networks using blocked operations
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube batch_ffbs_all_relations_blocked(const arma::mat& Z_4d,
                                           const arma::cube& M,
                                           const arma::cube& Aarray,
                                           const arma::cube& Barray,
                                           double sigma2,
                                           int m, int p, int Tt) {
    arma::cube Theta_4d_out(m, m, p * Tt);
    
    // process relations in blocks 
    const int block_size = std::max(1, std::min(4, p)); 
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for(int j_block = 0; j_block < p; j_block += block_size) {
        int j_end = std::min(j_block + block_size, p);
        
        // process block of relations
        for(int j = j_block; j < j_end; j++) {
            // preallocate workspace for this relation
            arma::cube Z_j(m, m, Tt);
            
            // vectorized extraction
            #pragma omp simd
            for(int t = 0; t < Tt; t++) {
                Z_j.slice(t) = reshape(Z_4d.col(j * Tt + t), m, m);
            }
            
            // ffbs
            arma::cube Theta_j = ffbs_theta_struct_5arg_cpp(Z_j, M.slice(j), Aarray, Barray, sigma2);
            
            // store
            #pragma omp simd
            for(int t = 0; t < Tt; t++) {
                Theta_4d_out.slice(j * Tt + t) = Theta_j.slice(t);
            }
        }
    }
    
    return Theta_4d_out;
}

// AB update for very large m (>100)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_AB_batch_large(const arma::mat& Theta_4d,
                          const arma::cube& Aarray_old,
                          const arma::cube& Barray_old,
                          double sigma2, double tauA2, double tauB2,
                          bool ar1, double rhoA, double rhoB,
                          int m, int p, int Tt) {
    arma::cube Aarray(m, m, Tt);
    arma::cube Barray(m, m, Tt);
    
    // precompute frequently used matrices
    arma::mat eye_m = eye(m, m);
    double inv_sigma2 = 1.0 / sigma2;
    double inv_tauA2 = 1.0 / tauA2;
    double inv_tauB2 = 1.0 / tauB2;
    
    // use tiled approach for very large m
    const int tile_size = std::min(32, m); // Smaller tiles 
    
    // update a with tiled computation
    #ifdef _OPENMP
    #pragma omp parallel
    {
        // thread-local workspace
        arma::mat F_local(p * m, m);
        arma::vec y_local(p * m);
        arma::mat V_local(m, m);
        
        #pragma omp for schedule(dynamic, 1)
        for(int i_tile = 0; i_tile < m; i_tile += tile_size) {
            int i_end = std::min(i_tile + tile_size, m);
            
            for(int i = i_tile; i < i_end; i++) {
                for(int t = 1; t < Tt; t++) {
                    // build design matrix
                    for(int j = 0; j < p; j++) {
                        int base_idx = j * m;
                        
                        arma::mat Theta_prev = reshape(Theta_4d.col(j * Tt + t - 1), m, m);
                        arma::mat Theta_curr = reshape(Theta_4d.col(j * Tt + t), m, m);
                        
                        F_local.rows(base_idx, base_idx + m - 1) = 
                            (Theta_prev * Barray_old.slice(t).t()).t();
                        y_local.subvec(base_idx, base_idx + m - 1) = Theta_curr.row(i).t();
                    }
                    
                    // compute posterior with Woodbury identity for large m
                    if (m > 200) {

                        // use Woodbury matrix identity ... thanks cohen
                        arma::mat FtF = F_local.t() * F_local;
                        arma::mat W = inv_tauA2 * eye_m + inv_sigma2 * FtF;
                        V_local = tauA2 * (eye_m - tauA2 * inv_sigma2 * 
                                          solve(W, FtF, solve_opts::likely_sympd));
                    } else {
                        arma::mat V_inv = inv_sigma2 * (F_local.t() * F_local) + inv_tauA2 * eye_m;
                        V_local = inv_sympd(V_inv);
                    }
                    
                    V_local = 0.5 * (V_local + V_local.t());
                    arma::vec m_post = V_local * (inv_sigma2 * (F_local.t() * y_local));
                    
                    if(ar1 && t > 1) {
                        m_post += (rhoA * inv_tauA2) * (V_local * Aarray.slice(t - 1).row(i).t());
                    }
                    
                    // sample new row
                    arma::vec a_new = mvnrnd(m_post, V_local);
                    Aarray.slice(t).row(i) = a_new.t();
                }
            }
        }
    }
    #else
    // serial fallback
    for(int i = 0; i < m; i++) {
        for(int t = 1; t < Tt; t++) {
            arma::mat F_it(p * m, m);
            arma::vec y_it(p * m);
            
            // build design matrix
            for(int j = 0; j < p; j++) {
                int base_idx = j * m;
                
                arma::mat Theta_prev = reshape(Theta_4d.col(j * Tt + t - 1), m, m);
                arma::mat Theta_curr = reshape(Theta_4d.col(j * Tt + t), m, m);
                
                F_it.rows(base_idx, base_idx + m - 1) = 
                    (Theta_prev * Barray_old.slice(t).t()).t();
                y_it.subvec(base_idx, base_idx + m - 1) = Theta_curr.row(i).t();
            }
            
            // compute posterior with Woodbury identity for large m
            arma::mat V;
            if (m > 200) {

                arma::mat FtF = F_it.t() * F_it;
                arma::mat W = inv_tauA2 * eye_m + inv_sigma2 * FtF;
                V = tauA2 * (eye_m - tauA2 * inv_sigma2 * 
                            solve(W, FtF, solve_opts::likely_sympd));
            } else {
                arma::mat V_inv = inv_sigma2 * (F_it.t() * F_it) + inv_tauA2 * eye_m;
                V = inv_sympd(V_inv);
            }
            
            V = 0.5 * (V + V.t());
            arma::vec m_post = V * (inv_sigma2 * (F_it.t() * y_it));
            
            if(ar1 && t > 1) {
                m_post += (rhoA * inv_tauA2) * (V * Aarray.slice(t - 1).row(i).t());
            }
            
            // sample new row
            arma::vec a_new = mvnrnd(m_post, V);
            Aarray.slice(t).row(i) = a_new.t();
        }
    }
    #endif
    
    // set A_1 = I
    Aarray.slice(0) = eye_m;
    
    // update B with same tiled approach
    #ifdef _OPENMP
    #pragma omp parallel
    {
        arma::mat F_local(p * m, m);
        arma::vec y_local(p * m);
        arma::mat V_local(m, m);
        
        #pragma omp for schedule(dynamic, 1)
        for(int k_tile = 0; k_tile < m; k_tile += tile_size) {
            int k_end = std::min(k_tile + tile_size, m);
            
            for(int k = k_tile; k < k_end; k++) {
                for(int t = 1; t < Tt; t++) {
                    for(int j = 0; j < p; j++) {
                        int base_idx = j * m;
                        
                        arma::mat Theta_prev = reshape(Theta_4d.col(j * Tt + t - 1), m, m);
                        arma::mat Theta_curr = reshape(Theta_4d.col(j * Tt + t), m, m);
                        
                        F_local.rows(base_idx, base_idx + m - 1) = Aarray.slice(t) * Theta_prev;
                        y_local.subvec(base_idx, base_idx + m - 1) = Theta_curr.col(k);
                    }
                    
                    // use Woodbury for large m
                    if (m > 200) {
                        arma::mat FtF = F_local.t() * F_local;
                        arma::mat W = inv_tauB2 * eye_m + inv_sigma2 * FtF;
                        V_local = tauB2 * (eye_m - tauB2 * inv_sigma2 * 
                                          solve(W, FtF, solve_opts::likely_sympd));
                    } else {
                        arma::mat V_inv = inv_sigma2 * (F_local.t() * F_local) + inv_tauB2 * eye_m;
                        V_local = inv_sympd(V_inv);
                    }
                    
                    V_local = 0.5 * (V_local + V_local.t());
                    arma::vec m_post = V_local * (inv_sigma2 * (F_local.t() * y_local));
                    
                    if(ar1 && t > 1) {
                        m_post += (rhoB * inv_tauB2) * (V_local * Barray.slice(t - 1).col(k));
                    }
                    
                    arma::vec b_new = mvnrnd(m_post, V_local);
                    Barray.slice(t).col(k) = b_new;
                }
            }
        }
    }
    #endif
    
    // set B_1 = I
    Barray.slice(0) = eye_m;
    
    return List::create(
        Named("Aarray") = Aarray,
        Named("Barray") = Barray
    );
}

// Variance computation for large networks using blocking
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_process_variance_blocked(const arma::mat& Theta_4d,
                                       const arma::cube& Aarray,
                                       const arma::cube& Barray,
                                       int m, int p, int Tt) {
    double proc_rss = 0.0;
    const int block_size = 64; 
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:proc_rss) schedule(static)
    #endif
    for(int j = 0; j < p; j++) {
        double local_rss = 0.0;
        
        // process time in blocks 
        for(int t_block = 1; t_block < Tt; t_block += block_size) {
            int t_end = std::min(t_block + block_size, Tt);
            
            for(int t = t_block; t < t_end; t++) {
                int idx_curr = j * Tt + t;
                int idx_prev = j * Tt + t - 1;
                
                arma::mat Theta_curr = reshape(Theta_4d.col(idx_curr), m, m);
                arma::mat Theta_prev = reshape(Theta_4d.col(idx_prev), m, m);
 
                // magic mult time ............. dont blow up
                arma::mat pred = Aarray.slice(t) * Theta_prev * Barray.slice(t).t();
                local_rss += accu(square(Theta_curr - pred));
            }
        }
        
        proc_rss += local_rss;
    }
    
    return proc_rss;
}