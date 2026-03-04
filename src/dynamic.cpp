#include <RcppArmadillo.h>
#include "dbn_stability.h"
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "thread_control.h"

using namespace Rcpp;
using namespace arma;

// Zero-copy view of a column of a 2D matrix as an n_row x n_col matrix.
// The returned mat shares memory with the source — no allocation or copy.
static inline arma::mat col_as_mat(const arma::mat& M, int col, int n_row, int n_col) {
    return arma::mat(const_cast<double*>(M.colptr(col)), n_row, n_col, false, true);
}

// Forward declarations for existing functions
arma::vec rz_fc_cpp(const arma::vec& R, const arma::vec& Z, const arma::vec& EZ, const List& iranks);
arma::cube ffbs_theta_struct_5arg_cpp(const arma::cube& Z, const arma::mat& mu,
                                     const arma::cube& A_array, const arma::cube& B_array,
                                     double sigma2);

// Safe symmetric positive definite inverse with regularization fallback
static arma::mat safe_inv_sympd(const arma::mat& M) {
    arma::mat M_sym = 0.5 * (M + M.t());
    arma::mat result;
    bool ok = arma::inv_sympd(result, M_sym);
    if (!ok) {
        double reg = 1e-6 * arma::norm(M_sym, "fro") + 1e-8;
        M_sym.diag() += reg;
        ok = arma::inv_sympd(result, M_sym);
        if (!ok) {
            result = arma::inv(M_sym);
        }
    }
    return result;
}

// Safe mvnrnd with regularization
static arma::vec safe_mvnrnd(const arma::vec& mu, const arma::mat& Sigma) {
    arma::mat S = 0.5 * (Sigma + Sigma.t());
    try {
        return arma::mvnrnd(mu, S);
    } catch (...) {
        double reg = 1e-6 * arma::norm(S, "fro") + 1e-8;
        S.diag() += reg;
        try {
            return arma::mvnrnd(mu, S);
        } catch (...) {
            return mu + arma::sqrt(arma::abs(S.diag())) % arma::randn(mu.n_elem);
        }
    }
}

// Batch update z for ordinal data
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat batch_update_Z_ordinal(const arma::mat& R_4d,
                                 const arma::mat& Z_4d,
                                 const arma::mat& Theta_4d,
                                 const arma::cube& M,
                                 const List& IR,
                                 int n_row, int n_col, int p, int Tt) {
    int nc = n_row * n_col;
    // create output array
    arma::mat Z_new(nc, p * Tt);

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
                    int time_idx = global_idx / nc;
                    if(time_idx == t) {
                        int spatial_idx = global_idx % nc;
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
                                       int n_row, int n_col, int p, int Tt) {
    set_dbn_threads(); // set threads from R options
    int nc = n_row * n_col;

    // allocate output once
    arma::mat Z_new(nc, p * Tt);

    // preallocate workspace for each thread
    #ifdef _OPENMP
    int n_threads = omp_get_max_threads();
    std::vector<arma::vec> workspace_R(n_threads, arma::vec(nc));
    std::vector<arma::vec> workspace_Z(n_threads, arma::vec(nc));
    std::vector<arma::vec> workspace_EZ(n_threads, arma::vec(nc));
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
        arma::vec local_R(nc);
        arma::vec local_Z(nc);
        arma::vec local_EZ(nc);
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
List precompute_time_indices(const List& IR, int n_row, int n_col, int p, int Tt) {
    int nc = n_row * n_col;
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
                arma::uvec time_mask = (arma::floor((idx_vec - 1) / nc) == t);
                arma::vec filtered_global = idx_vec.elem(find(time_mask));

                if(filtered_global.n_elem > 0) {
                    // convert to local indices
                    arma::uvec local_indices(filtered_global.n_elem);
                    for(unsigned int i = 0; i < filtered_global.n_elem; i++) {
                        int global_idx = filtered_global(i) - 1;
                        local_indices(i) = (global_idx % nc) + 1;
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

// mu update with blocked computation for large networks
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_mu_dynamic(const arma::mat& Z_4d,
                              const arma::mat& Theta_4d,
                              double g2,
                              double a_g, double b_g,
                              int n_row, int n_col, int p, int Tt) {
    int nc = n_row * n_col;
    double mu_var = 1.0 / (Tt + 1.0 / g2);
    arma::cube M(n_row, n_col, p);

    // process relations sequentially (noise.randn() is not thread-safe)
    for(int j = 0; j < p; j++) {
        // blocked summation for numerical stability with large Tt
        arma::mat sum_diff(n_row, n_col, fill::zeros);

        // sum over time
        for(int t = 0; t < Tt; t++) {
            int col_idx = j * Tt + t;
            sum_diff += col_as_mat(Z_4d, col_idx, n_row, n_col) - col_as_mat(Theta_4d, col_idx, n_row, n_col);
        }

        // compute posterior mean
        arma::mat mu_hat = mu_var * sum_diff;

        // sample from posterior
        arma::mat noise(n_row, n_col);
        noise.randn();
        M.slice(j) = mu_hat + sqrt(mu_var) * noise;
    }

    // update g2 with stable computation
    double M_sum_sq = 0.0;
    for(int j = 0; j < p; j++) {
        M_sum_sq += accu(square(M.slice(j)));
    }

    double shape = (a_g + nc * p) / 2.0;
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
                           int n_row, int n_col, int p, int Tt) {
    arma::cube Theta_4d_out(n_row, n_col, p * Tt);

    // sequential over relations (FFBS calls arma::randn which is not thread-safe)
    for(int j = 0; j < p; j++) {
        // extract data for this relation
        arma::cube Z_j(n_row, n_col, Tt);
        int nc = n_row * n_col;

        // vectorized extraction
        for(int t = 0; t < Tt; t++) {
            Z_j.slice(t) = col_as_mat(Z_4d, j * Tt + t, n_row, n_col);
        }

        // run ffbs — dimensions propagate automatically from Z_j
        arma::cube Theta_j = ffbs_theta_struct_5arg_cpp(Z_j, M.slice(j), Aarray, Barray, sigma2);

        // vectorized storage
        for(int t = 0; t < Tt; t++) {
            Theta_4d_out.slice(j * Tt + t) = Theta_j.slice(t);
        }
    }

    return Theta_4d_out;
}

// AB update for networks
// A is n_row x n_row (sender dynamics), B is n_col x n_col (receiver dynamics)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_AB_batch_extended(const arma::mat& Theta_4d,
                    const arma::cube& Aarray_old,
                    const arma::cube& Barray_old,
                    double sigma2, double tauA2, double tauB2,
                    bool ar1, double rhoA, double rhoB,
                    int n_row, int n_col, int p, int Tt) {
    int nc = n_row * n_col;
    arma::cube Aarray(n_row, n_row, Tt);
    arma::cube Barray(n_col, n_col, Tt);

    // precompute frequently used matrices
    arma::mat eye_nr = eye(n_row, n_row);
    arma::mat eye_nc = eye(n_col, n_col);
    double inv_sigma2 = 1.0 / sigma2;
    double inv_tauA2 = 1.0 / tauA2;
    double inv_tauB2 = 1.0 / tauB2;

    // --- Update A (sender dynamics, n_row x n_row) ---
    // For row i of A: Theta_t(i,:) = A(i,:) * Theta_{t-1} * B_t'
    // Response: row i of Theta_t -> length n_col
    // Design: (Theta_{t-1} * B_t')' -> n_col x n_row
    // Parameter: A(i,:) -> length n_row
    // sequential: safe_mvnrnd is not thread-safe
    for(int i = 0; i < n_row; i++) {
        // pre-allocate workspace outside t-loop
        arma::mat F_it(p * n_col, n_row);
        arma::vec y_it(p * n_col);
        for(int t = 1; t < Tt; t++) {

            // build design matrix with better memory access
            for(int j = 0; j < p; j++) {
                int base_idx = j * n_col;
                int col_idx_prev = j * Tt + t - 1;
                int col_idx_curr = j * Tt + t;

                // extract and reshape
                arma::mat Theta_prev = col_as_mat(Theta_4d, col_idx_prev, n_row, n_col);
                arma::mat Theta_curr = col_as_mat(Theta_4d, col_idx_curr, n_row, n_col);

                // Design: (Theta_{t-1} * B_t')' = B_t * Theta_{t-1}' -> n_col x n_row
                F_it.rows(base_idx, base_idx + n_col - 1) = (Theta_prev * Barray_old.slice(t).t()).t();
                // Response: row i of Theta_t -> n_col x 1
                y_it.subvec(base_idx, base_idx + n_col - 1) = Theta_curr.row(i).t();
            }

            // compute posterior with pre-computed constants
            // F_it' * F_it: (n_row x p*n_col) * (p*n_col x n_row) = n_row x n_row
            arma::mat V_inv = inv_sigma2 * (F_it.t() * F_it) + inv_tauA2 * eye_nr;
            arma::mat V = safe_inv_sympd(V_inv);
            V = 0.5 * (V + V.t());
            arma::vec m_post = V * (inv_sigma2 * (F_it.t() * y_it));

            if(ar1 && t > 1) {
                m_post += (rhoA * inv_tauA2) * (V * Aarray.slice(t - 1).row(i).t());
            }

            // sample new row
            arma::vec a_new = safe_mvnrnd(m_post, V);
            Aarray.slice(t).row(i) = a_new.t();
        }
    }

    // set a_1 = I
    Aarray.slice(0) = eye_nr;

    // --- Update B (receiver dynamics, n_col x n_col) ---
    // For column k of B: Theta_t(:,k) = A_t * Theta_{t-1} * B(k,:)'
    // Response: column k of Theta_t -> length n_row
    // Design: A_t * Theta_{t-1} -> n_row x n_col
    // Parameter: B(k,:)' -> length n_col
    // sequential: safe_mvnrnd is not thread-safe
    for(int k = 0; k < n_col; k++) {
        // pre-allocate workspace outside t-loop
        arma::mat F_kt(p * n_row, n_col);
        arma::vec y_kt(p * n_row);
        for(int t = 1; t < Tt; t++) {

            for(int j = 0; j < p; j++) {
                int base_idx = j * n_row;
                int col_idx_prev = j * Tt + t - 1;
                int col_idx_curr = j * Tt + t;

                arma::mat Theta_prev = col_as_mat(Theta_4d, col_idx_prev, n_row, n_col);
                arma::mat Theta_curr = col_as_mat(Theta_4d, col_idx_curr, n_row, n_col);

                // Design: A_t * Theta_{t-1} -> n_row x n_col
                F_kt.rows(base_idx, base_idx + n_row - 1) = Aarray.slice(t) * Theta_prev;
                // Response: column k of Theta_t -> n_row x 1
                y_kt.subvec(base_idx, base_idx + n_row - 1) = Theta_curr.col(k);
            }

            // F_kt' * F_kt: (n_col x p*n_row) * (p*n_row x n_col) = n_col x n_col
            arma::mat V_inv = inv_sigma2 * (F_kt.t() * F_kt) + inv_tauB2 * eye_nc;
            arma::mat V = safe_inv_sympd(V_inv);
            V = 0.5 * (V + V.t());
            arma::vec m_post = V * (inv_sigma2 * (F_kt.t() * y_kt));

            if(ar1 && t > 1) {
                m_post += (rhoB * inv_tauB2) * (V * Barray.slice(t - 1).col(k));
            }

            arma::vec b_new = safe_mvnrnd(m_post, V);
            Barray.slice(t).col(k) = b_new;
        }
    }

    // set b_1 = I
    Barray.slice(0) = eye_nc;

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
                           int n_row, int n_col, int p, int Tt,
                           bool is_gaussian = false) {
    int nc = n_row * n_col;
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
            arma::mat Theta_curr = col_as_mat(Theta_4d, idx_curr, n_row, n_col);
            arma::mat Theta_prev = col_as_mat(Theta_4d, idx_prev, n_row, n_col);

            // compute residual = Theta_curr - A_t * Theta_prev * B_t'
            arma::mat pred = Aarray.slice(t) * Theta_prev * Barray.slice(t).t();
            local_rss += accu(square(Theta_curr - pred));
        }

        proc_rss += local_rss;
    }

    double sigma2 = (b_sig + proc_rss / 2.0) / R::rgamma((a_sig + nc * (Tt - 1) * p) / 2.0, 1.0);

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

                arma::mat Z_jt = col_as_mat(Z_4d, idx, n_row, n_col);
                arma::mat Theta_jt = col_as_mat(Theta_4d, idx, n_row, n_col);

                local_rss += accu(square(Z_jt - (Theta_jt + M_j)));
            }

            obs_rss += local_rss;
        }

        sigma2_obs = (1.0 + obs_rss / 2.0) / R::rgamma((1.0 + nc * Tt * p) / 2.0, 1.0);
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
                                           int n_row, int n_col, int p, int Tt) {
    arma::cube Theta_4d_out(n_row, n_col, p * Tt);

    // process relations in blocks
    const int block_size = std::max(1, std::min(4, p));

    // sequential: FFBS calls arma::randn which is not thread-safe
    for(int j_block = 0; j_block < p; j_block += block_size) {
        int j_end = std::min(j_block + block_size, p);

        // process block of relations
        for(int j = j_block; j < j_end; j++) {
            // preallocate workspace for this relation
            arma::cube Z_j(n_row, n_col, Tt);

            // vectorized extraction
            for(int t = 0; t < Tt; t++) {
                Z_j.slice(t) = col_as_mat(Z_4d, j * Tt + t, n_row, n_col);
            }

            // ffbs
            arma::cube Theta_j = ffbs_theta_struct_5arg_cpp(Z_j, M.slice(j), Aarray, Barray, sigma2);

            // store
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
                          int n_row, int n_col, int p, int Tt) {
    int nc = n_row * n_col;
    arma::cube Aarray(n_row, n_row, Tt);
    arma::cube Barray(n_col, n_col, Tt);

    // precompute frequently used matrices
    arma::mat eye_nr = eye(n_row, n_row);
    arma::mat eye_nc = eye(n_col, n_col);
    double inv_sigma2 = 1.0 / sigma2;
    double inv_tauA2 = 1.0 / tauA2;
    double inv_tauB2 = 1.0 / tauB2;

    // use tiled approach for very large networks
    const int tile_size_A = std::min(32, n_row);
    const int tile_size_B = std::min(32, n_col);

    // update A with tiled computation
    // Note: parallel disabled because safe_mvnrnd is not thread-safe
    #ifdef _OPENMP
    {
        arma::mat F_local(p * n_col, n_row);
        arma::vec y_local(p * n_col);
        arma::mat V_local(n_row, n_row);

        for(int i_tile = 0; i_tile < n_row; i_tile += tile_size_A) {
            int i_end = std::min(i_tile + tile_size_A, n_row);

            for(int i = i_tile; i < i_end; i++) {
                for(int t = 1; t < Tt; t++) {
                    // build design matrix
                    for(int j = 0; j < p; j++) {
                        int base_idx = j * n_col;

                        arma::mat Theta_prev = col_as_mat(Theta_4d, j * Tt + t - 1, n_row, n_col);
                        arma::mat Theta_curr = col_as_mat(Theta_4d, j * Tt + t, n_row, n_col);

                        F_local.rows(base_idx, base_idx + n_col - 1) =
                            (Theta_prev * Barray_old.slice(t).t()).t();
                        y_local.subvec(base_idx, base_idx + n_col - 1) = Theta_curr.row(i).t();
                    }

                    // compute posterior with Woodbury identity for large n_row
                    if (n_row > 50) {
                        arma::mat FtF = F_local.t() * F_local;
                        arma::mat W = inv_tauA2 * eye_nr + inv_sigma2 * FtF;
                        V_local = tauA2 * (eye_nr - tauA2 * inv_sigma2 *
                                          solve(W, FtF, solve_opts::likely_sympd));
                    } else {
                        arma::mat V_inv = inv_sigma2 * (F_local.t() * F_local) + inv_tauA2 * eye_nr;
                        V_local = safe_inv_sympd(V_inv);
                    }

                    V_local = 0.5 * (V_local + V_local.t());
                    arma::vec m_post = V_local * (inv_sigma2 * (F_local.t() * y_local));

                    if(ar1 && t > 1) {
                        m_post += (rhoA * inv_tauA2) * (V_local * Aarray.slice(t - 1).row(i).t());
                    }

                    // sample new row
                    arma::vec a_new = safe_mvnrnd(m_post, V_local);
                    Aarray.slice(t).row(i) = a_new.t();
                }
            }
        }
    }
    #else
    // serial fallback
    for(int i = 0; i < n_row; i++) {
        arma::mat F_it(p * n_col, n_row);
        arma::vec y_it(p * n_col);
        for(int t = 1; t < Tt; t++) {

            // build design matrix
            for(int j = 0; j < p; j++) {
                int base_idx = j * n_col;

                arma::mat Theta_prev = col_as_mat(Theta_4d, j * Tt + t - 1, n_row, n_col);
                arma::mat Theta_curr = col_as_mat(Theta_4d, j * Tt + t, n_row, n_col);

                F_it.rows(base_idx, base_idx + n_col - 1) =
                    (Theta_prev * Barray_old.slice(t).t()).t();
                y_it.subvec(base_idx, base_idx + n_col - 1) = Theta_curr.row(i).t();
            }

            // compute posterior with Woodbury identity for large n_row
            arma::mat V;
            if (n_row > 50) {
                arma::mat FtF = F_it.t() * F_it;
                arma::mat W = inv_tauA2 * eye_nr + inv_sigma2 * FtF;
                V = tauA2 * (eye_nr - tauA2 * inv_sigma2 *
                            solve(W, FtF, solve_opts::likely_sympd));
            } else {
                arma::mat V_inv = inv_sigma2 * (F_it.t() * F_it) + inv_tauA2 * eye_nr;
                V = safe_inv_sympd(V_inv);
            }

            V = 0.5 * (V + V.t());
            arma::vec m_post = V * (inv_sigma2 * (F_it.t() * y_it));

            if(ar1 && t > 1) {
                m_post += (rhoA * inv_tauA2) * (V * Aarray.slice(t - 1).row(i).t());
            }

            // sample new row
            arma::vec a_new = safe_mvnrnd(m_post, V);
            Aarray.slice(t).row(i) = a_new.t();
        }
    }
    #endif

    // set A_1 = I
    Aarray.slice(0) = eye_nr;

    // update B with same tiled approach
    // Note: parallel disabled because safe_mvnrnd is not thread-safe
    #ifdef _OPENMP
    {
        arma::mat F_local(p * n_row, n_col);
        arma::vec y_local(p * n_row);
        arma::mat V_local(n_col, n_col);

        for(int k_tile = 0; k_tile < n_col; k_tile += tile_size_B) {
            int k_end = std::min(k_tile + tile_size_B, n_col);

            for(int k = k_tile; k < k_end; k++) {
                for(int t = 1; t < Tt; t++) {
                    for(int j = 0; j < p; j++) {
                        int base_idx = j * n_row;

                        arma::mat Theta_prev = col_as_mat(Theta_4d, j * Tt + t - 1, n_row, n_col);
                        arma::mat Theta_curr = col_as_mat(Theta_4d, j * Tt + t, n_row, n_col);

                        F_local.rows(base_idx, base_idx + n_row - 1) = Aarray.slice(t) * Theta_prev;
                        y_local.subvec(base_idx, base_idx + n_row - 1) = Theta_curr.col(k);
                    }

                    // use Woodbury for large n_col
                    if (n_col > 50) {
                        arma::mat FtF = F_local.t() * F_local;
                        arma::mat W = inv_tauB2 * eye_nc + inv_sigma2 * FtF;
                        V_local = tauB2 * (eye_nc - tauB2 * inv_sigma2 *
                                          solve(W, FtF, solve_opts::likely_sympd));
                    } else {
                        arma::mat V_inv = inv_sigma2 * (F_local.t() * F_local) + inv_tauB2 * eye_nc;
                        V_local = safe_inv_sympd(V_inv);
                    }

                    V_local = 0.5 * (V_local + V_local.t());
                    arma::vec m_post = V_local * (inv_sigma2 * (F_local.t() * y_local));

                    if(ar1 && t > 1) {
                        m_post += (rhoB * inv_tauB2) * (V_local * Barray.slice(t - 1).col(k));
                    }

                    arma::vec b_new = safe_mvnrnd(m_post, V_local);
                    Barray.slice(t).col(k) = b_new;
                }
            }
        }
    }
    #else
    // serial fallback for B
    for(int k = 0; k < n_col; k++) {
        arma::mat F_kt(p * n_row, n_col);
        arma::vec y_kt(p * n_row);
        for(int t = 1; t < Tt; t++) {

            for(int j = 0; j < p; j++) {
                int base_idx = j * n_row;

                arma::mat Theta_prev = col_as_mat(Theta_4d, j * Tt + t - 1, n_row, n_col);
                arma::mat Theta_curr = col_as_mat(Theta_4d, j * Tt + t, n_row, n_col);

                F_kt.rows(base_idx, base_idx + n_row - 1) = Aarray.slice(t) * Theta_prev;
                y_kt.subvec(base_idx, base_idx + n_row - 1) = Theta_curr.col(k);
            }

            arma::mat V;
            if (n_col > 50) {
                arma::mat FtF = F_kt.t() * F_kt;
                arma::mat W = inv_tauB2 * eye_nc + inv_sigma2 * FtF;
                V = tauB2 * (eye_nc - tauB2 * inv_sigma2 *
                            solve(W, FtF, solve_opts::likely_sympd));
            } else {
                arma::mat V_inv = inv_sigma2 * (F_kt.t() * F_kt) + inv_tauB2 * eye_nc;
                V = safe_inv_sympd(V_inv);
            }

            V = 0.5 * (V + V.t());
            arma::vec m_post = V * (inv_sigma2 * (F_kt.t() * y_kt));

            if(ar1 && t > 1) {
                m_post += (rhoB * inv_tauB2) * (V * Barray.slice(t - 1).col(k));
            }

            arma::vec b_new = safe_mvnrnd(m_post, V);
            Barray.slice(t).col(k) = b_new;
        }
    }
    #endif

    // set B_1 = I
    Barray.slice(0) = eye_nc;

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
                                       int n_row, int n_col, int p, int Tt) {
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

                arma::mat Theta_curr = col_as_mat(Theta_4d, idx_curr, n_row, n_col);
                arma::mat Theta_prev = col_as_mat(Theta_4d, idx_prev, n_row, n_col);

                arma::mat pred = Aarray.slice(t) * Theta_prev * Barray.slice(t).t();
                local_rss += accu(square(Theta_curr - pred));
            }
        }

        proc_rss += local_rss;
    }

    return proc_rss;
}
