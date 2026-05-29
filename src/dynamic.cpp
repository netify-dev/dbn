#include <RcppArmadillo.h>
#include "dbn_stability.h"
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include <random>
#include <vector>
#include "thread_control.h"

using namespace Rcpp;
using namespace arma;

// zero-copy view of a 2D matrix column as an n_row x n_col matrix
// shares memory with the source, no allocation or copy
static inline arma::mat col_as_mat(const arma::mat& M, int col, int n_row, int n_col) {
    return arma::mat(const_cast<double*>(M.colptr(col)), n_row, n_col, false, true);
}

// forward declarations
arma::vec rz_fc_cpp(const arma::vec& R, const arma::vec& Z, const arma::vec& EZ, const List& iranks);
arma::cube ffbs_theta_struct_5arg_cpp(const arma::cube& Z, const arma::mat& mu,
                                     const arma::cube& A_array, const arma::cube& B_array,
                                     double sigma2);

// safe symmetric PD inverse: delegates to shared helper in dbn_stability.h
static inline arma::mat safe_inv_sympd(const arma::mat& M) {
    return dbn_safe_inv_sympd(M);
}


// thread-safe multivariate normal sampling with thread-local RNG
// uses manual Cholesky + std::mt19937_64 instead of arma::mvnrnd
static arma::vec thread_safe_mvnrnd(const arma::vec& mu, const arma::mat& Sigma,
                                     std::mt19937_64& rng) {
    int d = mu.n_elem;
    arma::mat S = 0.5 * (Sigma + Sigma.t());

    // draw standard normals from thread-local RNG
    std::normal_distribution<double> norm(0.0, 1.0);
    arma::vec z(d);
    for (int i = 0; i < d; i++) {
        z(i) = norm(rng);
    }

    // try cholesky decomposition
    arma::mat L;
    bool ok = arma::chol(L, S, "lower");
    if (ok) {
        return mu + L * z;
    }

    // regularize and retry
    double reg = 1e-6 * arma::norm(S, "fro") + 1e-8;
    S.diag() += reg;
    ok = arma::chol(L, S, "lower");
    if (ok) {
        return mu + L * z;
    }

    // fallback to diagonal sampling
    return mu + arma::sqrt(arma::abs(S.diag())) % z;
}

// initialize n_units deterministic RNG engines derived from R's RNG.
// must be called from the main thread before any parallel region.
//
// IMPORTANT: this consumes R's RNG exactly ONCE (one R::runif() call) and then
// derives every per-unit RNG seed from that single base via SplitMix64. The
// caller should pass `n_units = number of independent draws to make` (e.g.
// n_row for an A-row update), NOT the OpenMP thread count -- then index into
// the returned vector with the LOOP variable (i, k, ...) so each independent
// draw uses its OWN deterministic RNG regardless of which thread happens to
// run it. Earlier versions keyed the RNG by `omp_get_thread_num()`, which
// silently rerouted the RNG stream when `set_dbn_threads()` changed.
static std::vector<std::mt19937_64> init_thread_rngs(int n_units) {
    std::vector<std::mt19937_64> rngs;
    rngs.reserve(n_units);
    uint64_t base = static_cast<uint64_t>(R::runif(0.0, 1.0) * 4294967296.0);
    for (int i = 0; i < n_units; i++) {
        // SplitMix64 mix of (base, unit_id): deterministic, well-distributed
        uint64_t z = base + static_cast<uint64_t>(i + 1) * 0x9E3779B97F4A7C15ULL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z = z ^ (z >> 31);
        rngs.emplace_back(z);
    }
    return rngs;
}

// batch update z for ordinal data
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
    // output array
    arma::mat Z_new(nc, p * Tt);

    // loop over relations
    for(int j = 0; j < p; j++) {
        // rank indices for this relation
        List IR_j = IR[j];

        // loop over time points
        for(int t = 0; t < Tt; t++) {
            // extract column vectors
            int col_idx = j * Tt + t;
            arma::vec R_vec = R_4d.col(col_idx);
            arma::vec Z_vec = Z_4d.col(col_idx);
            arma::vec Theta_vec = Theta_4d.col(col_idx);

            // expected z = Theta + M for this relation
            arma::mat M_j = M.slice(j);
            arma::vec M_vec = vectorise(M_j);
            arma::vec EZ_vec = Theta_vec + M_vec;

            // build rank index list for this time point
            List IR_t;
            CharacterVector names = IR_j.names();
            for(int k = 0; k < IR_j.size(); k++) {
                std::string rank_name = Rcpp::as<std::string>(names[k]);
                arma::vec idx_vec = Rcpp::as<arma::vec>(IR_j[k]);

                // keep indices belonging to this time point
                arma::uvec local_indices;
                for(unsigned int i = 0; i < idx_vec.n_elem; i++) {
                    int global_idx = idx_vec(i) - 1; // to 0-based
                    int time_idx = global_idx / nc;
                    if(time_idx == t) {
                        int spatial_idx = global_idx % nc;
                        local_indices.resize(local_indices.n_elem + 1);
                        local_indices(local_indices.n_elem - 1) = spatial_idx + 1; // back to 1-based
                    }
                }

                if(local_indices.n_elem > 0) {
                    IR_t[rank_name] = local_indices;
                } else {
                    IR_t[rank_name] = arma::uvec();
                }
            }

            // update z via rank likelihood
            arma::vec Z_updated = rz_fc_cpp(R_vec, Z_vec, EZ_vec, IR_t);

            // store updated column
            Z_new.col(col_idx) = Z_updated;
        }
    }

    return Z_new;
}

// z update with preallocated memory and better cache usage
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
    int nc = n_row * n_col;

    // allocate output once
    arma::mat Z_new(nc, p * Tt);

    // no OpenMP here -- rz_fc_cpp touches R SEXP objects
    // which aren't thread-safe (corrupts R's PROTECT stack)
    for(int j = 0; j < p; j++) {
        // precomputed time indices for this relation
        List IR_j_time = IR_time_indices[j];
        arma::mat M_j = M.slice(j);
        arma::vec M_vec = vectorise(M_j);

        // loop over time with cache-friendly access
        for(int t = 0; t < Tt; t++) {
            int col_idx = j * Tt + t;

            // precomputed indices for this time point
            List IR_t = IR_j_time[t];

            // expected z
            arma::vec EZ_vec = Theta_4d.col(col_idx) + M_vec;

            // rank likelihood z update
            Z_new.col(col_idx) = rz_fc_cpp(R_4d.col(col_idx), Z_4d.col(col_idx), EZ_vec, IR_t);
        }
    }

    return Z_new;
}

// precompute time-specific rank indices to avoid repeated work
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

        // split indices by time point
        for(int t = 0; t < Tt; t++) {
            List IR_t;

            for(int k = 0; k < IR_j.size(); k++) {
                std::string rank_name = Rcpp::as<std::string>(names[k]);
                arma::vec idx_vec = Rcpp::as<arma::vec>(IR_j[k]);

                // filter to indices at time t
                arma::uvec time_mask = (arma::floor((idx_vec - 1) / nc) == t);
                arma::vec filtered_global = idx_vec.elem(find(time_mask));

                if(filtered_global.n_elem > 0) {
                    // map global to local spatial indices
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

// mu (baseline mean) update
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_mu_dynamic(const arma::mat& Z_4d,
                              const arma::mat& Theta_4d,
                              const arma::cube& M_prev,
                              double g2,
                              double a_g, double b_g,
                              int n_row, int n_col, int p, int Tt) {
    int nc = n_row * n_col;
    double mu_var = 1.0 / (Tt + 1.0 / g2);
    arma::cube M(n_row, n_col, p);

    // sequential: noise.randn() is not thread-safe
    for(int j = 0; j < p; j++) {
        // accumulate z - theta across time
        arma::mat sum_diff(n_row, n_col, fill::zeros);

        // sum over time
        for(int t = 0; t < Tt; t++) {
            int col_idx = j * Tt + t;
            sum_diff += col_as_mat(Z_4d, col_idx, n_row, n_col) - col_as_mat(Theta_4d, col_idx, n_row, n_col);
        }
        // M is the long-run baseline of Theta in the model
        // Theta_t = M + A_t (Theta_{t-1} - M) B_t' + noise. Identifying M
        // separately from Theta requires the bilinear-residual conditional,
        // which is not implemented here -- instead Theta absorbs the baseline
        // and M shrinks toward the prior mean (0). See known limitations.

        // compute posterior mean
        arma::mat mu_hat = mu_var * sum_diff;

        // sample from posterior
        arma::mat noise(n_row, n_col);
        noise.randn();
        M.slice(j) = mu_hat + sqrt(mu_var) * noise;
    }

    // update g2 (prior variance on M)
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

// FFBS across all relations
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

    // sequential: FFBS uses arma::randn which isn't thread-safe
    for(int j = 0; j < p; j++) {
        // gather slices for this relation
        arma::cube Z_j(n_row, n_col, Tt);
        int nc = n_row * n_col;

        // copy into contiguous cube
        for(int t = 0; t < Tt; t++) {
            Z_j.slice(t) = col_as_mat(Z_4d, j * Tt + t, n_row, n_col);
        }

        // run FFBS for this relation
        arma::cube Theta_j = ffbs_theta_struct_5arg_cpp(Z_j, M.slice(j), Aarray, Barray, sigma2);

        // write back to output
        for(int t = 0; t < Tt; t++) {
            Theta_4d_out.slice(j * Tt + t) = Theta_j.slice(t);
        }
    }

    return Theta_4d_out;
}

// A/B update for networks (asymmetric / directed path).
//
// Three prior modes are supported, selected by `prior_kind`:
//   "rw"  (default): random-walk prior anchored at identity.
//        A_0 = I (pinned), A_t = A_{t-1} + eps_A_t, eps_A_t ~ N(0, tauA2 I).
//        Same for B. The anchor at t=0 resolves the (A, B) -> (cA, c^-1 B)
//        scale ambiguity that broke prior RW attempts -- without the anchor,
//        the chain drifts in scale and tauA2 blows up.
//   "iid" (legacy): iid zero-mean prior A_t ~ N(0, tauA2 I), independent across t.
//        Use only for backward compatibility; tau_A^2 is then a coarse
//        flexibility scale, not an innovation variance.
//   "ar1" (existing): AR(1) with persistence rho around zero.
//        A_t = rho * A_{t-1} + eps, eps ~ N(0, tauA2 I). Activated via the
//        legacy `ar1 = TRUE` boolean for backward compatibility.
//
// On top of any of the above, an optional contractivity (ridge) prior
// A_t ~ N(0, 1/kappaA_inv I) shrinks each operator slice toward zero,
// softly bounding its spectral radius. kappaA_inv = kappaB_inv = 0 (default)
// leaves behaviour unchanged. The per-time operator is under-identified
// (one transition feeds each A_t, B_t), so without this the random-walk
// prior lets the spectral radius drift into the explosive region; the ridge
// only adds to the posterior precision (mean target is zero), so it is a
// pure pull toward contractivity and changes no other code path.
//
// A is n_row x n_row (sender dynamics), B is n_col x n_col (receiver dynamics)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_AB_batch_extended(const arma::mat& Theta_4d,
                    const arma::cube& Aarray_old,
                    const arma::cube& Barray_old,
                    double sigma2, double tauA2, double tauB2,
                    bool ar1, double rhoA, double rhoB,
                    int n_row, int n_col, int p, int Tt,
                    std::string prior_kind = "rw",
                    double kappaA_inv = 0.0, double kappaB_inv = 0.0) {
    int nc = n_row * n_col;
    arma::cube Aarray(n_row, n_row, Tt);
    arma::cube Barray(n_col, n_col, Tt);

    // cache frequently used quantities
    arma::mat eye_nr = eye(n_row, n_row);
    arma::mat eye_nc = eye(n_col, n_col);
    double inv_sigma2 = 1.0 / sigma2;
    double inv_tauA2 = 1.0 / tauA2;
    double inv_tauB2 = 1.0 / tauB2;

    // --- update A (sender dynamics, n_row x n_row) ---
    // each row of A is independent: parallelize with thread-local RNG
    // private per-row storage avoids false sharing on column-major cubes
    set_dbn_threads();
#ifdef _OPENMP
    int n_threads_A = std::min(omp_get_max_threads(), n_row);
#else
    int n_threads_A = 1;
#endif
    auto rngs_A = init_thread_rngs(n_row);

    // A_private[i] holds row i's sampled values across time
    std::vector<arma::mat> A_private(n_row, arma::mat(n_row, Tt, arma::fill::zeros));

    // Prior precision multipliers on inv_tauA2 per t:
    //   iid / ar1: 1 for every t
    //   rw       : 2 for interior t (two RW neighbors), 1 for t = Tt-1 (last,
    //              only A_{T-1} -- A_T conditional reduces to 1-sided).
    //              t=1 is interior because A_0 = I is the other RW neighbor.
    bool use_rw = (prior_kind == "rw") && !ar1;
    // F_t = (Theta_{t-1} * B_t^T)^T stacked across relations is independent
    // of the row index i, so we precompute it (and V_t = (F'F/sigma2 + I/tau)^-1)
    // once per t and share across all i. F'F and the inverse were previously
    // recomputed n_row times. The sampler still calls thread_safe_mvnrnd with
    // the same Sigma argument so the RNG / Cholesky path is bit-identical.
    std::vector<arma::mat> F_t_cache(Tt);
    std::vector<arma::mat> V_t_cache(Tt);
    for(int t = 1; t < Tt; t++) {
        arma::mat F_t(p * n_col, n_row);
        for(int j = 0; j < p; j++) {
            int base_idx = j * n_col;
            int col_idx_prev = j * Tt + t - 1;
            arma::mat Theta_prev = col_as_mat(Theta_4d, col_idx_prev, n_row, n_col);
            F_t.rows(base_idx, base_idx + n_col - 1) =
                (Theta_prev * Barray_old.slice(t).t()).t();
        }
        // Prior precision multiplier: under RW, t < Tt-1 has 2 neighbors so
        // precision is 2/tau^2; t = Tt-1 has only 1 (the previous slice) so
        // precision is 1/tau^2. Under iid / ar1 (legacy), always 1/tau^2.
        double prior_prec_mult = 1.0;
        if (use_rw && t < Tt - 1) prior_prec_mult = 2.0;
        arma::mat V_inv = inv_sigma2 * (F_t.t() * F_t) + (prior_prec_mult * inv_tauA2 + kappaA_inv) * eye_nr;
        arma::mat V_t = safe_inv_sympd(V_inv);
        V_t = 0.5 * (V_t + V_t.t());
        F_t_cache[t] = std::move(F_t);
        V_t_cache[t] = std::move(V_t);
    }

    #pragma omp parallel for schedule(static) num_threads(n_threads_A)
    for(int i = 0; i < n_row; i++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        arma::vec y_it(p * n_col);
        for(int t = 1; t < Tt; t++) {

            for(int j = 0; j < p; j++) {
                int base_idx = j * n_col;
                int col_idx_curr = j * Tt + t;
                arma::mat Theta_curr = col_as_mat(Theta_4d, col_idx_curr, n_row, n_col);
                y_it.subvec(base_idx, base_idx + n_col - 1) = Theta_curr.row(i).t();
            }

            const arma::mat& F_t = F_t_cache[t];
            const arma::mat& V_t = V_t_cache[t];
            // Data contribution to posterior mean (common to all priors).
            arma::vec m_post = V_t * (inv_sigma2 * (F_t.t() * y_it));

            if (use_rw) {
                // RW prior contribution. Read previous neighbor:
                //   t == 1: A_0 = I, so neighbor_prev = row i of I (zero vector
                //           with a 1 in position i)
                //   t  > 1: A_{t-1, i, :} from this iteration's just-sampled
                //           A_private[i].col(t-1)
                arma::vec neighbor_prev(n_row, arma::fill::zeros);
                if (t == 1) {
                    neighbor_prev(i) = 1.0;  // row i of identity
                } else {
                    neighbor_prev = A_private[i].col(t - 1);
                }
                if (t < Tt - 1) {
                    // Interior: prior conditional has two neighbors. Read
                    // A_{t+1, i, :} from `Aarray_old` -- this is the previous-
                    // iteration snapshot (the local A_private buffer holds
                    // only this iteration's per-row results and hasn't been
                    // assembled into a full cube yet). Valid Gibbs ordering:
                    // backward neighbor is current-iter (Gauss-Seidel),
                    // forward neighbor is prev-iter (Jacobi).
                    arma::vec neighbor_next = Aarray_old.slice(t + 1).row(i).t();
                    m_post += V_t * (inv_tauA2 * (neighbor_prev + neighbor_next));
                } else {
                    // Last (t == Tt - 1): only previous neighbor.
                    m_post += V_t * (inv_tauA2 * neighbor_prev);
                }
            } else if (ar1 && t > 1) {
                m_post += (rhoA * inv_tauA2) * (V_t * A_private[i].col(t - 1));
            }
            // Legacy iid prior: zero-mean prior contribution (no addition to
            // m_post). Reached when prior_kind == "iid" and ar1 == FALSE.

            arma::vec a_new = thread_safe_mvnrnd(m_post, V_t, rngs_A[i]);
            A_private[i].col(t) = a_new;
        }
    }

    // copy private storage into output cube
    Aarray.slice(0) = eye_nr;
    for(int t = 1; t < Tt; t++) {
        for(int i = 0; i < n_row; i++) {
            Aarray.slice(t).row(i) = A_private[i].col(t).t();
        }
    }

    // --- update B (receiver dynamics, n_col x n_col) ---
    // private storage for AR(1) correctness and consistency with A
#ifdef _OPENMP
    int n_threads_B = std::min(omp_get_max_threads(), n_col);
#else
    int n_threads_B = 1;
#endif
    auto rngs_B = init_thread_rngs(n_col);

    std::vector<arma::mat> B_private(n_col, arma::mat(n_col, Tt, arma::fill::zeros));

    // Precompute F_kt and V_t once per t (independent of column index k).
    std::vector<arma::mat> FB_t_cache(Tt);
    std::vector<arma::mat> VB_t_cache(Tt);
    for(int t = 1; t < Tt; t++) {
        arma::mat F_t(p * n_row, n_col);
        for(int j = 0; j < p; j++) {
            int base_idx = j * n_row;
            int col_idx_prev = j * Tt + t - 1;
            arma::mat Theta_prev = col_as_mat(Theta_4d, col_idx_prev, n_row, n_col);
            F_t.rows(base_idx, base_idx + n_row - 1) = Aarray.slice(t) * Theta_prev;
        }
        double prior_prec_mult_B = 1.0;
        if (use_rw && t < Tt - 1) prior_prec_mult_B = 2.0;
        arma::mat V_inv = inv_sigma2 * (F_t.t() * F_t) + (prior_prec_mult_B * inv_tauB2 + kappaB_inv) * eye_nc;
        arma::mat V_t = safe_inv_sympd(V_inv);
        V_t = 0.5 * (V_t + V_t.t());
        FB_t_cache[t] = std::move(F_t);
        VB_t_cache[t] = std::move(V_t);
    }

    #pragma omp parallel for schedule(static) num_threads(n_threads_B)
    for(int k = 0; k < n_col; k++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        arma::vec y_kt(p * n_row);
        for(int t = 1; t < Tt; t++) {

            for(int j = 0; j < p; j++) {
                int base_idx = j * n_row;
                int col_idx_curr = j * Tt + t;
                arma::mat Theta_curr = col_as_mat(Theta_4d, col_idx_curr, n_row, n_col);
                y_kt.subvec(base_idx, base_idx + n_row - 1) = Theta_curr.col(k);
            }

            const arma::mat& F_t = FB_t_cache[t];
            const arma::mat& V_t = VB_t_cache[t];
            arma::vec m_post = V_t * (inv_sigma2 * (F_t.t() * y_kt));

            if (use_rw) {
                // B_0 = I anchor. For column k of B at t = 1, neighbor_prev =
                // col k of I (one-hot vector with 1 at position k).
                arma::vec neighbor_prev(n_col, arma::fill::zeros);
                if (t == 1) {
                    neighbor_prev(k) = 1.0;
                } else {
                    neighbor_prev = B_private[k].col(t - 1);
                }
                if (t < Tt - 1) {
                    arma::vec neighbor_next = Barray_old.slice(t + 1).col(k);
                    m_post += V_t * (inv_tauB2 * (neighbor_prev + neighbor_next));
                } else {
                    m_post += V_t * (inv_tauB2 * neighbor_prev);
                }
            } else if (ar1 && t > 1) {
                m_post += (rhoB * inv_tauB2) * (V_t * B_private[k].col(t - 1));
            }

            arma::vec b_new = thread_safe_mvnrnd(m_post, V_t, rngs_B[k]);
            B_private[k].col(t) = b_new;
        }
    }

    // copy private storage into output cube
    Barray.slice(0) = eye_nc;
    for(int t = 1; t < Tt; t++) {
        for(int k = 0; k < n_col; k++) {
            Barray.slice(t).col(k) = B_private[k].col(t);
        }
    }

    return List::create(
        Named("Aarray") = Aarray,
        Named("Barray") = Barray
    );
}

// sample process and observation variances
//
// when exclude_diagonal == true (used for symmetric / undirected
// networks where self-ties are structurally undefined), the residual
// sum of squares and the effective sample size both drop the
// diagonal entries i == j (for square networks only). this gives an
// unbiased posterior on sigma^2 / sigma^2_obs by treating the
// diagonal as not observed rather than as observed-to-be-zero.
//
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
                           bool is_gaussian = false,
                           bool exclude_diagonal = false) {
    int nc = n_row * n_col;
    // for unipartite networks, the diagonal has n_row entries per slice;
    // bipartite (n_row != n_col) has no meaningful diagonal so we keep
    // all entries even if exclude_diagonal is set.
    bool drop_diag = exclude_diagonal && (n_row == n_col);
    int n_off_per_slice = drop_diag ? (nc - n_row) : nc;
    // process variance RSS
    double proc_rss = 0.0;

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:proc_rss)
    #endif
    for(int j = 0; j < p; j++) {
        double local_rss = 0.0;

        for(int t = 1; t < Tt; t++) {
            int idx_curr = j * Tt + t;
            int idx_prev = j * Tt + t - 1;

            // current and previous theta slices
            arma::mat Theta_curr = col_as_mat(Theta_4d, idx_curr, n_row, n_col);
            arma::mat Theta_prev = col_as_mat(Theta_4d, idx_prev, n_row, n_col);

            // residual: theta_t - A_t * theta_{t-1} * B_t'
            arma::mat pred = Aarray.slice(t) * Theta_prev * Barray.slice(t).t();
            arma::mat resid = Theta_curr - pred;
            if (drop_diag) {
                local_rss += accu(square(resid)) -
                    accu(square(resid.diag()));
            } else {
                local_rss += accu(square(resid));
            }
        }

        proc_rss += local_rss;
    }

    double sigma2 = (b_sig + proc_rss / 2.0) /
        R::rgamma((a_sig + n_off_per_slice * (Tt - 1) * p) / 2.0, 1.0);

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

                arma::mat resid_obs = Z_jt - (Theta_jt + M_j);
                if (drop_diag) {
                    local_rss += accu(square(resid_obs)) -
                        accu(square(resid_obs.diag()));
                } else {
                    local_rss += accu(square(resid_obs));
                }
            }

            obs_rss += local_rss;
        }

        sigma2_obs = (1.0 + obs_rss / 2.0) /
            R::rgamma((1.0 + n_off_per_slice * Tt * p) / 2.0, 1.0);
    }

    return List::create(
        Named("sigma2") = sigma2,
        Named("sigma2_obs") = sigma2_obs
    );
}

// FFBS for large networks using blocked operations
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

    // block size for relation batching
    const int block_size = std::max(1, std::min(4, p));

    // sequential: FFBS uses arma::randn, not thread-safe
    for(int j_block = 0; j_block < p; j_block += block_size) {
        int j_end = std::min(j_block + block_size, p);

        // run FFBS for each relation in this block
        for(int j = j_block; j < j_end; j++) {
            // workspace for this relation
            arma::cube Z_j(n_row, n_col, Tt);

            // copy slices into contiguous cube
            for(int t = 0; t < Tt; t++) {
                Z_j.slice(t) = col_as_mat(Z_4d, j * Tt + t, n_row, n_col);
            }

            // run FFBS
            arma::cube Theta_j = ffbs_theta_struct_5arg_cpp(Z_j, M.slice(j), Aarray, Barray, sigma2);

            // write back
            for(int t = 0; t < Tt; t++) {
                Theta_4d_out.slice(j * Tt + t) = Theta_j.slice(t);
            }
        }
    }

    return Theta_4d_out;
}

// A/B update for large networks (m > 100)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_AB_batch_large(const arma::mat& Theta_4d,
                          const arma::cube& Aarray_old,
                          const arma::cube& Barray_old,
                          double sigma2, double tauA2, double tauB2,
                          bool ar1, double rhoA, double rhoB,
                          int n_row, int n_col, int p, int Tt,
                          std::string prior_kind = "rw",
                          double kappaA_inv = 0.0, double kappaB_inv = 0.0) {
    bool use_rw = (prior_kind == "rw") && !ar1;
    int nc = n_row * n_col;
    arma::cube Aarray(n_row, n_row, Tt);
    arma::cube Barray(n_col, n_col, Tt);

    // cache frequently used quantities
    arma::mat eye_nr = eye(n_row, n_row);
    arma::mat eye_nc = eye(n_col, n_col);
    double inv_sigma2 = 1.0 / sigma2;
    double inv_tauA2 = 1.0 / tauA2;
    double inv_tauB2 = 1.0 / tauB2;

    // --- update A (sender dynamics) ---
    // private per-row storage avoids false sharing
    set_dbn_threads();
#ifdef _OPENMP
    int n_threads_A = std::min(omp_get_max_threads(), n_row);
#else
    int n_threads_A = 1;
#endif
    auto rngs_A = init_thread_rngs(n_row);

    std::vector<arma::mat> A_private(n_row, arma::mat(n_row, Tt, arma::fill::zeros));

    #pragma omp parallel for schedule(static) num_threads(n_threads_A)
    for(int i = 0; i < n_row; i++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        arma::mat F_local(p * n_col, n_row);
        arma::vec y_local(p * n_col);
        arma::mat V_local(n_row, n_row);

        for(int t = 1; t < Tt; t++) {
            for(int j = 0; j < p; j++) {
                int base_idx = j * n_col;

                arma::mat Theta_prev = col_as_mat(Theta_4d, j * Tt + t - 1, n_row, n_col);
                arma::mat Theta_curr = col_as_mat(Theta_4d, j * Tt + t, n_row, n_col);

                F_local.rows(base_idx, base_idx + n_col - 1) =
                    (Theta_prev * Barray_old.slice(t).t()).t();
                y_local.subvec(base_idx, base_idx + n_col - 1) = Theta_curr.row(i).t();
            }

            // direct symmetric-PD inverse of the precision matrix.
            double prior_prec_mult = (use_rw && t < Tt - 1) ? 2.0 : 1.0;
            arma::mat V_inv = inv_sigma2 * (F_local.t() * F_local) + (prior_prec_mult * inv_tauA2 + kappaA_inv) * eye_nr;
            V_local = safe_inv_sympd(V_inv);

            V_local = 0.5 * (V_local + V_local.t());
            arma::vec m_post = V_local * (inv_sigma2 * (F_local.t() * y_local));

            if (use_rw) {
                arma::vec neighbor_prev(n_row, arma::fill::zeros);
                if (t == 1) {
                    neighbor_prev(i) = 1.0;
                } else {
                    neighbor_prev = A_private[i].col(t - 1);
                }
                if (t < Tt - 1) {
                    arma::vec neighbor_next = Aarray_old.slice(t + 1).row(i).t();
                    m_post += V_local * (inv_tauA2 * (neighbor_prev + neighbor_next));
                } else {
                    m_post += V_local * (inv_tauA2 * neighbor_prev);
                }
            } else if (ar1 && t > 1) {
                m_post += (rhoA * inv_tauA2) * (V_local * A_private[i].col(t - 1));
            }

            arma::vec a_new = thread_safe_mvnrnd(m_post, V_local, rngs_A[i]);
            A_private[i].col(t) = a_new;
        }
    }

    Aarray.slice(0) = eye_nr;
    for(int t = 1; t < Tt; t++) {
        for(int i = 0; i < n_row; i++) {
            Aarray.slice(t).row(i) = A_private[i].col(t).t();
        }
    }

    // --- update B (receiver dynamics) ---
#ifdef _OPENMP
    int n_threads_B = std::min(omp_get_max_threads(), n_col);
#else
    int n_threads_B = 1;
#endif
    auto rngs_B = init_thread_rngs(n_col);

    std::vector<arma::mat> B_private(n_col, arma::mat(n_col, Tt, arma::fill::zeros));

    #pragma omp parallel for schedule(static) num_threads(n_threads_B)
    for(int k = 0; k < n_col; k++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        arma::mat F_local(p * n_row, n_col);
        arma::vec y_local(p * n_row);
        arma::mat V_local(n_col, n_col);

        for(int t = 1; t < Tt; t++) {
            for(int j = 0; j < p; j++) {
                int base_idx = j * n_row;

                arma::mat Theta_prev = col_as_mat(Theta_4d, j * Tt + t - 1, n_row, n_col);
                arma::mat Theta_curr = col_as_mat(Theta_4d, j * Tt + t, n_row, n_col);

                F_local.rows(base_idx, base_idx + n_row - 1) = Aarray.slice(t) * Theta_prev;
                y_local.subvec(base_idx, base_idx + n_row - 1) = Theta_curr.col(k);
            }

            // direct symmetric-PD inverse (see comment in A update above).
            double prior_prec_mult_B = (use_rw && t < Tt - 1) ? 2.0 : 1.0;
            arma::mat V_inv = inv_sigma2 * (F_local.t() * F_local) + (prior_prec_mult_B * inv_tauB2 + kappaB_inv) * eye_nc;
            V_local = safe_inv_sympd(V_inv);

            V_local = 0.5 * (V_local + V_local.t());
            arma::vec m_post = V_local * (inv_sigma2 * (F_local.t() * y_local));

            if (use_rw) {
                arma::vec neighbor_prev(n_col, arma::fill::zeros);
                if (t == 1) {
                    neighbor_prev(k) = 1.0;
                } else {
                    neighbor_prev = B_private[k].col(t - 1);
                }
                if (t < Tt - 1) {
                    arma::vec neighbor_next = Barray_old.slice(t + 1).col(k);
                    m_post += V_local * (inv_tauB2 * (neighbor_prev + neighbor_next));
                } else {
                    m_post += V_local * (inv_tauB2 * neighbor_prev);
                }
            } else if (ar1 && t > 1) {
                m_post += (rhoB * inv_tauB2) * (V_local * B_private[k].col(t - 1));
            }

            arma::vec b_new = thread_safe_mvnrnd(m_post, V_local, rngs_B[k]);
            B_private[k].col(t) = b_new;
        }
    }

    Barray.slice(0) = eye_nc;
    for(int t = 1; t < Tt; t++) {
        for(int k = 0; k < n_col; k++) {
            Barray.slice(t).col(k) = B_private[k].col(t);
        }
    }

    return List::create(
        Named("Aarray") = Aarray,
        Named("Barray") = Barray
    );
}

// process variance RSS for large networks (blocked over time)
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

        // loop over time in blocks
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

// gaussian observation RSS for the dynamic model
// takes the flat matrix layout (nc x p*Tt) used by the MCMC loop
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_gaussian_obs_residuals_dynamic_cpp(const arma::mat& Z_4d,
                                                   const arma::mat& Theta_4d,
                                                   const arma::cube& M,
                                                   int n_row, int n_col, int p, int Tt) {
    set_dbn_threads();
    int n_total = p * Tt;
    double obs_rss = 0.0;

    #pragma omp parallel for reduction(+:obs_rss) schedule(static)
    for (int idx = 0; idx < n_total; idx++) {
        int j = idx / Tt;
        int t = idx % Tt;
        int col_idx = j * Tt + t;

        arma::mat Z_jt = col_as_mat(Z_4d, col_idx, n_row, n_col);
        arma::mat Theta_jt = col_as_mat(Theta_4d, col_idx, n_row, n_col);
        arma::mat diff = Z_jt - (Theta_jt + M.slice(j));
        obs_rss += arma::accu(diff % diff);
    }

    return obs_rss;
}

// AR(1) innovation sum of squares for tau update
// innov[t] = A[t] - rho * A[t-1] - (1 - rho) * I
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_ar1_innovation_ss_cpp(const arma::cube& ABarray, double rho, int n, int Tt) {
    set_dbn_threads();
    arma::mat I = arma::eye(n, n);
    double ss = 0.0;

    #pragma omp parallel for reduction(+:ss) schedule(static)
    for (int t = 1; t < Tt; t++) {
        arma::mat innov = ABarray.slice(t) - rho * ABarray.slice(t - 1) - (1.0 - rho) * I;
        ss += arma::accu(innov % innov);
    }

    return ss;
}

// numerator and denominator for AR(1) rho full conditional
// diffA_t   = A[t]   - I
// diffA_tm1 = A[t-1] - I
// num   = sum(diffA_t * diffA_tm1)
// denom = sum(diffA_tm1^2)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List compute_rho_update_cpp(const arma::cube& ABarray, int n, int Tt) {
    set_dbn_threads();
    arma::mat I = arma::eye(n, n);
    double num = 0.0;
    double denom = 0.0;

    #pragma omp parallel for reduction(+:num,denom) schedule(static)
    for (int t = 1; t < Tt; t++) {
        arma::mat diff_t   = ABarray.slice(t) - I;
        arma::mat diff_tm1 = ABarray.slice(t - 1) - I;
        num   += arma::accu(diff_t % diff_tm1);
        denom += arma::accu(diff_tm1 % diff_tm1);
    }

    return List::create(
        Named("num") = num,
        Named("denom") = denom
    );
}
