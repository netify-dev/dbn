#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// Fast log-sum-exp trick for numerical stability
inline double log_sum_exp(const arma::vec& log_vals) {
    double max_val = log_vals.max();
    if(!std::isfinite(max_val)) {
        return -datum::inf;
    }
    return max_val + log(sum(exp(log_vals - max_val)));
}

// Forward algorithm for HMM
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat forward_hmm(const arma::cube& Theta_avg,
                           const List& A_list,
                           const List& B_list,
                           const arma::mat& Pi,
                           double sigma2,
                           const arma::vec& pi0) {
    int Tt = Theta_avg.n_slices;
    int R = Pi.n_rows;
    int m = Theta_avg.n_rows;
    
    // precompute constants
    double log_norm_const = -0.5 * m * m * log(2 * M_PI * sigma2);
    double inv_2sigma2 = 0.5 / sigma2;
    
    // precompuite log transition matrix
    arma::mat log_Pi(R, R);
    for(int i = 0; i < R; i++) {
        for(int j = 0; j < R; j++) {
            log_Pi(i, j) = log(Pi(i, j) + 1e-300);
        }
    }
    
    // log forward probs
    arma::mat log_alpha(R, Tt);
    
    // prealloc workspace for parallel computation
    #ifdef _OPENMP
    int n_threads = omp_get_max_threads();
    std::vector<arma::mat> workspace(n_threads, arma::mat(m, m));
    #else
    arma::mat workspace(m, m);
    #endif
    
    // t = 1: initial state
    arma::vec log_p1(R);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int r = 0; r < R; r++) {
        #ifdef _OPENMP
        int tid = omp_get_thread_num();
        arma::mat& temp = workspace[tid];
        #else
        arma::mat& temp = workspace;
        #endif
        
        arma::mat A_r = as<arma::mat>(A_list[r]);
        arma::mat B_r = as<arma::mat>(B_list[r]);
        
        // compute mu_1 under regime r using temporary workspace
        temp = A_r * Theta_avg.slice(0);
        arma::mat mu_1 = temp * B_r.t();
        
        // log density
        mu_1 -= Theta_avg.slice(0);  // In-place subtraction
        double quad_form = arma::dot(arma::vectorise(mu_1), arma::vectorise(mu_1));
        
        log_p1(r) = log(pi0(r)) + log_norm_const - inv_2sigma2 * quad_form;
    }
    
    // normalize to prevent underflow
    double log_c1 = log_sum_exp(log_p1);
    log_alpha.col(0) = log_p1 - log_c1;
    
    // fwd recursion with improved cache usage
    for(int t = 1; t < Tt; t++) {
        arma::vec log_pt(R);
        
        // extract current and previous time slices for better cache usage
        arma::mat Theta_curr = Theta_avg.slice(t);
        arma::mat Theta_prev = Theta_avg.slice(t-1);
        
        // compute emission probs
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int r = 0; r < R; r++) {
            #ifdef _OPENMP
            int tid = omp_get_thread_num();
            arma::mat& temp = workspace[tid];
            #else
            arma::mat& temp = workspace;
            #endif
            
            arma::mat A_r = as<arma::mat>(A_list[r]);
            arma::mat B_r = as<arma::mat>(B_list[r]);
            
            // prediction under regime r with workspace reuse
            temp = A_r * Theta_prev;
            arma::mat pred = temp * B_r.t();
            pred -= Theta_curr;  // In-place residual computation
            
            double quad_form = arma::dot(arma::vectorise(pred), arma::vectorise(pred));
            log_pt(r) = -inv_2sigma2 * quad_form;
        }
        
        // transition step - more robust implementation
        arma::vec log_alpha_prev = log_alpha.col(t-1);
        arma::vec log_alpha_new(R);
        
        // compute transition probs
        for(int j = 0; j < R; j++) {
            arma::vec log_trans(R);
            for(int i = 0; i < R; i++) {
                log_trans(i) = log_alpha_prev(i) + log_Pi(i, j);
            }
            log_alpha_new(j) = log_sum_exp(log_trans) + log_pt(j);
        }
        
        // normlaize
        double log_ct = log_sum_exp(log_alpha_new);
        log_alpha.col(t) = log_alpha_new - log_ct;
    }
    
    return log_alpha;
}

// Forward algorithm with beam search for very long sequences
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat forward_hmm_fast(const arma::cube& Theta_avg,
                          const List& A_list,
                          const List& B_list,
                          const arma::mat& Pi,
                          double sigma2,
                          const arma::vec& pi0,
                          int beam_width = 0) {
    int Tt = Theta_avg.n_slices;
    int R = Pi.n_rows;
    int m = Theta_avg.n_rows;
    
    // use full forward if beam width not specified
    if(beam_width <= 0 || beam_width >= R) {
        return forward_hmm(Theta_avg, A_list, B_list, Pi, sigma2, pi0);
    }
    
    // log fwd probs - use sparse representation for beam
    arma::mat log_alpha(R, Tt);
    log_alpha.fill(-datum::inf);
    
    // precompute constants
    double log_norm_const = -0.5 * m * m * log(2 * M_PI * sigma2);
    double inv_sigma2 = 1.0 / sigma2;
    
    // cache A and B matrices for faster access
    std::vector<arma::mat> A_cache(R), B_cache(R);
    for(int r = 0; r < R; r++) {
        A_cache[r] = as<arma::mat>(A_list[r]);
        B_cache[r] = as<arma::mat>(B_list[r]);
    }
    
    // t = 1: initial state
    arma::vec log_p1(R);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int r = 0; r < R; r++) {
        arma::mat pred = A_cache[r] * Theta_avg.slice(0) * B_cache[r].t();
        arma::mat resid = Theta_avg.slice(0) - pred;
        double quad_form = accu(resid % resid);
        log_p1(r) = log(pi0(r)) + log_norm_const - 0.5 * quad_form * inv_sigma2;
    }
    
    // normalize
    double log_c1 = log_sum_exp(log_p1);
    log_alpha.col(0) = log_p1 - log_c1;
    
    // identify top beam_width states for first time
    arma::uvec active_states = arma::sort_index(log_alpha.col(0), "descend");
    active_states = active_states.head(beam_width);
    
    // fwd recursion with beam search
    for(int t = 1; t < Tt; t++) {
        arma::vec log_alpha_new(R);
        log_alpha_new.fill(-datum::inf);
        
        // only consider transitions from active states
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int j = 0; j < R; j++) {
            double max_trans = -datum::inf;
            for(size_t idx = 0; idx < active_states.n_elem; idx++) {
                int i = active_states(idx);
                double trans_prob = log_alpha(i, t-1) + log(Pi(i, j) + 1e-300);
                if(trans_prob > max_trans) {
                    max_trans = trans_prob;
                }
            }
            
            if(max_trans > -datum::inf) {
                // compute emission probability
                arma::mat pred = A_cache[j] * Theta_avg.slice(t-1) * B_cache[j].t();
                arma::mat resid = Theta_avg.slice(t) - pred;
                double quad_form = accu(resid % resid);
                log_alpha_new(j) = max_trans + log_norm_const - 0.5 * quad_form * inv_sigma2;
            }
        }
        
        // normalize and update active states
        double log_ct = log_sum_exp(log_alpha_new);
        log_alpha.col(t) = log_alpha_new - log_ct;
        
        // update beam
        active_states = arma::sort_index(log_alpha.col(t), "descend");
        active_states = active_states.head(beam_width);
    }
    
    return log_alpha;
}

// Backward sampling for HMM states
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerVector backward_sample(const arma::mat& log_alpha,
                                  const arma::mat& Pi) {
    int Tt = log_alpha.n_cols;
    int R = log_alpha.n_rows;
    
    IntegerVector S(Tt);
    
    // sample final state
    arma::vec probs = exp(log_alpha.col(Tt-1));
    probs = probs / sum(probs); // Ensure normalized
    
    double u = R::runif(0, 1);
    double cumsum = 0;
    for(int r = 0; r < R; r++) {
        cumsum += probs(r);
        if(u <= cumsum) {
            S[Tt-1] = r + 1; // R uses 1-based indexing
            break;
        }
    }
    
    // backward sampling
    for(int t = Tt-2; t >= 0; t--) {
        int s_next = S[t+1] - 1; // convert to 0-based
        
        // P(S_t = j | S_{t+1} = s_next, Y_{1:t})
        arma::vec log_probs(R);
        for(int j = 0; j < R; j++) {
            log_probs(j) = log_alpha(j, t) + log(Pi(j, s_next));
        }
        
        // normalize
        double log_norm = log_sum_exp(log_probs);
        arma::vec probs_t = exp(log_probs - log_norm);
        
        // sample
        u = R::runif(0, 1);
        cumsum = 0;
        for(int r = 0; r < R; r++) {
            cumsum += probs_t(r);
            if(u <= cumsum) {
                S[t] = r + 1;
                break;
            }
        }
    }
    
    return S;
}

// Fast backward sampling for beam search
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerVector backward_sample_fast(const arma::mat& log_alpha,
                                  const arma::mat& Pi) {
    int Tt = log_alpha.n_cols;
    int R = log_alpha.n_rows;
    
    IntegerVector S(Tt);
    
    // use Gumbel-max trick for faster sampling
    // sample final state
    arma::vec log_probs_T = log_alpha.col(Tt-1);
    arma::vec gumbel_noise(R);
    for(int r = 0; r < R; r++) {
        double u = R::runif(0, 1);
        gumbel_noise(r) = -log(-log(u));
    }
    arma::vec perturbed = log_probs_T + gumbel_noise;
    S[Tt-1] = perturbed.index_max() + 1;
    
    // backward sampling with caching
    arma::vec log_Pi_col(R);
    
    for(int t = Tt-2; t >= 0; t--) {
        int s_next = S[t+1] - 1;
        
        // precompute transition column
        for(int i = 0; i < R; i++) {
            log_Pi_col(i) = log(Pi(i, s_next) + 1e-300);
        }
        
        // compute posterior
        arma::vec log_probs = log_alpha.col(t) + log_Pi_col;
        
        // Gumbel-max sampling
        for(int r = 0; r < R; r++) {
            double u = R::runif(0, 1);
            gumbel_noise(r) = -log(-log(u));
        }
        perturbed = log_probs + gumbel_noise;
        S[t] = perturbed.index_max() + 1;
    }
    
    return S;
}

// Build time-varying arrays from regime assignments
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List build_regime_arrays(const IntegerVector& S,
                        const List& A_list,
                        const List& B_list,
                        int m, int Tt) {
    arma::cube Aarray(m, m, Tt);
    arma::cube Barray(m, m, Tt);
    
    // preextract matrices from lists to avoid repeated conversions
    int R = A_list.size();
    std::vector<arma::mat> A_vec(R), B_vec(R);
    for(int r = 0; r < R; r++) {
        A_vec[r] = as<arma::mat>(A_list[r]);
        B_vec[r] = as<arma::mat>(B_list[r]);
    }
    
    // process in blocks for better cache usage
    const int block_size = 32;
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int t_block = 0; t_block < Tt; t_block += block_size) {
        int t_end = std::min(t_block + block_size, Tt);
        
        for(int t = t_block; t < t_end; t++) {
            int regime = S[t] - 1; // Convert to 0-based
            
            // direct memory copy for efficiency
            std::memcpy(Aarray.slice_memptr(t), A_vec[regime].memptr(), m * m * sizeof(double));
            std::memcpy(Barray.slice_memptr(t), B_vec[regime].memptr(), m * m * sizeof(double));
        }
    }
    
    return List::create(
        Named("Aarray") = Aarray,
        Named("Barray") = Barray
    );
}

// Collect theta pairs for a specific regime
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List collect_regime_thetas(const arma::cube& Theta_avg,
                          const IntegerVector& S,
                          int regime, int m) {
    // find times where S_t = regime (for t >= 2)
    std::vector<int> idx;
    for(int t = 1; t < S.size(); t++) {
        if(S[t] == regime) {
            idx.push_back(t);
        }
    }
    
    int n_r = idx.size();
    
    if(n_r == 0) {
        return List::create(
            Named("Th_prev") = arma::cube(),
            Named("Th_curr") = arma::cube(),
            Named("n_obs") = 0
        );
    }
    
    arma::cube Th_prev(m, m, n_r);
    arma::cube Th_curr(m, m, n_r);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < n_r; i++) {
        int t = idx[i];
        Th_prev.slice(i) = Theta_avg.slice(t - 1);
        Th_curr.slice(i) = Theta_avg.slice(t);
    }
    
    return List::create(
        Named("Th_prev") = Th_prev,
        Named("Th_curr") = Th_curr,
        Named("n_obs") = n_r
    );
}

// Fast transition counting
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat count_transitions(const IntegerVector& S, int R) {
    arma::mat n_ij(R, R, fill::zeros);
    
    for(int t = 1; t < S.size(); t++) {
        int from = S[t-1] - 1; // Convert to 0-based
        int to = S[t] - 1;
        n_ij(from, to) += 1;
    }
    
    return n_ij;
}

// Compute residuals for all regimes at once
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_regime_residuals(const List& A_list,
                               const arma::mat& I_m,
                               int R, int m) {
    double rss = 0.0;
    
    for(int r = 0; r < R; r++) {
        arma::mat A_r = as<arma::mat>(A_list[r]);
        arma::mat diff = A_r - I_m;
        rss += accu(diff % diff);
    }
    
    return rss;
}

// State sequence initialization using spectral clustering
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerVector init_states_spectral(const arma::cube& Y,
                                  int R, int m, int p, int Tt) {
    IntegerVector S(Tt);
    
    // compute similarity matrix based on network snapshots
    // should incorporate this into netify
    arma::mat sim_mat(Tt, Tt, fill::zeros);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int t1 = 0; t1 < Tt; t1++) {
        for(int t2 = t1; t2 < Tt; t2++) {
            double sim = 0.0;
            int count = 0;
            
            // avg similarity across relations
            for(int rel = 0; rel < p; rel++) {
                // extract slices for this relation
                arma::mat Y1(m, m);
                arma::mat Y2(m, m);
                
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < m; j++) {
                        Y1(i,j) = Y(i + j*m + rel*m*m + t1*m*m*p);
                        Y2(i,j) = Y(i + j*m + rel*m*m + t2*m*m*p);
                    }
                }
                
                // count finite values
                arma::uvec finite1 = find_finite(Y1);
                arma::uvec finite2 = find_finite(Y2);
                arma::uvec common = intersect(finite1, finite2);
                
                if(common.n_elem > 0) {
                    arma::vec v1 = Y1.elem(common);
                    arma::vec v2 = Y2.elem(common);
                    sim += as_scalar(cor(v1, v2));
                    count++;
                }
            }
            
            if(count > 0) {
                sim_mat(t1, t2) = sim_mat(t2, t1) = sim / count;
            }
        }
    }
    
    // simple k-means on similarity features
    arma::mat features = sim_mat;
    arma::mat centroids;
    
    // params: means, data, k, seed_mode, n_iter, print_mode
    bool success = arma::kmeans(centroids, features.t(), R, 
                                random_subset, 10, false);
    
    if(success) {
        // assign each time point to nearest centroid
        for(int t = 0; t < Tt; t++) {
            arma::vec dists(R);
            for(int r = 0; r < R; r++) {
                dists(r) = norm(features.col(t) - centroids.col(r), 2);
            }
            S[t] = dists.index_min() + 1; // Convert to 1-based
        }
    } else {
        // fall back to random initialization
        for(int t = 0; t < Tt; t++) {
            S[t] = (int)(R::runif(0, 1) * R) + 1;
        }
    }
    
    return S;
}