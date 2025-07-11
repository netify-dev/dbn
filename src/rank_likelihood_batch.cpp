#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "thread_control.h"

using namespace Rcpp;
using namespace arma;

// Fast truncated normal sampler using accept-reject
inline double rtruncnorm(double mean, double sd, double lower, double upper) {
    // stdz bounds
    double a = (lower - mean) / sd;
    double b = (upper - mean) / sd;
    
    // use different algorithms based on bounds
    if (a > 0.0) {
        // both bounds above mean - use exponential proposal
        double alpha = (a + sqrt(a*a + 4.0)) / 2.0;
        double z, u, rho;
        do {
            z = -log(1.0 - R::runif(0.0, 1.0)) / alpha + a;
            rho = exp(-0.5 * pow(z - alpha, 2.0));
            u = R::runif(0.0, 1.0);
        } while (u > rho || z > b);
        return mean + sd * z;
    } else if (b < 0.0) {
        // both bounds below mean - flip and use exponential
        return -rtruncnorm(-mean, sd, -upper, -lower);
    } else {
        // bounds straddle mean - use normal rejection
        double z;
        do {
            z = R::rnorm(0.0, 1.0);
        } while (z < a || z > b);
        return mean + sd * z;
    }
}

// Rank likelihood sampler
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat rz_fc_matrix(const arma::mat& R, const arma::mat& Z_current,
                    const arma::mat& EZ, const List& iranks) {
    int n = R.n_elem;
    arma::mat Z_new = Z_current;
    
    // flatten matrices for easier indexing
    arma::vec R_vec = vectorise(R);
    arma::vec Z_vec = vectorise(Z_current);
    arma::vec EZ_vec = vectorise(EZ);
    
    // get unique ranks
    int n_ranks = iranks.size();
    
    // process each rank level
    CharacterVector rank_names = iranks.names();
    for(int r = 0; r < n_ranks; r++) {
        String rank_name = rank_names[r];
        
        if(rank_name == "NA") {
            // handle missing values - sample from prior
            IntegerVector idx = iranks[r];
            for(int i = 0; i < idx.size(); i++) {
                int pos = idx[i] - 1; // R uses 1-based indexing
                Z_vec(pos) = R::rnorm(EZ_vec(pos), 1.0);
            }
        } else {
            // get indices for this rank
            IntegerVector idx = iranks[r];
            int rank_val = std::atoi(rank_name.get_cstring());
            
            // find bounds for this rank
            double lower_bound = -datum::inf;
            double upper_bound = datum::inf;
            
            // check for rank below
            if(r > 0) {
                String prev_rank = rank_names[r-1];
                if(prev_rank != "NA") {
                    IntegerVector prev_idx = iranks[r-1];
                    // find maximum z value from lower rank
                    for(int i = 0; i < prev_idx.size(); i++) {
                        lower_bound = std::max(lower_bound, Z_vec(prev_idx[i] - 1));
                    }
                }
            }
            
            // check for rank above
            if(r < n_ranks - 1) {
                String next_rank = rank_names[r+1];
                if(next_rank != "NA") {
                    IntegerVector next_idx = iranks[r+1];
                    // find minimum z value from higher rank
                    for(int i = 0; i < next_idx.size(); i++) {
                        upper_bound = std::min(upper_bound, Z_vec(next_idx[i] - 1));
                    }
                }
            }
            
            // samp truncnorm for all indices at this rank
            for(int i = 0; i < idx.size(); i++) {
                int pos = idx[i] - 1;
                double mean = EZ_vec(pos);
                
                // ensure bounds are feasible
                if(lower_bound >= upper_bound) {
                    // bounds crossed - use small interval around mean
                    lower_bound = mean - 0.1;
                    upper_bound = mean + 0.1;
                }
                
                Z_vec(pos) = rtruncnorm(mean, 1.0, lower_bound, upper_bound);
            }
        }
    }
    
    // reshape back to matrix
    Z_new = reshape(Z_vec, R.n_rows, R.n_cols);
    return Z_new;
}

// Batch rank likelihood update for multiple relations
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube rz_fc_batch(const arma::cube& R, const arma::cube& Z_current,
                      const arma::cube& EZ, const List& IR_list,
                      int m, int p, int Tt) {
    set_dbn_threads(); // set threads from R options
    arma::cube Z_new = Z_current;
    
    // process each relation specified by user
    for(int j = 0; j < p; j++) {
        List IR_j = IR_list[j];
        
        // collect all slices for this relation
        arma::cube R_j(m, m, Tt);
        arma::cube Z_j(m, m, Tt);
        arma::cube EZ_j(m, m, Tt);
        
        for(int t = 0; t < Tt; t++) {
            int slice_idx = j * Tt + t;
            R_j.slice(t) = R.slice(slice_idx);
            Z_j.slice(t) = Z_current.slice(slice_idx);
            EZ_j.slice(t) = EZ.slice(slice_idx);
        }
        
        // flatten for processing
        arma::vec R_j_vec = vectorise(R_j);
        arma::vec Z_j_vec = vectorise(Z_j);
        arma::vec EZ_j_vec = vectorise(EZ_j);
        
        // process all ranks for this relation
        CharacterVector rank_names = IR_j.names();
        int n_ranks = IR_j.size();
        
        for(int r = 0; r < n_ranks; r++) {
            String rank_name = rank_names[r];
            
            if(rank_name == "NA") {
                // handle missing values
                IntegerVector idx = IR_j[r];
                for(int i = 0; i < idx.size(); i++) {
                    int pos = idx[i] - 1;
                    if(pos >= 0 && pos < Z_j_vec.n_elem) {
                        Z_j_vec(pos) = R::rnorm(EZ_j_vec(pos), 1.0);
                    }
                }
            } else {
                // get indices for this rank
                IntegerVector idx = IR_j[r];
                
                // find bounds
                double lower_bound = -datum::inf;
                double upper_bound = datum::inf;
                
                // check for rank below
                if(r > 0) {
                    String prev_rank = rank_names[r-1];
                    if(prev_rank != "NA") {
                        IntegerVector prev_idx = IR_j[r-1];
                        for(int i = 0; i < prev_idx.size(); i++) {
                            int pos = prev_idx[i] - 1;
                            if(pos >= 0 && pos < Z_j_vec.n_elem) {
                                lower_bound = std::max(lower_bound, Z_j_vec(pos));
                            }
                        }
                    }
                }
                
                // check for rank above
                if(r < n_ranks - 1) {
                    String next_rank = rank_names[r+1];
                    if(next_rank != "NA") {
                        IntegerVector next_idx = IR_j[r+1];
                        for(int i = 0; i < next_idx.size(); i++) {
                            int pos = next_idx[i] - 1;
                            if(pos >= 0 && pos < Z_j_vec.n_elem) {
                                upper_bound = std::min(upper_bound, Z_j_vec(pos));
                            }
                        }
                    }
                }
                
                // sample for all indices at this rank
                for(int i = 0; i < idx.size(); i++) {
                    int pos = idx[i] - 1;
                    if(pos >= 0 && pos < Z_j_vec.n_elem) {
                        double mean = EZ_j_vec(pos);
                        
                        // ensure bounds are feasible
                        if(lower_bound >= upper_bound) {
                            lower_bound = mean - 0.1;
                            upper_bound = mean + 0.1;
                        }
                        
                        Z_j_vec(pos) = rtruncnorm(mean, 1.0, lower_bound, upper_bound);
                    }
                }
            }
        }
        
        // reshape back and copy to output
        for(int t = 0; t < Tt; t++) {
            int slice_idx = j * Tt + t;
            // extract the portion of Z_j_vec for this time slice
            arma::vec Z_t_vec = Z_j_vec.subvec(t * m * m, (t + 1) * m * m - 1);
            arma::mat Z_t_mat = reshape(Z_t_vec, m, m);
            Z_new.slice(slice_idx) = Z_t_mat;
        }
    }
    
    return Z_new;
}

// Pre-compute rank structure for efficiency
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List precompute_rank_structure(const arma::cube& R, int m, int p, int Tt) {
    List IR_all(p);
    
    for(int j = 0; j < p; j++) {
        // collect all data for this relation across all time points
        arma::vec R_j_vec;
        R_j_vec.set_size(m * m * Tt);
        
        // extract data for relation j across all time points
        for(int t = 0; t < Tt; t++) {
            int slice_idx = j * Tt + t;
            arma::vec slice_vec = vectorise(R.slice(slice_idx));
            R_j_vec.subvec(t * m * m, (t + 1) * m * m - 1) = slice_vec;
        }
        
        // first find finite values to avoid NaN issues with unique()
        arma::uvec finite_idx = find_finite(R_j_vec);
        arma::vec unique_ranks;
        
        if(finite_idx.n_elem > 0) {
            // get unique ranks only from finite values
            unique_ranks = unique(R_j_vec.elem(finite_idx));
            unique_ranks = sort(unique_ranks);
        }
        
        // build index structure
        List IR_j;
        CharacterVector rank_names;
        
        // add NA category if there are missing values
        arma::uvec na_idx = find_nonfinite(R_j_vec);
        if(na_idx.n_elem > 0) {
            IntegerVector na_vec(na_idx.n_elem);
            for(size_t i = 0; i < na_idx.n_elem; i++) {
                na_vec[i] = na_idx[i] + 1; // Convert to 1-based
            }
            IR_j.push_back(na_vec);
            rank_names.push_back("NA");
        }
        
        // add each rank level ... might need to be rewritten for efficiency later
        if(finite_idx.n_elem > 0) {
            for(size_t r = 0; r < unique_ranks.n_elem; r++) {
                arma::uvec rank_idx = find(R_j_vec == unique_ranks[r]);
                IntegerVector rank_vec(rank_idx.n_elem);
                for(size_t i = 0; i < rank_idx.n_elem; i++) {
                    rank_vec[i] = rank_idx[i] + 1; // Convert to 1-based
                }
                IR_j.push_back(rank_vec);
                rank_names.push_back(std::to_string(static_cast<int>(unique_ranks[r])));
            }
        }
        
        IR_j.names() = rank_names;
        IR_all[j] = IR_j;
    }
    
    return IR_all;
}