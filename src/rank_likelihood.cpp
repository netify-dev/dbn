#include <RcppArmadillo.h>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

//' Fast rank-likelihood sampler in C++
//' 
//' @description Port of rz_fc() to C++ for ~8-15x speedup on large ordinal networks
//' @param R Observed ranks (vectorized matrix)
//' @param Z Current latent values (vectorized)
//' @param EZ Expected values (vectorized)
//' @param iranks List of indices for each rank level
//' @return Updated Z vector
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec rz_fc_cpp(const arma::vec& R,
                    const arma::vec& Z, 
                    const arma::vec& EZ,
                    const Rcpp::List& iranks) {
  
  const int n = Z.n_elem;
  arma::vec Z_new = Z;  // copy to preserve input
  
  // get number of rank levels (excluding na)
  int n_ranks = iranks.size();
  if (iranks.containsElementNamed("NA")) {
    n_ranks--;
  }
  
  // get max rank value
  int max_rank = 0;
  Rcpp::CharacterVector names = iranks.names();
  for (int k = 0; k < iranks.size(); k++) {
    if (Rcpp::as<std::string>(names[k]) != "NA") {
      int rank_val = std::atoi(Rcpp::as<std::string>(names[k]).c_str());
      if (rank_val > max_rank) max_rank = rank_val;
    }
  }
  
  // process each rank level
  for (int k = 1; k <= max_rank; k++) {
    // check if this rank exists in the data
    std::string rank_key = std::to_string(k);
    bool rank_exists = false;
    int rank_idx = -1;
    
    for (int j = 0; j < iranks.size(); j++) {
      if (Rcpp::as<std::string>(names[j]) == rank_key) {
        rank_exists = true;
        rank_idx = j;
        break;
      }
    }
    
    if (!rank_exists) continue;
    
    // get indices for this rank
    arma::uvec idx = Rcpp::as<arma::uvec>(iranks[rank_idx]) - 1;  // R to C++ indexing
    
    if (idx.n_elem == 0) continue;
    
    // determine bounds
    double lb, ub;
    
    // lower bound
    if (k == 1) {
      lb = -arma::datum::inf;
    } else {
      // find max Z value from rank k-1
      lb = -arma::datum::inf;
      std::string prev_rank_key = std::to_string(k-1);
      for (int j = 0; j < iranks.size(); j++) {
        if (Rcpp::as<std::string>(names[j]) == prev_rank_key) {
          arma::uvec prev_idx = Rcpp::as<arma::uvec>(iranks[j]) - 1;
          if (prev_idx.n_elem > 0) {
            lb = Z_new.elem(prev_idx).max();
          }
          break;
        }
      }
    }
    
    // upper bound
    if (k == max_rank) {
      ub = arma::datum::inf;
    } else {
      // find min Z value from rank k+1
      ub = arma::datum::inf;
      std::string next_rank_key = std::to_string(k+1);
      for (int j = 0; j < iranks.size(); j++) {
        if (Rcpp::as<std::string>(names[j]) == next_rank_key) {
          arma::uvec next_idx = Rcpp::as<arma::uvec>(iranks[j]) - 1;
          if (next_idx.n_elem > 0) {
            ub = Z_new.elem(next_idx).min();
          }
          break;
        }
      }
    }
    
    // sample from truncated normal for this rank
    arma::vec mu_k = EZ.elem(idx);
    arma::vec Z_k(idx.n_elem);
    
    for (arma::uword i = 0; i < idx.n_elem; i++) {
      double mu_i = mu_k(i);
      
      // handle non-finite mu_i
      if (!std::isfinite(mu_i)) {
        mu_i = 0.0;  // use zero as default
      }
      
      // compute cdf bounds
      double p_lo = R::pnorm(lb, mu_i, 1.0, true, false);
      double p_hi = R::pnorm(ub, mu_i, 1.0, true, false);
      
      // hande edge cases 
      const double eps = std::numeric_limits<double>::epsilon();
      p_lo = std::max(p_lo, eps);
      p_hi = std::min(p_hi, 1.0 - eps);
      
      // samp uniform and transform
      double u;
      if (p_lo < p_hi) {
        u = R::runif(p_lo, p_hi);
      } else {
        // bounds are inverted or equal, use the mean
        u = (p_lo + p_hi) / 2.0;
      }
      double z_val = R::qnorm(u, mu_i, 1.0, true, false);
      
      // check for non-finite values and handle edge cases
      if (!std::isfinite(z_val)) {
        if (p_hi - p_lo < eps) {
          // very tight bounds, use midpoint
          if (std::isfinite(lb) && std::isfinite(ub)) {
            z_val = (lb + ub) / 2.0;
          } else if (std::isfinite(lb)) {
            z_val = lb + 1.0;
          } else if (std::isfinite(ub)) {
            z_val = ub - 1.0;
          } else {
            z_val = mu_i;
          }
        } else {
          // fallback to mean
          z_val = mu_i;
        }
      }
      
      Z_k(i) = z_val;
    }
    
    // update Z_new
    Z_new.elem(idx) = Z_k;
  }
  
  // handle NA values
  if (iranks.containsElementNamed("NA")) {
    // find index of "NA" in the list
    int na_idx = -1;
    for (int j = 0; j < iranks.size(); j++) {
      if (Rcpp::as<std::string>(names[j]) == "NA") {
        na_idx = j;
        break;
      }
    }
    
    if (na_idx >= 0) {
      arma::uvec idx_na = Rcpp::as<arma::uvec>(iranks[na_idx]) - 1;  
      if (idx_na.n_elem > 0) {
        arma::vec mu_na = EZ.elem(idx_na);
        arma::vec Z_na(idx_na.n_elem);
        
        for (arma::uword i = 0; i < idx_na.n_elem; i++) {
          // for NA vals, samp from norm distr without truncating
          if (std::isfinite(mu_na(i))) {
            Z_na(i) = R::rnorm(mu_na(i), 1.0);
          } else {
            Z_na(i) = R::rnorm(0.0, 1.0);
          }
        }
        
        Z_new.elem(idx_na) = Z_na;
      }
    }
  }
  
  return Z_new;
}

//' Build rank indices from rank matrix
//' 
//' @description Helper to convert rank matrix to iranks format for C++
//' @param R Rank matrix
//' @return List of indices for each rank
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List build_rank_indices(const arma::mat& R) {
  arma::vec R_vec = arma::vectorise(R);
  
  // find finite values first
  arma::uvec finite_idx = arma::find_finite(R_vec);
  
  // find unique ranks only among finite values
  arma::vec unique_ranks;
  if (finite_idx.n_elem > 0) {
    unique_ranks = arma::unique(R_vec(finite_idx));
  }
  
  Rcpp::List iranks;
  Rcpp::CharacterVector names;
  
  // build index list for each rank
  for (arma::uword k = 0; k < unique_ranks.n_elem; k++) {
    double rank_val = unique_ranks(k);
    arma::uvec idx = arma::find(R_vec == rank_val) + 1;  // C++ to R indexing
    
    iranks.push_back(idx);
    names.push_back(std::to_string(static_cast<int>(rank_val)));
  }
  
  // add NA indices
  arma::uvec na_idx = arma::find_nonfinite(R_vec) + 1;
  if (na_idx.n_elem > 0) {
    iranks.push_back(na_idx);
    names.push_back("NA");
  }
  
  iranks.attr("names") = names;
  return iranks;
}