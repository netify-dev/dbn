#include <RcppArmadillo.h>
#include "dbn_stability.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Validate transition matrix proposals
//' 
//' @param A_prop Proposed actor transition matrix
//' @param B_prop Proposed latent transition matrix
//' @param p Number of actors
//' @param q Latent dimension
//' @param check_stationary Check stationarity condition
//' @return Validity indicator
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
bool validate_transition_proposal(const arma::mat& A_prop,
                                const arma::mat& B_prop,
                                int p, int q,
                                bool check_stationary = true) {
    // check dimensions
    if (A_prop.n_rows != p || A_prop.n_cols != p) {
        return false;
    }
    if (B_prop.n_rows != q || B_prop.n_cols != q) {
        return false;
    }
    
    // check for nan or inf
    if (!A_prop.is_finite() || !B_prop.is_finite()) {
        return false;
    }
    
    // check spectral radius
    arma::cx_vec eigenvalues_A, eigenvalues_B;
    arma::eig_gen(eigenvalues_A, A_prop);
    arma::eig_gen(eigenvalues_B, B_prop);
    
    double max_abs_eig_A = 0.0, max_abs_eig_B = 0.0;
    for (size_t i = 0; i < eigenvalues_A.n_elem; i++) {
        max_abs_eig_A = std::max(max_abs_eig_A, std::abs(eigenvalues_A(i)));
    }
    for (size_t i = 0; i < eigenvalues_B.n_elem; i++) {
        max_abs_eig_B = std::max(max_abs_eig_B, std::abs(eigenvalues_B(i)));
    }
    
    if (max_abs_eig_A >= SPECTRAL_RADIUS_THRESHOLD || 
        max_abs_eig_B >= SPECTRAL_RADIUS_THRESHOLD) {
        return false;
    }
    
    // check stationarity of combined system
    if (check_stationary && !is_stationary(A_prop, B_prop, p, q)) {
        return false;
    }
    
    return true;
}

//' Adaptive random walk proposal for transition matrices
//' 
//' @param A_current Current actor transition matrix
//' @param B_current Current latent transition matrix
//' @param step_size_A Step size for A
//' @param step_size_B Step size for B
//' @param max_attempts Maximum proposal attempts
//' @return List with proposed matrices and validity indicator
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List propose_transition_adaptive(const arma::mat& A_current,
                                     const arma::mat& B_current,
                                     double step_size_A = 0.1,
                                     double step_size_B = 0.1,
                                     int max_attempts = 10) {
    int p = A_current.n_rows;
    int q = B_current.n_rows;
    
    for (int attempt = 0; attempt < max_attempts; attempt++) {
        // random walk proposals
        arma::mat eps_A = arma::randn(p, p) * step_size_A;
        arma::mat eps_B = arma::randn(q, q) * step_size_B;
        
        arma::mat A_prop = A_current + eps_A;
        arma::mat B_prop = B_current + eps_B;
        
        // stabilize if needed ... or try to 
        A_prop = stabilize_spectral_radius(A_prop);
        B_prop = stabilize_spectral_radius(B_prop);
        
        // validate proposal
        if (validate_transition_proposal(A_prop, B_prop, p, q)) {
            return Rcpp::List::create(
                Rcpp::Named("A") = A_prop,
                Rcpp::Named("B") = B_prop,
                Rcpp::Named("valid") = true,
                Rcpp::Named("attempts") = attempt + 1
            );
        }
        
        // reduce step size for next attempt
        step_size_A *= 0.8;
        step_size_B *= 0.8;
    }
    
    // failed to find valid proposal ... return current values
    return Rcpp::List::create(
        Rcpp::Named("A") = A_current,
        Rcpp::Named("B") = B_current,
        Rcpp::Named("valid") = false,
        Rcpp::Named("attempts") = max_attempts
    );
}

//' Check for edge cases in the data
//' 
//' @param Y Data cube (p x p x T)
//' @return List with edge case indicators
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List check_edge_cases(const arma::cube& Y) {
    int p = Y.n_rows;
    int T = Y.n_slices;
    
    bool has_missing = false;
    bool has_infinite = false;
    bool is_constant = true;
    double first_val = Y(0, 0, 0);
    
    int n_missing = 0;
    int n_infinite = 0;
    
    for (int t = 0; t < T; t++) {
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                double val = Y(i, j, t);
                
                if (std::isnan(val)) {
                    has_missing = true;
                    n_missing++;
                } else if (std::isinf(val)) {
                    has_infinite = true;
                    n_infinite++;
                } else if (std::abs(val - first_val) > 1e-10) {
                    is_constant = false;
                }
            }
        }
    }
    
    // check for zero variance
    arma::vec Y_vec = arma::vectorise(Y);
    // Remove NaN/Inf
    Y_vec = Y_vec(arma::find_finite(Y_vec)); 
    double variance = 0.0;
    if (Y_vec.n_elem > 1) {
        variance = arma::var(Y_vec);
    }
    
    return Rcpp::List::create(
        Rcpp::Named("has_missing") = has_missing,
        Rcpp::Named("has_infinite") = has_infinite,
        Rcpp::Named("is_constant") = is_constant,
        Rcpp::Named("n_missing") = n_missing,
        Rcpp::Named("n_infinite") = n_infinite,
        Rcpp::Named("variance") = variance,
        Rcpp::Named("n_obs") = p * p * T
    );
}