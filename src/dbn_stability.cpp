#include "dbn_stability.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Stabilize spectral radius of a matrix
//' 
//' @param M Matrix to stabilize
//' @param threshold Maximum allowed spectral radius
//' @return Stabilized matrix
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat stabilize_spectral_radius(const arma::mat& M, double threshold) {
    arma::cx_vec eigenvalues;
    arma::cx_mat eigenvectors;
    
    arma::eig_gen(eigenvalues, eigenvectors, M);
    
    double max_abs_eigenvalue = 0.0;
    for (size_t i = 0; i < eigenvalues.n_elem; i++) {
        double abs_val = std::abs(eigenvalues(i));
        if (abs_val > max_abs_eigenvalue) {
            max_abs_eigenvalue = abs_val;
        }
    }
    
    if (max_abs_eigenvalue > threshold) {
        // rescale to enforce spectral radius bound
        return M * (threshold / max_abs_eigenvalue);
    }
    
    return M;
}

//' Safe Cholesky decomposition with regularization
//' 
//' @param L Output lower triangular matrix
//' @param A Input matrix
//' @param reg Regularization parameter
//' @return Success indicator
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
bool safe_cholesky(arma::mat& L, const arma::mat& A, double reg) {
    for (int attempt = 0; attempt < MAX_CHOLESKY_ATTEMPTS; attempt++) {
        // fresh copy each attempt to avoid accumulating regularization
        arma::mat A_reg = 0.5 * (A + A.t());  // enforce symmetry

        if (attempt > 0) {
            // add regularization to diagonal elements
            double reg_amount = reg * std::pow(10, attempt - 1);
            A_reg.diag() += reg_amount;
        }

        bool success = arma::chol(L, A_reg, "lower");
        if (success) {
            return true;
        }
    }
    
    // last resort: force positive definiteness
    arma::mat A_final = ensure_positive_definite(A, reg);
    A_final = 0.5 * (A_final + A_final.t());  // enforce symmetry
    bool ok = arma::chol(L, A_final, "lower");
    if (!ok) {
        // absolute last resort: identity
        L = arma::eye(arma::size(A));
    }
    return ok;
}

//' Ensure a matrix is positive definite
//' 
//' @param M Input matrix
//' @param min_eigenvalue Minimum eigenvalue threshold
//' @return Positive definite matrix
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat ensure_positive_definite(const arma::mat& M, double min_eigenvalue) {
    // enforce symmetry
    arma::mat M_sym = 0.5 * (M + M.t());
    
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    
    arma::eig_sym(eigenvalues, eigenvectors, M_sym);
    
    // clamp eigenvalues to minimum threshold
    for (size_t i = 0; i < eigenvalues.n_elem; i++) {
        if (eigenvalues(i) < min_eigenvalue) {
            eigenvalues(i) = min_eigenvalue;
        }
    }
    
    // reconstruct from corrected eigenvalues
    return eigenvectors * arma::diagmat(eigenvalues) * eigenvectors.t();
}

//' Check if VAR process is stationary
//' 
//' @param A Autoregressive coefficient matrix
//' @param B Bilinear coefficient matrix  
//' @param p Number of actors
//' @param q Latent dimension
//' @return Stationarity indicator
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
bool is_stationary(const arma::mat& A, const arma::mat& B, int p, int q) {
    // p, q unused but kept in signature for API consistency
    (void)p; (void)q;
    
    // bilinear stationarity check: rho(A)*rho(B) < 1
    // avoids forming the full Kronecker product
    
    // spectral radius of A
    arma::cx_vec eigenvalues_A;
    arma::eig_gen(eigenvalues_A, A);
    double rho_A = 0.0;
    for (size_t i = 0; i < eigenvalues_A.n_elem; i++) {
        double abs_val = std::abs(eigenvalues_A(i));
        if (abs_val > rho_A) rho_A = abs_val;
    }
    
    // spectral radius of B
    arma::cx_vec eigenvalues_B;
    arma::eig_gen(eigenvalues_B, B);
    double rho_B = 0.0;
    for (size_t i = 0; i < eigenvalues_B.n_elem; i++) {
        double abs_val = std::abs(eigenvalues_B(i));
        if (abs_val > rho_B) rho_B = abs_val;
    }
    
    // stationary if rho(A)*rho(B) < 1
    return (rho_A * rho_B) < SPECTRAL_RADIUS_THRESHOLD;
}

// cached Kronecker product computation
arma::mat KroneckerCache::get_or_compute(const arma::mat& A, const arma::mat& B) {
    // cache key from memory addresses and dimensions
    std::uintptr_t idA = reinterpret_cast<std::uintptr_t>(A.memptr());
    std::uintptr_t idB = reinterpret_cast<std::uintptr_t>(B.memptr());
    CacheKey key = std::make_tuple(idA, idB, A.n_rows, A.n_cols, B.n_rows, B.n_cols);
    
    auto it = cache.find(key);
    if (it != cache.end()) {
        return it->second;
    }
    
    // compute and cache
    arma::mat result = arma::kron(A, B);
    cache[key] = result;
    return result;
}