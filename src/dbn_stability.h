#ifndef DBN_STABILITY_H
#define DBN_STABILITY_H

#include <RcppArmadillo.h>
#include <map>

// stability constants
constexpr double STABILITY_TOL = 1e-10;
constexpr double SPECTRAL_RADIUS_THRESHOLD = 0.995;
constexpr double CHOLESKY_REGULARIZATION = 1e-6;
constexpr int MAX_CHOLESKY_ATTEMPTS = 5;

// force symmetry
inline arma::mat force_sym(const arma::mat& M) {
    return 0.5 * (M + M.t());
}

// stability helpers
arma::mat stabilize_spectral_radius(const arma::mat& M, double threshold = SPECTRAL_RADIUS_THRESHOLD);
bool safe_cholesky(arma::mat& L, const arma::mat& A, double reg = CHOLESKY_REGULARIZATION);
arma::mat ensure_positive_definite(const arma::mat& M, double min_eigenvalue = STABILITY_TOL);
bool is_stationary(const arma::mat& A, const arma::mat& B, int p, int q);

// memoized kronecker product
struct KroneckerCache {
    // key: memory addresses + dimensions
    using CacheKey = std::tuple<std::uintptr_t, std::uintptr_t, int, int, int, int>;
    std::map<CacheKey, arma::mat> cache;
    
    arma::mat get_or_compute(const arma::mat& A, const arma::mat& B);
    void clear() { cache.clear(); }
};

#endif // DBN_STABILITY_H