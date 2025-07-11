#ifndef DBN_STABILITY_H
#define DBN_STABILITY_H

#include <RcppArmadillo.h>
#include <map>

// constants for keeping calculations numerically stable
constexpr double STABILITY_TOL = 1e-10;
constexpr double SPECTRAL_RADIUS_THRESHOLD = 0.995;
constexpr double CHOLESKY_REGULARIZATION = 1e-6;
constexpr int MAX_CHOLESKY_ATTEMPTS = 5;

// function declarations for stability utility functions
arma::mat stabilize_spectral_radius(const arma::mat& M, double threshold = SPECTRAL_RADIUS_THRESHOLD);
bool safe_cholesky(arma::mat& L, const arma::mat& A, double reg = CHOLESKY_REGULARIZATION);
arma::mat ensure_positive_definite(const arma::mat& M, double min_eigenvalue = STABILITY_TOL);
bool is_stationary(const arma::mat& A, const arma::mat& B, int p, int q);

// cache structures for storing kronecker products
struct KroneckerCache {
    // key includes dimensions and memory addresses to prevent hash collisions
    using CacheKey = std::tuple<std::uintptr_t, std::uintptr_t, int, int, int, int>;
    std::map<CacheKey, arma::mat> cache;
    
    arma::mat get_or_compute(const arma::mat& A, const arma::mat& B);
    void clear() { cache.clear(); }
};

#endif // DBN_STABILITY_H