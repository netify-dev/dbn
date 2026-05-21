#ifndef DBN_STABILITY_H
#define DBN_STABILITY_H

#include <RcppArmadillo.h>
#include <map>
#include <random>
#include <cmath>

// stability constants
constexpr double STABILITY_TOL = 1e-10;
constexpr double SPECTRAL_RADIUS_THRESHOLD = 0.995;
constexpr double CHOLESKY_REGULARIZATION = 1e-6;
constexpr int MAX_CHOLESKY_ATTEMPTS = 5;

// force symmetry
inline arma::mat force_sym(const arma::mat& M) {
    return 0.5 * (M + M.t());
}

// shared safe symmetric-PD inverse with regularization fallback.
// This consolidates duplicate copies that previously lived in
// dynamic.cpp / update_AB_batch.cpp / lowrank_parallel.cpp.
//
// This helper never throws. The previous final fallback used the
// throwing form arma::inv(M_sym); when M_sym is non-finite (e.g. the
// asymmetric sampler has drifted along its scale orbit until A_t/B_t
// overflow) that throw escaped an OpenMP parallel region and triggered
// std::terminate -- a silent process death with no R-level error. All
// paths below use bool-return Armadillo calls and a finite diagonal
// last resort, so callers always get a finite-dimensioned matrix back
// and any non-finiteness is surfaced by the R-side divergence guard.
inline arma::mat dbn_safe_inv_sympd(const arma::mat& M) {
    arma::mat M_sym = 0.5 * (M + M.t());
    arma::mat result;
    bool ok = arma::inv_sympd(result, M_sym);
    if (!ok) {
        double fro = arma::norm(M_sym, "fro");
        double reg = (std::isfinite(fro) ? 1e-6 * fro : 1e-6) + 1e-8;
        M_sym.diag() += reg;
        ok = arma::inv_sympd(result, M_sym);
        if (!ok) {
            ok = arma::inv(result, M_sym);   // bool form: does not throw
        }
        if (!ok) {
            // last resort: invert the (regularized) diagonal only; this
            // never throws and never returns wrong dimensions.
            arma::vec d = M_sym.diag();
            d.transform([](double v) {
                return (std::isfinite(v) && std::abs(v) > 1e-12) ? 1.0 / v : 0.0;
            });
            result = arma::diagmat(d);
        }
    }
    return result;
}

// shared safe multivariate-normal draw (uses R's RNG via arma::mvnrnd).
// NOT thread safe. For thread-safe sampling use the explicit RNG-passing
// version in dynamic.cpp (thread_safe_mvnrnd).
inline arma::vec dbn_safe_mvnrnd(const arma::vec& mu, const arma::mat& Sigma) {
    arma::mat S = 0.5 * (Sigma + Sigma.t());
    if (S.n_rows != mu.n_elem || S.n_cols != mu.n_elem) {
        return mu + 1e-4 * arma::randn(mu.n_elem);
    }
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