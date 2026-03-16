#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include <random>
#include "thread_control.h"

using namespace Rcpp;
using namespace arma;

//' Fast parallel Gaussian approximation for ordinal data
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube rz_gaussian_approx_cpp(const arma::cube& R, const arma::cube& Z,
                                  const arma::cube& EZ, double sigma = 1.0) {
    set_dbn_threads();

    arma::cube Z_new = Z; // copy
    int n_slices = R.n_slices;

    // initialize per-thread RNGs from R's RNG (main thread only)
#ifdef _OPENMP
    int n_threads = std::min(omp_get_max_threads(), n_slices);
#else
    int n_threads = 1;
#endif
    std::vector<std::mt19937_64> rngs;
    rngs.reserve(n_threads);
    for (int i = 0; i < n_threads; i++) {
        uint64_t seed = static_cast<uint64_t>(R::runif(0.0, 1.0) * 4294967296.0);
        seed ^= static_cast<uint64_t>(i + 1) * 2654435761ULL;
        rngs.emplace_back(seed);
    }

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
    #endif
    for (int s = 0; s < n_slices; s++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        arma::mat R_slice = R.slice(s);
        arma::mat Z_slice = Z.slice(s);
        arma::mat EZ_slice = EZ.slice(s);

        // non-NA entries
        arma::uvec valid_idx = find_finite(R_slice);

        if (valid_idx.n_elem > 0) {
            // valid values
            arma::vec R_valid = R_slice.elem(valid_idx);
            arma::vec EZ_valid = EZ_slice.elem(valid_idx);

            // normalize ranks
            double R_min = R_valid.min();
            double R_max = R_valid.max();

            arma::vec Z_updates(valid_idx.n_elem);

            // thread-safe normal sampling
            std::normal_distribution<double> norm_dist(0.0, sigma);
            arma::vec noise(valid_idx.n_elem);
            for (arma::uword i = 0; i < valid_idx.n_elem; i++) {
                noise(i) = norm_dist(rngs[tid]);
            }

            if (R_max > R_min) {
                arma::vec R_norm = (R_valid - R_min) / (R_max - R_min);

                // bias toward rank ordering
                arma::vec means = EZ_valid + 0.1 * (R_norm - 0.5);
                Z_updates = means + noise;
            } else {
                // constant rank: sample from prior
                Z_updates = EZ_valid + noise;
            }

            // update Z_new
            Z_slice.elem(valid_idx) = Z_updates;
            Z_new.slice(s) = Z_slice;
        }
    }

    return Z_new;
}
