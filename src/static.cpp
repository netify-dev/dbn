#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "thread_control.h"

using namespace Rcpp;
using namespace arma;

//' Reshape 4D array to 3D for C++ processing
//' @description Efficiently reshape Z from n_row x n_col x p x Tt to n_row x n_col x (p*Tt)
//' @param Z_4d Input 4D array as R array
//' @param n_row Number of row nodes (senders)
//' @param n_col Number of column nodes (receivers)
//' @param p Number of relations
//' @param Tt Number of time points
//' @return 3D cube for C++ processing
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube reshape_Z_to_cube(const NumericVector& Z_4d, int n_row, int n_col, int p, int Tt) {
    int nc = n_row * n_col;
    arma::cube Z_cube(n_row, n_col, p * Tt);

    const double* Z_ptr = Z_4d.begin();

    for (int slice = 0; slice < p * Tt; slice++) {
        int j = slice / Tt;
        int t = slice % Tt;

        double* cube_slice = Z_cube.slice_memptr(slice);
        const double* z_offset = Z_ptr + (t * nc * p + j * nc);

        std::memcpy(cube_slice, z_offset, nc * sizeof(double));
    }

    return Z_cube;
}

//' Compute SSE for diagonal elements efficiently
//' @description Compute sum of squared errors for diagonal elements of B matrices
//' @param B_list List of B matrices
//' @param K Number of matrices
//' @return SSE value
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_diagonal_sse(const List& B_list, int K) {
    double sse = 0.0;

    for (int k = 0; k < K; k++) {
        arma::mat B = as<arma::mat>(B_list[k]);
        arma::vec diag_B = B.diag();
        arma::vec diff = diag_B - 1.0;
        sse += arma::dot(diff, diff);
    }

    return sse;
}


//' Compute sum of squared deviations from identity for A or B arrays
//' @description Efficiently compute sum of squared deviations from the identity matrix
//'   for each matrix slice of A or B from time t = 2 to Tt
//' @param ABarray 3D array of A or B matrices (m x m x Tt)
//' @param m Number of nodes (dimension of the square matrix)
//' @param Tt Number of time points
//' @return Sum of squared deviations
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_deviation_sum(const arma::cube& ABarray, int m, int Tt) {
    double sum_sq = 0.0;
    arma::mat I = arma::eye(m, m);

    // start from t=2 (index 1), matching the R code
    for (int t = 1; t < Tt; t++) {
        arma::mat diff = ABarray.slice(t) - I;
        sum_sq += arma::accu(diff % diff);
    }

    return sum_sq;
}

//' Compute mean M for static model
//' @description Efficiently compute mean of Z across time for each relation
//' @param Z_flat Flattened Z array (n_row*n_col x p*Tt)
//' @param n_row Number of row nodes
//' @param n_col Number of column nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return M array (n_row x n_col x p)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_M_static(const arma::mat& Z_flat, int n_row, int n_col, int p, int Tt) {
    set_dbn_threads();
    int nc = n_row * n_col;
    arma::cube M(n_row, n_col, p);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int j = 0; j < p; j++) {

        arma::vec ones_vec(Tt, arma::fill::ones);
        arma::mat Z_j_cols(nc, Tt);

        for (int t = 0; t < Tt; t++) {
            Z_j_cols.col(t) = Z_flat.col(j * Tt + t);
        }

        arma::vec sum_vec = Z_j_cols * ones_vec;
        M.slice(j) = arma::reshape(sum_vec / Tt, n_row, n_col);
    }

    return M;
}

//' Compute residual sum of squares for static model
//' @description Compute sum((Z - M)^2) across all times and relations
//' @param Z_4d Original Z array (n_row x n_col x p x Tt) as flattened vector
//' @param M Mean array (n_row x n_col x p)
//' @param n_row Number of row nodes
//' @param n_col Number of column nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return Sum of squared residuals
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_rss_static(const NumericVector& Z_4d, const arma::cube& M,
                         int n_row, int n_col, int p, int Tt) {
    double rss = 0.0;
    int nc = n_row * n_col;

    // wrap raw R vector for zero-copy column-major views
    const double* z_ptr = REAL(Z_4d);

    for (int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);

        for (int t = 0; t < Tt; t++) {
            int offset = j * nc + t * nc * p;
            // zero-copy view of this n_row x n_col slice
            arma::mat Z_jt(const_cast<double*>(z_ptr + offset), n_row, n_col, false, true);
            arma::mat diff = Z_jt - M_j;
            rss += arma::accu(diff % diff);
        }
    }

    return rss;
}

//' Reshape for large networks with parallelization
//' @description Parallel version of reshape_Z_to_cube for large networks
//' @param Z_4d Input 4D array as R array
//' @param n_row Number of row nodes
//' @param n_col Number of column nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return 3D cube for C++ processing
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube reshape_Z_to_cube_parallel(const NumericVector& Z_4d, int n_row, int n_col, int p, int Tt) {
    set_dbn_threads();
    int nc = n_row * n_col;
    arma::cube Z_cube(n_row, n_col, p * Tt);

    const double* Z_ptr = Z_4d.begin();

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 4)
    #endif
    for (int slice = 0; slice < p * Tt; slice++) {
        int j = slice / Tt;
        int t = slice % Tt;

        double* cube_slice = Z_cube.slice_memptr(slice);
        const double* z_offset = Z_ptr + (t * nc * p + j * nc);

        std::memcpy(cube_slice, z_offset, nc * sizeof(double));
    }

    return Z_cube;
}

//' M computation with blocked summation for numerical stability
//' @description Compute mean with blocked summation for large networks
//' @param Z_flat Flattened Z array (n_row*n_col x p*Tt)
//' @param n_row Number of row nodes
//' @param n_col Number of column nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return M array (n_row x n_col x p)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube compute_M_static_blocked(const arma::mat& Z_flat, int n_row, int n_col, int p, int Tt) {
    set_dbn_threads();
    int nc = n_row * n_col;
    arma::cube M(n_row, n_col, p, arma::fill::zeros);
    const int block_size = 64;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int j = 0; j < p; j++) {
        arma::mat Z_j_sum(n_row, n_col, arma::fill::zeros);

        for (int t_block = 0; t_block < Tt; t_block += block_size) {
            int t_end = std::min(t_block + block_size, Tt);

            for (int t = t_block; t < t_end; t++) {
                int col_idx = j * Tt + t;
                arma::vec z_vec = Z_flat.col(col_idx);
                Z_j_sum += arma::reshape(z_vec, n_row, n_col);
            }
        }

        M.slice(j) = Z_j_sum / Tt;
    }

    return M;
}

//' RSS computation with parallel reduction
//' @description Parallel computation of residual sum of squares
//' @param Z_cube Z array as cube (n_row x n_col x p*Tt)
//' @param M Mean array (n_row x n_col x p)
//' @param n_row Number of row nodes
//' @param n_col Number of column nodes
//' @param p Number of relations
//' @param Tt Number of time points
//' @return Sum of squared residuals
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_rss_static_parallel(const arma::cube& Z_cube, const arma::cube& M,
                                  int n_row, int n_col, int p, int Tt) {
    double rss = 0.0;

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:rss)
    #endif
    for (int j = 0; j < p; j++) {
        double local_rss = 0.0;
        arma::mat M_j = M.slice(j);

        for (int t = 0; t < Tt; t++) {
            int idx = j * Tt + t;
            arma::mat diff = Z_cube.slice(idx) - M_j;
            local_rss += arma::accu(diff % diff);
        }
        rss += local_rss;
    }

    return rss;
}

//' Cache-efficient B update for static model (sender effects)
//' @description B1 (sender effects, n_row x n_row) update using tiled matrix operations.
//'   Model: Z_jt = B * M_j + E, so regression is column-wise:
//'   z_c = B * m_c, giving XtX = sum M_j * M_j' and XtY = sum Z_jt * M_j'
//' @param Z_cube Z array as cube (n_row x n_col x p*Tt)
//' @param M Mean array (n_row x n_col x p)
//' @param s2 Observation variance
//' @param t2 Prior variance
//' @param n_row Number of row nodes (senders)
//' @param n_col Number of column nodes (receivers)
//' @param p Number of relations
//' @param Tt Number of time points
//' @return Updated B matrix (n_row x n_row)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat update_B_static_tiled(const arma::cube& Z_cube, const arma::cube& M,
                                double s2, double t2, int n_row, int n_col, int p, int Tt) {
    set_dbn_threads();

    // sender effects: n_row x n_row
    arma::mat XtX(n_row, n_row, arma::fill::zeros);
    arma::mat XtY(n_row, n_row, arma::fill::zeros);

    // XtX = Tt * sum_j M_j * M_j', XtY = sum_{j,t} Z_jt * M_j'
    #ifdef _OPENMP
    #pragma omp parallel
    {
        arma::mat local_XtX(n_row, n_row, arma::fill::zeros);
        arma::mat local_XtY(n_row, n_row, arma::fill::zeros);

        #pragma omp for schedule(static)
        for (int j = 0; j < p; j++) {
            arma::mat M_j = M.slice(j);
            // M_j * M_j': (n_row x n_col) * (n_col x n_row) = n_row x n_row
            local_XtX += Tt * (M_j * M_j.t());

            arma::mat sum_Z(n_row, n_col, arma::fill::zeros);
            for (int t = 0; t < Tt; t++) {
                sum_Z += Z_cube.slice(j * Tt + t);
            }
            // (n_row x n_col) * (n_col x n_row) = n_row x n_row
            local_XtY += sum_Z * M_j.t();
        }

        #pragma omp critical
        {
            XtX += local_XtX;
            XtY += local_XtY;
        }
    }
    #else
    for (int j = 0; j < p; j++) {
        arma::mat M_j = M.slice(j);
        XtX += Tt * (M_j * M_j.t());

        arma::mat sum_Z(n_row, n_col, arma::fill::zeros);
        for (int t = 0; t < Tt; t++) {
            sum_Z += Z_cube.slice(j * Tt + t);
        }
        XtY += sum_Z * M_j.t();
    }
    #endif

    // posterior precision and mean
    double lambda = 1.0 / t2;
    double gamma = 1.0 / s2;

    XtX.diag() += (lambda / gamma) * arma::ones(n_row);

    arma::mat gXtX = gamma * XtX;
    gXtX = 0.5 * (gXtX + gXtX.t());
    arma::mat L_post;
    bool chol_success = arma::chol(L_post, gXtX, "lower");

    if (!chol_success) {
        arma::mat U, V;
        arma::vec s;
        bool svd_ok = arma::svd_econ(U, s, V, XtX);

        if (!svd_ok || s.n_elem == 0) {
            // degenerate: return identity + noise
            return arma::eye(n_row, n_row) + 0.01 * arma::randn(n_row, n_row);
        }

        double eps = 1e-8 * s.max();
        s = arma::max(s, eps * arma::ones(s.n_elem));

        arma::mat post_cov = V * arma::diagmat(1.0 / (gamma * s + lambda)) * V.t();
        arma::mat post_mean = post_cov * (gamma * XtY + lambda * arma::eye(n_row, n_row));

        arma::mat noise = arma::randn(n_row, n_row);
        return post_mean + V * arma::diagmat(arma::sqrt(1.0 / (gamma * s + lambda))) * V.t() * noise;
    }

    arma::mat post_mean_part = arma::solve(arma::trimatl(L_post), gamma * XtY + lambda * arma::eye(n_row, n_row));
    arma::mat post_mean = arma::solve(arma::trimatu(L_post.t()), post_mean_part);

    arma::mat Z_sample = arma::randn(n_row, n_row);
    arma::mat B_new = post_mean + arma::solve(arma::trimatu(L_post.t()), Z_sample) / sqrt(gamma);

    return B_new;
}
