#include <RcppArmadillo.h>
#include "dbn_stability.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::field<arma::cube> piecewise_theta_update_cpp(
    const arma::field<arma::cube>& Z_block,
    const arma::cube& Theta_prev,
    const arma::cube& M,
    const arma::mat& A,
    const arma::mat& B,
    double s2,
    bool is_first_block) {

  int n_row = M.n_rows;
  int n_col = M.n_cols;
  int p = M.n_slices;
  int block_len = Z_block.n_elem;

  // floor variance
  s2 = std::max(s2, 1e-6);

  // output: field of cubes, one per time point
  arma::field<arma::cube> Theta_block(block_len);
  for (int t = 0; t < block_len; ++t) {
    Theta_block(t) = arma::cube(n_row, n_col, p);
  }

  // posterior variance for Theta given Z
  double post_var = s2 / (1.0 + s2);
  double post_sd = std::sqrt(post_var);

  // update each time point
  for (int t = 0; t < block_len; ++t) {
    for (int j = 0; j < p; ++j) {
      arma::mat E_Theta(n_row, n_col);

      if (t == 0 && is_first_block) {
        // first time point of first block: prior centered at M
        E_Theta = M.slice(j);
      } else {
        // AR update: E[Theta_t] = A * (Theta_{t-1} - M) * B' + M
        arma::mat Theta_prev_j;
        if (t == 0) {
          // first time point but not first block: use previous block's last Theta
          Theta_prev_j = Theta_prev.slice(j);
        } else {
          Theta_prev_j = Theta_block(t - 1).slice(j);
        }
        E_Theta = A * (Theta_prev_j - M.slice(j)) * B.t() + M.slice(j);
      }

      // sample from posterior: Theta_t | Z_t ~ N(post_mean, post_var)
      arma::mat Z_tj = Z_block(t).slice(j);
      arma::mat resid = Z_tj - E_Theta;
      arma::mat post_mean = E_Theta + post_var * resid;

      // add noise
      Theta_block(t).slice(j) = post_mean + post_sd * arma::randn(n_row, n_col);
    }
  }

  return Theta_block;
}


//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List piecewise_suff_stats_A_cpp(
    const arma::field<arma::cube>& Theta_block,
    const arma::cube& Theta_prev,
    const arma::cube& M,
    const arma::mat& B,
    bool is_first_block) {

  int n_row = M.n_rows;
  int p = M.n_slices;
  int block_len = Theta_block.n_elem;

  arma::mat sum_YX(n_row, n_row, arma::fill::zeros);
  arma::mat sum_XX(n_row, n_row, arma::fill::zeros);

  // start from t=1 (need t-1 for lagged value)
  int t_start = is_first_block ? 1 : 0;

  for (int t = t_start; t < block_len; ++t) {
    for (int j = 0; j < p; ++j) {
      // Y_t = Theta_t - M
      arma::mat Y_tj = Theta_block(t).slice(j) - M.slice(j);

      // X_{t-1} = Theta_{t-1} - M
      arma::mat X_tj;
      if (t == 0) {
        // use previous block's last Theta
        X_tj = Theta_prev.slice(j) - M.slice(j);
      } else {
        X_tj = Theta_block(t - 1).slice(j) - M.slice(j);
      }

      // X_t * B'
      arma::mat X_tj_B = X_tj * B.t();

      // accumulate sufficient statistics
      sum_YX += Y_tj * X_tj_B.t();
      sum_XX += X_tj_B * X_tj_B.t();
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("sum_YX") = sum_YX,
    Rcpp::Named("sum_XX") = sum_XX
  );
}


//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List piecewise_suff_stats_B_cpp(
    const arma::field<arma::cube>& Theta_block,
    const arma::cube& Theta_prev,
    const arma::cube& M,
    const arma::mat& A,
    bool is_first_block) {

  int n_col = M.n_cols;
  int p = M.n_slices;
  int block_len = Theta_block.n_elem;

  arma::mat sum_YX(n_col, n_col, arma::fill::zeros);
  arma::mat sum_XX(n_col, n_col, arma::fill::zeros);

  // start from t=1 (need t-1 for lagged value)
  int t_start = is_first_block ? 1 : 0;

  for (int t = t_start; t < block_len; ++t) {
    for (int j = 0; j < p; ++j) {
      // Y_t' = (Theta_t - M)'
      arma::mat Y_tj = (Theta_block(t).slice(j) - M.slice(j)).t();

      // X_{t-1}' = (Theta_{t-1} - M)'
      arma::mat X_tj;
      if (t == 0) {
        X_tj = (Theta_prev.slice(j) - M.slice(j)).t();
      } else {
        X_tj = (Theta_block(t - 1).slice(j) - M.slice(j)).t();
      }

      // X_t' * A'
      arma::mat X_tj_A = X_tj * A.t();

      // accumulate sufficient statistics
      sum_YX += Y_tj * X_tj_A.t();
      sum_XX += X_tj_A * X_tj_A.t();
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("sum_YX") = sum_YX,
    Rcpp::Named("sum_XX") = sum_XX
  );
}


//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List piecewise_block_update_cpp(
    const arma::cube& Z_4d,
    arma::cube& Theta_4d,
    const arma::mat& M,
    Rcpp::List A_blocks,
    Rcpp::List B_blocks,
    const arma::ivec& block_starts,
    const arma::ivec& block_ends,
    double s2,
    double t2,
    bool symmetric) {

  int n_row = Z_4d.n_rows;
  int n_col = Z_4d.n_cols;
  int K = block_starts.n_elem;

  // floor variances
  s2 = std::max(s2, 1e-6);
  t2 = std::max(t2, 1e-6);

  double post_var = s2 / (1.0 + s2);
  double post_sd = std::sqrt(post_var);

  double sse_total = 0.0;

  // process each block
  for (int k = 0; k < K; ++k) {
    int t_start = block_starts(k);
    int t_end = block_ends(k);
    int block_len = t_end - t_start + 1;

    arma::mat A_k = Rcpp::as<arma::mat>(A_blocks[k]);
    arma::mat B_k = Rcpp::as<arma::mat>(B_blocks[k]);

    // update Theta within block
    for (int t = t_start; t <= t_end; ++t) {
      arma::mat E_Theta(n_row, n_col);

      if (t == 0) {
        // first time point: prior centered at M
        E_Theta = M;
      } else {
        // AR update: E[Theta_t] = A * (Theta_{t-1} - M) * B' + M
        arma::mat Theta_prev = Theta_4d.slice(t - 1);
        E_Theta = A_k * (Theta_prev - M) * B_k.t() + M;
      }

      // sample from posterior
      arma::mat Z_t = Z_4d.slice(t);
      arma::mat resid = Z_t - E_Theta;
      arma::mat post_mean = E_Theta + post_var * resid;
      Theta_4d.slice(t) = post_mean + post_sd * arma::randn(n_row, n_col);
    }

    // update A_k if block has dynamics (more than 1 time point)
    if (block_len > 1) {
      arma::mat sum_YX(n_row, n_row, arma::fill::zeros);
      arma::mat sum_XX(n_row, n_row, arma::fill::zeros);

      for (int t = t_start + 1; t <= t_end; ++t) {
        arma::mat Y_t = Theta_4d.slice(t) - M;
        arma::mat X_t = Theta_4d.slice(t - 1) - M;
        arma::mat X_t_B = X_t * B_k.t();

        sum_YX += Y_t * X_t_B.t();
        sum_XX += X_t_B * X_t_B.t();
      }

      // posterior for A_k
      arma::mat precision = sum_XX / s2 + arma::eye(n_row, n_row) / t2;
      arma::mat V_A;
      bool solve_ok = arma::inv_sympd(V_A, force_sym(precision));
      if (!solve_ok) {
        precision.diag() += 1e-6;
        solve_ok = arma::inv_sympd(V_A, force_sym(precision));
        if (!solve_ok) {
          V_A = arma::inv(precision);
        }
      }

      // prior mean centered at scaled identity
      double diag_mean = arma::mean(A_k.diag());
      double diag_var = 1.0 / (n_row / t2 + 1.0);
      double E_diag = R::rnorm(diag_mean, std::sqrt(diag_var));
      arma::mat E_A = E_diag * arma::eye(n_row, n_row);

      // sample A from posterior
      arma::mat A_mean = (sum_YX / s2 + E_A / t2) * V_A;

      // cholesky for sampling
      arma::mat L_A;
      bool chol_ok = arma::chol(L_A, force_sym(V_A), "lower");
      arma::mat A_new;
      if (chol_ok) {
        A_new = A_mean + arma::randn(n_row, n_row) * L_A.t();
      } else {
        // fallback: diagonal noise
        arma::vec diag_sd = arma::sqrt(arma::abs(V_A.diag()));
        A_new = A_mean + arma::randn(n_row, n_row) * arma::diagmat(diag_sd);
      }

      A_blocks[k] = A_new;
      sse_total += arma::accu(arma::square(A_new - E_A));

      // update B_k (unless symmetric)
      if (!symmetric) {
        arma::mat sum_YX_B(n_col, n_col, arma::fill::zeros);
        arma::mat sum_XX_B(n_col, n_col, arma::fill::zeros);

        for (int t = t_start + 1; t <= t_end; ++t) {
          arma::mat Y_t = (Theta_4d.slice(t) - M).t();
          arma::mat X_t = (Theta_4d.slice(t - 1) - M).t();
          arma::mat X_t_A = X_t * A_new.t();

          sum_YX_B += Y_t * X_t_A.t();
          sum_XX_B += X_t_A * X_t_A.t();
        }

        // posterior for B_k
        arma::mat precision_B = sum_XX_B / s2 + arma::eye(n_col, n_col) / t2;
        arma::mat V_B;
        solve_ok = arma::inv_sympd(V_B, force_sym(precision_B));
        if (!solve_ok) {
          precision_B.diag() += 1e-6;
          solve_ok = arma::inv_sympd(V_B, force_sym(precision_B));
          if (!solve_ok) {
            V_B = arma::inv(precision_B);
          }
        }

        // prior mean
        double diag_mean_B = arma::mean(B_k.diag());
        double E_diag_B = R::rnorm(diag_mean_B, std::sqrt(diag_var));
        arma::mat E_B = E_diag_B * arma::eye(n_col, n_col);

        // sample B
        arma::mat B_mean = (sum_YX_B / s2 + E_B / t2) * V_B;
        arma::mat L_B;
        chol_ok = arma::chol(L_B, force_sym(V_B), "lower");
        arma::mat B_new;
        if (chol_ok) {
          B_new = B_mean + arma::randn(n_col, n_col) * L_B.t();
        } else {
          arma::vec diag_sd_B = arma::sqrt(arma::abs(V_B.diag()));
          B_new = B_mean + arma::randn(n_col, n_col) * arma::diagmat(diag_sd_B);
        }

        B_blocks[k] = B_new;
        sse_total += arma::accu(arma::square(B_new - E_B));
      } else {
        // symmetric: B = A
        B_blocks[k] = A_new;
        sse_total += arma::accu(arma::square(A_new - E_A));
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("Theta") = Theta_4d,
    Rcpp::Named("A_blocks") = A_blocks,
    Rcpp::Named("B_blocks") = B_blocks,
    Rcpp::Named("sse") = sse_total
  );
}


//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List piecewise_mcmc_iter_cpp(
    arma::cube Z,
    arma::cube Theta,
    arma::mat M,
    Rcpp::List A_blocks,
    Rcpp::List B_blocks,
    const arma::ivec& block_starts,
    const arma::ivec& block_ends,
    double s2,
    double t2,
    double g2,
    bool symmetric) {

  int n_row = Z.n_rows;
  int n_col = Z.n_cols;
  int T = Z.n_slices;
  int nc = n_row * n_col;

  // Update M
  arma::mat M_sum(n_row, n_col, arma::fill::zeros);
  for (int t = 0; t < T; ++t) {
    M_sum += Z.slice(t) - Theta.slice(t);
  }
  M = M_sum / (T + 1.0 / g2) + arma::randn(n_row, n_col) / std::sqrt(T + 1.0 / g2);

  // Update g2
  double g2_shape = (1.0 + nc) / 2.0;
  double g2_rate = (1.0 + arma::accu(arma::square(M))) / 2.0;
  g2 = 1.0 / R::rgamma(g2_shape, 1.0 / g2_rate);

  // Block updates (Theta, A, B)
  Rcpp::List block_result = piecewise_block_update_cpp(
    Z, Theta, M, A_blocks, B_blocks,
    block_starts, block_ends, s2, t2, symmetric
  );

  Theta = Rcpp::as<arma::cube>(block_result["Theta"]);
  A_blocks = block_result["A_blocks"];
  B_blocks = block_result["B_blocks"];
  double sse_total = Rcpp::as<double>(block_result["sse"]);

  // Update s2
  double sse_theta = 0.0;
  for (int t = 1; t < T; ++t) {
    sse_theta += arma::accu(arma::square(Theta.slice(t) - Theta.slice(t - 1)));
  }
  double s2_shape = (1.0 + nc * (T - 1)) / 2.0;
  double s2_rate = (1.0 + sse_theta) / 2.0;
  s2 = 1.0 / R::rgamma(s2_shape, 1.0 / s2_rate);

  // Update t2
  int K = block_starts.n_elem;
  int n_AB_params = K * (n_row * n_row + n_col * n_col);
  if (symmetric) n_AB_params = K * n_row * n_row;
  double t2_shape = (1.0 + n_AB_params) / 2.0;
  double t2_rate = (1.0 + sse_total) / 2.0;
  t2 = 1.0 / R::rgamma(t2_shape, 1.0 / t2_rate);

  return Rcpp::List::create(
    Rcpp::Named("Z") = Z,
    Rcpp::Named("Theta") = Theta,
    Rcpp::Named("M") = M,
    Rcpp::Named("A_blocks") = A_blocks,
    Rcpp::Named("B_blocks") = B_blocks,
    Rcpp::Named("s2") = s2,
    Rcpp::Named("t2") = t2,
    Rcpp::Named("g2") = g2
  );
}
