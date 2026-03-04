#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' Update A and B matrices for static model (C++ version)
//'
//' @description Fast C++ implementation of update_AB_static.
//'   Supports bipartite networks where Theta is n_row x n_col (rectangular).
//'   A is n_row x n_row (sender dynamics), B is n_col x n_col (receiver dynamics).
//' @param Theta_prev Previous Theta values (n_row x n_col x n_times)
//' @param Theta_curr Current Theta values (n_row x n_col x n_times)
//' @param B_init Initial B matrix (n_col x n_col)
//' @param tau_A2 Prior variance for A
//' @param tau_B2 Prior variance for B
//' @param sigma2 Innovation variance
//' @return List with updated A and B matrices
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List update_AB_static_cpp(const arma::cube& Theta_prev,
                                const arma::cube& Theta_curr,
                                const arma::mat& B_init,
                                double tau_A2,
                                double tau_B2,
                                double sigma2) {

  int n_row = Theta_prev.n_rows;
  int n_col = Theta_prev.n_cols;
  int n_times = Theta_prev.n_slices;

  // initialize return matrices: A is n_row x n_row, B is n_col x n_col
  arma::mat A(n_row, n_row, arma::fill::zeros);
  arma::mat B = B_init;

  // handle empty case
  if (n_times == 0) {
    A.eye();
    B.eye();
    return Rcpp::List::create(
      Rcpp::Named("A") = A,
      Rcpp::Named("B") = B
    );
  }

  // regularization parameter
  double reg_A = sigma2 / tau_A2;
  double reg_B = sigma2 / tau_B2;

  // update A row by row: row i of A has n_row parameters
  // model: Theta_curr[i,:] = A[i,:] * Theta_prev * B'
  // each time step gives n_col observations
  for (int i = 0; i < n_row; i++) {
    arma::mat Xi(n_times * n_col, n_row);
    arma::vec yi(n_times * n_col);

    for (int t = 0; t < n_times; t++) {
      // Theta_B = Theta_prev * B' is (n_row x n_col)
      arma::mat Theta_B = Theta_prev.slice(t) * B.t();

      // design: Theta_B.t() is (n_col x n_row) — n_col observations per time step
      Xi.rows(t * n_col, (t + 1) * n_col - 1) = Theta_B.t();

      // response: row i of Theta_curr (n_col elements)
      yi.subvec(t * n_col, (t + 1) * n_col - 1) = Theta_curr.slice(t).row(i).t();
    }

    // ridge regression: n_row parameters
    arma::mat XtX = Xi.t() * Xi + reg_A * arma::eye(n_row, n_row);
    arma::vec XtY = Xi.t() * yi;

    arma::vec ai;
    bool solved = arma::solve(ai, XtX, XtY, arma::solve_opts::likely_sympd);
    if (!solved) {
      ai = arma::pinv(XtX) * XtY;
    }

    A.row(i) = ai.t();
  }

  // update B row by row: row j of B has n_col parameters
  // model: Theta_curr[:,j] = A * Theta_prev * B[j,:]'
  // each time step gives n_row observations
  for (int j = 0; j < n_col; j++) {
    arma::mat Xj(n_times * n_row, n_col);
    arma::vec yj(n_times * n_row);

    for (int t = 0; t < n_times; t++) {
      // A_Theta = A * Theta_prev is (n_row x n_col)
      arma::mat A_Theta = A * Theta_prev.slice(t);

      // design: A_Theta is (n_row x n_col) — n_row observations per time step
      Xj.rows(t * n_row, (t + 1) * n_row - 1) = A_Theta;

      // response: column j of Theta_curr (n_row elements)
      yj.subvec(t * n_row, (t + 1) * n_row - 1) = Theta_curr.slice(t).col(j);
    }

    // ridge regression: n_col parameters
    arma::mat XtX = Xj.t() * Xj + reg_B * arma::eye(n_col, n_col);
    arma::vec XtY = Xj.t() * yj;

    arma::vec bj;
    bool solved = arma::solve(bj, XtX, XtY, arma::solve_opts::likely_sympd);
    if (!solved) {
      bj = arma::pinv(XtX) * XtY;
    }

    B.col(j) = bj;
  }

  return Rcpp::List::create(
    Rcpp::Named("A") = A,
    Rcpp::Named("B") = B
  );
}

//' Build design matrix F for alpha updates (C++ version)
//' 
//' @description Fast C++ implementation of build_F_alpha
//' @param U Orthonormal matrix (m x r)
//' @param Theta_prev Previous Theta values (m x m x T)
//' @param B_array B matrices over time (m x m x T)
//' @param compute_all Compute for all time points at once
//' @return Design matrix F (m^2*T x r) if compute_all=true, else (m^2 x r)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat build_F_alpha_cpp(const arma::mat& U,
                            const arma::cube& Theta_prev,
                            const arma::cube& B_array,
                            bool compute_all = false) {
  
  int m = U.n_rows;
  int r = U.n_cols;
  int T = Theta_prev.n_slices;
  
  if (!compute_all) {
    // single time point (use first slice)
    arma::mat F(m * m, r);
    arma::mat Theta_B = Theta_prev.slice(0) * B_array.slice(0).t();
    
    // compute each column: vec(S_k * Theta_prev * B^T)
    for (int k = 0; k < r; k++) {
      arma::vec u_k = U.col(k);
      arma::mat S_k = u_k * u_k.t();
      F.col(k) = arma::vectorise(S_k * Theta_B);
    }
    
    return F;
  } else {
    // all time points
    arma::mat F(m * m * T, r);
    
    for (int t = 0; t < T; t++) {
      int row_start = t * m * m;
      arma::mat Theta_B = Theta_prev.slice(t) * B_array.slice(t).t();
      
      for (int k = 0; k < r; k++) {
        arma::vec u_k = U.col(k);
        arma::mat S_k = u_k * u_k.t();
        F(arma::span(row_start, row_start + m*m - 1), k) = 
          arma::vectorise(S_k * Theta_B);
      }
    }
    
    return F;
  }
}