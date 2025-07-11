#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' Update A and B matrices for static model (C++ version)
//' 
//' @description Fast C++ implementation of update_AB_static
//' @param Theta_prev Previous Theta values (m x m x n_times)
//' @param Theta_curr Current Theta values (m x m x n_times) 
//' @param B_init Initial B matrix (m x m)
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
  
  int m = Theta_prev.n_rows;
  int n_times = Theta_prev.n_slices;
  
  // initialize return matrices
  arma::mat A(m, m, arma::fill::zeros);
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
  
  // update A row by row to avoid large matrices
  for (int i = 0; i < m; i++) {
    // build design matrix for row i of A
    arma::mat Xi(n_times * m, m);
    arma::vec yi(n_times * m);
    
    for (int t = 0; t < n_times; t++) {
      // for row i of A: Y_i = A_i * (Theta_prev * B^T)
      arma::mat Theta_B = Theta_prev.slice(t) * B.t();
      
      // design matrix: each row is a row of Theta_B
      Xi.rows(t * m, (t + 1) * m - 1) = Theta_B.t();
      
      // response: column i of current Theta
      yi.subvec(t * m, (t + 1) * m - 1) = Theta_curr.slice(t).col(i);
    }
    
    // ridge regression for row i
    arma::mat XtX = Xi.t() * Xi + reg_A * arma::eye(m, m);
    arma::vec XtY = Xi.t() * yi;
    
    // solve using Cholesky decomposition
    arma::vec ai;
    bool solved = arma::solve(ai, XtX, XtY, arma::solve_opts::likely_sympd);
    if (!solved) {
      // fallback to SVD if Cholesky fails (make sure to implement in batch)
      ai = arma::pinv(XtX) * XtY;
    }
    
    A.row(i) = ai.t();
  }
  
  // update B column by column using the new A
  for (int j = 0; j < m; j++) {
    // build design matrix for column j of B
    arma::mat Xj(n_times * m, m);
    arma::vec yj(n_times * m);
    
    for (int t = 0; t < n_times; t++) {
      // for column j of B: Y_j = (A * Theta_prev) * B_j
      arma::mat A_Theta = A * Theta_prev.slice(t);
      
      // design matrix: each row is a row of A_Theta
      Xj.rows(t * m, (t + 1) * m - 1) = A_Theta.t();
      
      // response: column j of current Theta
      yj.subvec(t * m, (t + 1) * m - 1) = Theta_curr.slice(t).col(j);
    }
    
    // ridge regression for column j
    arma::mat XtX = Xj.t() * Xj + reg_B * arma::eye(m, m);
    arma::vec XtY = Xj.t() * yj;
    
    // solve using Cholesky decomposition
    arma::vec bj;
    bool solved = arma::solve(bj, XtX, XtY, arma::solve_opts::likely_sympd);
    if (!solved) {
      // fallback to SVD if Cholesky fails (make sure to implement in batch)
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