#include <RcppArmadillo.h>
#include "dbn_stability.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Sample from inverse Wishart distribution with stability checks
//' 
//' @param nu Degrees of freedom
//' @param S Scale matrix
//' @return Sample from inverse Wishart distribution
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat rinvwishart_stable(int nu, const arma::mat& S) {
    int p = S.n_rows;
    
    // make sure nu > p + 1 for proper distribution
    if (nu <= p + 1) {
        Rcpp::stop("Degrees of freedom must be greater than dimension + 1");
    }
    
    // make sure s is positive definite
    arma::mat S_pd = ensure_positive_definite(S);
    
    // cholesky decomposition of s
    arma::mat L_S;
    if (!safe_cholesky(L_S, S_pd)) {
        Rcpp::stop("Failed to compute Cholesky decomposition of scale matrix");
    }
    
    // generate wishart(nu, i) variate
    arma::mat A(p, p, arma::fill::zeros);
    
    // fill diagonal with chi-squared variates
    for (int i = 0; i < p; i++) {
        A(i, i) = std::sqrt(R::rchisq(nu - i));
    }
    
    // fill lower triangle with standard normals
    for (int i = 1; i < p; i++) {
        for (int j = 0; j < i; j++) {
            A(i, j) = R::rnorm(0, 1);
        }
    }
    
    // w = a * a'
    arma::mat W = A * A.t();
    
    // transform to inverse wishart
    arma::mat W_inv = arma::inv_sympd(W);
    arma::mat result = L_S * W_inv * L_S.t();
    
    // make sure result is symmetric and positive definite
    result = 0.5 * (result + result.t());
    return ensure_positive_definite(result);
}

//' Check if proposed innovation covariance is valid
//' 
//' @param Sigma_e Proposed innovation covariance
//' @param min_eigenvalue Minimum allowed eigenvalue
//' @param max_condition Maximum allowed condition number
//' @return Validity indicator
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
bool is_valid_innovation_cov(const arma::mat& Sigma_e, 
                           double min_eigenvalue = 1e-6,
                           double max_condition = 1e6) {
    arma::vec eigenvalues;
    arma::eig_sym(eigenvalues, Sigma_e);
    
    // check minimum eigenvalue
    double min_eig = eigenvalues.min();
    if (min_eig < min_eigenvalue) {
        return false;
    }
    
    // check condition number
    double max_eig = eigenvalues.max();
    double condition_number = max_eig / min_eig;
    if (condition_number > max_condition) {
        return false;
    }
    
    return true;
}

//' Adaptive Metropolis-Hastings update for innovation covariance
//' 
//' @param Sigma_e_current Current innovation covariance
//' @param nu Prior degrees of freedom
//' @param S Prior scale matrix
//' @param X Latent states (p x q x T)
//' @param A Actor transition matrix
//' @param B Latent transition matrix
//' @param adapt_scale Adaptation scale parameter
//' @return List with new covariance and acceptance indicator
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List update_innovation_cov_adaptive(const arma::mat& Sigma_e_current,
                                        int nu,
                                        const arma::mat& S,
                                        const arma::cube& X,
                                        const arma::mat& A,
                                        const arma::mat& B,
                                        double adapt_scale = 0.1) {
    int T = X.n_slices;
    int p = X.n_rows;
    int q = X.n_cols;
    
    // compute sum of squared innovations
    arma::mat SS(q, q, arma::fill::zeros);
    for (int t = 1; t < T; t++) {
        arma::mat innovation = X.slice(t) - A * X.slice(t-1) * B.t();
        SS += innovation.t() * innovation;
    }
    
    // posterior parameters
    int nu_post = nu + (T - 1) * p;
    arma::mat S_post = S + SS;
    
    // make sure s_post is positive definite
    S_post = ensure_positive_definite(S_post);
    
    // sample proposal from posterior
    arma::mat Sigma_e_prop = rinvwishart_stable(nu_post, S_post);
    
    // early rejection based on validity
    if (!is_valid_innovation_cov(Sigma_e_prop)) {
        return Rcpp::List::create(
            Rcpp::Named("Sigma_e") = Sigma_e_current,
            Rcpp::Named("accepted") = false
        );
    }
    
    // for conjugate update, always accept
    return Rcpp::List::create(
        Rcpp::Named("Sigma_e") = Sigma_e_prop,
        Rcpp::Named("accepted") = true
    );
}