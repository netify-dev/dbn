#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute impulse response for constant A,B matrices
//' 
//' @param A Transition matrix A (m x m)
//' @param B Transition matrix B (m x m)
//' @param S Shock matrix (m x m)
//' @param H Number of horizons to compute
//' @return Cube of impulse responses (m x m x H+1)
//' @export
// [[Rcpp::export]]
arma::cube impulse_response_const(const arma::mat& A,
                                  const arma::mat& B,
                                  const arma::mat& S,
                                  int H) {
    int m = A.n_rows;
    
    // Input validation
    if (A.n_rows != A.n_cols || B.n_rows != B.n_cols || 
        A.n_rows != B.n_rows || A.n_rows != S.n_rows || 
        S.n_rows != S.n_cols) {
        Rcpp::stop("Matrix dimensions must be consistent and square");
    }
    
    if (H < 0) {
        Rcpp::stop("H must be non-negative");
    }
    
    arma::cube Delta(m, m, H + 1);
    Delta.slice(0) = S;

    arma::mat A_pow = A;           // A^1
    arma::mat B_powT = B.t();      // (B^1)^T

    for (int h = 1; h <= H; ++h) {
        Delta.slice(h) = A_pow * S * B_powT;

        // update powers for next horizon
        if (h < H) {
            A_pow = A_pow * A;          // A^{h+1}
            B_powT = B_powT * B.t();    // (B^{h+1})^T
        }
    }
    return Delta;           // m × m × (H+1)
}

//' Compute impulse response for time-varying A,B matrices
//' 
//' @param Aarray Cube of A matrices over time (m x m x T)
//' @param Barray Cube of B matrices over time (m x m x T)
//' @param S Shock matrix (m x m)
//' @param t0 Time index of shock (0-based)
//' @param H Number of horizons to compute
//' @return Cube of impulse responses (m x m x H+1)
//' @export
// [[Rcpp::export]]
arma::cube impulse_response_dynamic(const arma::cube& Aarray,
                                    const arma::cube& Barray,
                                    const arma::mat& S,
                                    int t0,
                                    int H) {
    int m = Aarray.n_rows;
    int T = Aarray.n_slices;
    
    // Input validation
    if (Aarray.n_rows != Aarray.n_cols || Barray.n_rows != Barray.n_cols || 
        Aarray.n_rows != Barray.n_rows || Aarray.n_rows != S.n_rows || 
        S.n_rows != S.n_cols) {
        Rcpp::stop("Matrix dimensions must be consistent and square");
    }
    
    if (Aarray.n_slices != Barray.n_slices) {
        Rcpp::stop("Aarray and Barray must have same number of time slices");
    }
    
    if (t0 < 0 || t0 >= T) {
        Rcpp::stop("t0 must be between 0 and T-1");
    }
    
    if (H < 0) {
        Rcpp::stop("H must be non-negative");
    }
    
    if (t0 + H >= T) {
        Rcpp::stop("t0 + H must be less than T (number of time slices)");
    }
    
    arma::cube Delta(m, m, H + 1, arma::fill::none);
    Delta.slice(0) = S;

    for (int h = 1; h <= H; ++h) {
        int tp = t0 + h;                    // absolute time index
        Delta.slice(h) = Aarray.slice(tp) *
                         Delta.slice(h-1) *
                         Barray.slice(tp).t();
    }
    return Delta;       // m × m × (H+1)
}