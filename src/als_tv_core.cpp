// RW(1)-smoothed time-varying ALS core solver (Stage 1A)
//
// The directed-Gaussian solver alternates between solving the entire
// A-trajectory conditional on B (block-tridiagonal system per row, Thomas
// algorithm) and the entire B-trajectory conditional on A.
//
// Conventions:
//   Phi:    [n_row x n_col x Tt] cube of centered states (Z - M)
//   A_arr:  [n_row x n_row x Tt] cube; slice 0 (t=1) is the initial state
//           and is not updated by these routines (the optimization is
//           over A_arr[, , 1..Tt-1] which correspond to model times 2..T)
//   B_arr:  [n_col x n_col x Tt] cube; same convention
//   Omega:  [n_row x n_col x Tt] cube of response masks (0/1)
//
// Time indexing here: we use 0-based C++ indexing internally; the caller
// (R wrapper) is responsible for the model's 1-based time convention.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Solve a tridiagonal system of n x n matrix blocks (Thomas algorithm).
//
// Stores its arguments as length-(T-1) vectors of matrices. We solve
//   Diag[s] * X[s] - Lower * X[s-1] - Upper * X[s+1] = RHS[s]
// where Lower and Upper are scalar lambda (the RW innovation precision).
//
// Returns the solved X[s] for s = 0..T-2 corresponding to model t = 2..T.
//
// Inputs:
//   diag_mats:  length-S list of (n x n) symmetric PD matrices
//               diag_mats[s] = M_ti D_ti M_ti^T + (mu + lambda * d_s) * I
//   rhs_vecs:   length-S vectors of length n (each is a row of RHS)
//   lambda:     RW penalty strength
//   S:          number of time blocks (T - 1)
//   n:          block dimension (n_row for A, n_col for B)
// Output:
//   length-S list of n-vectors (solved row trajectory)
//
// This routine is called once per row of A_t (or per row of B_t). The
// block-tridiagonal structure means we Cholesky-factor in a forward pass,
// then back-substitute. Each block is dimension n.
//
// For numerical stability we use Cholesky factorization (the diagonal
// blocks are symmetric PD: data-fit Gram + ridge I).
static void thomas_tridiag_solve(const std::vector<arma::mat>& diag_mats,
                                  const std::vector<arma::vec>& rhs_vecs,
                                  double lambda,
                                  int S,
                                  int n,
                                  std::vector<arma::vec>& out_x) {
    // Forward elimination: reduce to upper-triangular.
    //   D'_0 = D_0
    //   c'_0 = rhs_0
    //   D'_s = D_s - lambda^2 * (D'_{s-1})^{-1}     (s >= 1)
    //   c'_s = rhs_s + lambda * (D'_{s-1})^{-1} c'_{s-1}
    // Backward substitution:
    //   x_{S-1} = (D'_{S-1})^{-1} c'_{S-1}
    //   x_s = (D'_s)^{-1} (c'_s + lambda * x_{s+1})
    //
    // We store Cholesky factors L_s such that D'_s = L_s L_s^T.

    std::vector<arma::mat> Lchol(S);
    std::vector<arma::vec> cprime(S);

    // s = 0
    Lchol[0] = arma::chol(diag_mats[0], "lower");
    cprime[0] = rhs_vecs[0];

    for (int s = 1; s < S; ++s) {
        // Compute D'_s = D_s - lambda^2 * (D'_{s-1})^{-1}
        // We need (D'_{s-1})^{-1} as a dense matrix to subtract from D_s.
        // Use the Cholesky to solve I -> M = L L^T \ I
        arma::mat Mprev_inv = arma::solve(arma::trimatl(Lchol[s - 1]),
                                          arma::eye(n, n),
                                          arma::solve_opts::fast);
        Mprev_inv = arma::solve(arma::trimatu(Lchol[s - 1].t()), Mprev_inv,
                                arma::solve_opts::fast);
        arma::mat Dprime_s = diag_mats[s] - (lambda * lambda) * Mprev_inv;
        // Re-symmetrize for numerical safety
        Dprime_s = 0.5 * (Dprime_s + Dprime_s.t());
        Lchol[s] = arma::chol(Dprime_s, "lower");

        // Compute c'_s = rhs_s + lambda * (D'_{s-1})^{-1} c'_{s-1}
        //   = rhs_s + lambda * solve(L_{s-1} L_{s-1}^T, c'_{s-1})
        arma::vec y = arma::solve(arma::trimatl(Lchol[s - 1]), cprime[s - 1],
                                  arma::solve_opts::fast);
        arma::vec Minv_c = arma::solve(arma::trimatu(Lchol[s - 1].t()), y,
                                       arma::solve_opts::fast);
        cprime[s] = rhs_vecs[s] + lambda * Minv_c;
    }

    // Backward substitution
    out_x.resize(S);

    arma::vec y_last = arma::solve(arma::trimatl(Lchol[S - 1]), cprime[S - 1],
                                   arma::solve_opts::fast);
    out_x[S - 1] = arma::solve(arma::trimatu(Lchol[S - 1].t()), y_last,
                               arma::solve_opts::fast);

    for (int s = S - 2; s >= 0; --s) {
        arma::vec rhs_s = cprime[s] + lambda * out_x[s + 1];
        arma::vec y = arma::solve(arma::trimatl(Lchol[s]), rhs_s,
                                  arma::solve_opts::fast);
        out_x[s] = arma::solve(arma::trimatu(Lchol[s].t()), y,
                               arma::solve_opts::fast);
    }
}

//' Trajectory-wise A-update for directed Gaussian RW-smoothed ALS
//'
//' Solves the block-tridiagonal system for the entire A trajectory
//' conditional on the B trajectory.
//'
//' @param Phi cube [n_row x n_col x Tt] of centered states
//' @param B_arr cube [n_col x n_col x Tt] of B operators (slices 1..Tt-1 used)
//' @param Omega cube [n_row x n_col x Tt] of 0/1 response masks
//' @param lambda RW penalty
//' @param mu level penalty
//' @return cube [n_row x n_row x Tt] of updated A operators
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube solve_A_trajectory_cpp(const arma::cube& Phi,
                                   const arma::cube& B_arr,
                                   const arma::cube& Omega,
                                   double lambda,
                                   double mu) {
    int n_row = Phi.n_rows;
    int n_col = Phi.n_cols;
    int Tt = Phi.n_slices;
    int S = Tt - 1;  // number of operator time points (t=2..T)
    if (S < 1) {
        return arma::cube(n_row, n_row, Tt, arma::fill::zeros);
    }

    // Build M_s = Phi_{s-1+1} B_s^T for s = 0..S-1  (s indexes t=2..T)
    // i.e. for model t = 2..T, M_t = Phi_{t-1} B_t^T
    // In 0-based: time index s = 0..S-1 corresponds to model t = s+2,
    // and M[s] = Phi.slice(s) * B_arr.slice(s+1).t()
    // (Phi.slice(s) is Phi_{t-1}, B_arr.slice(s+1) is B_t)
    std::vector<arma::mat> M_list(S);
    for (int s = 0; s < S; ++s) {
        M_list[s] = Phi.slice(s) * B_arr.slice(s + 1).t();
    }

    // Right-hand response: Y_s = Phi.slice(s+1) (i.e. Phi_t for t = s+2)
    // Mask: Omega.slice(s+1)
    // boundary multiplicity d_s: 1 at s=0 and s=S-1, 2 interior (when S >= 3)
    // (translating to model: t=2 and t=T are boundaries.)
    // When S == 1, d_0 = 0 (no RW penalty contribution at all)
    // When S == 2, d_0 = 1, d_1 = 1 (both are boundaries, each has one neighbor)

    arma::cube A_out(n_row, n_row, Tt, arma::fill::zeros);
    A_out.slice(0) = arma::eye(n_row, n_row);  // unchanged initial slice

    if (S == 1) {
        // No RW penalty; just a ridge least-squares solve per row.
        const arma::mat& M0 = M_list[0];
        const arma::mat& Y0 = Phi.slice(1);
        const arma::mat& Omg0 = Omega.slice(1);
        arma::mat A_t = arma::zeros<arma::mat>(n_row, n_row);
        for (int i = 0; i < n_row; ++i) {
            arma::rowvec omega_i = Omg0.row(i);
            arma::mat D = arma::diagmat(omega_i);
            arma::mat Gram = M0 * D * M0.t() + mu * arma::eye(n_row, n_row);
            arma::vec rhs = M0 * (omega_i.t() % Y0.row(i).t());
            arma::vec a_i = arma::solve(Gram, rhs);
            A_t.row(i) = a_i.t();
        }
        A_out.slice(1) = A_t;
        return A_out;
    }

    // For S >= 2: build the tridiagonal systems per row.
    // d_s is 1 at boundaries (s=0, s=S-1) and 2 in interior
    arma::vec d_vec(S);
    d_vec.fill(2.0);
    d_vec(0) = 1.0;
    d_vec(S - 1) = 1.0;

    // Build the row-i system, solve, store in A_out.
    // For row i, we need:
    //   diag_s = M_s * D_{s,i} * M_s^T + (mu + lambda * d_s) * I
    //   rhs_s  = M_s * (omega_{s,i} .* y_{s,i})
    // where omega_{s,i} = Omega.slice(s+1).row(i)  (as col vec)
    //       y_{s,i}     = Phi.slice(s+1).row(i)    (as col vec)
    //       D_{s,i}     = diag(omega_{s,i})

    std::vector<arma::mat> diag_mats(S);
    std::vector<arma::vec> rhs_vecs(S);
    std::vector<arma::vec> sol(S);

    arma::mat I_n = arma::eye(n_row, n_row);

    for (int i = 0; i < n_row; ++i) {
        for (int s = 0; s < S; ++s) {
            const arma::mat& M_s = M_list[s];
            const arma::mat& Omg_s = Omega.slice(s + 1);
            const arma::mat& Y_s = Phi.slice(s + 1);
            arma::rowvec omega_i = Omg_s.row(i);
            arma::mat D = arma::diagmat(omega_i);
            arma::mat Gram = M_s * D * M_s.t() + (mu + lambda * d_vec(s)) * I_n;
            // Symmetrize for numerical robustness
            Gram = 0.5 * (Gram + Gram.t());
            arma::vec rhs = M_s * (omega_i.t() % Y_s.row(i).t());
            diag_mats[s] = Gram;
            rhs_vecs[s] = rhs;
        }
        thomas_tridiag_solve(diag_mats, rhs_vecs, lambda, S, n_row, sol);
        for (int s = 0; s < S; ++s) {
            A_out.slice(s + 1).row(i) = sol[s].t();
        }
    }

    return A_out;
}

//' Trajectory-wise B-update for directed Gaussian RW-smoothed ALS
//'
//' Analogous to solve_A_trajectory_cpp, with N_t = X_t^T A_t^T as the
//' design matrix. Note the corrected matrix order per reviewer 2.
//'
//' @param Phi cube [n_row x n_col x Tt] of centered states
//' @param A_arr cube [n_row x n_row x Tt] of A operators (slices 1..Tt-1 used)
//' @param Omega cube [n_row x n_col x Tt] of 0/1 response masks
//' @param lambda RW penalty
//' @param mu level penalty
//' @return cube [n_col x n_col x Tt] of updated B operators
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube solve_B_trajectory_cpp(const arma::cube& Phi,
                                   const arma::cube& A_arr,
                                   const arma::cube& Omega,
                                   double lambda,
                                   double mu) {
    int n_row = Phi.n_rows;
    int n_col = Phi.n_cols;
    int Tt = Phi.n_slices;
    int S = Tt - 1;
    if (S < 1) {
        return arma::cube(n_col, n_col, Tt, arma::fill::zeros);
    }

    // N_t = X_t^T A_t^T where X_t = Phi_{t-1}.
    // In 0-based: N[s] = Phi.slice(s).t() * A_arr.slice(s+1).t()
    // shape: [n_col x n_row]
    std::vector<arma::mat> N_list(S);
    for (int s = 0; s < S; ++s) {
        N_list[s] = Phi.slice(s).t() * A_arr.slice(s + 1).t();
    }

    arma::cube B_out(n_col, n_col, Tt, arma::fill::zeros);
    B_out.slice(0) = arma::eye(n_col, n_col);

    if (S == 1) {
        const arma::mat& N0 = N_list[0];
        const arma::mat& Y0 = Phi.slice(1);            // [n_row x n_col]
        const arma::mat& Omg0 = Omega.slice(1);        // [n_row x n_col]
        arma::mat B_t = arma::zeros<arma::mat>(n_col, n_col);
        for (int j = 0; j < n_col; ++j) {
            // Row j of B operates on Y^T (column j of Y): we need column j
            arma::colvec omega_j = Omg0.col(j);       // [n_row]
            arma::mat D = arma::diagmat(omega_j);
            arma::mat Gram = N0 * D * N0.t() + mu * arma::eye(n_col, n_col);
            arma::vec rhs = N0 * (omega_j % Y0.col(j));
            arma::vec b_j = arma::solve(Gram, rhs);
            B_t.row(j) = b_j.t();
        }
        B_out.slice(1) = B_t;
        return B_out;
    }

    arma::vec d_vec(S);
    d_vec.fill(2.0);
    d_vec(0) = 1.0;
    d_vec(S - 1) = 1.0;

    std::vector<arma::mat> diag_mats(S);
    std::vector<arma::vec> rhs_vecs(S);
    std::vector<arma::vec> sol(S);

    arma::mat I_n = arma::eye(n_col, n_col);

    for (int j = 0; j < n_col; ++j) {
        for (int s = 0; s < S; ++s) {
            const arma::mat& N_s = N_list[s];
            const arma::mat& Omg_s = Omega.slice(s + 1);
            const arma::mat& Y_s = Phi.slice(s + 1);
            arma::colvec omega_j = Omg_s.col(j);
            arma::mat D = arma::diagmat(omega_j);
            arma::mat Gram = N_s * D * N_s.t() + (mu + lambda * d_vec(s)) * I_n;
            Gram = 0.5 * (Gram + Gram.t());
            arma::vec rhs = N_s * (omega_j % Y_s.col(j));
            diag_mats[s] = Gram;
            rhs_vecs[s] = rhs;
        }
        thomas_tridiag_solve(diag_mats, rhs_vecs, lambda, S, n_col, sol);
        for (int s = 0; s < S; ++s) {
            B_out.slice(s + 1).row(j) = sol[s].t();
        }
    }

    return B_out;
}

//' Compute the RW-smoothed ALS objective value
//'
//' @param Phi cube [n_row x n_col x Tt]
//' @param A_arr cube [n_row x n_row x Tt]
//' @param B_arr cube [n_col x n_col x Tt]
//' @param Omega cube [n_row x n_col x Tt]
//' @param lambda_A RW penalty on A
//' @param lambda_B RW penalty on B
//' @param mu_A level penalty on A
//' @param mu_B level penalty on B
//' @return scalar objective value
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_tv_als_objective_cpp(const arma::cube& Phi,
                                     const arma::cube& A_arr,
                                     const arma::cube& B_arr,
                                     const arma::cube& Omega,
                                     double lambda_A,
                                     double lambda_B,
                                     double mu_A,
                                     double mu_B) {
    int Tt = Phi.n_slices;
    int S = Tt - 1;
    double obj = 0.0;

    // Data fit: sum over t=2..T of 0.5 * ||Omega_t .* (Phi_t - A_t Phi_{t-1} B_t^T)||_F^2
    for (int s = 0; s < S; ++s) {
        arma::mat pred = A_arr.slice(s + 1) * Phi.slice(s) * B_arr.slice(s + 1).t();
        arma::mat resid = Omega.slice(s + 1) % (Phi.slice(s + 1) - pred);
        obj += 0.5 * arma::accu(resid % resid);
    }

    // RW penalty: 0.5 * lambda * sum over t=3..T of ||A_t - A_{t-1}||_F^2
    if (S >= 2) {
        for (int s = 1; s < S; ++s) {
            arma::mat diff_A = A_arr.slice(s + 1) - A_arr.slice(s);
            obj += 0.5 * lambda_A * arma::accu(diff_A % diff_A);
            arma::mat diff_B = B_arr.slice(s + 1) - B_arr.slice(s);
            obj += 0.5 * lambda_B * arma::accu(diff_B % diff_B);
        }
    }

    // Level penalty: 0.5 * mu * sum over t=2..T of ||A_t||_F^2
    for (int s = 0; s < S; ++s) {
        obj += 0.5 * mu_A * arma::accu(A_arr.slice(s + 1) % A_arr.slice(s + 1));
        obj += 0.5 * mu_B * arma::accu(B_arr.slice(s + 1) % B_arr.slice(s + 1));
    }

    return obj;
}
