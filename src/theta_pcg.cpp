// matrix-free preconditioned conjugate gradient sampler for the
// dyad-vector latent state Theta_{1:T} under the corrected
// single-operator model. avoids forming the m x m transition matrix
// G_t = G_t(A_t) where m = n(n-1)/2; instead applies G_t implicitly
// via embed_0 -> A_t X A_t' -> vech_0, costing O(n^3) per G action
// instead of O(m^2) = O(n^4) for the matrix form.
//
// the dyad state x_t = vech_0(Theta_t) has block-tridiagonal precision
//     Q[t, t]   = (1/sigma_obs^2) I + (1/sigma^2) I + 1[t<T] (1/sigma^2) G_{t+1}' G_{t+1}
//                                                   + 1[t==1] (1/sigma_init^2) I
//     Q[t, t-1] = -(1/sigma^2) G_t
// (with symmetric off-diagonal blocks). the right-hand side is
//     b[t] = (y_t - mu) / sigma_obs^2 + noise contributions for sampling.
// the sampler uses perturb-and-solve: draw noise z ~ N(0, Q) by adding
// one independent gaussian per factor (init prior, T-1 process factors,
// T observation factors), then solve Q x = b + z by CG.

#include <RcppArmadillo.h>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// embed an m-vector of upper-triangle dyads into an n x n symmetric
// matrix with zero diagonal. pairs is m x 2 with 0-based (i, j),
// i < j.
static inline arma::mat embed_0_cpp(const arma::vec& x,
	int n, const arma::umat& pairs) {
	int m = pairs.n_rows;
	arma::mat M(n, n, arma::fill::zeros);
	for (int d = 0; d < m; d++) {
		arma::uword i = pairs(d, 0);
		arma::uword j = pairs(d, 1);
		M(i, j) = x(d);
		M(j, i) = x(d);
	}
	return M;
}

// inverse: extract upper-triangle dyads from a symmetric matrix
static inline arma::vec vech_0_cpp(const arma::mat& M,
	const arma::umat& pairs) {
	int m = pairs.n_rows;
	arma::vec x(m);
	for (int d = 0; d < m; d++) {
		x(d) = M(pairs(d, 0), pairs(d, 1));
	}
	return x;
}

// G_t(x) = vech_0(A_t * embed_0(x) * A_t').   cost: O(n^3).
static inline arma::vec apply_G_cpp(const arma::vec& x,
	const arma::mat& A_t, int n, const arma::umat& pairs) {
	arma::mat X = embed_0_cpp(x, n, pairs);
	arma::mat Y = A_t * X * A_t.t();
	return vech_0_cpp(Y, pairs);
}

// G_t^*(z) = vech_0(A_t' * embed_0(z) * A_t).  cost: O(n^3).
static inline arma::vec apply_Gt_cpp(const arma::vec& z,
	const arma::mat& A_t, int n, const arma::umat& pairs) {
	arma::mat Z = embed_0_cpp(z, n, pairs);
	arma::mat Y = A_t.t() * Z * A_t;
	return vech_0_cpp(Y, pairs);
}

// precision-vector multiply Q v for the block-tridiagonal path
// precision. v is m x T; output is m x T. inputs:
//   A_arr: n x n x T cube of operators
//   sigma2: process variance
//   sigma2_obs: observation variance
//   sigma2_init: variance of the (weak) prior on x_1
//   pairs: dyad index table
static arma::mat Q_apply_cpp(const arma::mat& v,
	const arma::cube& A_arr,
	double sigma2, double sigma2_obs, double sigma2_init,
	int n, const arma::umat& pairs) {
	int m = v.n_rows;
	int Tt = v.n_cols;
	arma::mat out(m, Tt, arma::fill::zeros);
	double inv_s2 = 1.0 / sigma2;
	double inv_s2_obs = 1.0 / sigma2_obs;
	double inv_s2_init = 1.0 / sigma2_init;
	out += inv_s2_obs * v;
	out.col(0) += inv_s2_init * v.col(0);
	for (int s = 1; s < Tt; s++) {
		arma::vec v_prev = v.col(s - 1);
		arma::vec G_v_prev = apply_G_cpp(v_prev,
			A_arr.slice(s), n, pairs);
		out.col(s) += inv_s2 * v.col(s);
		out.col(s) -= inv_s2 * G_v_prev;
		arma::vec Gt_term = apply_Gt_cpp(G_v_prev - v.col(s),
			A_arr.slice(s), n, pairs);
		out.col(s - 1) += inv_s2 * Gt_term;
	}
	return out;
}

// build the RHS b + noise for perturb-and-solve. when sample_noise is
// true, the result is b' = b + z with z ~ N(0, Q) so solving Q x = b'
// gives x ~ N(Q^{-1} b, Q^{-1}). when false, returns b alone (used to
// get the posterior mean for sanity checks).
static arma::mat build_rhs_cpp(const arma::cube& Y_arr_3d,
	const arma::mat& M_mat,
	const arma::cube& A_arr,
	double sigma2, double sigma2_obs, double sigma2_init,
	int n, const arma::umat& pairs, bool sample_noise) {
	int m = pairs.n_rows;
	int Tt = A_arr.n_slices;
	arma::vec mu_d = vech_0_cpp(M_mat, pairs);
	arma::mat b(m, Tt);
	for (int t = 0; t < Tt; t++) {
		arma::vec y_d = vech_0_cpp(Y_arr_3d.slice(t), pairs);
		b.col(t) = (y_d - mu_d) / sigma2_obs;
	}
	if (sample_noise) {
		// init noise on x_1 prior
		arma::vec u_init = arma::randn(m);
		b.col(0) += u_init / std::sqrt(sigma2_init);
		// process noise at each t = 1..Tt-1 (0-based)
		for (int t = 1; t < Tt; t++) {
			arma::vec eps = arma::randn(m);
			b.col(t)     += eps / std::sqrt(sigma2);
			b.col(t - 1) -= apply_Gt_cpp(eps, A_arr.slice(t),
				n, pairs) / std::sqrt(sigma2);
		}
		// obs noise at each t = 0..Tt-1
		for (int t = 0; t < Tt; t++) {
			arma::vec eta = arma::randn(m);
			b.col(t) += eta / std::sqrt(sigma2_obs);
		}
	}
	return b;
}

// per-dyad diagonal D_t of the operator G_t, computed directly from A_t.
// for dyad d = (i, j),
//   [G_t]_{dd} = A_t(i,i) A_t(j,j) + A_t(i,j) A_t(j,i),
// the exact diagonal of the symmetric-congruence operator. replacing G_t
// by diag(D_t) (rather than by a single shared scalar) is the
// Frobenius-optimal *diagonal* surrogate and is exact whenever A_t is
// itself diagonal. for conflict networks the per-actor persistence
// A_t(i,i) varies widely, so [G_t]_{dd} varies strongly across dyads and
// the per-dyad diagonal is a far better surrogate than a shared scalar.
// returned as an m x Tt matrix, one row per dyad.
static arma::mat compute_D_diag_cpp(const arma::cube& A_arr,
	const arma::umat& pairs) {
	int Tt = A_arr.n_slices;
	int m = pairs.n_rows;
	arma::mat D(m, Tt, arma::fill::zeros);
	for (int t = 0; t < Tt; t++) {
		const arma::mat& A = A_arr.slice(t);
		for (int d = 0; d < m; d++) {
			arma::uword i = pairs(d, 0);
			arma::uword j = pairs(d, 1);
			D(d, t) = A(i, i) * A(j, j) + A(i, j) * A(j, i);
		}
	}
	return D;
}

// per-dyad block-tridiagonal-in-time preconditioner: the path precision Q
// with every operator G_t replaced by the diagonal matrix diag(D_t).
// because diag(D_t) acts on each dyad coordinate independently, the
// resulting precision M is block-diagonal across the m dyads, each block a
// Tt x Tt symmetric tridiagonal. M^{-1} is therefore m independent
// Tt-point tridiagonal solves, each with its own LDL^T factorization:
// O(m Tt) total, negligible against the O(n^3 Tt) Q-apply. M is a valid
// Gaussian precision (SPD), hence a valid CG preconditioner; preconditioning
// changes only the CG iteration count, never the converged solution. an
// earlier surrogate used a single shared scalar g_t in place of diag(D_t);
// that washed out the wide cross-dyad spread of [G_t]_{dd} and left CG
// essentially unpreconditioned once A_t departed from a scalar multiple of
// the identity.
struct TimeTridiagPrecond {
	arma::mat piv;    // m x Tt LDL^T diagonal, one row per dyad
	arma::mat Lsub;   // m x Tt LDL^T sub-diagonal; Lsub.col(t) couples t, t-1

	void build(const arma::mat& D, double sigma2, double sigma2_obs,
		double sigma2_init) {
		int m = D.n_rows;
		int Tt = D.n_cols;
		double inv_s2 = 1.0 / sigma2;
		double inv_s2_obs = 1.0 / sigma2_obs;
		double inv_s2_init = 1.0 / sigma2_init;
		arma::mat d(m, Tt), e(m, Tt, arma::fill::zeros);
		for (int t = 0; t < Tt; t++) {
			arma::vec diag(m);
			diag.fill(inv_s2_obs);
			if (t == 0) diag += inv_s2_init;
			if (t >= 1) diag += inv_s2;
			if (t < Tt - 1) diag += inv_s2 * arma::square(D.col(t + 1));
			d.col(t) = diag;
		}
		for (int t = 1; t < Tt; t++) e.col(t) = -inv_s2 * D.col(t);
		piv.set_size(m, Tt);
		Lsub.set_size(m, Tt);
		Lsub.zeros();
		piv.col(0) = d.col(0);
		for (int t = 1; t < Tt; t++) {
			Lsub.col(t) = e.col(t) / piv.col(t - 1);
			piv.col(t) = d.col(t) - Lsub.col(t) % e.col(t);
			piv.col(t) = arma::clamp(piv.col(t), 1e-300, arma::datum::inf);
		}
	}

	// apply M^{-1} to an m x Tt residual (columns indexed by time)
	arma::mat apply(const arma::mat& r) const {
		int Tt = r.n_cols;
		arma::mat y(r.n_rows, Tt);
		y.col(0) = r.col(0);
		for (int t = 1; t < Tt; t++) {
			y.col(t) = r.col(t) - Lsub.col(t) % y.col(t - 1);
		}
		arma::mat v(r.n_rows, Tt);
		v.col(Tt - 1) = y.col(Tt - 1) / piv.col(Tt - 1);
		for (int t = Tt - 2; t >= 0; t--) {
			v.col(t) = y.col(t) / piv.col(t)
				- Lsub.col(t + 1) % v.col(t + 1);
		}
		return v;
	}
};

// simpler block-diagonal preconditioner: ignores time-coupling in M,
// using only process/prior/observation variance structure. cheaper to
// apply (just per-dyad-per-time division) but may help conditioning by
// focusing on dominant variance structure rather than trying to capture
// the operator G_t's off-diagonal structure.
struct BlockDiagonalPrecond {
	arma::mat inv_diag;  // m x Tt: reciprocal of diagonal, one per dyad-time

	void build(const arma::mat& D, double sigma2, double sigma2_obs,
		double sigma2_init) {
		int m = D.n_rows;
		int Tt = D.n_cols;
		double inv_s2 = 1.0 / sigma2;
		double inv_s2_obs = 1.0 / sigma2_obs;
		double inv_s2_init = 1.0 / sigma2_init;
		inv_diag.set_size(m, Tt);

		for (int t = 0; t < Tt; t++) {
			arma::vec diag_val(m);
			diag_val.fill(inv_s2_obs);
			if (t == 0) diag_val += inv_s2_init;
			if (t >= 1) diag_val += inv_s2;
			// omit G_t coupling term; focus on variance structure only
			inv_diag.col(t) = 1.0 / diag_val;
		}
	}

	arma::mat apply(const arma::mat& r) const {
		// just per-element multiplication: M^{-1} is diagonal
		return r % inv_diag;
	}
};

// preconditioned CG with time-tridiagonal preconditioner
static arma::mat cg_solve_tridiag(const arma::mat& b,
	const arma::cube& A_arr,
	double sigma2, double sigma2_obs, double sigma2_init,
	int n, const arma::umat& pairs,
	double tol = 1e-8, int max_iter = 200, int* iters = nullptr) {
	int m = b.n_rows;
	int Tt = b.n_cols;
	TimeTridiagPrecond precond;
	precond.build(compute_D_diag_cpp(A_arr, pairs),
		sigma2, sigma2_obs, sigma2_init);
	arma::mat x(m, Tt, arma::fill::zeros);
	arma::mat r = b;   // r = b - Q x  with x = 0
	arma::mat z = precond.apply(r);
	arma::mat p = z;
	double rz_old = arma::accu(r % z);
	double b_norm = std::sqrt(arma::accu(b % b));
	if (b_norm < 1e-15) b_norm = 1.0;
	int it = 0;
	for (it = 0; it < max_iter; it++) {
		arma::mat Qp = Q_apply_cpp(p, A_arr, sigma2, sigma2_obs,
			sigma2_init, n, pairs);
		double alpha = rz_old / arma::accu(p % Qp);
		x += alpha * p;
		r -= alpha * Qp;
		double r_norm = std::sqrt(arma::accu(r % r));
		if (r_norm / b_norm < tol) { it++; break; }
		z = precond.apply(r);
		double rz_new = arma::accu(r % z);
		double beta = rz_new / rz_old;
		p = z + beta * p;
		rz_old = rz_new;
	}
	if (iters != nullptr) *iters = it;
	return x;
}

// preconditioned CG with block-diagonal preconditioner (process/prior/obs structure only)
static arma::mat cg_solve_blockdiag(const arma::mat& b,
	const arma::cube& A_arr,
	double sigma2, double sigma2_obs, double sigma2_init,
	int n, const arma::umat& pairs,
	double tol = 1e-8, int max_iter = 200, int* iters = nullptr) {
	int m = b.n_rows;
	int Tt = b.n_cols;
	BlockDiagonalPrecond precond;
	precond.build(compute_D_diag_cpp(A_arr, pairs),
		sigma2, sigma2_obs, sigma2_init);
	arma::mat x(m, Tt, arma::fill::zeros);
	arma::mat r = b;   // r = b - Q x  with x = 0
	arma::mat z = precond.apply(r);
	arma::mat p = z;
	double rz_old = arma::accu(r % z);
	double b_norm = std::sqrt(arma::accu(b % b));
	if (b_norm < 1e-15) b_norm = 1.0;
	int it = 0;
	for (it = 0; it < max_iter; it++) {
		arma::mat Qp = Q_apply_cpp(p, A_arr, sigma2, sigma2_obs,
			sigma2_init, n, pairs);
		double alpha = rz_old / arma::accu(p % Qp);
		x += alpha * p;
		r -= alpha * Qp;
		double r_norm = std::sqrt(arma::accu(r % r));
		if (r_norm / b_norm < tol) { it++; break; }
		z = precond.apply(r);
		double rz_new = arma::accu(r % z);
		double beta = rz_new / rz_old;
		p = z + beta * p;
		rz_old = rz_new;
	}
	if (iters != nullptr) *iters = it;
	return x;
}

// main entry: preconditioned CG for the block-tridiagonal path precision. tol is a
// relative residual stop; the converged solution is independent of the
// preconditioner. iters, if non-null, receives the iteration count.
// uses the block-diagonal preconditioner (BlockDiagonalPrecond), which focuses on
// process/prior/observation variance structure and omits off-diagonal G_t coupling.
// empirically, this simpler preconditioner is more effective than attempting to
// capture the rough diagonal approximation of G_t.
static arma::mat cg_solve_cpp(const arma::mat& b,
	const arma::cube& A_arr,
	double sigma2, double sigma2_obs, double sigma2_init,
	int n, const arma::umat& pairs,
	double tol = 1e-8, int max_iter = 200, int* iters = nullptr) {
	return cg_solve_blockdiag(b, A_arr, sigma2, sigma2_obs, sigma2_init,
		n, pairs, tol, max_iter, iters);
}

// diagnostic: benchmark CG with time-tridiagonal preconditioner
// [[Rcpp::export]]
Rcpp::List cg_benchmark_tridiag_cpp(const arma::cube& Y_arr_3d,
	const arma::mat& M_mat,
	const arma::cube& A_arr,
	double sigma2, double sigma2_obs,
	const arma::umat& pairs) {
	int n = A_arr.n_rows;
	int Tt = A_arr.n_slices;
	double sigma2_init = sigma2 * 10.0;
	arma::mat b = build_rhs_cpp(Y_arr_3d, M_mat, A_arr,
		sigma2, sigma2_obs, sigma2_init,
		n, pairs, /*sample_noise=*/ false);  // no noise for deterministic comparison

	auto start = std::chrono::high_resolution_clock::now();
	int cg_iters = 0;
	arma::mat x = cg_solve_tridiag(b, A_arr, sigma2, sigma2_obs,
		sigma2_init, n, pairs, 1e-8, 200, &cg_iters);
	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

	return Rcpp::List::create(
		Rcpp::Named("cg_iters") = cg_iters,
		Rcpp::Named("elapsed_ms") = elapsed_ms
	);
}

// diagnostic: benchmark CG with block-diagonal preconditioner
// [[Rcpp::export]]
Rcpp::List cg_benchmark_blockdiag_cpp(const arma::cube& Y_arr_3d,
	const arma::mat& M_mat,
	const arma::cube& A_arr,
	double sigma2, double sigma2_obs,
	const arma::umat& pairs) {
	int n = A_arr.n_rows;
	int Tt = A_arr.n_slices;
	double sigma2_init = sigma2 * 10.0;
	arma::mat b = build_rhs_cpp(Y_arr_3d, M_mat, A_arr,
		sigma2, sigma2_obs, sigma2_init,
		n, pairs, /*sample_noise=*/ false);  // no noise for deterministic comparison

	auto start = std::chrono::high_resolution_clock::now();
	int cg_iters = 0;
	arma::mat x = cg_solve_blockdiag(b, A_arr, sigma2, sigma2_obs,
		sigma2_init, n, pairs, 1e-8, 200, &cg_iters);
	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

	return Rcpp::List::create(
		Rcpp::Named("cg_iters") = cg_iters,
		Rcpp::Named("elapsed_ms") = elapsed_ms
	);
}

// public sampler: returns a draw of Theta_{1:T} as an n x n x T cube
// from the gaussian full conditional implied by the observation model
// y_t = M + Theta_t + obs noise, the dynamics Theta_t = A_t Theta_{t-1}
// A_t' + process noise (off-diagonal entries only, zero diagonal), and
// a weakly informative prior on Theta_1.
//
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::cube theta_pcg_sampler_cpp(const arma::cube& Y_arr_3d,
	const arma::mat& M_mat,
	const arma::cube& A_arr,
	double sigma2, double sigma2_obs,
	const arma::umat& pairs) {
	int n = A_arr.n_rows;
	int Tt = A_arr.n_slices;
	int m = pairs.n_rows;
	if ((int)Y_arr_3d.n_slices != Tt) {
		stop("theta_pcg_sampler_cpp: Y_arr and A_arr must have "
			"matching time dimensions");
	}
	if ((int)Y_arr_3d.n_rows != n || (int)Y_arr_3d.n_cols != n) {
		stop("theta_pcg_sampler_cpp: Y_arr slices must be square "
			"and match A dimensions");
	}
	double sigma2_init = sigma2 * 10.0;
	arma::mat b = build_rhs_cpp(Y_arr_3d, M_mat, A_arr,
		sigma2, sigma2_obs, sigma2_init,
		n, pairs, /*sample_noise=*/ true);
	arma::mat x = cg_solve_cpp(b, A_arr, sigma2, sigma2_obs,
		sigma2_init, n, pairs);
	arma::cube Theta_out(n, n, Tt);
	for (int t = 0; t < Tt; t++) {
		Theta_out.slice(t) = embed_0_cpp(x.col(t), n, pairs);
	}
	return Theta_out;
}

// ============================================================================
// ASYMMETRIC PCG FOR BIPARTITE/DIRECTED NETWORKS
// ============================================================================
// The asymmetric model has state x_t = vec(Theta_t) where Theta_t is n_row x n_col.
// The transition is Theta_t = A_t * Theta_{t-1} * B_t^T + noise.
// The precision is block-tridiagonal in the (n_row * n_col)-dimensional space.
// The Kronecker product structure G_t(x) = vec(A_t * mat(x) * B_t^T)
// is applied implicitly to avoid forming the huge operator matrix.

// Apply G_t to a vector x (asymmetric case).
// x is n_row*n_col, A_t is n_row x n_row, B_t is n_col x n_col.
// Returns vec(A_t * mat(x, n_row, n_col) * B_t^T).
static inline arma::vec apply_G_asym(const arma::vec& x,
	const arma::mat& A_t, const arma::mat& B_t,
	int n_row, int n_col) {
	arma::mat X = arma::reshape(x, n_row, n_col);
	arma::mat Y = A_t * X * B_t.t();
	return arma::vectorise(Y);
}

// Apply G_t^* (adjoint) to a vector z.
// Returns vec(A_t^T * mat(z, n_row, n_col) * B_t).
static inline arma::vec apply_Gt_asym(const arma::vec& z,
	const arma::mat& A_t, const arma::mat& B_t,
	int n_row, int n_col) {
	arma::mat Z = arma::reshape(z, n_row, n_col);
	arma::mat Y = A_t.t() * Z * B_t;
	return arma::vectorise(Y);
}

// Q_apply_asym: precision-vector multiply for block-tridiagonal asymmetric path precision.
static arma::mat Q_apply_asym(const arma::mat& v,
	const arma::cube& A_arr, const arma::cube& B_arr,
	double sigma2, double sigma2_obs, double sigma2_init,
	int n_row, int n_col) {
	int m = v.n_rows;  // m = n_row * n_col
	int Tt = v.n_cols;
	arma::mat out(m, Tt, arma::fill::zeros);
	double inv_s2 = 1.0 / sigma2;
	double inv_s2_obs = 1.0 / sigma2_obs;
	double inv_s2_init = 1.0 / sigma2_init;

	out += inv_s2_obs * v;
	out.col(0) += inv_s2_init * v.col(0);

	for (int s = 1; s < Tt; s++) {
		arma::vec v_prev = v.col(s - 1);
		arma::vec G_v_prev = apply_G_asym(v_prev,
			A_arr.slice(s), B_arr.slice(s), n_row, n_col);
		out.col(s) += inv_s2 * v.col(s);
		out.col(s) -= inv_s2 * G_v_prev;
		arma::vec Gt_term = apply_Gt_asym(G_v_prev - v.col(s),
			A_arr.slice(s), B_arr.slice(s), n_row, n_col);
		out.col(s - 1) += inv_s2 * Gt_term;
	}
	return out;
}

// Build RHS b + noise for asymmetric case.
static arma::mat build_rhs_asym(const arma::cube& Y_arr_3d,
	const arma::mat& M_mat,
	const arma::cube& A_arr, const arma::cube& B_arr,
	double sigma2, double sigma2_obs, double sigma2_init,
	int n_row, int n_col, bool sample_noise) {
	int m = n_row * n_col;
	int Tt = A_arr.n_slices;
	arma::vec mu_d = arma::vectorise(M_mat);
	arma::mat b(m, Tt);

	for (int t = 0; t < Tt; t++) {
		arma::vec y_d = arma::vectorise(Y_arr_3d.slice(t));
		b.col(t) = (y_d - mu_d) / sigma2_obs;
	}

	if (sample_noise) {
		// init noise on x_1 prior
		arma::vec u_init = arma::randn(m);
		b.col(0) += u_init / std::sqrt(sigma2_init);
		// process noise at each t = 1..Tt-1 (0-based)
		for (int t = 1; t < Tt; t++) {
			arma::vec eps = arma::randn(m);
			b.col(t) += eps / std::sqrt(sigma2);
			b.col(t - 1) -= apply_Gt_asym(eps,
				A_arr.slice(t), B_arr.slice(t), n_row, n_col) / std::sqrt(sigma2);
		}
		// obs noise at each t = 0..Tt-1
		for (int t = 0; t < Tt; t++) {
			arma::vec eta = arma::randn(m);
			b.col(t) += eta / std::sqrt(sigma2_obs);
		}
	}
	return b;
}

// Block-diagonal preconditioner for asymmetric case: D(i,j,t) = A(i,i) * B(j,j).
struct BlockDiagonalPrecond_asym {
	arma::mat inv_diag;  // (n_row * n_col) x Tt

	void build(const arma::cube& A_arr, const arma::cube& B_arr,
		double sigma2, double sigma2_obs, double sigma2_init,
		int n_row, int n_col) {
		int m = n_row * n_col;
		int Tt = A_arr.n_slices;
		double inv_s2 = 1.0 / sigma2;
		double inv_s2_obs = 1.0 / sigma2_obs;
		double inv_s2_init = 1.0 / sigma2_init;
		inv_diag.set_size(m, Tt);

		for (int t = 0; t < Tt; t++) {
			const arma::mat& A_t = A_arr.slice(t);
			const arma::mat& B_t = B_arr.slice(t);
			arma::vec diag_val(m);
			int idx = 0;
			for (int i = 0; i < n_row; i++) {
				for (int j = 0; j < n_col; j++) {
					diag_val(idx) = inv_s2_obs;
					if (t == 0) diag_val(idx) += inv_s2_init;
					if (t >= 1) diag_val(idx) += inv_s2 * A_t(i, i) * B_t(j, j);
					idx++;
				}
			}
			inv_diag.col(t) = 1.0 / diag_val;
		}
	}

	arma::mat apply(const arma::mat& r) const {
		return r % inv_diag;
	}
};

// Preconditioned CG solver for asymmetric case.
static arma::mat cg_solve_asym(const arma::mat& b,
	const arma::cube& A_arr, const arma::cube& B_arr,
	double sigma2, double sigma2_obs, double sigma2_init,
	int n_row, int n_col,
	double tol = 1e-8, int max_iter = 200, int* iters = nullptr) {
	int m = b.n_rows;
	int Tt = b.n_cols;
	BlockDiagonalPrecond_asym precond;
	precond.build(A_arr, B_arr, sigma2, sigma2_obs, sigma2_init, n_row, n_col);

	arma::mat x(m, Tt, arma::fill::zeros);
	arma::mat r = b;
	arma::mat z = precond.apply(r);
	arma::mat p = z;
	double rz_old = arma::accu(r % z);
	double b_norm = std::sqrt(arma::accu(b % b));
	if (b_norm < 1e-15) b_norm = 1.0;

	int it = 0;
	for (it = 0; it < max_iter; it++) {
		arma::mat Qp = Q_apply_asym(p, A_arr, B_arr, sigma2, sigma2_obs,
			sigma2_init, n_row, n_col);
		double alpha = rz_old / arma::accu(p % Qp);
		x += alpha * p;
		r -= alpha * Qp;
		double r_norm = std::sqrt(arma::accu(r % r));
		if (r_norm / b_norm < tol) { it++; break; }
		z = precond.apply(r);
		double rz_new = arma::accu(r % z);
		double beta = rz_new / rz_old;
		p = z + beta * p;
		rz_old = rz_new;
	}
	if (iters != nullptr) *iters = it;
	return x;
}

// Public asymmetric PCG sampler: returns Theta_{1:T} as an n_row x n_col x T cube.
// [[Rcpp::export]]
arma::cube theta_pcg_asym_sampler_cpp(const arma::cube& Y_arr_3d,
	const arma::mat& M_mat,
	const arma::cube& A_arr,
	const arma::cube& B_arr,
	double sigma2, double sigma2_obs) {
	int n_row = A_arr.n_rows;
	int n_col = B_arr.n_cols;
	int Tt = A_arr.n_slices;

	if ((int)Y_arr_3d.n_slices != Tt) {
		stop("theta_pcg_asym_sampler_cpp: Y_arr and A_arr must have "
			"matching time dimensions");
	}
	if ((int)Y_arr_3d.n_rows != n_row || (int)Y_arr_3d.n_cols != n_col) {
		stop("theta_pcg_asym_sampler_cpp: Y_arr dimensions must match "
			"A_arr rows and B_arr columns");
	}
	if ((int)B_arr.n_slices != Tt) {
		stop("theta_pcg_asym_sampler_cpp: B_arr must have Tt time slices");
	}

	double sigma2_init = sigma2 * 10.0;
	arma::mat b = build_rhs_asym(Y_arr_3d, M_mat, A_arr, B_arr,
		sigma2, sigma2_obs, sigma2_init,
		n_row, n_col, /*sample_noise=*/ true);
	arma::mat x = cg_solve_asym(b, A_arr, B_arr, sigma2, sigma2_obs,
		sigma2_init, n_row, n_col);

	arma::cube Theta_out(n_row, n_col, Tt);
	for (int t = 0; t < Tt; t++) {
		Theta_out.slice(t) = arma::reshape(x.col(t), n_row, n_col);
	}
	return Theta_out;
}
