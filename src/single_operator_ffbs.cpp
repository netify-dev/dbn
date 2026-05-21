// row-trajectory FFBS for the undirected single-operator model
// (B = A, arbitrary A, dyad-level likelihood with structural-no-loop
// diagonal). under B = A and the unipartite undirected convention,
// the row-r conditional on a_{r, 1:T} given the other rows of A and
// the latent state Theta is linear-Gaussian, so exact Gibbs draws
// are available via standard FFBS without any M-H step.
//
// per row r, time t, the dyad-level likelihood includes one
// observation per dyad incident to r:
//     y_d = Theta_{i, j, t}   for d = (i, j) with i < j and r in {i, j}
// the design row at dyad d is
//     h_d = C_t * a_{partner, t}'    where partner is the other actor
//                                    in dyad d and C_t = Theta_{t-1}
// the orientation-specific formulation keeps both branches of the
// dyad design explicit even when Theta is symmetric, so the
// likelihood remains correct under tiny numerical asymmetry.
//
// random-walk prior on a_{r, 1:T}: a_{r, t} = a_{r, t-1} + omega,
// omega entries iid N(0, tauA2).
//
// signature mirrors update_A_batch (src/update_AB_batch.cpp) but
// operates row-by-row using the latest A_t for the design (sequential
// Gauss-Seidel scan), which is exact Gibbs for the B = A target.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// safe symmetric PD inverse with regularization fallback
static arma::mat safe_inv_sympd_so(const arma::mat& M) {
	arma::mat M_sym = 0.5 * (M + M.t());
	arma::mat result;
	bool ok = arma::inv_sympd(result, M_sym);
	if (!ok) {
		double reg = 1e-6 * arma::norm(M_sym, "fro") + 1e-8;
		M_sym.diag() += reg;
		ok = arma::inv_sympd(result, M_sym);
		if (!ok) {
			result = arma::inv(M_sym);
		}
	}
	return result;
}

// safe multivariate normal draw with regularization
static arma::vec safe_mvnrnd_so(const arma::vec& mu, const arma::mat& Sigma) {
	arma::mat S = 0.5 * (Sigma + Sigma.t());
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

// row-conditional math used by the suff-stat row-FFBS below. given
// row r of A and the other rows fixed, the row conditional likelihood
// for a_r at time t is
//    sum_{j != r} (Theta_{r, j, t} - h_{j, t}' a_r)^2 / sigma^2
// with h_{j, t} = (Theta_{t-1} a_j')_{n x 1}. let
//    W_t = Theta_{t-1} A_t'       (n x n: column j is Theta_{t-1} a_j')
//    S_t = W_t W_t'                (n x n: sum over all j of h_j h_j')
//    B_t = W_t Theta_t            (n x n: column j is sum_k W_t[k,] Theta_{kj,t})
// then for row r,
//    H^T H = S_t - w_r w_r'        (sum over j != r)
//    H^T y = B_t[, r]
// (the j = r contribution to H^T y is zero because Theta has
// structural zero diagonal).

// helper: precompute W, S, B from current A and Theta.
static void compute_AB_suff(const arma::cube& A_arr,
	const arma::cube& Theta_arr,
	arma::cube& W, arma::cube& S, arma::cube& B) {
	int n  = (int)A_arr.n_rows;
	int Tt = (int)A_arr.n_slices;
	W.set_size(n, n, Tt); S.set_size(n, n, Tt); B.set_size(n, n, Tt);
	W.zeros(); S.zeros(); B.zeros();
	for (int t = 1; t < Tt; t++) {
		W.slice(t) = Theta_arr.slice(t - 1) * A_arr.slice(t).t();
		S.slice(t) = W.slice(t) * W.slice(t).t();
		B.slice(t) = W.slice(t) * Theta_arr.slice(t);
	}
}

// incremental update: row r changed in A_arr.slice(t). updates W, S, B
// to reflect the new row. only column r of W changes; S and B get
// rank-one updates plus a correction from the Theta_t row.
static void update_AB_suff_after_row(arma::cube& W, arma::cube& S,
	arma::cube& B,
	const arma::cube& A_arr, const arma::cube& Theta_arr,
	int r, int t) {
	arma::vec new_col = Theta_arr.slice(t - 1) * A_arr.slice(t).row(r).t();
	arma::vec old_col = W.slice(t).col(r);
	W.slice(t).col(r) = new_col;
	S.slice(t) = S.slice(t) - (old_col * old_col.t()) +
		(new_col * new_col.t());
	arma::rowvec theta_r = Theta_arr.slice(t).row(r);
	B.slice(t) = B.slice(t) + (new_col - old_col) * theta_r;
}

// one row's trajectory FFBS using suff stats. modifies Aarray and
// the W/S/B cubes in place. uses the latest other-rows of A for the
// design (sequential Gauss-Seidel sweep).
//
// tau2_init_r is the prior variance on a_{r, 1}. allowing this to
// vary per row (e.g., from a horseshoe-style hierarchy on row L2
// norms) breaks the row-collapse degeneracy that arises when all
// rows have the same A_1 prior scale.
static void row_ffbs_suffstat(int r,
	arma::cube& Aarray,
	const arma::cube& Theta_arr,
	arma::cube& W, arma::cube& S, arma::cube& B,
	double sigma2, double tau2,
	double tau2_init_r) {
	int n  = (int)Aarray.n_rows;
	int Tt = (int)Aarray.n_slices;
	arma::mat W_proc = tau2 * arma::eye(n, n);
	// forward filter init: a_{r, 1} ~ N(0, tau2_init_r I)
	arma::vec m_cur = arma::zeros(n);
	arma::mat C_cur = tau2_init_r * arma::eye(n, n);
	arma::mat m_fwd(n, Tt, arma::fill::zeros);
	arma::cube C_fwd(n, n, Tt, arma::fill::zeros);
	m_fwd.col(0) = m_cur; C_fwd.slice(0) = C_cur;
	for (int t = 1; t < Tt; t++) {
		arma::mat C_pred = C_cur + W_proc;
		arma::mat C_pred_inv;
		bool ok = arma::inv_sympd(C_pred_inv,
			0.5 * (C_pred + C_pred.t()));
		if (!ok) C_pred_inv = arma::inv(C_pred);
		arma::vec w_r = W.slice(t).col(r);
		arma::mat HtH = S.slice(t) - (w_r * w_r.t());
		arma::vec Hty = B.slice(t).col(r);
		arma::mat post_prec = C_pred_inv + HtH / sigma2;
		arma::mat post_var;
		ok = arma::inv_sympd(post_var,
			0.5 * (post_prec + post_prec.t()));
		if (!ok) post_var = arma::inv(post_prec);
		post_var = 0.5 * (post_var + post_var.t());
		m_cur = post_var * (C_pred_inv * m_cur + Hty / sigma2);
		C_cur = post_var;
		m_fwd.col(t) = m_cur; C_fwd.slice(t) = C_cur;
	}
	// backward sample
	arma::vec a_now;
	{
		arma::mat C_T = C_fwd.slice(Tt - 1);
		C_T = 0.5 * (C_T + C_T.t()) + 1e-10 * arma::eye(n, n);
		try {
			a_now = arma::mvnrnd(m_fwd.col(Tt - 1), C_T);
		} catch (...) {
			a_now = m_fwd.col(Tt - 1);
		}
	}
	Aarray.slice(Tt - 1).row(r) = a_now.t();
	update_AB_suff_after_row(W, S, B, Aarray, Theta_arr, r, Tt - 1);
	for (int t = Tt - 2; t >= 0; t--) {
		arma::mat C_t = C_fwd.slice(t);
		arma::mat C_pred = C_t + W_proc;
		arma::mat C_pred_inv;
		bool ok = arma::inv_sympd(C_pred_inv,
			0.5 * (C_pred + C_pred.t()));
		if (!ok) C_pred_inv = arma::inv(C_pred);
		arma::mat B_t = C_t * C_pred_inv;
		arma::vec h_t = m_fwd.col(t) + B_t * (a_now - m_fwd.col(t));
		arma::mat H_t = C_t - B_t * C_pred * B_t.t();
		H_t = 0.5 * (H_t + H_t.t()) + 1e-10 * arma::eye(n, n);
		try {
			a_now = arma::mvnrnd(h_t, H_t);
		} catch (...) {
			a_now = h_t;
		}
		Aarray.slice(t).row(r) = a_now.t();
		if (t >= 1) update_AB_suff_after_row(W, S, B, Aarray,
			Theta_arr, r, t);
	}
}

// largest absolute eigenvalue of a real square matrix.
static double spectral_radius(const arma::mat& A) {
	arma::cx_vec ev = arma::eig_gen(A);
	double sr = 0.0;
	for (arma::uword i = 0; i < ev.n_elem; i++) {
		double m = std::abs(ev(i));
		if (m > sr) sr = m;
	}
	return sr;
}

// max over t of spectral_radius(Aarray.slice(t)). early exits if any
// slice already exceeds `bound`.
static double max_spectral_radius(const arma::cube& Aarray, double bound) {
	int Tt = (int)Aarray.n_slices;
	double max_sr = 0.0;
	for (int t = 0; t < Tt; t++) {
		double sr = spectral_radius(Aarray.slice(t));
		if (sr > max_sr) max_sr = sr;
		if (max_sr > bound) return max_sr;
	}
	return max_sr;
}

// row-FFBS with rejection on max_t rho(A_t) > rho_max. saves the
// row-r entries of Aarray and the W/S/B cubes before each attempt;
// on rejection, restores them and resamples (different RNG draw).
// returns the number of rejections before acceptance, or -max_rejects
// if we gave up (in which case the row r values are restored to their
// pre-call state, i.e., no update this sweep for row r).
static int row_ffbs_suffstat_rejecting(int r,
	arma::cube& Aarray,
	const arma::cube& Theta_arr,
	arma::cube& W, arma::cube& S, arma::cube& B,
	double sigma2, double tau2,
	double tau2_init_r,
	double rho_max,
	int max_rejects) {
	int n  = (int)Aarray.n_rows;
	int Tt = (int)Aarray.n_slices;

	// snapshot what the row update would modify
	arma::mat A_row_save(Tt, n);
	for (int t = 0; t < Tt; t++) {
		A_row_save.row(t) = Aarray.slice(t).row(r);
	}
	arma::cube W_save = W;
	arma::cube S_save = S;
	arma::cube B_save = B;

	int rejects = 0;
	while (true) {
		row_ffbs_suffstat(r, Aarray, Theta_arr, W, S, B,
			sigma2, tau2, tau2_init_r);
		double sr = max_spectral_radius(Aarray, rho_max);
		if (sr <= rho_max) return rejects;

		rejects++;
		// restore state for the next attempt
		for (int t = 0; t < Tt; t++) {
			Aarray.slice(t).row(r) = A_row_save.row(t);
		}
		W = W_save; S = S_save; B = B_save;

		if (rejects >= max_rejects) {
			return -max_rejects;  // gave up; row r left unchanged
		}
	}
}

// public wrapper: updates each row of A_{1:T} sequentially using the
// suff-stat row-FFBS. Theta is held fixed. returns a list with the
// updated A cube plus rejection diagnostics.
//
// tau2_init_per_row is a length-n vector of per-row A_1 prior
// variances. pass an empty / zero-length vector to tie each row's
// A_1 prior variance to tau2 (the default).
//
// rho_max: if > 0, enable rejection sampling on the row trajectory
// to keep max_t spectral_radius(A_t) <= rho_max. set to 0 to disable.
// max_rejects: number of resampling attempts before giving up on a
// row and leaving it unchanged for this sweep.
//
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List update_A_suffstat_cpp(const arma::cube& Theta_arr,
	const arma::cube& Aarray_old,
	double sigma2, double tau2,
	const arma::vec& tau2_init_per_row,
	double rho_max = 0.0,
	int max_rejects = 50) {
	int n  = (int)Aarray_old.n_rows;
	int Tt = (int)Aarray_old.n_slices;
	if ((int)Theta_arr.n_slices != Tt) {
		stop("update_A_suffstat_cpp: Theta and A must have matching "
			"time dimensions");
	}
	arma::vec tau2_init_use(n);
	if (tau2_init_per_row.n_elem == 0) {
		tau2_init_use.fill(tau2);
	} else {
		if ((int)tau2_init_per_row.n_elem != n) {
			stop("update_A_suffstat_cpp: tau2_init_per_row must have "
				"length n (or be empty to tie to tau2)");
		}
		tau2_init_use = tau2_init_per_row;
	}

	arma::cube Aarray = Aarray_old;
	arma::cube W, S, B;
	compute_AB_suff(Aarray, Theta_arr, W, S, B);

	arma::ivec rejects_per_row(n, arma::fill::zeros);
	int gave_up_rows = 0;

	for (int r = 0; r < n; r++) {
		if (rho_max > 0.0) {
			int rj = row_ffbs_suffstat_rejecting(r, Aarray, Theta_arr,
				W, S, B, sigma2, tau2, tau2_init_use(r),
				rho_max, max_rejects);
			if (rj < 0) {
				gave_up_rows++;
				rejects_per_row(r) = -rj;
			} else {
				rejects_per_row(r) = rj;
			}
		} else {
			row_ffbs_suffstat(r, Aarray, Theta_arr, W, S, B,
				sigma2, tau2, tau2_init_use(r));
		}
	}

	return List::create(
		Named("Aarray")          = Aarray,
		Named("rejects_per_row") = rejects_per_row,
		Named("gave_up_rows")    = gave_up_rows
	);
}

// ----------------------------------------------------------------------
// per-entry Metropolis-Gibbs scan for symmetric A_t (canonicalization
// of the orthogonal congruence orbit). bilinear model is
//   Z_t = A_t Z_{t-1} A_t + sigma E_t
// with A_t = A_t^T enforced and Z_t symmetric with zero diagonal.
//
// uses incremental residual updates for the proposal evaluation: a
// scalar perturbation at A_t entry (i, j) (with the symmetric pair
// at (j, i) when i != j) affects only O(n) entries of the bilinear
// product P_t = A_t Z_{t-1} A_t, so each M-H step costs O(n) instead
// of O(n^3).
//
// inputs:
//   Z_arr        observed/oracle state cube, n x n x Tt
//   Aarray_init  symmetric initial cube, n x n x Tt
//   sigma2       process noise variance (fixed during this call)
//   tau_A2       random-walk innovation variance (fixed)
//   kappa_A2     A_1 prior variance (per entry, fixed)
//   proposal_sd  scalar M-H proposal scale (adapted during burn)
//   n_sweeps     total sweeps to run
//   burn         burn-in count; samples returned from sweep burn+1
//   odens        thinning factor for retained samples
//   adapt        whether to adapt proposal_sd during burn
// ----------------------------------------------------------------------

// [[Rcpp::export]]
List update_A_symmetric_mh_cpp(const arma::cube& Z_arr,
	const arma::cube& Aarray_init,
	double sigma2, double tau_A2, double kappa_A2,
	double proposal_sd,
	int n_sweeps, int burn = 0, int odens = 1,
	bool adapt = true) {
	int n  = (int)Z_arr.n_rows;
	int Tt = (int)Z_arr.n_slices;
	if ((int)Aarray_init.n_rows != n ||
		(int)Aarray_init.n_cols != n ||
		(int)Aarray_init.n_slices != Tt) {
		stop("update_A_symmetric_mh_cpp: Aarray_init dimensions must "
			"match Z_arr");
	}
	arma::cube A = Aarray_init;

	// precompute P_t = A_t Z_{t-1} A_t and G_t = Z_{t-1} A_t for t = 1..Tt-1
	arma::cube P(n, n, Tt, arma::fill::zeros);
	arma::cube G(n, n, Tt, arma::fill::zeros);
	for (int t = 1; t < Tt; t++) {
		G.slice(t) = Z_arr.slice(t - 1) * A.slice(t);
		P.slice(t) = A.slice(t) * G.slice(t);
	}

	// storage for retained draws
	int n_keep = std::max(0, (n_sweeps - burn) / odens);
	arma::field<arma::cube> A_store(n_keep);
	arma::vec accept_rates(n_sweeps, arma::fill::zeros);

	double cur_sd = proposal_sd;
	int keep_id = 0;

	for (int sweep = 0; sweep < n_sweeps; sweep++) {
		int sweep_accepts = 0;
		long sweep_proposals = 0;

		for (int t = 0; t < Tt; t++) {
			for (int i = 0; i < n; i++) {
				for (int j = i; j < n; j++) {
					sweep_proposals++;
					double a_cur = A(i, j, t);
					double delta = cur_sd * R::rnorm(0.0, 1.0);
					double a_prop = a_cur + delta;

					double ll_diff = 0.0;
					double drss = 0.0;

					if (t >= 1) {
						const arma::mat& X = Z_arr.slice(t - 1);
						const arma::mat& Z_t = Z_arr.slice(t);
						const arma::mat& P_t = P.slice(t);
						const arma::mat& G_t = G.slice(t);

						if (i == j) {
							// H = E_{ii}: only entry (i, i) of A_t perturbed
							// affected upper-triangle off-diagonal entries:
							//   (i, l) for l > i: dP = delta * G[i, l]
							//   (k, i) for k < i: dP = delta * G[i, k]
							for (int l = i + 1; l < n; l++) {
								double dp = delta * G_t(i, l);
								double r_old = Z_t(i, l) - P_t(i, l);
								drss += -2.0 * r_old * dp + dp * dp;
							}
							for (int k = 0; k < i; k++) {
								double dp = delta * G_t(i, k);
								double r_old = Z_t(k, i) - P_t(k, i);
								drss += -2.0 * r_old * dp + dp * dp;
							}
						} else {
							// i < j; H = E_{ij} + E_{ji}: entries (i, j) and (j, i)
							// linear dP[k, l] = delta * (1{k=i} G[j, l] + 1{k=j} G[i, l]
							//                          + 1{l=j} G[i, k] + 1{l=i} G[j, k])
							// quadratic only at (i, j), (i, i), (j, j), (j, i)
							// in upper off-diagonal, only (i, j) gets quadratic contribution.

							// row i, l > i (l != i is automatic):
							for (int l = i + 1; l < n; l++) {
								double dp_lin = delta * G_t(j, l);
								if (l == j) {
									// from 1{l=j} G[i, k=i]
									dp_lin += delta * G_t(i, i);
								}
								double dp = dp_lin;
								if (l == j) {
									// quadratic at (i, j)
									dp += delta * delta * X(j, i);
								}
								double r_old = Z_t(i, l) - P_t(i, l);
								drss += -2.0 * r_old * dp + dp * dp;
							}
							// row j, l > j:
							for (int l = j + 1; l < n; l++) {
								double dp = delta * G_t(i, l);  // from 1{k=j}
								double r_old = Z_t(j, l) - P_t(j, l);
								drss += -2.0 * r_old * dp + dp * dp;
							}
							// col i, k < i:
							for (int k = 0; k < i; k++) {
								double dp = delta * G_t(j, k);  // from 1{l=i}
								double r_old = Z_t(k, i) - P_t(k, i);
								drss += -2.0 * r_old * dp + dp * dp;
							}
							// col j, k < j, k != i (since (i,j) already in row-i case):
							for (int k = 0; k < j; k++) {
								if (k == i) continue;
								double dp = delta * G_t(i, k);  // from 1{l=j}
								double r_old = Z_t(k, j) - P_t(k, j);
								drss += -2.0 * r_old * dp + dp * dp;
							}
						}

						ll_diff = -drss / (2.0 * sigma2);
					}

					// prior change (random walk + initial)
					double prior_diff = 0.0;
					if (t == 0) {
						prior_diff = (-0.5 * a_prop * a_prop / kappa_A2) -
							(-0.5 * a_cur * a_cur / kappa_A2);
					} else {
						double a_prev = A(i, j, t - 1);
						prior_diff += (-0.5 * (a_prop - a_prev) * (a_prop - a_prev) / tau_A2) -
							(-0.5 * (a_cur - a_prev) * (a_cur - a_prev) / tau_A2);
					}
					if (t < Tt - 1) {
						double a_next = A(i, j, t + 1);
						prior_diff += (-0.5 * (a_next - a_prop) * (a_next - a_prop) / tau_A2) -
							(-0.5 * (a_next - a_cur) * (a_next - a_cur) / tau_A2);
					}

					double log_alpha = ll_diff + prior_diff;
					if (std::isfinite(log_alpha) &&
						std::log(R::runif(0.0, 1.0)) < log_alpha) {
						// accept
						sweep_accepts++;

						// O(n) rank-2 update to G_t and P_t. We MUST apply the P
						// update first, while G_t_mut still holds G_old (since the
						// rank-2 P update reads G_old rows). Then update G, then A.
						if (t >= 1) {
							const arma::mat& X = Z_arr.slice(t - 1);
							arma::mat& G_t_mut = G.slice(t);
							arma::mat& P_t_mut = P.slice(t);
							if (i == j) {
								// P update (uses G_old): row i and column i, plus
								// the quadratic (i,i) correction.
								for (int p = 0; p < n; p++) {
									double Gi_p = G_t_mut(i, p); // G_old(i, p)
									P_t_mut(i, p) += delta * Gi_p; // row i
									P_t_mut(p, i) += delta * Gi_p; // col i: P[k,i] += delta * G_old(i, k)
								}
								P_t_mut(i, i) += delta * delta * X(i, i);
								// now update G: G[:, i] += delta * X[:, i]
								for (int p = 0; p < n; p++) {
									G_t_mut(p, i) += delta * X(p, i);
								}
							} else {
								// P update (uses G_old): rows i,j and cols i,j,
								// plus four quadratic corrections.
								for (int p = 0; p < n; p++) {
									double Gi_p = G_t_mut(i, p); // G_old(i, p)
									double Gj_p = G_t_mut(j, p); // G_old(j, p)
									// row i: P[i, l] += delta * G_old(j, l)
									P_t_mut(i, p) += delta * Gj_p;
									// row j: P[j, l] += delta * G_old(i, l)
									P_t_mut(j, p) += delta * Gi_p;
									// col j: P[k, j] += delta * G_old(i, k)
									P_t_mut(p, j) += delta * Gi_p;
									// col i: P[k, i] += delta * G_old(j, k)
									P_t_mut(p, i) += delta * Gj_p;
								}
								// quadratic corrections (X is symmetric)
								double dd = delta * delta;
								P_t_mut(i, j) += dd * X(i, j);
								P_t_mut(j, i) += dd * X(i, j);
								P_t_mut(i, i) += dd * X(j, j);
								P_t_mut(j, j) += dd * X(i, i);
								// now update G: G[:, j] += delta * X[:, i],
								//               G[:, i] += delta * X[:, j]
								for (int p = 0; p < n; p++) {
									G_t_mut(p, j) += delta * X(p, i);
									G_t_mut(p, i) += delta * X(p, j);
								}
							}
						}
						// finally update A (after P/G; A itself is not read by the
						// rank-2 update above)
						A(i, j, t) = a_prop;
						if (i != j) A(j, i, t) = a_prop;
					}
				}
			}
		}

		accept_rates(sweep) = (double)sweep_accepts /
			(double)sweep_proposals;

		// adapt during burn
		if (adapt && sweep > 25 && sweep < burn && (sweep % 25) == 0) {
			double recent_ar = arma::mean(
				accept_rates.subvec(sweep - 25, sweep - 1));
			if (recent_ar < 0.2) cur_sd *= 0.85;
			else if (recent_ar > 0.5) cur_sd *= 1.15;
		}

		// store retained sample
		if (sweep >= burn && ((sweep - burn) % odens) == 0 &&
			keep_id < n_keep) {
			A_store(keep_id) = A;
			keep_id++;
		}
	}

	return List::create(
		Named("A_store")           = A_store,
		Named("A_final")           = A,
		Named("accept_rates")      = accept_rates,
		Named("final_proposal_sd") = cur_sd
	);
}

// ----------------------------------------------------------------------
// diagonal-included variant of update_A_ar1_mh_cpp. additionally
// penalizes the diagonal residual of A_t Z_{t-1} A_t, treating
// Z_t[k,k] = 0 as a structural-zero observation. this adds n
// equations per time step which the off-diagonal-only likelihood
// was missing (the gap that left an n-dimensional flat per t).
//
// sigma2_diag is the noise variance for the diagonal residual; can
// be set equal to sigma2 (treating self-loops as observations on
// same scale) or larger (treating the structural zero as a softer
// constraint).
// ----------------------------------------------------------------------

// [[Rcpp::export]]
List update_A_ar1_mh_diag_cpp(const arma::cube& Z_arr,
	const arma::cube& Aarray_init,
	const arma::mat& Abar,
	double phi,
	double sigma2, double sigma2_diag, double tau_Delta2,
	double proposal_sd,
	int n_sweeps, int burn = 0, int odens = 1,
	bool adapt = true) {
	int n  = (int)Z_arr.n_rows;
	int Tt = (int)Z_arr.n_slices;
	arma::cube A = Aarray_init;

	double sigma_init2 = tau_Delta2 / (1.0 - phi * phi);

	arma::cube P(n, n, Tt, arma::fill::zeros);
	arma::cube G(n, n, Tt, arma::fill::zeros);
	for (int t = 1; t < Tt; t++) {
		G.slice(t) = Z_arr.slice(t - 1) * A.slice(t);
		P.slice(t) = A.slice(t) * G.slice(t);
	}

	int n_keep = std::max(0, (n_sweeps - burn) / odens);
	arma::field<arma::cube> A_store(n_keep);
	arma::vec accept_rates(n_sweeps, arma::fill::zeros);
	double cur_sd = proposal_sd;
	int keep_id = 0;

	for (int sweep = 0; sweep < n_sweeps; sweep++) {
		int sweep_accepts = 0;
		long sweep_proposals = 0;

		for (int t = 0; t < Tt; t++) {
			for (int i = 0; i < n; i++) {
				for (int j = i; j < n; j++) {
					sweep_proposals++;
					double a_cur = A(i, j, t);
					double delta = cur_sd * R::rnorm(0.0, 1.0);
					double a_prop = a_cur + delta;
					double abar_ij = Abar(i, j);

					double ll_diff = 0.0;
					double drss_off = 0.0;
					double drss_diag = 0.0;

					if (t >= 1) {
						const arma::mat& X = Z_arr.slice(t - 1);
						const arma::mat& Z_t = Z_arr.slice(t);
						const arma::mat& P_t = P.slice(t);
						const arma::mat& G_t = G.slice(t);

						if (i == j) {
							// off-diagonal residuals
							for (int l = i + 1; l < n; l++) {
								double dp = delta * G_t(i, l);
								double r_old = Z_t(i, l) - P_t(i, l);
								drss_off += -2.0 * r_old * dp + dp * dp;
							}
							for (int k = 0; k < i; k++) {
								double dp = delta * G_t(i, k);
								double r_old = Z_t(k, i) - P_t(k, i);
								drss_off += -2.0 * r_old * dp + dp * dp;
							}
							// diagonal residual at (i, i)
							// dP[i,i] = 2*delta*G[i,i] + delta^2 * X[i,i]
							double dp_diag = 2.0 * delta * G_t(i, i) +
								delta * delta * X(i, i);
							// r_old (i, i) = 0 - P_t(i, i)
							double r_old_diag = -P_t(i, i);
							drss_diag += -2.0 * r_old_diag * dp_diag +
								dp_diag * dp_diag;
						} else {
							// off-diagonal residuals
							for (int l = i + 1; l < n; l++) {
								double dp_lin = delta * G_t(j, l);
								if (l == j) dp_lin += delta * G_t(i, i);
								double dp = dp_lin;
								if (l == j) dp += delta * delta * X(j, i);
								double r_old = Z_t(i, l) - P_t(i, l);
								drss_off += -2.0 * r_old * dp + dp * dp;
							}
							for (int l = j + 1; l < n; l++) {
								double dp = delta * G_t(i, l);
								double r_old = Z_t(j, l) - P_t(j, l);
								drss_off += -2.0 * r_old * dp + dp * dp;
							}
							for (int k = 0; k < i; k++) {
								double dp = delta * G_t(j, k);
								double r_old = Z_t(k, i) - P_t(k, i);
								drss_off += -2.0 * r_old * dp + dp * dp;
							}
							for (int k = 0; k < j; k++) {
								if (k == i) continue;
								double dp = delta * G_t(i, k);
								double r_old = Z_t(k, j) - P_t(k, j);
								drss_off += -2.0 * r_old * dp + dp * dp;
							}
							// diagonal contributions:
							//  dP[i,i] = 2*delta*G[j,i] + delta^2 * X[j,j]
							//  dP[j,j] = 2*delta*G[i,j] + delta^2 * X[i,i]
							double dp_ii = 2.0 * delta * G_t(j, i) +
								delta * delta * X(j, j);
							double r_old_ii = -P_t(i, i);
							drss_diag += -2.0 * r_old_ii * dp_ii + dp_ii * dp_ii;
							double dp_jj = 2.0 * delta * G_t(i, j) +
								delta * delta * X(i, i);
							double r_old_jj = -P_t(j, j);
							drss_diag += -2.0 * r_old_jj * dp_jj + dp_jj * dp_jj;
						}
						ll_diff = -drss_off / (2.0 * sigma2) -
							drss_diag / (2.0 * sigma2_diag);
					}

					double prior_diff = 0.0;
					if (t == 0) {
						double d_prop = a_prop - abar_ij;
						double d_cur  = a_cur  - abar_ij;
						prior_diff += -(d_prop * d_prop - d_cur * d_cur) /
							(2.0 * sigma_init2);
					} else {
						double a_prev = A(i, j, t - 1);
						double m = (1.0 - phi) * abar_ij + phi * a_prev;
						prior_diff += -((a_prop - m) * (a_prop - m) -
							(a_cur - m) * (a_cur - m)) / (2.0 * tau_Delta2);
					}
					if (t < Tt - 1) {
						double a_next = A(i, j, t + 1);
						double m_next_prop = (1.0 - phi) * abar_ij + phi * a_prop;
						double m_next_cur  = (1.0 - phi) * abar_ij + phi * a_cur;
						prior_diff += -((a_next - m_next_prop) * (a_next - m_next_prop) -
							(a_next - m_next_cur) * (a_next - m_next_cur)) /
							(2.0 * tau_Delta2);
					}

					double log_alpha = ll_diff + prior_diff;
					if (std::isfinite(log_alpha) &&
						std::log(R::runif(0.0, 1.0)) < log_alpha) {
						sweep_accepts++;
						// O(n) rank-2 update to G_t and P_t. P update first
						// (reads G_old); then G update; then A update.
						if (t >= 1) {
							const arma::mat& X = Z_arr.slice(t - 1);
							arma::mat& G_t_mut = G.slice(t);
							arma::mat& P_t_mut = P.slice(t);
							if (i == j) {
								for (int p = 0; p < n; p++) {
									double Gi_p = G_t_mut(i, p);
									P_t_mut(i, p) += delta * Gi_p;
									P_t_mut(p, i) += delta * Gi_p;
								}
								P_t_mut(i, i) += delta * delta * X(i, i);
								for (int p = 0; p < n; p++) {
									G_t_mut(p, i) += delta * X(p, i);
								}
							} else {
								for (int p = 0; p < n; p++) {
									double Gi_p = G_t_mut(i, p);
									double Gj_p = G_t_mut(j, p);
									P_t_mut(i, p) += delta * Gj_p;
									P_t_mut(j, p) += delta * Gi_p;
									P_t_mut(p, j) += delta * Gi_p;
									P_t_mut(p, i) += delta * Gj_p;
								}
								double dd = delta * delta;
								P_t_mut(i, j) += dd * X(i, j);
								P_t_mut(j, i) += dd * X(i, j);
								P_t_mut(i, i) += dd * X(j, j);
								P_t_mut(j, j) += dd * X(i, i);
								for (int p = 0; p < n; p++) {
									G_t_mut(p, j) += delta * X(p, i);
									G_t_mut(p, i) += delta * X(p, j);
								}
							}
						}
						A(i, j, t) = a_prop;
						if (i != j) A(j, i, t) = a_prop;
					}
				}
			}
		}

		accept_rates(sweep) = (double)sweep_accepts /
			(double)sweep_proposals;
		if (adapt && sweep > 25 && sweep < burn && (sweep % 25) == 0) {
			double recent_ar = arma::mean(
				accept_rates.subvec(sweep - 25, sweep - 1));
			if (recent_ar < 0.2) cur_sd *= 0.85;
			else if (recent_ar > 0.5) cur_sd *= 1.15;
		}
		if (sweep >= burn && ((sweep - burn) % odens) == 0 &&
			keep_id < n_keep) {
			A_store(keep_id) = A;
			keep_id++;
		}
	}

	return List::create(
		Named("A_store")           = A_store,
		Named("A_final")           = A,
		Named("accept_rates")      = accept_rates,
		Named("final_proposal_sd") = cur_sd
	);
}

// ----------------------------------------------------------------------
// per-entry M-H Gibbs for symmetric A_t = Abar + Delta_t with stationary
// AR(1) on Delta_t. model:
//   A_t = Abar + Delta_t,  Abar = Abar^T,  Delta_t = Delta_t^T
//   Delta_1 ~ N(0, tau_Delta^2 / (1 - phi^2))   (stationary marginal)
//   Delta_t = phi Delta_{t-1} + eta_t           (AR(1) with eta_t ~ N(0, tau_Delta^2))
// equivalently on A_t directly:
//   A_t | A_{t-1}, Abar ~ N((1-phi) Abar + phi A_{t-1}, tau_Delta^2 I)
// likelihood: Z_t = A_t Z_{t-1} A_t + sigma E_t (triangle-only RSS)
//
// this kernel updates A_{1:T} conditional on Abar and phi. Abar should
// be updated separately in R via Gaussian conjugate given A_{1:T}.
// ----------------------------------------------------------------------

// [[Rcpp::export]]
List update_A_ar1_mh_cpp(const arma::cube& Z_arr,
	const arma::cube& Aarray_init,
	const arma::mat& Abar,
	double phi,
	double sigma2, double tau_Delta2,
	double proposal_sd,
	int n_sweeps, int burn = 0, int odens = 1,
	bool adapt = true) {
	int n  = (int)Z_arr.n_rows;
	int Tt = (int)Z_arr.n_slices;
	arma::cube A = Aarray_init;

	// stationary marginal variance for Delta_1
	double sigma_init2 = tau_Delta2 / (1.0 - phi * phi);

	// precompute P_t and G_t for the bilinear at each t
	arma::cube P(n, n, Tt, arma::fill::zeros);
	arma::cube G(n, n, Tt, arma::fill::zeros);
	for (int t = 1; t < Tt; t++) {
		G.slice(t) = Z_arr.slice(t - 1) * A.slice(t);
		P.slice(t) = A.slice(t) * G.slice(t);
	}

	int n_keep = std::max(0, (n_sweeps - burn) / odens);
	arma::field<arma::cube> A_store(n_keep);
	arma::vec accept_rates(n_sweeps, arma::fill::zeros);
	double cur_sd = proposal_sd;
	int keep_id = 0;

	for (int sweep = 0; sweep < n_sweeps; sweep++) {
		int sweep_accepts = 0;
		long sweep_proposals = 0;

		for (int t = 0; t < Tt; t++) {
			for (int i = 0; i < n; i++) {
				for (int j = i; j < n; j++) {
					sweep_proposals++;
					double a_cur = A(i, j, t);
					double delta = cur_sd * R::rnorm(0.0, 1.0);
					double a_prop = a_cur + delta;
					double abar_ij = Abar(i, j);

					double ll_diff = 0.0;
					double drss = 0.0;

					if (t >= 1) {
						const arma::mat& X = Z_arr.slice(t - 1);
						const arma::mat& Z_t = Z_arr.slice(t);
						const arma::mat& P_t = P.slice(t);
						const arma::mat& G_t = G.slice(t);

						if (i == j) {
							for (int l = i + 1; l < n; l++) {
								double dp = delta * G_t(i, l);
								double r_old = Z_t(i, l) - P_t(i, l);
								drss += -2.0 * r_old * dp + dp * dp;
							}
							for (int k = 0; k < i; k++) {
								double dp = delta * G_t(i, k);
								double r_old = Z_t(k, i) - P_t(k, i);
								drss += -2.0 * r_old * dp + dp * dp;
							}
						} else {
							for (int l = i + 1; l < n; l++) {
								double dp_lin = delta * G_t(j, l);
								if (l == j) dp_lin += delta * G_t(i, i);
								double dp = dp_lin;
								if (l == j) dp += delta * delta * X(j, i);
								double r_old = Z_t(i, l) - P_t(i, l);
								drss += -2.0 * r_old * dp + dp * dp;
							}
							for (int l = j + 1; l < n; l++) {
								double dp = delta * G_t(i, l);
								double r_old = Z_t(j, l) - P_t(j, l);
								drss += -2.0 * r_old * dp + dp * dp;
							}
							for (int k = 0; k < i; k++) {
								double dp = delta * G_t(j, k);
								double r_old = Z_t(k, i) - P_t(k, i);
								drss += -2.0 * r_old * dp + dp * dp;
							}
							for (int k = 0; k < j; k++) {
								if (k == i) continue;
								double dp = delta * G_t(i, k);
								double r_old = Z_t(k, j) - P_t(k, j);
								drss += -2.0 * r_old * dp + dp * dp;
							}
						}
						ll_diff = -drss / (2.0 * sigma2);
					}

					// AR(1) around Abar prior for A_t[i,j]
					double prior_diff = 0.0;
					if (t == 0) {
						// Delta_1 = A_1 - Abar; stationary marginal variance
						double d_prop = a_prop - abar_ij;
						double d_cur  = a_cur  - abar_ij;
						prior_diff += -(d_prop * d_prop - d_cur * d_cur) /
							(2.0 * sigma_init2);
					} else {
						// transition from t-1 to t:
						//   A_t = (1-phi) Abar + phi A_{t-1} + eta_t
						double a_prev = A(i, j, t - 1);
						double m = (1.0 - phi) * abar_ij + phi * a_prev;
						prior_diff += -((a_prop - m) * (a_prop - m) -
							(a_cur - m) * (a_cur - m)) /
							(2.0 * tau_Delta2);
					}
					if (t < Tt - 1) {
						// transition from t to t+1:
						//   A_{t+1} = (1-phi) Abar + phi A_t + eta_{t+1}
						double a_next = A(i, j, t + 1);
						double m_next_prop = (1.0 - phi) * abar_ij + phi * a_prop;
						double m_next_cur  = (1.0 - phi) * abar_ij + phi * a_cur;
						prior_diff += -((a_next - m_next_prop) * (a_next - m_next_prop) -
							(a_next - m_next_cur) * (a_next - m_next_cur)) /
							(2.0 * tau_Delta2);
					}

					double log_alpha = ll_diff + prior_diff;
					if (std::isfinite(log_alpha) &&
						std::log(R::runif(0.0, 1.0)) < log_alpha) {
						sweep_accepts++;
						// O(n) rank-2 update to G_t and P_t. P update first
						// (reads G_old); then G update; then A update.
						if (t >= 1) {
							const arma::mat& X = Z_arr.slice(t - 1);
							arma::mat& G_t_mut = G.slice(t);
							arma::mat& P_t_mut = P.slice(t);
							if (i == j) {
								for (int p = 0; p < n; p++) {
									double Gi_p = G_t_mut(i, p);
									P_t_mut(i, p) += delta * Gi_p;
									P_t_mut(p, i) += delta * Gi_p;
								}
								P_t_mut(i, i) += delta * delta * X(i, i);
								for (int p = 0; p < n; p++) {
									G_t_mut(p, i) += delta * X(p, i);
								}
							} else {
								for (int p = 0; p < n; p++) {
									double Gi_p = G_t_mut(i, p);
									double Gj_p = G_t_mut(j, p);
									P_t_mut(i, p) += delta * Gj_p;
									P_t_mut(j, p) += delta * Gi_p;
									P_t_mut(p, j) += delta * Gi_p;
									P_t_mut(p, i) += delta * Gj_p;
								}
								double dd = delta * delta;
								P_t_mut(i, j) += dd * X(i, j);
								P_t_mut(j, i) += dd * X(i, j);
								P_t_mut(i, i) += dd * X(j, j);
								P_t_mut(j, j) += dd * X(i, i);
								for (int p = 0; p < n; p++) {
									G_t_mut(p, j) += delta * X(p, i);
									G_t_mut(p, i) += delta * X(p, j);
								}
							}
						}
						A(i, j, t) = a_prop;
						if (i != j) A(j, i, t) = a_prop;
					}
				}
			}
		}

		accept_rates(sweep) = (double)sweep_accepts /
			(double)sweep_proposals;
		if (adapt && sweep > 25 && sweep < burn && (sweep % 25) == 0) {
			double recent_ar = arma::mean(
				accept_rates.subvec(sweep - 25, sweep - 1));
			if (recent_ar < 0.2) cur_sd *= 0.85;
			else if (recent_ar > 0.5) cur_sd *= 1.15;
		}
		if (sweep >= burn && ((sweep - burn) % odens) == 0 &&
			keep_id < n_keep) {
			A_store(keep_id) = A;
			keep_id++;
		}
	}

	return List::create(
		Named("A_store")           = A_store,
		Named("A_final")           = A,
		Named("accept_rates")      = accept_rates,
		Named("final_proposal_sd") = cur_sd
	);
}
