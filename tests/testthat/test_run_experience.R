# Regression tests for run-experience fixes:
#   Issue 1 - progress-bar percentage / ETA
#   Issue 2 - asymmetric-fit stability (no silent process death)
#   Issue 3 - set_dbn_threads honesty
#   Issue 4 - keep argument controlling object size

# ---------------------------------------------------------------------------
# Issue 1: progress bar must reflect true completion of (burn + nscan).
#
# The bug: cli_progress_bar() was created with total = n_iter, but
# cli_progress_update() (a relative +1) was called only every `verbose`
# iterations. With the default verbose = 100 the bar received n_iter / 100
# updates and so could only ever reach 1% of n_iter. The fix calls
# cli_progress_update(set = iter), an absolute set, so the bar reflects the
# true iteration count regardless of update frequency.
#
# cli progress bars render only to an interactive terminal, so we cannot
# read the displayed percentage back. Instead we exercise the update path
# (numeric verbose => the "every k iterations" branch) and assert the fit
# completes cleanly and reports the correct total iteration count.
test_that("Issue 1: dynamic fit with periodic progress updates completes", {
	sim <- simulate_dynamic_dbn(n = 5, time = 6, seed = 1)
	expect_no_error(
		fit <- dbn_dynamic(sim$Y, family = "gaussian",
						   nscan = 40, burn = 20, odens = 2,
						   verbose = 5L, seed = 1)
	)
	# n_iter is burn + nscan; the absolute-set progress path must agree
	expect_equal(fit$n_iter, 60)
})

# ---------------------------------------------------------------------------
# Issue 2: the asymmetric (symmetric = FALSE) dynamic sampler must not die
# silently. The crash was a C++ exception thrown inside an OpenMP parallel
# region (the throwing arma::inv fallback in dbn_safe_inv_sympd), which
# escapes the region and triggers std::terminate. dbn_safe_inv_sympd is now
# non-throwing, and an R-side divergence guard converts any non-finite
# influence matrices into an informative R error.
test_that("Issue 2: asymmetric dynamic fit completes with finite output", {
	sim <- simulate_dynamic_dbn(n = 5, time = 6, seed = 2)
	expect_no_error(
		fit <- dbn_dynamic(sim$Y, family = "gaussian", symmetric = FALSE,
						   nscan = 40, burn = 20, odens = 2, seed = 2)
	)
	# the divergence guard must not have fired on benign data, and the
	# influence-matrix draws must be finite
	expect_true(all(vapply(fit$A, function(a) all(is.finite(a)), logical(1))))
	expect_true(all(vapply(fit$B, function(b) all(is.finite(b)), logical(1))))
})

# ---------------------------------------------------------------------------
# Issue 3: set_dbn_threads must round-trip honestly and emit a one-time note
# explaining the limited within-fit speedup of the dynamic sampler.
test_that("Issue 3: set_dbn_threads round-trips and notes its scope", {
	old <- get_dbn_threads()
	on.exit(set_dbn_threads(old, quiet = TRUE), add = TRUE)

	set_dbn_threads(1L, quiet = TRUE)
	expect_equal(get_dbn_threads(), 1L)

	# the honesty note fires once for n_threads > 1
	options(dbn.threads_note_shown = NULL)
	expect_message(set_dbn_threads(2L), "limited speedup")
	expect_equal(get_dbn_threads(), 2L)

	# quiet = TRUE suppresses the note
	options(dbn.threads_note_shown = NULL)
	expect_no_message(set_dbn_threads(2L, quiet = TRUE))

	# NULL resets to single-threaded
	set_dbn_threads(NULL, quiet = TRUE)
	expect_equal(get_dbn_threads(), 1L)
})

# ---------------------------------------------------------------------------
# Issue 4: the keep argument must drop the requested draw-indexed components
# and yield a materially smaller object.
test_that("Issue 4: keep argument drops Theta and shrinks the fit object", {
	sim <- simulate_dynamic_dbn(n = 6, time = 8, seed = 4)

	full <- dbn_dynamic(sim$Y, family = "gaussian",
					   nscan = 60, burn = 20, odens = 2, verbose = FALSE, seed = 4)
	lean <- dbn_dynamic(sim$Y, family = "gaussian",
					   nscan = 60, burn = 20, odens = 2, verbose = FALSE, seed = 4,
					   keep = c("A", "B", "M"))

	# full object retains the latent-state draws; lean object drops them
	expect_false(is.null(full$Theta))
	expect_true(is.null(lean$Theta))
	expect_true(is.null(lean$draws$theta))

	# A/B/M are retained in the lean object
	expect_false(is.null(lean$A))
	expect_false(is.null(lean$B))
	expect_false(is.null(lean$M))

	# the lean object is smaller
	expect_lt(as.numeric(object.size(lean)), as.numeric(object.size(full)))

	# an unrecognized keep value is rejected
	expect_error(
		dbn_dynamic(sim$Y, family = "gaussian", nscan = 20, burn = 10,
				   verbose = FALSE, keep = c("A", "Bogus")),
		"Unrecognized"
	)
})

# ---------------------------------------------------------------------------
# Proactive run-experience features (v1.1.1):
#   - asymmetric scale renormalization (no scale drift on directed fits)
#   - estimate_memory() accuracy (accounts for double-stored components)
#   - log-safe heartbeat, pre-flight panel, closing summary
#   - dbn_fit_many() multi-fit wrapper

test_that(".dbn_fmt_duration formats h/m/s compactly", {
	expect_equal(dbn:::.dbn_fmt_duration(45), "45s")
	expect_equal(dbn:::.dbn_fmt_duration(200), "3m20s")
	expect_equal(dbn:::.dbn_fmt_duration(3661), "1h01m")
	expect_equal(dbn:::.dbn_fmt_duration(NA_real_), "--")
	expect_equal(dbn:::.dbn_fmt_duration(-5), "--")
})

test_that(".dbn_data_looks_undirected detects symmetric adjacency", {
	sym <- array(0, c(5, 5, 3))
	for (t in 1:3) {
		M <- matrix(stats::rnorm(25), 5, 5)
		sym[, , t] <- M + t(M)
	}
	expect_true(dbn:::.dbn_data_looks_undirected(sym, 5, 5))
	expect_false(dbn:::.dbn_data_looks_undirected(array(stats::rnorm(75), c(5, 5, 3)), 5, 5))
	# a non-square (bipartite) network is never flagged
	expect_false(dbn:::.dbn_data_looks_undirected(array(stats::rnorm(30), c(5, 6, 1)), 5, 6))
})

test_that(".dbn_mixing_note guards trivial traces", {
	expect_null(dbn:::.dbn_mixing_note(rep(1, 50)))
	expect_null(dbn:::.dbn_mixing_note(c(1, 2, 3)))
	note <- dbn:::.dbn_mixing_note(stats::rnorm(100), "sigma2")
	expect_true(is.character(note) && grepl("autocorrelation", note))
})

# ---------------------------------------------------------------------------
# estimate_memory must predict the true fitted-object size. The fix accounts
# for Theta and M being held in two representations; the earlier formula
# counted them once and so under-reported.
test_that("estimate_memory predicts the fitted object size", {
	sim <- simulate_dynamic_dbn(n = 10, time = 12, seed = 7)
	fit <- dbn_dynamic(sim$Y, family = "gaussian",
					   nscan = 200, burn = 50, odens = 2,
					   verbose = FALSE, seed = 7)
	# estimate with the dimensions the fit actually used
	pred <- estimate_memory(n_row = fit$dims$n_row, n_col = fit$dims$n_col,
							p = fit$dims$p, Tt = fit$dims$Tt,
							nscan = 200, burn = 50, odens = 2,
							family = "gaussian", quiet = TRUE)
	expect_named(pred, c("object_gb", "peak_gb"))
	actual_gb <- length(serialize(fit, NULL)) / 1024^3
	# a genuine relative check: expect_equal()'s tolerance acts as an
	# absolute bound at this object scale and would pass for any value.
	ratio <- pred[["object_gb"]] / actual_gb
	expect_gt(ratio, 0.85)
	expect_lt(ratio, 1.15)
})

test_that("estimate_memory reflects the keep argument", {
	full <- estimate_memory(n_row = 25, Tt = 88, nscan = 1000, odens = 10,
							family = "gaussian", quiet = TRUE)
	lean <- estimate_memory(n_row = 25, Tt = 88, nscan = 1000, odens = 10,
							family = "gaussian", keep = c("A", "B", "M"),
							quiet = TRUE)
	# dropping the Theta draws is a substantial reduction
	expect_lt(lean[["object_gb"]], 0.8 * full[["object_gb"]])
})

# ---------------------------------------------------------------------------
# The heartbeat is a plain cli_inform line, so it survives non-interactive
# and redirected output where an animated progress bar does not render.
test_that("dynamic fit emits a log-safe heartbeat line", {
	sim <- simulate_dynamic_dbn(n = 5, time = 6, seed = 11)
	expect_message(
		dbn_dynamic(sim$Y, family = "gaussian", nscan = 40, burn = 20,
					odens = 2, verbose = 10L, seed = 11),
		"iter"
	)
})

# ---------------------------------------------------------------------------
# The asymmetric model is scale-ambiguous; without renormalization the A/B
# scale random-walks until a matrix overflows. Renormalization rebalances A
# and B to a common RMS each iteration, so every stored draw has matching
# A and B scales and the fit stays finite.
test_that("asymmetric dynamic fit stays scale-balanced", {
	sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 21)
	fit <- dbn_dynamic(sim$Y, family = "gaussian", symmetric = FALSE,
					   nscan = 200, burn = 100, odens = 2,
					   verbose = FALSE, seed = 21)
	expect_true(all(vapply(fit$A, function(a) all(is.finite(a)), logical(1))))
	expect_true(all(vapply(fit$B, function(b) all(is.finite(b)), logical(1))))
	rms <- function(x) sqrt(mean(x^2))
	ratios <- mapply(function(a, b) rms(a) / rms(b), fit$A, fit$B)
	expect_equal(ratios, rep(1, length(ratios)), tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# dbn_fit_many runs independent fits and isolates a failed fit instead of
# aborting the whole batch.
test_that("dbn_fit_many fits a list and isolates failures", {
	sims <- list(
		ok_a = simulate_dynamic_dbn(n = 5, time = 6, seed = 31)$Y,
		ok_b = simulate_dynamic_dbn(n = 5, time = 6, seed = 32)$Y
	)
	fits <- dbn_fit_many(sims, model = "dynamic", family = "gaussian",
						 nscan = 30, burn = 10, odens = 2,
						 mc.cores = 1, verbose = FALSE)
	expect_length(fits, 2)
	expect_named(fits, c("ok_a", "ok_b"))
	expect_s3_class(fits$ok_a, "dbn")
	expect_s3_class(fits$ok_b, "dbn")

	mixed <- list(good = sims$ok_a, bad = array(1, dim = c(2, 2)))
	res <- suppressWarnings(
		dbn_fit_many(mixed, model = "dynamic", family = "gaussian",
					 nscan = 30, burn = 10, mc.cores = 1, verbose = FALSE)
	)
	expect_s3_class(res$good, "dbn")
	expect_s3_class(res$bad, "dbn_fit_error")

	# seed must be passed via `seeds`, not through `...`
	expect_error(
		dbn_fit_many(sims, model = "dynamic", family = "gaussian",
					 nscan = 30, burn = 10, seed = 1, verbose = FALSE),
		"seed"
	)
})
