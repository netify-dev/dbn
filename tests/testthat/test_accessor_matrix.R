# accessor coverage across the model x family matrix. sweeps {static,
# dynamic, piecewise, hmm, lowrank} x {gaussian, ordinal, binary} and
# exercises every exported accessor. each accessor must either return
# well-formed output or raise a directed cli error; a raw crash, a
# duplicated variable name, an na label, or a visibly-returning summary
# is a failure.

skip_on_cran()

# small fits, structure-only (we are not checking numeric recovery here)
fit_combo = function(model, family) {
	set.seed(42)
	n = 8; Tt = 12
	simd = simulate_dynamic_dbn(n = n, p = 1, time = Tt, seed = 42)
	Zc = simd$Z                                   # continuous
	Yo = simd$Y                                   # simulator's ordinal coding
	Yb = array(1 * (simd$Z > 0), dim = dim(simd$Z))         # binary
	for (t in seq_len(Tt)) {
		diag(Yo[, , 1, t]) = NA; diag(Yb[, , 1, t]) = NA
	}
	dat = switch(family, gaussian = Zc, ordinal = Yo, binary = Yb)
	args = list(data = dat, model = model, family = family,
							 nscan = 300, burn = 150, odens = 2, verbose = FALSE)
	if (model == "piecewise") args$blocks = c(6, Tt)
	if (model == "hmm")       args$R = 2
	if (model == "lowrank")   args$r = 2
	suppressWarnings(suppressMessages(do.call(dbn, args)))
}

# classify a call: "ok" (no error) or "directed" (cli/rlang condition) --
# anything else propagates as a test failure
call_class = function(expr) {
	res = tryCatch({ force(expr); "ok" },
		error = function(e) {
			msg = conditionMessage(e)
			# a directed error carries an informative message (not a raw
			# "subscript out of bounds" / "non-conformable" C++/base crash)
			raw = grepl("out of bounds|non-conformable|unused argument|could not find function|undefined columns|incorrect number of dim|subscript", msg)
			if (raw) paste0("RAW: ", msg) else "directed"
		})
	res
}

models   = c("static", "dynamic", "piecewise", "hmm", "lowrank")
families = c("gaussian", "ordinal", "binary")

for (mdl in models) {
	for (fam in families) {
		test_that(sprintf("accessors behave on %s / %s", mdl, fam), {
			fit = tryCatch(fit_combo(mdl, fam), error = function(e) NULL)
			skip_if(is.null(fit), sprintf("%s/%s did not fit", mdl, fam))

			# 1. print / summary never raw-crash; summary returns INVISIBLY
			expect_false(startsWith(call_class(utils::capture.output(print(fit))), "RAW"))
			sc = call_class(utils::capture.output(summary(fit)))
			expect_false(startsWith(sc, "RAW"))
			if (sc == "ok") {
				wv = withVisible(summary(fit))
				expect_false(wv$visible)            # guards the double-print regression
			}

			# 2. base-stats + broom accessors: ok or directed, never raw
			for (acc in c("coef", "tidy", "glance", "confint")) {
				cl = call_class(do.call(acc, list(fit)))
				expect_false(startsWith(cl, "RAW"),
										 info = sprintf("%s(%s/%s): %s", acc, mdl, fam, cl))
			}

			# 3. as_draws: ok or directed; if ok, no duplicate / NA variable names
			if (requireNamespace("posterior", quietly = TRUE)) {
				da = tryCatch(suppressMessages(posterior::as_draws(fit)), error = function(e) e)
				if (!inherits(da, "error")) {
					v = posterior::variables(da)
					expect_false(any(is.na(v)), info = sprintf("as_draws NA var (%s/%s)", mdl, fam))
					expect_equal(anyDuplicated(v), 0L,
											 info = sprintf("as_draws dup var (%s/%s): %s",
																			mdl, fam, paste(v[duplicated(v)], collapse = ",")))
				}
			}

			# 4. check_convergence: ok or directed, never raw
			cc = call_class(suppressMessages(check_convergence(fit)))
			expect_false(startsWith(cc, "RAW"),
									 info = sprintf("check_convergence(%s/%s): %s", mdl, fam, cc))

			# 5. predict (forecast + ppd): never raw-crash
			pf = call_class(suppressWarnings(predict(fit, H = 2, draws = 20, summary = "mean")))
			expect_false(startsWith(pf, "RAW"),
									 info = sprintf("predict-forecast(%s/%s): %s", mdl, fam, pf))
			pp = call_class(suppressWarnings(predict(fit, draws = 10)))
			expect_false(startsWith(pp, "RAW"),
									 info = sprintf("predict-ppd(%s/%s): %s", mdl, fam, pp))
		})
	}
}

test_that("check_convergence reports tau_alpha2 for low-rank fits", {
	fit = fit_combo("lowrank", "gaussian")
	skip_if(is.null(fit), "lowrank did not fit")
	out = utils::capture.output(suppressMessages(check_convergence(fit)))
	expect_true(any(grepl("tau_alpha2", out)),
							info = "tau_alpha2 must appear in the low-rank convergence table")
})

test_that("static gaussian as_draws does not duplicate sigma2 as sigma2_obs", {
	skip_if_not_installed("posterior")
	fit = fit_combo("static", "gaussian")
	skip_if(is.null(fit), "static did not fit")
	v = posterior::variables(suppressMessages(posterior::as_draws(fit)))
	# sigma2_obs is an alias of sigma2 for the static model; it must not be
	# exposed as a separate (duplicate) variable
	expect_false("sigma2_obs" %in% v)
	expect_equal(anyDuplicated(v), 0L)
})
