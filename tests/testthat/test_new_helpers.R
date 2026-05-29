# build_shock accepts character actor names via a fit
test_that("build_shock accepts character actor names via a fit", {
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 6, seed = 1)
	dimnames(sim$Z) = list(
		c("a", "b", "c", "d", "e"),
		c("a", "b", "c", "d", "e"),
		"tie",
		as.character(seq_len(dim(sim$Z)[4]))
	)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	S = build_shock(fit, type = "unit_edge", i = "a", j = "c", magnitude = 1)
	expect_equal(dim(S), c(5, 5))
	expect_equal(S[1, 3], 1)
	expect_error(
		build_shock(fit, type = "unit_edge", i = "ZZ", j = "c"),
		regexp = "not found"
	)
})

# dyad_path title includes character actor labels
test_that("dyad_path title includes character actor labels", {
	skip_if_not_installed("ggplot2")
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 6, seed = 1)
	dimnames(sim$Z) = list(
		c("a", "b", "c", "d", "e"),
		c("a", "b", "c", "d", "e"),
		"tie",
		as.character(seq_len(dim(sim$Z)[4]))
	)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	p = dyad_path(fit, i = "a", j = "c")
	expect_true(grepl("a.*c", p$labels$title))
})

# as_dbn_array_edgelist builds the right 4D shape
test_that("as_dbn_array_edgelist builds the right 4D shape", {
	edges = data.frame(
		from = c("a", "b", "c"),
		to = c("b", "c", "a"),
		time = c(1, 1, 2),
		w = c(1, 1, 1)
	)
	Y = as_dbn_array_edgelist(edges, weight = "w")
	expect_equal(dim(Y), c(3, 3, 1, 2))
	expect_true(all(is.na(diag(Y[, , 1, 1]))))
})

# node_embedding returns n x r matrix with row names
test_that("node_embedding returns n x r matrix with row names", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 6, p = 1, time = 8, r = 2, seed = 1)
	dimnames(sim$Z) = list(
		letters[1:6],
		letters[1:6],
		"tie",
		as.character(seq_len(dim(sim$Z)[4]))
	)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "lowrank",
			family = "gaussian",
			r = 2,
			nscan = 100,
			burn = 50,
			odens = 2,
			verbose = FALSE
		)
	))
	U = node_embedding(fit)
	expect_equal(dim(U), c(6, 2))
	expect_equal(rownames(U), letters[1:6])
})

# posterior_predict.dbn returns a draws x obs matrix
test_that("posterior_predict.dbn returns a draws x obs matrix", {
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 6, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	yrep = posterior_predict.dbn(fit, ndraws = 4)
	expect_true(is.matrix(yrep))
	expect_equal(nrow(yrep), 4)
})

# log_lik.dbn returns n_draws x n_obs for gaussian dynamic fit
test_that("log_lik.dbn returns n_draws x n_obs for gaussian dynamic fit", {
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 6, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	ll = log_lik.dbn(fit)
	expect_true(is.matrix(ll))
	expect_true(all(is.finite(ll)))
})

# coupling_trajectory returns tidy frame with ribbon columns
test_that("coupling_trajectory returns tidy frame with ribbon columns", {
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 6, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	ct = coupling_trajectory(fit)
	expect_true(all(c("actor", "t", "mean", "lo", "hi") %in% names(ct)))
	expect_equal(nrow(ct), 5 * dim(fit$Theta)[4])
})

# posterior_communities returns named membership
test_that("posterior_communities returns named membership", {
	skip_on_cran()
	sim = simulate_lowrank_dbn(n = 8, p = 1, time = 8, r = 2, seed = 1)
	dimnames(sim$Z) = list(
		letters[1:8],
		letters[1:8],
		"tie",
		as.character(seq_len(dim(sim$Z)[4]))
	)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "lowrank",
			family = "gaussian",
			r = 2,
			nscan = 100,
			burn = 50,
			odens = 2,
			verbose = FALSE
		)
	))
	pc = posterior_communities(fit, k = 3)
	expect_length(pc$membership, 8)
	expect_equal(names(pc$membership), letters[1:8])
})

# format.dbn produces a one-line description
test_that("format.dbn produces a one-line description", {
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 6, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	s = format(fit)
	expect_true(grepl("dbn", s) && grepl("dynamic", s) && grepl("gaussian", s))
})

# dbn_operator returns W, W_abs_mean, and W_frob
test_that("dbn_operator returns W, W_abs_mean, and W_frob", {
	sim = simulate_dynamic_dbn(n = 5, p = 1, time = 6, seed = 1)
	fit = suppressWarnings(suppressMessages(
		dbn(
			sim$Z,
			model = "dynamic",
			family = "gaussian",
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE
		)
	))
	op = dbn_operator(fit)
	expect_true(all(c("W", "W_abs_mean", "W_frob") %in% names(op)))
	expect_true(all(op$W_abs_mean >= abs(op$W)))
})

# stat_reciprocity returns NA on symmetric input
test_that("stat_reciprocity returns NA on symmetric input", {
	X_sym = matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3)
	expect_true(is.na(stat_reciprocity(X_sym)))
	expect_true(is.na(stat_reciprocity(matrix(0, 4, 4))))
})
