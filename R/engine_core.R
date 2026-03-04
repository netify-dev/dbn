####
# shared preprocessing
####

#' Generic MCMC Engine for DBN Models
#'
#' @description Core Gibbs sampler that drives all DBN model variants
#' @name engine_core
#' @keywords internal
NULL

#' Shared Preprocessing
#'
#' @description Common preprocessing for all DBN models
#' @param Y Data array
#' @param family Distribution family
#' @return List with preprocessed components
#' @keywords internal
shared_preprocess <- function(Y, family = "ordinal") {
	# accept sparse input
	if (inherits(Y, "dgCMatrix")) {
		Y <- as(Y, "TsparseMatrix")
	}

	n_row <- dim(Y)[1]
	n_col <- dim(Y)[2]
	is_bipartite <- (n_row != n_col)

	dims <- list(
		m = n_row,
		n_row = n_row,
		n_col = n_col,
		p = dim(Y)[3],
		Tt = dim(Y)[4],
		is_bipartite = is_bipartite
	)

	# family-specific data validation
	if (family == "binary") {
		Y_vals <- Y[!is.na(Y)]
		if (length(Y_vals) > 0 && !all(Y_vals %in% c(0, 1))) {
			cli::cli_abort(c(
				"Binary family requires 0/1 data.",
				"x" = "Found values outside {{0, 1}}.",
				"i" = "Use {.code family = \"gaussian\"} for continuous data or {.code family = \"ordinal\"} for integer data."
			))
		}
	} else if (family == "ordinal") {
		Y_vals <- Y[!is.na(Y)]
		if (length(Y_vals) > 0 && !all(Y_vals == floor(Y_vals))) {
			cli::cli_warn(c(
				"Ordinal family expects integer-valued data.",
				"i" = "Non-integer values detected; data will be treated as continuous ranks."
			))
		}
	}

	# precompute ranks
	IR <- precompute_ranks(Y)

	# initialize z based on family
	Z <- Y
	if (family == "ordinal") {
		means <- numeric(dims$p)
		sds <- numeric(dims$p)
		for (j in 1:dims$p) {
			Y_j <- Y[, , j, ]
			means[j] <- mean(Y_j, na.rm = TRUE)
			sds[j] <- sd(Y_j, na.rm = TRUE)
		}
		Y_flat <- array(Y, dim = c(n_row, n_col, dims$p * dims$Tt))
		Z_flat <- compute_zscores_batch(Y_flat, means, sds, n_row, n_col, dims$p, dims$Tt)
		Z <- array(Z_flat, dim = c(n_row, n_col, dims$p, dims$Tt))
		Z[is.na(Z)] <- 0
	} else {
		Z <- Y
		if (family == "gaussian") {
			Z[!is.finite(Z)] <- 0
		}
	}

	# initialize mean
	M <- array(apply(Z, c(1, 2, 3), mean, na.rm = TRUE), dim = c(n_row, n_col, dims$p))

	# initialize theta
	Theta <- sweep(Z, c(1, 2, 3), M, "-") + rsan(dim(Z))

	list(
		Y = Y,
		Z = Z,
		M = M,
		Theta = Theta,
		R = Y,
		IR = IR,
		dims = dims
	)
}

####
# z and mu update
####

#' Generic Update for Z and Mu
#'
#' @description Updates latent Z values and baseline mu
#' @param pre Current preprocessing object
#' @param g2 Variance for mu
#' @param FAM Family object
#' @param Aarray Optional A array for binary family
#' @param Barray Optional B array for binary family
#' @return Updated preprocessing object
#' @keywords internal
update_Z_mu <- function(pre, g2, FAM, Aarray = NULL, Barray = NULL) {
	# update mu
	mu_var <- 1 / (pre$dims$Tt + 1 / g2)
	mu_hat <- mu_var * apply(pre$Z - pre$Theta, c(1, 2, 3), sum)
	pre$M <- mu_hat + sqrt(mu_var) * rsan(dim(pre$M))

	# family-specific z refresh
	pre <- FAM$draw_latent(pre, Aarray = Aarray, Barray = Barray)
	pre
}

####
# model function registry
####

#' Create Model Functions List
#'
#' @description Creates list of model-specific update functions
#' @param model Model type
#' @return List of update functions
#' @keywords internal
create_model_funs <- function(model) {
	switch(model,
		static = list(
			init = init_static,
			update_Z = update_Z_static,
			update_Theta = update_Theta_static,
			update_AB = update_AB_static,
			update_hyper = update_hyper_static,
			collect = collect_static
		),
		dynamic = list(
			init = init_dynamic,
			update_Z = update_Z_dynamic,
			update_Theta = update_Theta_dynamic,
			update_AB = update_AB_dynamic,
			update_hyper = update_hyper_dynamic,
			collect = collect_dynamic
		),
		lowrank = list(
			init = init_lowrank,
			update_Z = update_Z_lowrank,
			update_Theta = update_Theta_lowrank,
			update_factor = update_factor_lowrank,
			update_hyper = update_hyper_lowrank,
			collect = collect_lowrank
		),
		hmm = list(
			init = init_hmm,
			update_Z = update_Z_hmm,
			update_Theta = update_Theta_hmm,
			update_state = update_state_hmm,
			update_regime = update_regime_hmm,
			update_hyper = update_hyper_hmm,
			collect = collect_hmm
		)
	)
}

#' Get package environment
#' @keywords internal
get_pkg_env <- function() {
	ns <- getNamespace("dbn")
	if (!exists(".pkg_env", envir = ns)) {
		assign(".pkg_env", new.env(parent = emptyenv()), envir = ns)
		assign("models", list(), envir = get(".pkg_env", envir = ns))
	}
	get(".pkg_env", envir = ns)
}

#' Register DBN Model
#'
#' @description Register a new DBN model variant
#' @param name Model name
#' @param fun_list List of model functions
#' @examples
#' \dontrun{
#' custom_model <- list(
#'     init = function(data) list(theta = rnorm(10)),
#'     update_Z = function(state, data) state,
#'     update_Theta = function(state, data) {
#'         state$theta <- state$theta + rnorm(10, 0, 0.1)
#'         state
#'     },
#'     collect = function(state) list(theta = state$theta)
#' )
#' register_dbn_model("custom", custom_model)
#' }
#' @keywords internal
register_dbn_model <- function(name, fun_list) {
	env <- get_pkg_env()
	env$models[[name]] <- fun_list
}
