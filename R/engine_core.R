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
#' @return List with preprocessed components
#' @keywords internal
shared_preprocess <- function(Y, family = "ordinal") {
    # accept sparse input
    if (inherits(Y, "dgCMatrix")) {
        Y <- as(Y, "TsparseMatrix") # uniform sparse layout
    }

    dims <- list(
        m = dim(Y)[1],
        p = dim(Y)[3],
        Tt = dim(Y)[4]
    )

    # precompute ranks (always compute for ordinal compatibility)
    IR <- precompute_ranks(Y)

    # initialize z based on family
    Z <- Y
    if (family == "ordinal") {
        # rank-likelihood transformation - compute means and sds for vectorized operation
        means <- numeric(dims$p)
        sds <- numeric(dims$p)
        for (j in 1:dims$p) {
            Y_j <- Y[, , j, ]
            means[j] <- mean(Y_j, na.rm = TRUE)
            sds[j] <- sd(Y_j, na.rm = TRUE)
        }
        # 
        Y_flat <- array(Y, dim = c(dims$m, dims$m, dims$p * dims$Tt))
        Z_flat <- compute_zscores_batch(Y_flat, means, sds, dims$m, dims$p, dims$Tt)
        Z <- array(Z_flat, dim = c(dims$m, dims$m, dims$p, dims$Tt))
        Z[is.na(Z)] <- 0
    } else {
        # gaussian and binary use y directly
        Z <- Y
        if (family == "gaussian") {
            # any na / inf / -inf would crash the state-space sampler - coerce to 0
            Z[!is.finite(Z)] <- 0
        }
    }

    # initialize mean - ensure proper dimensions
    M <- array(apply(Z, c(1, 2, 3), mean, na.rm = TRUE), dim = c(dims$m, dims$m, dims$p))

    # initialize theta
    Theta <- sweep(Z, c(1, 2, 3), M, "-") + rsan(dim(Z))

    list(
        Y = Y,
        Z = Z,
        M = M,
        Theta = Theta,
        R = Y, # keep original data for all families
        IR = IR,
        dims = dims
    )
}

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
    # update mu (identical for every family)
    mu_var <- 1 / (pre$dims$Tt + 1 / g2)
    mu_hat <- mu_var * apply(pre$Z - pre$Theta, c(1, 2, 3), sum)
    pre$M <- mu_hat + sqrt(mu_var) * rsan(dim(pre$M))

    # family-specific z refresh
    pre <- FAM$draw_latent(pre, Aarray = Aarray, Barray = Barray)

    pre
}

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
    # use a cached environment stored in the package namespace
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
#' @param fun_list List of model functions including: init, update_Z,
#'   update_Theta, and collect functions
#' @examples
#' \dontrun{
#' # Register a custom DBN model variant
#' custom_model <- list(
#'     init = function(data) {
#'         # Initialize model state
#'         list(theta = rnorm(10))
#'     },
#'     update_Z = function(state, data) {
#'         # Update latent positions
#'         state
#'     },
#'     update_Theta = function(state, data) {
#'         # Update parameters
#'         state$theta <- state$theta + rnorm(10, 0, 0.1)
#'         state
#'     },
#'     collect = function(state) {
#'         # Collect samples for storage
#'         list(theta = state$theta)
#'     }
#' )
#' register_dbn_model("custom", custom_model)
#' }
#' @keywords internal
register_dbn_model <- function(name, fun_list) {
    env <- get_pkg_env()
    env$models[[name]] <- fun_list
}
