#' Family Objects for DBN Models
#'
#' @description Internal functions to create likelihood family objects for different outcome types
#' @name family
#' @keywords internal
NULL

#' Internal constructor for likelihood "family" objects
#' @keywords internal
dbn_make_family <- function(name,
                            draw_latent, # function(pre, ...) -> updated z
                            ffbs_wrapper, # function(z, mu, ...) -> theta
                            loglik, # function(a,b,theta_{t-1},theta_t, fam_pars) -> l
                            linkinv, # function(theta, misc, ...) -> expected value
                            rgen_obs, # function(theta, misc, ...) -> simulated observations
                            init_pars = list()) {
    structure(
        list(
            name = name,
            draw_latent = draw_latent,
            ffbs_wrapper = ffbs_wrapper,
            loglik = loglik,
            linkinv = linkinv,
            rgen_obs = rgen_obs,
            init_pars = init_pars
        ),
        class = "dbn_family"
    )
}

#' Ordinal Family
#'
#' @description Family object for ordinal outcomes using rank likelihood
#' @return A dbn_family object
#' @keywords internal
family_ordinal <- function() {
    dbn_make_family(
        name = "ordinal",
        draw_latent = function(pre, ...) {
            # 
            for (j in seq_len(pre$dims$p)) {
                # create ez = theta + m, broadcasting m across time
                EZ <- pre$Theta[, , j, ]
                for (t in 1:pre$dims$Tt) {
                    EZ[, , t] <- EZ[, , t] + pre$M[, , j]
                }
                pre$Z[, , j, ] <- rz_fc(pre$R[, , j, ], pre$Z[, , j, ], EZ, pre$IR[[j]])
            }
            pre
        },
        ffbs_wrapper = function(Z, mu, Aarray, Barray, sigma2_proc, sigma2_obs, ...) {
            # for ordinal outcomes, observation variance is fixed at 1 for identifiability
            # the rank likelihood model assumes z ~ n(theta + m, i), not n(theta + m, sigma2_obs * i)
            # ... just a constraint of probit-style rank models
            if (!missing(sigma2_obs) && sigma2_obs != 1) {
                warning("observation variance sigma2_obs is fixed at 1 for ordinal family (identifiability constraint)")
            }
            ffbs_theta(Z, mu, Aarray, Barray, sigma2_proc) # obs var fixed = 1
        },
        loglik = step_ll, # unchanged
        linkinv = function(theta, misc) {
            # for ordinal: return expected rank given latent value and cutpoints
            if (is.null(misc$cuts)) {
                # if no cutpoints, just return theta (latent scale)
                return(theta)
            }

            # compute probabilities for each category
            K <- length(misc$cuts) + 1 # number of cats
            cts <- seq_len(K)

            # for each element of theta, compute expected category
            expected <- array(NA, dim = dim(theta))

            for (i in seq_len(length(theta))) {
                # cumulative probabilities
                cum_probs <- c(0, pnorm(misc$cuts - theta[i]), 1)
                # category probabilities
                probs <- diff(cum_probs)
                # expected value
                expected[i] <- sum(cts * probs)
            }

            expected
        },
        rgen_obs = function(theta, misc) {
            # generate ordinal observations given latent values
            # add baseline mean if provided
            if (!is.null(misc$M)) {
                # Check dimensions
                dim_theta <- dim(theta)
                dim_M <- dim(misc$M)
                
                # m may have fewer dimensions than theta
                if (length(dim_M) == 3 && length(dim_theta) == 4) {
                    # broadcast m across time
                    theta_mean <- sweep(theta, 1:3, misc$M, "+")
                } else if (length(dim_M) == length(dim_theta) && all(dim_M == dim_theta)) {
                    # Same dimensions, can add directly
                    theta_mean <- theta + misc$M
                } else {
                    stop(sprintf("Incompatible dimensions: theta is %s, M is %s",
                                paste(dim_theta, collapse="x"),
                                paste(dim_M, collapse="x")))
                }
            } else {
                theta_mean <- theta
            }

            if (is.null(misc$cuts)) {
                # default to 5 categories with standard normal cutpoints
                misc$cuts <- qnorm(seq(0.2, 0.8, by = 0.2))
            }

            # for each latent value, determine which category it falls into
            u <- array(runif(length(theta_mean)), dim = dim(theta_mean))

            # compute cumulative probabilities and find intervals
            ranks <- array(NA, dim = dim(theta_mean))

            for (i in seq_len(length(theta_mean))) {
                cum_probs <- pnorm(misc$cuts - theta_mean[i])
                ranks[i] <- findInterval(u[i], c(0, cum_probs, 1)) + 1L
            }

            # ensure ranks are in valid range (1 to number of categories)
            n_cats <- length(misc$cuts) + 1
            ranks <- pmin(pmax(ranks, 1L), n_cats)

            ranks
        }
    )
}

#' Gaussian Family
#'
#' @description Family object for Gaussian outcomes
#' @return A dbn_family object
#' @keywords internal
family_gaussian <- function() {
    dbn_make_family(
        name = "gaussian",
        draw_latent = function(pre, ...) pre, # z ≡ y (no latent augmentation)
        ffbs_wrapper = function(Z, mu, Aarray, Barray, sigma2_proc, sigma2_obs, ...) {
            ffbs_theta_struct_cpp(Z, mu, Aarray, Barray, sigma2_proc, sigma2_obs)
        },
        loglik = function(A, B, Theta_prev, Theta_curr, fam_pars) {
            sigma2_obs <- fam_pars$sigma2_obs
            resid <- Theta_curr - A %*% Theta_prev %*% t(B)
            -0.5 * sum(resid^2) / sigma2_obs -
                0.5 * length(resid) * log(2 * pi * sigma2_obs)
        },
        linkinv = function(theta, misc, ...) {
            # for gaussian: identity link, expected value is theta itself
            theta
        },
        rgen_obs = function(theta, misc, sigma2_obs = 1) {
            # generate gaussian observations with observation variance
            # add baseline mean if provided
            if (!is.null(misc$M)) {
                # Check dimensions
                dim_theta <- dim(theta)
                dim_M <- dim(misc$M)
                
                # m may have fewer dimensions than theta
                if (length(dim_M) == 3 && length(dim_theta) == 4) {
                    # broadcast m across time
                    theta_mean <- sweep(theta, 1:3, misc$M, "+")
                } else if (length(dim_M) == length(dim_theta) && all(dim_M == dim_theta)) {
                    # Same dimensions, can add directly
                    theta_mean <- theta + misc$M
                } else {
                    stop(sprintf("Incompatible dimensions: theta is %s, M is %s",
                                paste(dim_theta, collapse="x"),
                                paste(dim_M, collapse="x")))
                }
            } else {
                theta_mean <- theta
            }

            if (!is.null(misc$sigma2_obs)) {
                sigma2_obs <- misc$sigma2_obs
            }
            # add gaussian noise
            theta_mean + array(rnorm(length(theta_mean), sd = sqrt(sigma2_obs)), dim(theta_mean))
        },
        init_pars = list(sigma2_obs = 1) # will be overwritten by sampler
    )
}

#' Binary Family
#'
#' @description Family object for binary outcomes using probit link
#' @return A dbn_family object
#' @keywords internal
family_binary <- function() {
    dbn_make_family(
        name = "binary",
        draw_latent = function(pre, Aarray, Barray, ...) {
            # check if truncnorm is available
            if (!requireNamespace("truncnorm", quietly = TRUE)) {
                stop("Package 'truncnorm' is required for binary outcomes. Please install it.")
            }

            # sample z_{ijt} | y, theta
            for (t in 1:pre$dims$Tt) {
                for (rel in 1:pre$dims$p) {
                    # use current latent mean: theta_t + m
                    eta_t <- pre$Theta[, , rel, t] + pre$M[, , rel]

                    Y_t <- pre$R[, , rel, t]
                    Z_t <- pre$Z[, , rel, t]

                    # probit link: sample from truncated normal
                    pos <- which(Y_t == 1)
                    neg <- which(Y_t == 0)

                    if (length(pos) > 0) {
                        Z_t[pos] <- truncnorm::rtruncnorm(length(pos),
                            a = 0, b = Inf,
                            mean = eta_t[pos], sd = 1
                        )
                    }
                    if (length(neg) > 0) {
                        Z_t[neg] <- truncnorm::rtruncnorm(length(neg),
                            a = -Inf, b = 0,
                            mean = eta_t[neg], sd = 1
                        )
                    }

                    pre$Z[, , rel, t] <- Z_t
                }
            }
            pre
        },
        ffbs_wrapper = function(Z, mu, Aarray, Barray, sigma2_proc, sigma2_obs = 1, ...) {
            # for binary outcomes, observation variance is fixed at 1 for identifiability
            # the probit model assumes z ~ n(theta + m, i) with truncation
            # ... standard constraint in probit mods to get at the scale
            if (!missing(sigma2_obs) && sigma2_obs != 1) {
                warning("observation variance sigma2_obs is fixed at 1 for binary family (probit identifiability)")
            }
            ffbs_theta(Z, mu, Aarray, Barray, sigma2_proc) # obs var fixed to 1
        },
        loglik = step_ll, # same latent-Gaussian form
        linkinv = function(theta, misc) {
            # for binary probit: return probability
            pnorm(theta)
        },
        rgen_obs = function(theta, misc) {
            # generate binary observations from probit probabilities
            # add baseline mean if provided
            if (!is.null(misc$M)) {
                # Check dimensions
                dim_theta <- dim(theta)
                dim_M <- dim(misc$M)
                
                # m may have fewer dimensions than theta
                if (length(dim_M) == 3 && length(dim_theta) == 4) {
                    # broadcast m across time
                    theta_mean <- sweep(theta, 1:3, misc$M, "+")
                } else if (length(dim_M) == length(dim_theta) && all(dim_M == dim_theta)) {
                    # same dimensions, can add directly
                    theta_mean <- theta + misc$M
                } else {
                    stop(sprintf("Incompatible dimensions: theta is %s, M is %s",
                                paste(dim_theta, collapse="x"),
                                paste(dim_M, collapse="x")))
                }
            } else {
                theta_mean <- theta
            }
            probs <- pnorm(theta_mean)
            array(rbinom(length(theta_mean), 1, probs), dim = dim(theta_mean))
        },
        init_pars = list() 
    )
}
