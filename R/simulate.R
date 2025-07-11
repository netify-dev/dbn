#' Simulate Data from DBN Models
#'
#' @description Functions to simulate data from various DBN model types
#' @name simulate
NULL

#' Simulate from Static DBN Model
#'
#' @description Generates data from a static DBN model with fixed A and B matrices
#' @param n Number of actors (default: 30)
#' @param p Number of relation types (default: 2)
#' @param time Number of time points (default: 50)
#' @param sigma2 Innovation variance (default: 0.5)
#' @param tau2 Variance for A/B deviations from identity (default: 0.1)
#' @param K Number of ordinal categories (default: 5)
#' @param return_truth Whether to return true parameters in a 'truth' sub-list (default: TRUE)
#' @param seed Random seed for reproducibility
#' @return List containing:
#'   \item{Y}{Observed ordinal data array}
#'   \item{Z}{Latent continuous data array}
#'   \item{Theta}{True latent mean array}
#'   \item{A}{True A matrix}
#'   \item{B}{True B matrix}
#'   \item{M}{True baseline means}
#'   \item{cuts}{List of cutpoints used for each relation}
#'   \item{truth}{If return_truth=TRUE, contains A, B, Theta, and cuts}
#' @export
#' @examples
#' \dontrun{
#' # Generate static model data
#' data <- simulate_static_dbn(n = 30, p = 2, time = 50, seed = 123)
#' }
simulate_static_dbn <- function(n = 30, p = 2, time = 50,
                                sigma2 = 0.5, tau2 = 0.1, K = 5,
                                return_truth = TRUE, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # generate true parameters
    # a and b are close to identity with small deviations
    A <- diag(n) + matrix(rnorm(n^2, 0, sqrt(tau2)), n, n)
    B <- diag(n) + matrix(rnorm(n^2, 0, sqrt(tau2)), n, n)

    # normalize to ensure stability (target spectral radius 0.99 < 0.995 threshold)
    eig_A <- eigen(A)$values
    eig_B <- eigen(B)$values
    # safeguard against zero eigenvalues
    max_eig_A <- max(abs(eig_A))
    max_eig_B <- max(abs(eig_B))
    if (max_eig_A > 1e-10) A <- A / (max_eig_A * 1.01) else A <- diag(n)
    if (max_eig_B > 1e-10) B <- B / (max_eig_B * 1.01) else B <- diag(n)

    # baseline means for each dyad and relation
    M <- array(rnorm(n * n * p, 0, 1), dim = c(n, n, p))

    # initialize latent z
    Z <- array(NA, dim = c(n, n, p, time))

    # first time point: just baseline + noise
    for (r in 1:p) {
        Z[, , r, 1] <- M[, , r] + matrix(rnorm(n^2, 0, sqrt(sigma2)), n, n)
    }

    # generate subsequent time points
    if (time > 1) {
        for (t in 2:time) {
            for (r in 1:p) {
                # mean for time t: mu + a * (z[t-1] - mu) * b^t
                deviation <- Z[, , r, t - 1] - M[, , r]
                mean_t <- M[, , r] + A %*% deviation %*% t(B)
                Z[, , r, t] <- mean_t + matrix(rnorm(n^2, 0, sqrt(sigma2)), n, n)
            }
        }
    }

    # compute theta (latent mean structure)
    Theta <- array(NA, dim = c(n, n, p, time))
    Theta[, , , 1] <- M # first time point is just baseline
    if (time > 1) {
        for (t in 2:time) {
            for (r in 1:p) {
                deviation <- Z[, , r, t - 1] - M[, , r]
                Theta[, , r, t] <- M[, , r] + A %*% deviation %*% t(B)
            }
        }
    }

    # convert to ordinal y
    Y <- array(NA, dim = c(n, n, p, time))
    cuts <- vector("list", p)

    # use relation-specific thresholds
    for (r in 1:p) {
        # flatten z for this relation
        z_flat <- c(Z[, , r, ])

        # create thresholds based on quantiles
        probs <- seq(0, 1, length.out = K + 1)
        cuts[[r]] <- quantile(z_flat, probs = probs, na.rm = TRUE)
        cuts[[r]][1] <- -Inf
        cuts[[r]][K + 1] <- Inf

        # discretize
        for (t in 1:time) {
            Y[, , r, t] <- cut(Z[, , r, t], breaks = cuts[[r]], labels = 1:K)
        }
    }

    # add self-loops as missing
    for (t in 1:time) {
        for (r in 1:p) {
            diag(Y[, , r, t]) <- NA
        }
    }

    out <- list(
        Y = Y,
        Z = Z,
        Theta = Theta,
        A = A,
        B = B,
        M = M,
        cuts = cuts,
        sigma2 = sigma2,
        tau2 = tau2,
        K = K
    )

    if (return_truth) {
        out$truth <- list(
            A = A,
            B = B,
            Theta = Theta,
            cuts = cuts
        )
    }

    return(out)
}

#' Simulate from Dynamic DBN Model
#'
#' @description Generates data from a dynamic DBN model with time-varying A and B matrices
#' @param n Number of actors (default: 30)
#' @param p Number of relation types (default: 2)
#' @param time Number of time points (default: 50)
#' @param sigma2 Innovation variance (default: 0.5)
#' @param tauA2 Variance for A innovations (default: 0.05)
#' @param tauB2 Variance for B innovations (default: 0.05)
#' @param ar1 Use AR(1) dynamics instead of random walk (default: FALSE)
#' @param rhoA AR coefficient for A (default: 0.9)
#' @param rhoB AR coefficient for B (default: 0.9)
#' @param K Number of ordinal categories (default: 5)
#' @param return_truth Whether to return true parameters in a 'truth' sub-list (default: TRUE)
#' @param seed Random seed
#' @return List containing simulated data and true parameters
#' @export
simulate_dynamic_dbn <- function(n = 30, p = 2, time = 50,
                                 sigma2 = 0.5, tauA2 = 0.05, tauB2 = 0.05,
                                 ar1 = FALSE, rhoA = 0.9, rhoB = 0.9,
                                 K = 5, return_truth = TRUE, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # initialize time-varying a and b
    Aarray <- array(0, dim = c(n, n, time))
    Barray <- array(0, dim = c(n, n, time))

    # start from identity
    Aarray[, , 1] <- diag(n)
    Barray[, , 1] <- diag(n)

    # generate time-varying matrices
    if (time > 1) {
        for (t in 2:time) {
            if (ar1) {
                # ar(1) dynamics
                innovA <- matrix(rnorm(n^2, 0, sqrt(tauA2)), n, n)
                innovB <- matrix(rnorm(n^2, 0, sqrt(tauB2)), n, n)

                Aarray[, , t] <- rhoA * Aarray[, , t - 1] + (1 - rhoA) * diag(n) + innovA
                Barray[, , t] <- rhoB * Barray[, , t - 1] + (1 - rhoB) * diag(n) + innovB
            } else {
                # random walk
                Aarray[, , t] <- Aarray[, , t - 1] + matrix(rnorm(n^2, 0, sqrt(tauA2)), n, n)
                Barray[, , t] <- Barray[, , t - 1] + matrix(rnorm(n^2, 0, sqrt(tauB2)), n, n)
            }

            # ensure stability by normalizing (target spectral radius 0.99 < 0.995 threshold)
            eig_A <- eigen(Aarray[, , t])$values
            eig_B <- eigen(Barray[, , t])$values
            max_eig_A <- max(abs(eig_A))
            max_eig_B <- max(abs(eig_B))
            if (max_eig_A > 1e-10) Aarray[, , t] <- Aarray[, , t] / (max_eig_A * 1.01)
            if (max_eig_B > 1e-10) Barray[, , t] <- Barray[, , t] / (max_eig_B * 1.01)
        }
    }

    # baseline means
    M <- array(rnorm(n * n * p, 0, 1), dim = c(n, n, p))

    # generate latent z
    Z <- array(NA, dim = c(n, n, p, time))

    # first time point
    for (r in 1:p) {
        Z[, , r, 1] <- M[, , r] + matrix(rnorm(n^2, 0, sqrt(sigma2)), n, n)
    }

    # subsequent time points
    if (time > 1) {
        for (t in 2:time) {
            for (r in 1:p) {
                # bilinear form: a_t * (z_{t-1} - mu) * b_t^t + mu
                deviation <- Z[, , r, t - 1] - M[, , r]
                mean_t <- M[, , r] + Aarray[, , t] %*% deviation %*% t(Barray[, , t])
                Z[, , r, t] <- mean_t + matrix(rnorm(n^2, 0, sqrt(sigma2)), n, n)
            }
        }
    }

    # compute theta (latent mean structure)
    Theta <- array(NA, dim = c(n, n, p, time))
    Theta[, , , 1] <- M # first time point
    if (time > 1) {
        for (t in 2:time) {
            for (r in 1:p) {
                deviation <- Z[, , r, t - 1] - M[, , r]
                Theta[, , r, t] <- M[, , r] + Aarray[, , t] %*% deviation %*% t(Barray[, , t])
            }
        }
    }

    # convert to ordinal y
    Y <- array(NA, dim = c(n, n, p, time))
    cuts <- vector("list", p)

    for (r in 1:p) {
        z_flat <- c(Z[, , r, ])
        probs <- seq(0, 1, length.out = K + 1)
        cuts[[r]] <- quantile(z_flat, probs = probs, na.rm = TRUE)
        cuts[[r]][1] <- -Inf
        cuts[[r]][K + 1] <- Inf

        for (t in 1:time) {
            Y[, , r, t] <- cut(Z[, , r, t], breaks = cuts[[r]], labels = 1:K)
        }
    }

    # add self-loops as missing
    for (t in 1:time) {
        for (r in 1:p) {
            diag(Y[, , r, t]) <- NA
        }
    }

    out <- list(
        Y = Y,
        Z = Z,
        Theta = Theta,
        A = Aarray,
        B = Barray,
        M = M,
        cuts = cuts,
        sigma2 = sigma2,
        tauA2 = tauA2,
        tauB2 = tauB2,
        ar1 = ar1,
        rhoA = if (ar1) rhoA else NULL,
        rhoB = if (ar1) rhoB else NULL,
        K = K
    )

    if (return_truth) {
        out$truth <- list(
            A = Aarray,
            B = Barray,
            Theta = Theta,
            cuts = cuts
        )
    }

    return(out)
}

#' Simulate from Low-rank DBN Model
#'
#' @description Generates data from a low-rank DBN model
#' @param n Number of actors (default: 30)
#' @param p Number of relation types (default: 2)
#' @param time Number of time points (default: 50)
#' @param r Rank (default: 3)
#' @param sigma2 Innovation variance (default: 0.5)
#' @param tau_alpha2 Variance for alpha innovations (default: 0.1)
#' @param tauB2 Variance for B innovations (default: 0.05)
#' @param ar1_alpha Use AR(1) for alpha dynamics (default: TRUE)
#' @param rho_alpha AR coefficient for alpha (default: 0.9)
#' @param seed Random seed
#' @param return_truth Whether to return true latent factors and parameters (default: TRUE)
#' @return List containing simulated data and true parameters
#' @export
#' @examples
#' \dontrun{
#' # Generate low-rank model data
#' data <- simulate_lowrank_dbn(n = 30, p = 2, time = 50, r = 3, seed = 123)
#'
#' # Example from patch: stable low-rank simulation that converges quickly
#' sim <- simulate_lowrank_dbn(
#'     n = 25,
#'     p = 1, # single relation - most stable
#'     time = 12,
#'     r = 1, # matches the advice above
#'     sigma2 = 0.2, # high SNR
#'     tau_alpha2 = 0.02, # slow latent drift
#'     seed = 42
#' )
#'
#' fit <- dbn_lowrank(sim$Y,
#'     r          = 1,
#'     n_iter     = 600,
#'     burn       = 200,
#'     thin       = 2,
#'     ar1_alpha  = FALSE, # fewer parameters
#'     verbose    = 100
#' )
#' # returns ESS > 200 for sigma2 and tau2 in ~20 s on a laptop
#' }
simulate_lowrank_dbn <- function(
    n = 30,
    p = 2,
    time = 50,
    r = 3,
    sigma2 = 0.5,
    tau_alpha2 = 0.1,
    tauB2 = 0.05,
    ar1_alpha = TRUE,
    rho_alpha = 0.9,
    seed = NULL,
    return_truth = TRUE) {
    if (!is.null(seed)) set.seed(seed)

    ## -----  latent roles --------------------------------------------------
    U_init <- matrix(rnorm(n * r), n, r)
    U <- qr.Q(qr(U_init))

    ## -----  factor paths --------------------------------------------------
    alpha <- matrix(0, r, time)
    alpha[, 1] <- rnorm(r, 0, sqrt(tau_alpha2))
    if (time > 1) {
        for (t in 2:time) {
            innov <- rnorm(r, 0, sqrt(tau_alpha2))
            alpha[, t] <- if (ar1_alpha) {
                rho_alpha * alpha[, t - 1] + innov
            } else {
                alpha[, t - 1] + innov
            }
        }
    }

    ## -----  build a_t, b_t, theta_t ---------------------------------------
    Aarray <- array(0, c(n, n, time))
    for (t in 1:time) {
        # handle r=1 case where alpha[,t] is a scalar
        if (r == 1) {
            Aarray[, , t] <- as.numeric(alpha[, t]) * (U %*% t(U))
        } else {
            Aarray[, , t] <- U %*% diag(alpha[, t]) %*% t(U)
        }
    }

    Barray <- array(0, c(n, n, time))
    Barray[, , 1] <- diag(n)
    if (time > 1) {
        for (t in 2:time) {
            Barray[, , t] <- Barray[, , t - 1] +
                matrix(rnorm(n * n, 0, sqrt(tauB2)), n, n)
            # normalize to spectral radius 0.99
            Barray[, , t] <- Barray[, , t] /
                (max(abs(eigen(Barray[, , t])$values)) * 1.01)
        }
    }

    Theta <- array(0, c(n, n, p, time))
    for (t in 1:time) {
        for (rel in 1:p) {
            Theta[, , rel, t] <- Aarray[, , t] %*% t(Barray[, , t])
        }
    }

    ## -----  latent z and ordinal y ----------------------------------------
    Z <- Theta + array(rnorm(n * n * p * time, 0, sqrt(sigma2)),
        dim = dim(Theta)
    )

    # relation-specific cut-points (5-category default kept)
    cuts <- vector("list", p)
    Y <- array(NA_integer_, dim = c(n, n, p, time))
    for (rel in 1:p) {
        cuts[[rel]] <- quantile(c(Z[, , rel, ]), probs = seq(0, 1, length = 6))
        cuts[[rel]][1] <- -Inf
        cuts[[rel]][length(cuts[[rel]])] <- Inf
        for (t in 1:time) {
            Y[, , rel, t] <- findInterval(Z[, , rel, t], cuts[[rel]][-1]) + 1
        }
    }

    diag_inds <- rep(seq_len(n) + (seq_len(n) - 1L) * n, p * time)
    Y[diag_inds] <- NA # remove self-loops

    ## -----  output --------------------------------------------------------
    out <- list(
        Y = Y,
        Z = Z,
        Theta = Theta, # << new
        U = U,
        alpha = alpha,
        A = Aarray,
        B = Barray,
        cuts = cuts, # << new – list of length p
        r = r,
        sigma2 = sigma2,
        tau_alpha2 = tau_alpha2,
        tauB2 = tauB2,
        ar1_alpha = ar1_alpha,
        rho_alpha = if (ar1_alpha) rho_alpha else NULL
    )

    if (return_truth) {
        out$truth <- list(
            U = U, alpha = alpha, A = Aarray,
            B = Barray, Theta = Theta, cuts = cuts
        )
    }

    out
}

#' Simulate from HMM DBN Model
#'
#' @description Generates data from a regime-switching HMM DBN model
#' @param n Number of actors (default: 30)
#' @param p Number of relation types (default: 2)
#' @param time Number of time points (default: 50)
#' @param R Number of regimes (default: 3)
#' @param sigma2 Innovation variance (default: 0.5)
#' @param tau_A2 Prior variance for regime A matrices (default: 0.2)
#' @param tau_B2 Prior variance for regime B matrices (default: 0.2)
#' @param transition_prob Diagonal transition probability (default: 0.8)
#' @param stickiness Deprecated alias for transition_prob
#' @param seed Random seed
#' @param return_truth Whether to return true latent states and parameters (default: TRUE)
#' @return List containing simulated data and true parameters
#' @export
#' @examples
#' \dontrun{
#' # Generate HMM model data
#' data <- simulate_hmm_dbn(n = 20, p = 1, time = 30, R = 2, seed = 123)
#'
#' # Example from patch: stable HMM simulation that converges quickly
#' sim_hmm <- simulate_hmm_dbn(
#'     n = 20,
#'     p = 1,
#'     time = 30,
#'     R = 2, # two regimes
#'     transition_prob = 0.9, # sticky
#'     sigma2 = 0.1, # very clean data
#'     seed = 123
#' )
#'
#' fit_hmm <- dbn_hmm(sim_hmm$Y,
#'     R        = 2,
#'     n_iter   = 800,
#'     burn     = 200,
#'     thin     = 2,
#'     family   = "ordinal",
#'     delta    = c(10, 1), # sticky prior
#'     verbose  = 100
#' )
#' # Converges (R-hat < 1.05) for sigma2, tau_A2, tau_B2 in < 2 min
#' }
simulate_hmm_dbn <- function(
    n = 30,
    p = 2,
    time = 50,
    R = 3,
    sigma2 = 0.5,
    tau_A2 = 0.2,
    tau_B2 = 0.2,
    transition_prob = 0.8,
    stickiness = NULL,
    seed = NULL,
    return_truth = TRUE) {
    if (!is.null(stickiness)) {
        transition_prob <- stickiness
    } # backward-compat.

    if (!is.null(seed)) set.seed(seed)

    ## -----  regime-specific parameters ------------------------------------
    A_list <- B_list <- vector("list", R)
    for (r in 1:R) {
        A_list[[r]] <- diag(n) + matrix(rnorm(n^2, 0, sqrt(tau_A2)), n)
        B_list[[r]] <- diag(n) + matrix(rnorm(n^2, 0, sqrt(tau_B2)), n)

        # accentuate differences for pedagogy
        if (r == 1) {
            A_list[[r]] <- 1.2 * A_list[[r]]
        } else if (r == 2) B_list[[r]] <- 1.2 * B_list[[r]]

        # normalize to spectral radius 0.99
        A_list[[r]] <- A_list[[r]] / (max(abs(eigen(A_list[[r]])$values)) * 1.01)
        B_list[[r]] <- B_list[[r]] / (max(abs(eigen(B_list[[r]])$values)) * 1.01)
    }

    ## -----  markov chain ---------------------------------------------------
    Pi <- matrix((1 - transition_prob) / (R - 1), R, R)
    diag(Pi) <- transition_prob

    S <- integer(time)
    S[1] <- sample.int(R, 1)
    if (time > 1) {
        for (t in 2:time) {
            S[t] <- sample.int(R, 1, prob = Pi[S[t - 1], ])
        }
    }

    ## -----  latent mean theta ---------------------------------------------
    Theta <- array(0, c(n, n, p, time))
    for (t in 1:time) {
        for (rel in 1:p) {
            Theta[, , rel, t] <- A_list[[S[t]]] %*% t(B_list[[S[t]]])
        }
    }

    ## -----  latent z and ordinal y ----------------------------------------
    Z <- Theta + array(rnorm(n * n * p * time, 0, sqrt(sigma2)),
        dim = dim(Theta)
    )

    cuts <- vector("list", p)
    Y <- array(NA_integer_, c(n, n, p, time))
    for (rel in 1:p) {
        cuts[[rel]] <- quantile(c(Z[, , rel, ]), probs = seq(0, 1, length = 6))
        cuts[[rel]][1] <- -Inf
        cuts[[rel]][length(cuts[[rel]])] <- Inf
        for (t in 1:time) {
            Y[, , rel, t] <- findInterval(Z[, , rel, t], cuts[[rel]][-1]) + 1
        }
    }

    diag_inds <- rep(seq_len(n) + (seq_len(n) - 1L) * n, p * time)
    Y[diag_inds] <- NA

    ## -----  output ---------------------------------------------------------
    out <- list(
        Y = Y,
        Z = Z,
        Theta = Theta, # << new
        S = S,
        A_list = A_list,
        B_list = B_list,
        Pi = Pi,
        cuts = cuts, # << new
        R = R,
        sigma2 = sigma2,
        tau_A2 = tau_A2,
        tau_B2 = tau_B2,
        transition_prob = transition_prob
    )

    if (return_truth) {
        out$truth <- list(
            S = S, A = A_list, B = B_list,
            Theta = Theta, cuts = cuts, Pi = Pi
        )
    }

    out
}

#' Simulate Simple Test Data
#'
#' @description Generates simple test data for quick testing
#' @param n Number of actors
#' @param p Number of relations
#' @param time Number of time points
#' @param seed Random seed
#' @return Array of ordinal network data
#' @export
simulate_test_data <- function(n = 10, p = 2, time = 20, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    Y <- array(NA, dim = c(n, n, p, time))

    # simple structure: different baseline for each relation
    for (t in 1:time) {
        for (r in 1:p) {
            if (r == 1) {
                # relation 1: higher values
                Y[, , r, t] <- matrix(sample(3:5, n * n, replace = TRUE, prob = c(0.3, 0.5, 0.2)), n, n)
            } else {
                # relation 2: lower values
                Y[, , r, t] <- matrix(sample(1:3, n * n, replace = TRUE, prob = c(0.5, 0.3, 0.2)), n, n)
            }
        }
    }

    # add temporal correlation
    if (time > 1) {
        for (t in 2:time) {
            for (r in 1:p) {
                # 70% of edges persist
                persist_idx <- sample(1:(n^2), round(0.7 * n^2))
                Y[, , r, t][persist_idx] <- Y[, , r, t - 1][persist_idx]
            }
        }
    }

    # remove self-loops
    for (t in 1:time) {
        for (r in 1:p) {
            diag(Y[, , r, t]) <- NA
        }
    }

    return(Y)
}