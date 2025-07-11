#' Low-rank helper functions
#' @name lowrank-helpers
#' @keywords internal
NULL

#' Initialize low-rank model
#' @param Y Data array
#' @param r Rank
#' @return List with U and alpha
#' @keywords internal
init_lowrank <- function(Y, r) {
    m <- dim(Y)[1]
    Tt <- dim(Y)[4]
    
    # svd of time-averaged adjacency for starting directions
    Ybar <- apply(Y, c(1, 2), mean, na.rm = TRUE)
    Ybar[is.na(Ybar)] <- 0
    sv <- svd(Ybar, nu = r, nv = r)
    U <- sv$u[, 1:r, drop = FALSE]
    alpha <- matrix(0, r, Tt)
    
    list(U = U, alpha = alpha)
}

#' Log-likelihood for U
#' @keywords internal
loglik_U <- function(U, alpha, Theta, Barray, sigma2) {
    m <- nrow(U)
    Tt <- ncol(alpha)
    ll <- 0
    
    for (t in 2:Tt) {
        alpha_t <- alpha[, t, drop = TRUE]
        if (length(alpha_t) == 1) {
            A_t <- alpha_t * U %*% t(U)
        } else {
            A_t <- U %*% diag(alpha_t) %*% t(U)
        }
        pred <- A_t %*% Theta[, , t - 1] %*% t(Barray[, , t])
        resid <- Theta[, , t] - pred
        ll <- ll - 0.5 * sum(resid^2) / sigma2
    }
    ll
}



#' DBN Low-rank Model (More accurate version ... but slow as heck)
#'
#' @description Fits DBN with low-rank sender effects
#' @param Y Data array (nodes x nodes x relations x time)
#' @param r Rank for low-rank factorization
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param ar1_alpha Use AR(1) for alpha dynamics
#' @param update_rho_alpha Update AR coefficient for alpha
#' @param ar1_B Use AR(1) for B dynamics
#' @param update_rho_B Update AR coefficient for B
#' @param seed Random seed
#' @param verbose Print progress
#' @param time_thin Save every nth time point
#' @param family Data family
#' @param previous Previous results to continue from
#' @param init Initial values
#' @param ... Additional arguments (currently unused)
#' @return List containing MCMC results
#' @export
dbn_lowrank_accurate <- function(Y, 
                        family = c("ordinal", "gaussian", "binary"),
                        r = 2,
                        nscan = 5000,
                        burn = 1000,
                        odens = 1,
                        ar1_alpha = TRUE,
                        update_rho_alpha = FALSE,
                        ar1_B = FALSE,
                        update_rho_B = FALSE,
                        seed = 6886,
                        verbose = 100,
                        time_thin = 1,
                        previous = NULL,
                        init = NULL,
                        ...) {

    # family setup
    family <- match.arg(family)
    FAM <- switch(family,
        ordinal = family_ordinal(),
        gaussian = family_gaussian(),
        binary = family_binary()
    )

    set.seed(seed)

    # shared preprocessing
    pre <- shared_preprocess(Y, family = family)
    dims <- pre$dims
    m <- dims$m
    p <- dims$p
    Tt <- dims$Tt
    
    # no field conversion needed - work with arrays directly

    # initialize or continue from previous
    if (!is.null(previous)) {
        # extract from previous run
        last_idx <- length(previous$sigma2_proc)
        sigma2_proc <- previous$sigma2_proc[last_idx]
        tau_alpha2 <- previous$tau_alpha2[last_idx]
        tau_B2 <- previous$tau_B2[last_idx]
        g2 <- previous$g2[last_idx]
        U <- previous$U[[last_idx]]
        alpha <- previous$alpha[[last_idx]]
        
        # expand alpha if needed
        if (ncol(alpha) < Tt) {
            alpha_full <- matrix(0, r, Tt)
            prev_times <- seq(1, Tt, by = previous$settings$time_thin)
            for (t in 1:Tt) {
                nearest <- which.min(abs(prev_times - t))
                alpha_full[, t] <- alpha[, nearest]
            }
            alpha <- alpha_full
        }
        
        Barray <- previous$B[[last_idx]]
        if (dim(Barray)[3] < Tt) {
            Barray_full <- array(0, dim = c(m, m, Tt))
            prev_times <- seq(1, Tt, by = previous$settings$time_thin)
            for (t in 1:Tt) {
                nearest <- which.min(abs(prev_times - t))
                Barray_full[, , t] <- Barray[, , nearest]
            }
            Barray <- Barray_full
        }
        
        rho_alpha <- if (ar1_alpha && !is.null(previous$rho_alpha)) previous$rho_alpha[last_idx] else if (ar1_alpha) 0.9 else 0
        rho_B <- if (ar1_B && !is.null(previous$rho_B)) previous$rho_B[last_idx] else if (ar1_B) 0.9 else 0
        sigma2_obs <- if (FAM$name == "gaussian" && !is.null(previous$sigma2_obs)) previous$sigma2_obs[last_idx] else FAM$init_pars$sigma2_obs %||% 1
        
        if (verbose) cli::cli_inform("Continuing from previous run")
    } else if (!is.null(init)) {
        # use provided initial values
        U <- init$U %||% init_lowrank(Y, r)$U
        alpha <- init$alpha %||% init_lowrank(Y, r)$alpha
        Barray <- init$B %||% {
            B <- array(0, dim = c(m, m, Tt))
            for (t in 1:Tt) B[, , t] <- diag(m)
            B
        }
        sigma2_proc <- init$sigma2_proc %||% 0.5
        sigma2_obs <- init$sigma2_obs %||% (FAM$init_pars$sigma2_obs %||% 1)
        tau_alpha2 <- init$tau_alpha2 %||% 1
        tau_B2 <- init$tau_B2 %||% 1
        g2 <- init$g2 %||% 1
        rho_alpha <- if (ar1_alpha) (init$rho_alpha %||% 0.9) else 0
        rho_B <- if (ar1_B) (init$rho_B %||% 0.9) else 0
        
        if (verbose) cli::cli_inform("Using provided initial values")
    } else {
        # default initialization
        lr <- init_lowrank(Y, r)
        U <- lr$U
        alpha <- lr$alpha
        tau_alpha2 <- 1
        rho_alpha <- if (ar1_alpha) 0.9 else 0
        sigma2_proc <- 0.5
        sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
        g2 <- 1
        
        Barray <- array(0, dim = c(m, m, Tt))
        for (t in 1:Tt) Barray[, , t] <- diag(m)
        tau_B2 <- 1
        rho_B <- if (ar1_B) 0.9 else 0
    }

    # mcmc setup
    n_iter <- burn + nscan
    keep_idx <- ((burn + 1):n_iter)[((burn + 1):n_iter) %% odens == 0]
    S <- length(keep_idx)
    
    # storage
    Usave <- vector("list", S)
    alphasave <- vector("list", S)
    Bsave <- vector("list", S)
    sigmasave <- numeric(S)
    tau_alpha_save <- numeric(S)
    tau_B_save <- numeric(S)
    g2save <- numeric(S)
    if (FAM$name == "gaussian") sigma2_obs_save <- numeric(S)
    if (ar1_alpha && update_rho_alpha) rho_alpha_save <- numeric(S)
    if (ar1_B && update_rho_B) rho_B_save <- numeric(S)
    
    time_keep <- seq(1, Tt, by = time_thin)
    
    # progress
    if (verbose) {
        cli::cli_progress_step("Running low-rank DBN MCMC")
        cli::cli_progress_bar("MCMC iterations", total = n_iter)
    }
    
    s <- 0
    epsilon <- 0.005
    accept_U <- 0
    
    # hyperparameter settings
    a_tau <- 10  # shape for tau priors
    b_tau <- 10  # rate for tau priors

    # pre-compute initial aarray
    Aarray <- compute_all_A_lowrank(U, alpha, Tt)
    
    # pre-flatten arrays to avoid repeated conversions
    total_slices <- p * Tt
    Z_flat <- array(0, dim = c(m, m, total_slices))
    Theta_flat <- array(0, dim = c(m, m, total_slices))
    M_flat <- array(0, dim = c(m, m, p))
    
    # vectorized copy
    for (j in 1:p) {
        M_flat[, , j] <- pre$M[, , j]
        slice_idx <- ((j-1)*Tt + 1):(j*Tt)
        for (i in seq_along(slice_idx)) {
            Z_flat[, , slice_idx[i]] <- pre$Z[, , j, i]
            Theta_flat[, , slice_idx[i]] <- pre$Theta[, , j, i]
        }
    }
    
    for (g in seq_len(n_iter)) {
        # 1. update z for ordinal/binary - temporarily use original method
        if (FAM$name == "ordinal") {
            # use z update
            pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "ordinal")
            # flatten back
            for (j in 1:p) {
                for (t in 1:Tt) {
                    Z_flat[, , (j-1)*Tt + t] <- pre$Z[, , j, t]
                }
            }
        } else if (FAM$name == "binary") {
            # use binary update
            pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "binary")
            # flatten back
            for (j in 1:p) {
                for (t in 1:Tt) {
                    Z_flat[, , (j-1)*Tt + t] <- pre$Z[, , j, t]
                }
            }
        }
        
        # 2. vectorized mu update directly on flattened arrays
        for (j in 1:p) {
            offset <- (j-1) * Tt
            mu_sum <- matrix(0, m, m)
            for (t in 1:Tt) {
                mu_sum <- mu_sum + Z_flat[, , offset + t] - Theta_flat[, , offset + t]
            }
            mu_var <- 1 / (Tt + 1 / g2)
            M_flat[, , j] <- mu_sum * mu_var / Tt + matrix(sqrt(mu_var) * rnorm(m*m), m, m)
        }
        M_ss <- sum(M_flat^2)
        g2 <- safe_rinv_gamma(10 + length(M_flat)/2, 10 + M_ss/2)
        
        # 3. batch theta update - already have flattened arrays
        Theta_flat <- ffbs_theta_all_relations(Z_flat, M_flat, Aarray, Barray,
                                              sigma2_proc, sigma2_obs, m, p, Tt)
        
        # 4. batch alpha update using block operations
        alpha <- update_alpha_batch(Theta_flat, U, Barray, 
                                   sigma2_proc, 
                                   tau_alpha2,
                                   ar1_alpha, rho_alpha, m, p, Tt, r)
        
        # update aarray with new alpha
        Aarray <- compute_all_A_lowrank(U, alpha, Tt)
        
        # 5. parallel b update
        Barray <- update_B_parallel(Theta_flat, Aarray,
                                         sigma2_proc,
                                         tau_B2,
                                         ar1_B, rho_B, m, p, Tt)
        
        # 6. update tau_alpha2
        if (ar1_alpha) {
            innov <- alpha[, -1, drop = FALSE] - rho_alpha * alpha[, -Tt, drop = FALSE]
        } else {
            innov <- alpha[, -1, drop = FALSE] - alpha[, -Tt, drop = FALSE]
        }
        innov_ss <- sum(innov^2)
        tau_alpha2 <- safe_rinv_gamma(a_tau + 0.5 * r * (Tt - 1), b_tau + 0.5 * innov_ss)
        
        # 7. update rho_alpha (optional)
        if (ar1_alpha && update_rho_alpha) {
            rho_prop <- rho_alpha + rnorm(1, 0, 0.05)
            if (abs(rho_prop) < 1) {
                ll_prop <- -sum((alpha[, -1] - rho_prop * alpha[, -Tt])^2) / (2 * tau_alpha2)
                ll_curr <- -sum((alpha[, -1] - rho_alpha * alpha[, -Tt])^2) / (2 * tau_alpha2)
                if (log(runif(1)) < ll_prop - ll_curr) {
                    rho_alpha <- rho_prop
                }
            }
        }
        
        # 8. update u using metropolis on stiefel manifold
        # adaptive update frequency based on network size
        u_update_freq <- max(5, floor(m / 10))
        if (g %% u_update_freq == 0) {
            # generate skew-symmetric proposal
            W_vec <- rnorm(m * (m - 1) / 2)
            W <- matrix(0, m, m)
            W[lower.tri(W)] <- W_vec
            W <- W - t(W)
            
            # adaptive norm capping
            norm_cap <- sqrt(m) * epsilon
            W_norm <- sqrt(sum(W^2))
            if (W_norm > norm_cap) {
                W <- W * (norm_cap / W_norm)
            }
            
            # fast cayley transform using woodbury identity
            I_m <- diag(m)
            half_W <- epsilon * W / 2
            U_prop <- U + 2 * half_W %*% solve(I_m - half_W %*% half_W) %*% U
            
            # compute theta average
            Theta_avg <- array(0, dim = c(m, m, Tt))
            if (p == 1) {
                Theta_avg <- Theta_flat[, , 1:Tt]
            } else {
                for (t in 1:Tt) {
                    for (j in 1:p) {
                        Theta_avg[, , t] <- Theta_avg[, , t] + Theta_flat[, , (j-1)*Tt + t]
                    }
                    Theta_avg[, , t] <- Theta_avg[, , t] / p
                }
            }
            
            log_acc <- loglik_U(U_prop, alpha, Theta_avg, Barray, sigma2_proc) -
                       loglik_U(U, alpha, Theta_avg, Barray, sigma2_proc)
            
            if (log(runif(1)) < log_acc) {
                U <- U_prop
                accept_U <- accept_U + 1
                # update aarray with new u
                Aarray <- compute_all_A_lowrank(U, alpha, Tt)
            }
        }
        
        # 9. update tau_b2
        if (ar1_B) {
            innovB <- Barray[, , 2:Tt] - rho_B * Barray[, , 1:(Tt - 1)]
        } else {
            innovB <- Barray[, , 2:Tt] - Barray[, , 1:(Tt - 1)]
        }
        innovB_ss <- sum(innovB^2)
        tau_B2 <- safe_rinv_gamma(a_tau + m^2 * (Tt - 1) / 2, b_tau + innovB_ss / 2)
        
        # 10. update rho_b (optional)
        if (ar1_B && update_rho_B) {
            rho_prop <- rho_B + rnorm(1, 0, 0.05)
            if (abs(rho_prop) < 1) {
                innovB_prop <- Barray[, , 2:Tt] - rho_prop * Barray[, , 1:(Tt - 1)]
                ll_prop <- -sum(innovB_prop^2) / (2 * tau_B2)
                innovB_curr <- Barray[, , 2:Tt] - rho_B * Barray[, , 1:(Tt - 1)]
                ll_curr <- -sum(innovB_curr^2) / (2 * tau_B2)
                if (log(runif(1)) < ll_prop - ll_curr) {
                    rho_B <- rho_prop
                }
            }
        }
        
        # 11. update sigma2 - batch computation using c++
        resid_ss <- compute_sigma2_lowrank_batch(Theta_flat, Aarray, Barray, m, p, Tt)
        sigma2_proc <- safe_rinv_gamma(a_tau + m^2 * p * (Tt - 1) / 2, b_tau + resid_ss / 2)
        
        # update observation variance for gaussian
        if (FAM$name == "gaussian") {
            resid_obs <- compute_gaussian_obs_residuals_flat(Z_flat, Theta_flat, M_flat, m, p, Tt)
            sigma2_obs <- safe_rinv_gamma(a_tau + m^2 * p * Tt / 2, b_tau + resid_obs / 2)
        }
        
        # 12. save
        if (g %in% keep_idx) {
            s <- s + 1
            Usave[[s]] <- U
            alphasave[[s]] <- alpha[, time_keep, drop = FALSE]
            Bsave[[s]] <- Barray[, , time_keep, drop = FALSE]
            sigmasave[s] <- sigma2_proc
            tau_alpha_save[s] <- tau_alpha2
            tau_B_save[s] <- tau_B2
            g2save[s] <- g2
            
            if (FAM$name == "gaussian") sigma2_obs_save[s] <- sigma2_obs
            if (ar1_alpha && update_rho_alpha) rho_alpha_save[s] <- rho_alpha
            if (ar1_B && update_rho_B) rho_B_save[s] <- rho_B
        }
        
        # adapt epsilon based on u update frequency
        if (g %% (u_update_freq * 20) == 0) {
            accept_rate <- accept_U / 20
            if (accept_rate < 0.2) {
                epsilon <- max(epsilon * 0.9, 1e-4)
            } else if (accept_rate > 0.3) {
                epsilon <- min(epsilon * 1.1, 0.1)
            }
            accept_U <- 0
        }
        
        if (verbose) cli::cli_progress_update()
        
        if (is.numeric(verbose) && g %% verbose == 0) {
            msg <- "iter {g}: sigma^2={round(sigma2_proc, 3)} tauA^2={round(tau_alpha2, 3)} tauB^2={round(tau_B2, 3)}"
            if (ar1_alpha && update_rho_alpha) msg <- paste0(msg, " rhoA={round(rho_alpha, 3)}")
            if (ar1_B && update_rho_B) msg <- paste0(msg, " rhoB={round(rho_B, 3)}")
            cli::cli_inform(msg)
        }
    }
    
    if (verbose) cli::cli_progress_done()
    
    # copy back to original structure only at the end - vectorized
    for (j in 1:p) {
        pre$M[, , j] <- M_flat[, , j]
        slice_idx <- ((j-1)*Tt + 1):(j*Tt)
        for (i in seq_along(slice_idx)) {
            pre$Z[, , j, i] <- Z_flat[, , slice_idx[i]]
            pre$Theta[, , j, i] <- Theta_flat[, , slice_idx[i]]
        }
    }
    
    # construct a arrays for output
    Asave <- vector("list", S)
    for (s in 1:S) {
        A_s <- array(0, dim = c(m, m, length(time_keep)))
        for (i in seq_along(time_keep)) {
            alpha_i <- alphasave[[s]][, i, drop = TRUE]
            if (length(alpha_i) == 1) {
                A_s[, , i] <- alpha_i * Usave[[s]] %*% t(Usave[[s]])
            } else {
                A_s[, , i] <- Usave[[s]] %*% diag(alpha_i) %*% t(Usave[[s]])
            }
        }
        Asave[[s]] <- A_s
    }
    
    # prepare output
    pars_df <- data.frame(
        sigma2_proc = sigmasave,
        tau_alpha2 = tau_alpha_save,
        tau_B2 = tau_B_save,
        g2 = g2save
    )
    
    if (FAM$name == "gaussian") pars_df$sigma2_obs <- sigma2_obs_save
    if (ar1_alpha && update_rho_alpha) pars_df$rho_alpha <- rho_alpha_save
    if (ar1_B && update_rho_B) pars_df$rho_B <- rho_B_save
    
    # create output structure
    out <- structure(
        list(
            draws = list(
                theta = NULL,
                z = NULL,
                pars = pars_df,
                misc = list(
                    U = Usave,
                    alpha = alphasave,
                    A = Asave,
                    B = Bsave
                )
            ),
            meta = list(
                dims = dims,
                draws = S,
                thin = odens,
                time_thin = time_thin,
                model = "lowrank"
            ),
            family = FAM,
            model = "lowrank",
            U = Usave,
            alpha = alphasave,
            A = Asave,
            B = Bsave,
            sigma2_proc = sigmasave,
            tau_alpha2 = tau_alpha_save,
            tau_B2 = tau_B_save,
            g2 = g2save,
            Y = Y,
            dims = dims,
            settings = list(
                r = r,
                nscan = nscan,
                burn = burn,
                odens = odens,
                ar1_alpha = ar1_alpha,
                ar1_B = ar1_B,
                time_thin = time_thin,
                family = family
            )
        ),
        class = "dbn"
    )
    
    # add sigma2_obs for gaussian family
    if (FAM$name == "gaussian") {
        out$sigma2_obs <- sigma2_obs_save
    }
    
    # add rho parameters if they were updated
    if (ar1_alpha && update_rho_alpha) {
        out$rho_alpha <- rho_alpha_save
    }
    if (ar1_B && update_rho_B) {
        out$rho_B <- rho_B_save
    }
    
    out
}

#' DBN Low-rank Model
#'
#' @description Fits DBN with low-rank sender effects 
#' @inheritParams dbn_lowrank_accurate
#' @param ... Additional arguments (currently unused)
#' @export
dbn_lowrank <- function(Y, 
                        family = c("ordinal", "gaussian", "binary"),
                        r = 2,
                        nscan = 5000,
                        burn = 1000,
                        odens = 1,
                        ar1_alpha = TRUE,
                        update_rho_alpha = FALSE,
                        ar1_B = FALSE,
                        update_rho_B = FALSE,
                        seed = 6886,
                        verbose = 100,
                        time_thin = 1,
                        previous = NULL,
                        init = NULL,
                        ...) {

    # family setup
    family <- match.arg(family)
    FAM <- switch(family,
        ordinal = family_ordinal(),
        gaussian = family_gaussian(),
        binary = family_binary()
    )

    set.seed(seed)

    # shared preprocessing
    pre <- shared_preprocess(Y, family = family)
    dims <- pre$dims
    m <- dims$m
    p <- dims$p
    Tt <- dims$Tt
    
    # initialize
    if (!is.null(previous)) {
        # extract from previous run
        last_idx <- length(previous$sigma2_proc)
        sigma2_proc <- previous$sigma2_proc[last_idx]
        tau_alpha2 <- previous$tau_alpha2[last_idx]
        tau_B2 <- previous$tau_B2[last_idx]
        g2 <- previous$g2[last_idx]
        U <- previous$U[[last_idx]]
        alpha <- previous$alpha[[last_idx]]
        
        # expand alpha if needed
        if (ncol(alpha) < Tt) {
            alpha_full <- matrix(0, r, Tt)
            prev_times <- seq(1, Tt, by = previous$settings$time_thin)
            for (t in 1:Tt) {
                nearest <- which.min(abs(prev_times - t))
                alpha_full[, t] <- alpha[, nearest]
            }
            alpha <- alpha_full
        }
        
        Barray <- previous$B[[last_idx]]
        if (dim(Barray)[3] < Tt) {
            Barray_full <- array(0, dim = c(m, m, Tt))
            prev_times <- seq(1, Tt, by = previous$settings$time_thin)
            for (t in 1:Tt) {
                nearest <- which.min(abs(prev_times - t))
                Barray_full[, , t] <- Barray[, , nearest]
            }
            Barray <- Barray_full
        }
        
        rho_alpha <- if (ar1_alpha && !is.null(previous$rho_alpha)) previous$rho_alpha[last_idx] else if (ar1_alpha) 0.9 else 0
        rho_B <- if (ar1_B && !is.null(previous$rho_B)) previous$rho_B[last_idx] else if (ar1_B) 0.9 else 0
        sigma2_obs <- if (FAM$name == "gaussian" && !is.null(previous$sigma2_obs)) previous$sigma2_obs[last_idx] else FAM$init_pars$sigma2_obs %||% 1
        
        if (verbose) cli::cli_inform("Continuing from previous run")
    } else if (!is.null(init)) {
        # use provided initial values
        U <- init$U %||% init_lowrank(Y, r)$U
        alpha <- init$alpha %||% init_lowrank(Y, r)$alpha
        Barray <- init$B %||% {
            B <- array(0, dim = c(m, m, Tt))
            for (t in 1:Tt) B[, , t] <- diag(m)
            B
        }
        sigma2_proc <- init$sigma2_proc %||% 0.5
        sigma2_obs <- init$sigma2_obs %||% (FAM$init_pars$sigma2_obs %||% 1)
        tau_alpha2 <- init$tau_alpha2 %||% 1
        tau_B2 <- init$tau_B2 %||% 1
        g2 <- init$g2 %||% 1
        rho_alpha <- if (ar1_alpha) (init$rho_alpha %||% 0.9) else 0
        rho_B <- if (ar1_B) (init$rho_B %||% 0.9) else 0
        
        if (verbose) cli::cli_inform("Using provided initial values")
    } else {
        # default initialization
        lr <- init_lowrank(Y, r)
        U <- lr$U
        alpha <- lr$alpha
        tau_alpha2 <- 1
        rho_alpha <- if (ar1_alpha) 0.9 else 0
        sigma2_proc <- 0.5
        sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
        g2 <- 1
        
        Barray <- array(0, dim = c(m, m, Tt))
        for (t in 1:Tt) Barray[, , t] <- diag(m)
        tau_B2 <- 1
        rho_B <- if (ar1_B) 0.9 else 0
    }

    # mcmc setup
    n_iter <- burn + nscan
    keep_idx <- ((burn + 1):n_iter)[((burn + 1):n_iter) %% odens == 0]
    S <- length(keep_idx)
    
    # storage
    Usave <- vector("list", S)
    alphasave <- vector("list", S)
    Bsave <- vector("list", S)
    sigmasave <- numeric(S)
    tau_alpha_save <- numeric(S)
    tau_B_save <- numeric(S)
    g2save <- numeric(S)
    if (FAM$name == "gaussian") sigma2_obs_save <- numeric(S)
    if (ar1_alpha && update_rho_alpha) rho_alpha_save <- numeric(S)
    if (ar1_B && update_rho_B) rho_B_save <- numeric(S)
    
    time_keep <- seq(1, Tt, by = time_thin)
    
    # progress
    if (verbose) {
        cli::cli_progress_step("Running low-rank DBN MCMC")
        cli::cli_progress_bar("MCMC iterations", total = n_iter)
    }
    
    s <- 0
    epsilon <- 0.005
    accept_U <- 0
    
    # hyperparameter settings
    a_tau <- 10
    b_tau <- 10
    
    # pre-compute initial A array
    Aarray <- compute_all_A_lowrank(U, alpha, Tt)
    
    # pre-flatten arrays ONCE to avoid repeated conversions
    total_slices <- p * Tt
    Z_flat <- array(0, dim = c(m, m, total_slices))
    Theta_flat <- array(0, dim = c(m, m, total_slices))
    M_flat <- array(0, dim = c(m, m, p))
    
    # vectorized copy
    for (j in 1:p) {
        M_flat[, , j] <- pre$M[, , j]
        for (t in 1:Tt) {
            idx <- (j-1)*Tt + t
            Z_flat[, , idx] <- pre$Z[, , j, t]
            Theta_flat[, , idx] <- pre$Theta[, , j, t]
        }
    }
    
    for (g in seq_len(n_iter)) {
        # 1. update z for ordinal/binary
        if (FAM$name == "ordinal") {
            # use z update
            pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "ordinal")
            # flatten back
            for (j in 1:p) {
                for (t in 1:Tt) {
                    Z_flat[, , (j-1)*Tt + t] <- pre$Z[, , j, t]
                }
            }
        } else if (FAM$name == "binary") {
            # use binary update
            pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "binary")
            # flatten back
            for (j in 1:p) {
                for (t in 1:Tt) {
                    Z_flat[, , (j-1)*Tt + t] <- pre$Z[, , j, t]
                }
            }
        }
        
        # 2. vectorized mu update directly on flattened arrays
        for (j in 1:p) {
            offset <- (j-1) * Tt
            mu_sum <- matrix(0, m, m)
            for (t in 1:Tt) {
                mu_sum <- mu_sum + Z_flat[, , offset + t] - Theta_flat[, , offset + t]
            }
            mu_var <- 1 / (Tt + 1 / g2)
            M_flat[, , j] <- mu_sum * mu_var / Tt + matrix(sqrt(mu_var) * rnorm(m*m), m, m)
        }
        M_ss <- sum(M_flat^2)
        g2 <- safe_rinv_gamma(10 + length(M_flat)/2, 10 + M_ss/2)
        
        # 3. theta update using blocked C++ function
        Theta_flat <- ffbs_theta_blocked(Z_flat, M_flat, Aarray, Barray,
                                        sigma2_proc, sigma2_obs, m, p, Tt)
        
        # 4. alpha update
        alpha <- update_alpha_optimized(Theta_flat, U, Barray, 
                                       sigma2_proc, 
                                       tau_alpha2,
                                       ar1_alpha, rho_alpha, m, p, Tt, r)
        
        # update Aarray with new alpha
        Aarray <- compute_all_A_lowrank(U, alpha, Tt)
        
        # 5. B update
        Barray <- update_B_parallel(Theta_flat, Aarray,
                                    sigma2_proc,
                                    tau_B2,
                                    ar1_B, rho_B, m, p, Tt)
        
        # 6. update tau_alpha2
        if (ar1_alpha) {
            innov <- alpha[, -1, drop = FALSE] - rho_alpha * alpha[, -Tt, drop = FALSE]
        } else {
            innov <- alpha[, -1, drop = FALSE] - alpha[, -Tt, drop = FALSE]
        }
        innov_ss <- sum(innov^2)
        tau_alpha2 <- safe_rinv_gamma(a_tau + 0.5 * r * (Tt - 1), b_tau + 0.5 * innov_ss)
        
        # 7. update rho_alpha (optional)
        if (ar1_alpha && update_rho_alpha) {
            rho_prop <- rho_alpha + rnorm(1, 0, 0.05)
            if (abs(rho_prop) < 1) {
                ll_prop <- -sum((alpha[, -1] - rho_prop * alpha[, -Tt])^2) / (2 * tau_alpha2)
                ll_curr <- -sum((alpha[, -1] - rho_alpha * alpha[, -Tt])^2) / (2 * tau_alpha2)
                if (log(runif(1)) < ll_prop - ll_curr) {
                    rho_alpha <- rho_prop
                }
            }
        }
        
        # 8. update U using Stiefel manifold update
        u_update_freq <- max(5, floor(m / 10))  # adaptive frequency
        if (g %% u_update_freq == 0) {
            # generate skew-symmetric proposal 
            W <- generate_skew_proposal(m, sqrt(m) * epsilon)
            
            # Cayley transform for Stiefel manifold
            U_prop <- cayley_transform(U, W, epsilon)
            
            # compute theta average for likelihood
            Theta_avg <- array(0, dim = c(m, m, Tt))
            if (p == 1) {
                Theta_avg <- Theta_flat[, , 1:Tt]
            } else {
                for (t in 1:Tt) {
                    for (j in 1:p) {
                        Theta_avg[, , t] <- Theta_avg[, , t] + Theta_flat[, , (j-1)*Tt + t]
                    }
                    Theta_avg[, , t] <- Theta_avg[, , t] / p
                }
            }
            
            # compute log acceptance ratio
            log_acc <- loglik_U(U_prop, alpha, Theta_avg, Barray, sigma2_proc) -
                       loglik_U(U, alpha, Theta_avg, Barray, sigma2_proc)
            
            if (log(runif(1)) < log_acc) {
                U <- U_prop
                accept_U <- accept_U + 1
                # update Aarray with new U
                Aarray <- compute_all_A_lowrank(U, alpha, Tt)
            }
        }
        
        # 9. update tau_B2
        if (ar1_B) {
            innovB <- Barray[, , 2:Tt] - rho_B * Barray[, , 1:(Tt - 1)]
        } else {
            innovB <- Barray[, , 2:Tt] - Barray[, , 1:(Tt - 1)]
        }
        innovB_ss <- sum(innovB^2)
        tau_B2 <- safe_rinv_gamma(a_tau + m^2 * (Tt - 1) / 2, b_tau + innovB_ss / 2)
        
        # 10. update rho_B (optional)
        if (ar1_B && update_rho_B) {
            rho_prop <- rho_B + rnorm(1, 0, 0.05)
            if (abs(rho_prop) < 1) {
                innovB_prop <- Barray[, , 2:Tt] - rho_prop * Barray[, , 1:(Tt - 1)]
                ll_prop <- -sum(innovB_prop^2) / (2 * tau_B2)
                innovB_curr <- Barray[, , 2:Tt] - rho_B * Barray[, , 1:(Tt - 1)]
                ll_curr <- -sum(innovB_curr^2) / (2 * tau_B2)
                if (log(runif(1)) < ll_prop - ll_curr) {
                    rho_B <- rho_prop
                }
            }
        }
        
        # 11. update sigma2 
        resid_ss <- compute_sigma2_simd(Theta_flat, U, alpha, Barray, m, p, Tt, r)
        sigma2_proc <- safe_rinv_gamma(a_tau + m^2 * p * (Tt - 1) / 2, b_tau + resid_ss / 2)
        
        # update observation variance for gaussian
        if (FAM$name == "gaussian") {
            resid_obs <- compute_gaussian_obs_residuals_flat(Z_flat, Theta_flat, M_flat, m, p, Tt)
            sigma2_obs <- safe_rinv_gamma(a_tau + m^2 * p * Tt / 2, b_tau + resid_obs / 2)
        }
        
        # 12. save
        if (g %in% keep_idx) {
            s <- s + 1
            Usave[[s]] <- U
            alphasave[[s]] <- alpha[, time_keep, drop = FALSE]
            Bsave[[s]] <- Barray[, , time_keep, drop = FALSE]
            sigmasave[s] <- sigma2_proc
            tau_alpha_save[s] <- tau_alpha2
            tau_B_save[s] <- tau_B2
            g2save[s] <- g2
            
            if (FAM$name == "gaussian") sigma2_obs_save[s] <- sigma2_obs
            if (ar1_alpha && update_rho_alpha) rho_alpha_save[s] <- rho_alpha
            if (ar1_B && update_rho_B) rho_B_save[s] <- rho_B
        }
        
        # adapt epsilon based on U update frequency
        if (g %% (u_update_freq * 20) == 0) {
            accept_rate <- accept_U / 20
            if (accept_rate < 0.2) {
                epsilon <- max(epsilon * 0.9, 1e-4)
            } else if (accept_rate > 0.3) {
                epsilon <- min(epsilon * 1.1, 0.1)
            }
            accept_U <- 0
        }
        
        if (verbose) cli::cli_progress_update()
        
        if (is.numeric(verbose) && g %% verbose == 0) {
            msg <- "iter {g}: sigma^2={round(sigma2_proc, 3)} tauA^2={round(tau_alpha2, 3)} tauB^2={round(tau_B2, 3)}"
            if (ar1_alpha && update_rho_alpha) msg <- paste0(msg, " rhoA={round(rho_alpha, 3)}")
            if (ar1_B && update_rho_B) msg <- paste0(msg, " rhoB={round(rho_B, 3)}")
            cli::cli_inform(msg)
        }
    }
    
    if (verbose) cli::cli_progress_done()
    
    # copy back to original structure only at the end - vectorized
    for (j in 1:p) {
        pre$M[, , j] <- M_flat[, , j]
        for (t in 1:Tt) {
            idx <- (j-1)*Tt + t
            pre$Z[, , j, t] <- Z_flat[, , idx]
            pre$Theta[, , j, t] <- Theta_flat[, , idx]
        }
    }
    
    # construct a arrays for output
    Asave <- vector("list", S)
    for (s in 1:S) {
        A_s <- array(0, dim = c(m, m, length(time_keep)))
        for (i in seq_along(time_keep)) {
            alpha_i <- alphasave[[s]][, i, drop = TRUE]
            A_s[, , i] <- Usave[[s]] %*% (alpha_i * t(Usave[[s]]))
        }
        Asave[[s]] <- A_s
    }
    
    # prepare output
    pars_df <- data.frame(
        sigma2_proc = sigmasave,
        tau_alpha2 = tau_alpha_save,
        tau_B2 = tau_B_save,
        g2 = g2save
    )
    
    if (FAM$name == "gaussian") pars_df$sigma2_obs <- sigma2_obs_save
    if (ar1_alpha && update_rho_alpha) pars_df$rho_alpha <- rho_alpha_save
    if (ar1_B && update_rho_B) pars_df$rho_B <- rho_B_save
    
    # create output structure
    out <- structure(
        list(
            draws = list(
                theta = NULL,
                z = NULL,
                pars = pars_df,
                misc = list(
                    U = Usave,
                    alpha = alphasave,
                    A = Asave,
                    B = Bsave
                )
            ),
            meta = list(
                dims = dims,
                draws = S,
                thin = odens,
                time_thin = time_thin,
                model = "lowrank"
            ),
            family = FAM,
            model = "lowrank",
            U = Usave,
            alpha = alphasave,
            A = Asave,
            B = Bsave,
            sigma2_proc = sigmasave,
            tau_alpha2 = tau_alpha_save,
            tau_B2 = tau_B_save,
            g2 = g2save,
            Y = Y,
            dims = dims,
            settings = list(
                r = r,
                nscan = nscan,
                burn = burn,
                odens = odens,
                ar1_alpha = ar1_alpha,
                ar1_B = ar1_B,
                time_thin = time_thin,
                family = family
            )
        ),
        class = "dbn"
    )
    
    # add sigma2_obs for gaussian family
    if (FAM$name == "gaussian") {
        out$sigma2_obs <- sigma2_obs_save
    }
    
    # add rho parameters if they were updated
    if (ar1_alpha && update_rho_alpha) {
        out$rho_alpha <- rho_alpha_save
    }
    if (ar1_B && update_rho_B) {
        out$rho_B <- rho_B_save
    }
    
    out
}
