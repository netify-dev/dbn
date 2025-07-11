#' DBN HMM Model
#'
#' @description Regime-switching DBN for large-scale networks (100+ actors, 50+ time steps)
#' @param Y Data array (nodes x nodes x relations x time)
#' @param R Number of regimes
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param delta Dirichlet prior for transition matrix
#' @param seed Random seed
#' @param verbose Print progress
#' @param time_thin Save every nth time point, when 1, save all time points
#' @param family Data family (ordinal, gaussian, or binary)
#' @param previous Previous dbn_hmm results to continue from (optional)
#' @param init List of initial values: S, A_list, B_list, Pi, sigma2_proc, tau_A2, tau_B2, g2, pi0 (optional)
#' @param ... Additional arguments
#' @return List containing MCMC results
#' @export
dbn_hmm <- function(Y, 
                    family = c("ordinal", "gaussian", "binary"),
                    R = 3,
                    nscan = 5000,
                    burn = 1000,
                    odens = 1,
                    delta = rep(1, R),
                    seed = 6886,
                    verbose = 100,
                    time_thin = 1,
                    previous = NULL,
                    init = NULL,
                    ...) {
    
    # meet the fam
    family <- match.arg(family)
    FAM <- switch(family,
        ordinal = family_ordinal(),
        gaussian = family_gaussian(),
        binary = family_binary()
    )
    
    # check for edge cases in data ... because of course ...
    check_edge_cases_hmm <- function(Y) {
        has_zero_var <- TRUE
        has_inf <- FALSE
        n_inf <- 0
        
        for (rel in 1:dim(Y)[3]) {
            Y_rel <- Y[, , rel, , drop = FALSE]
            if (length(dim(Y_rel)) == 4) {
                Y_rel <- array(Y_rel, c(dim(Y)[1], dim(Y)[2], dim(Y)[4]))
            }
            edge_cases <- check_edge_cases(Y_rel)
            if (edge_cases$variance > 0) has_zero_var <- FALSE
            if (edge_cases$has_infinite) {
                has_inf <- TRUE
                n_inf <- n_inf + edge_cases$n_infinite
            }
        }
        
        list(
            variance = if (has_zero_var) 0 else 1,
            has_infinite = has_inf,
            n_infinite = n_inf,
            n_obs = prod(dim(Y))
        )
    }
    
    edge_cases <- check_edge_cases_hmm(Y)
    if (edge_cases$variance == 0) {
        cli::cli_abort(c(
            "Input data has zero variance",
            "i" = "All values in the data array are constant",
            "x" = "DBN models require variation in the data to estimate parameters"
        ))
    }
    if (edge_cases$has_infinite) {
        cli::cli_abort(c(
            "Input data contains infinite values",
            "i" = "Found {edge_cases$n_infinite} infinite values out of {edge_cases$n_obs} total observations",
            "x" = "Please check your data for invalid entries"
        ))
    }
    
    set.seed(seed)
    
    # shared preprocessing
    pre <- shared_preprocess(Y, family = family)
    dims <- pre$dims
    m <- dims$m
    p <- dims$p
    Tt <- dims$Tt
    
    # optim params based on problem size
    large_scale <- (m > 50 || Tt > 100)
    use_beam_search <- (Tt > 200 && R > 5)
    beam_width <- if(use_beam_search) min(R, max(3, R/2)) else 0
    
    # prealloc work arrays
    m2 <- m * m
    
    # init from previous or defaults
    if (!is.null(previous)) {
        # Validate previous results
        if (!inherits(previous, "dbn") || previous$model != "hmm") {
            cli::cli_abort("previous must be results from dbn_hmm()")
        }
        
        # extract from previous run
        last_idx <- length(previous$sigma2_proc)
        sigma2_proc <- previous$sigma2_proc[last_idx]
        tau_A2 <- previous$tau_A2[last_idx]
        tau_B2 <- previous$tau_B2[last_idx]
        g2 <- previous$g2[last_idx]
        S <- previous$S[[last_idx]]
        
        # expand states if needed
        if (length(S) < Tt) {
            S_full <- integer(Tt)
            prev_times <- seq(1, Tt, by = previous$settings$time_thin)
            for (t in 1:Tt) {
                nearest <- which.min(abs(prev_times - t))
                S_full[t] <- S[nearest]
            }
            S <- S_full
        }
        
        A_list <- previous$A[[last_idx]]
        B_list <- previous$B[[last_idx]]
        Pi <- previous$Pi[[last_idx]]
        pi0 <- if (!is.null(previous$pi0)) previous$pi0 else rep(1 / R, R)
        sigma2_obs <- if (FAM$name == "gaussian" && !is.null(previous$sigma2_obs)) {
            previous$sigma2_obs[last_idx]
        } else {
            FAM$init_pars$sigma2_obs %||% 1
        }
        
        if (verbose) cli::cli_inform("Continuing from previous run with sigma2_proc={round(sigma2_proc,3)}, tau_A2={round(tau_A2,3)}, tau_B2={round(tau_B2,3)}")
    } else if (!is.null(init)) {
        # use provided initial values
        S <- init$S %||% sample(R, Tt, replace = TRUE)
        Pi <- init$Pi %||% matrix(1/R, R, R)
        diag(Pi) <- diag(Pi) + 0.5
        Pi <- Pi / rowSums(Pi)
        pi0 <- init$pi0 %||% rep(1 / R, R)
        A_list <- init$A_list %||% replicate(R, diag(m) + matrix(rnorm(m2, 0, 0.05), m, m), simplify = FALSE)
        B_list <- init$B_list %||% replicate(R, diag(m) + matrix(rnorm(m2, 0, 0.05), m, m), simplify = FALSE)
        sigma2_proc <- init$sigma2_proc %||% 1
        sigma2_obs <- init$sigma2_obs %||% (FAM$init_pars$sigma2_obs %||% 1)
        tau_A2 <- init$tau_A2 %||% 1
        tau_B2 <- init$tau_B2 %||% 1
        g2 <- init$g2 %||% 1
        
        if (verbose) cli::cli_inform("Using initial values: sigma2_proc={round(sigma2_proc,3)}, tau_A2={round(tau_A2,3)}, tau_B2={round(tau_B2,3)}")
    } else {
        # smart initialization based on problem size
        pi0 <- rep(1 / R, R)
        
        if (Tt >= R * 2) {
            # k means for the win ... in future iterations try out alternative initialization algos
            n_init <- min(20, Tt)
            Y_early <- apply(Y[, , , 1:n_init, drop = FALSE], 4, function(x) c(x))
            Y_early_clean <- Y_early
            
            # replace NA with row means
            for (i in 1:nrow(Y_early_clean)) {
                na_idx <- is.na(Y_early_clean[i, ])
                if (any(na_idx)) {
                    row_mean <- mean(Y_early_clean[i, !na_idx])
                    Y_early_clean[i, na_idx] <- if (is.finite(row_mean)) row_mean else 0
                }
            }
            
            if (any(is.finite(Y_early_clean))) {
                km <- tryCatch({
                    kmeans(t(Y_early_clean), centers = R, nstart = 10, iter.max = 50)
                }, error = function(e) {
                    list(cluster = sample(R, ncol(Y_early_clean), replace = TRUE))
                })
                
                S <- integer(Tt)
                S[1:n_init] <- km$cluster
                
                # propagate with persistence
                if (Tt > n_init) {
                    for (t in (n_init + 1):Tt) {
                        if (runif(1) < 0.85) {
                            S[t] <- S[t-1]
                        } else {
                            S[t] <- sample(R, 1)
                        }
                    }
                }
            } else {
                S <- sample(R, Tt, replace = TRUE)
            }
        } else {
            S <- sample(R, Tt, replace = TRUE)
        }
        
        # initialize transition matrix with persistence
        Pi <- matrix(delta[1] / sum(delta), R, R)
        diag(Pi) <- diag(Pi) + 2
        Pi <- Pi / rowSums(Pi)
        
        # initialize A and B matrices
        A_list <- vector("list", R)
        B_list <- vector("list", R)
        
        # data-driven initialization scale
        Y_vals <- as.vector(Y[!is.na(Y)])
        Y_scale <- if(length(Y_vals) > 0) sd(Y_vals) else 1
        noise_scale <- min(0.1, Y_scale / (5 * sqrt(m)))
        
        for (r in 1:R) {
            A_list[[r]] <- diag(m) + matrix(rnorm(m2, 0, noise_scale), m, m)
            B_list[[r]] <- diag(m) + matrix(rnorm(m2, 0, noise_scale), m, m)
            
            # ensure stability
            A_list[[r]] <- stabilize_spectral_radius(A_list[[r]], 0.95)
            B_list[[r]] <- stabilize_spectral_radius(B_list[[r]], 0.95)
        }
        
        # variance initialization based on data scale
        tau_A2 <- tau_B2 <- 1
        sigma2_proc <- max(0.01, Y_scale^2 / 100)
        sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
        g2 <- 1
    }
    
    # MCMC setup
    n_iter <- burn + nscan
    keep_idx <- ((burn + 1):n_iter)[((burn + 1):n_iter) %% odens == 0]
    n_keep <- length(keep_idx)
    
    # set up trashcans
    Ssave <- vector("list", n_keep)
    Asave <- vector("list", n_keep)
    Bsave <- vector("list", n_keep)
    Pisave <- vector("list", n_keep)
    sigmasave <- numeric(n_keep)
    tau_A_save <- numeric(n_keep)
    tau_B_save <- numeric(n_keep)
    g2save <- numeric(n_keep)
    if (FAM$name == "gaussian") sigma2_obs_save <- numeric(n_keep)
    
    # 
    Msave <- vector("list", n_keep)
    Thetasave <- if(time_thin > 1 || m > 100) NULL else vector("list", n_keep)
    if (FAM$name %in% c("ordinal", "binary")) {
        Zsave <- if(time_thin > 1 || m > 100) NULL else vector("list", n_keep)
    }
    
    # pick what time points to keep
    time_keep <- seq(1, Tt, by = time_thin)
    
    # fancy progress bar
    if (verbose) {
        cli::cli_progress_step("Running HMM DBN MCMC")
        cli::cli_progress_bar("MCMC iterations", total = n_iter)
    }
    
    s <- 0
    
    # precompute constants
    I_m <- diag(m)
    
    # preallocate work arrays to avoid repeated allocation
    Theta_avg <- array(0, dim = c(m, m, Tt))
    regime_counts <- matrix(0, n_iter, R)
    
    # cache for frequently used computations
    log_2pi_sigma2 <- -0.5 * m2 * log(2 * pi)
    
    # 
    for (g in seq_len(n_iter)) {
        # 1. update latent variables Z and M
        regime_arrays <- build_regime_arrays(S, A_list, B_list, m, Tt)
        Aarray <- regime_arrays$Aarray
        Barray <- regime_arrays$Barray
        
        # update M
        mu_var <- 1 / (pre$dims$Tt + 1 / g2)
        mu_hat <- mu_var * apply(pre$Z - pre$Theta, c(1, 2, 3), sum)
        pre$M <- mu_hat + sqrt(mu_var) * rsan(dim(pre$M))
        
        # update Z using
        if (FAM$name == "ordinal") {
            pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "ordinal")
        } else if (FAM$name == "binary") {
            # use binary update
            pre$Z <- update_Z_optimized(pre$R, pre$Z, pre$Theta, pre$M, pre$IR, family = "binary")
        }
        
        # 2. update theta using FFBS
        for (rel in 1:p) {
            pre$Theta[, , rel, ] <- FAM$ffbs_wrapper(
                pre$Z[, , rel, ], pre$M[, , rel],
                Aarray, Barray, 
                sigma2_proc = sigma2_proc,
                sigma2_obs = sigma2_obs
            )
        }
        
        # average theta over relations for state updates
        Theta_avg <- apply(pre$Theta, c(1, 2, 4), mean)
        
        # 3. forward-backward sample S_{1:T}
        if (use_beam_search) {
            # use beam search for very long sequences
            log_alpha <- forward_hmm_fast(Theta_avg, A_list, B_list, Pi, sigma2_proc, pi0, beam_width)
            S <- backward_sample_fast(log_alpha, Pi)
        } else {
            # standard forward-backward
            log_alpha <- forward_hmm(Theta_avg, A_list, B_list, Pi, sigma2_proc, pi0)
            S <- backward_sample(log_alpha, Pi)
        }
        
        # track regime occupancy
        for (r in 1:R) {
            regime_counts[g, r] <- sum(S == r)
        }
        
        # update initial state probabilities (optional)
        if (g > burn) {
            pi0 <- (pi0 * (g - burn - 1) + as.numeric(S[1] == 1:R)) / (g - burn)
        }
        
        # 4. update transition matrix Pi
        n_ij <- count_transitions(S, R)
        # add diagonal pseudo-counts to prevent empty rows
        n_ij <- n_ij + diag(R)
        for (i in 1:R) {
            Pi[i, ] <- rdirichlet(delta + n_ij[i, ])
        }
        
        # 5. update regime-specific (A, B) parameters
        for (r in 1:R) {
            regime_data <- collect_regime_thetas(Theta_avg, S, r, m)
            if (regime_data$n_obs > 0) {
                # use C++ update and hope ...
                upd <- update_AB_static_cpp(
                    regime_data$Th_prev, regime_data$Th_curr,
                    B_list[[r]], tau_A2, tau_B2, sigma2_proc
                )
                A_list[[r]] <- upd$A
                B_list[[r]] <- upd$B
            } else {
                # no observations for this regime - use prior
                A_list[[r]] <- diag(m) + matrix(rnorm(m2, 0, sqrt(tau_A2 / 10)), m, m)
                B_list[[r]] <- diag(m) + matrix(rnorm(m2, 0, sqrt(tau_B2 / 10)), m, m)
            }
        }
        
        # 6. update hyperparameters (pooled across regimes)
        # tau_A2
        A_resid <- compute_regime_residuals(A_list, I_m, R, m)
        n_active <- sum(regime_counts[g, ] > 0)
        shape0 <- 10 + m2
        rate0 <- 10
        tau_A2 <- safe_rinv_gamma(shape0 + m2 * n_active / 2, rate0 + A_resid / 2)
        
        # tau_B2
        B_resid <- compute_regime_residuals(B_list, I_m, R, m)
        tau_B2 <- safe_rinv_gamma(shape0 + m2 * n_active / 2, rate0 + B_resid / 2)
        
        # 7. update sigma2_proc using residual computation
        if (large_scale) {
            # use bil resids for large problems
            resid <- compute_bilinear_residuals_fast(pre$Theta, Aarray, Barray, m, p, Tt)
        } else {
            # flatten theta for standard computation
            Theta_flat <- array(pre$Theta, dim = c(m, m, p * Tt))
            resid <- compute_bilinear_residuals(Theta_flat, Aarray, Barray, m, p, Tt)
        }
        sigma2_proc <- safe_rinv_gamma(shape0 + m2 * p * (Tt - 1) / 2, rate0 + resid / 2)
        
        # update observation variance for gaussian family
        if (FAM$name == "gaussian") {
            resid_obs <- compute_gaussian_obs_residuals(pre$Z, pre$Theta, pre$M)
            sigma2_obs <- safe_rinv_gamma(shape0 + length(pre$Z) / 2, rate0 + resid_obs / 2)
        }
        
        # 8. Update g2
        g2 <- safe_rinv_gamma(shape0 + prod(dim(pre$M)) / 2, (rate0 + sum(pre$M^2)) / 2)
        
        # 9. store
        if (g %in% keep_idx) {
            s <- s + 1
            Ssave[[s]] <- S[time_keep]
            Asave[[s]] <- A_list
            Bsave[[s]] <- B_list
            Pisave[[s]] <- Pi
            sigmasave[s] <- sigma2_proc
            tau_A_save[s] <- tau_A2
            tau_B_save[s] <- tau_B2
            g2save[s] <- g2
            
            if (FAM$name == "gaussian") {
                sigma2_obs_save[s] <- sigma2_obs
            }
            
            # save M and optionally Theta/Z
            Msave[[s]] <- pre$M
            if (!is.null(Thetasave)) {
                if (time_thin > 1) {
                    Thetasave[[s]] <- pre$Theta[, , , time_keep, drop = FALSE]
                } else {
                    Thetasave[[s]] <- pre$Theta
                }
            }
            
            if (FAM$name %in% c("ordinal", "binary") && !is.null(Zsave)) {
                if (time_thin > 1) {
                    Zsave[[s]] <- pre$Z[, , , time_keep, drop = FALSE]
                } else {
                    Zsave[[s]] <- pre$Z
                }
            }
        }
        
        if (verbose) cli::cli_progress_update()
        
        # 
        if (is.numeric(verbose) && g %% verbose == 0) {
            regime_table <- table(factor(S, levels = 1:R))
            msg <- "iter {g}: sigma^2={round(sigma2_proc, 3)} tauA^2={round(tau_A2, 3)} tauB^2={round(tau_B2, 3)}"
            msg <- paste0(msg, " regimes=", paste(regime_table, collapse = "/"))
            cli::cli_inform(msg)
        }
    }
    
    if (verbose) cli::cli_progress_done()
    
    # report regime occupancy diagnostics if verbose
    if (verbose) {
        cli::cli_h3("Regime occupancy summary")
        regime_props <- colMeans(regime_counts[(burn + 1):n_iter, ]) / Tt
        for (r in 1:R) {
            cli::cli_inform("  Regime {r}: {round(regime_props[r] * 100, 1)}% ({round(regime_props[r] * Tt)} time points on average)")
        }
        
        # warn if any regime is rarely occupied
        if (any(regime_props < 0.05)) {
            cli::cli_warn("Some regimes are rarely occupied (<5%). Consider reducing R.")
        }
    }
    
    # convert A and B lists to arrays for compatibility
    Asave_array <- vector("list", n_keep)
    Bsave_array <- vector("list", n_keep)
    
    for (s in 1:n_keep) {
        A_s <- array(0, dim = c(m, m, R))
        B_s <- array(0, dim = c(m, m, R))
        for (r in 1:R) {
            A_s[, , r] <- Asave[[s]][[r]]
            B_s[, , r] <- Bsave[[s]][[r]]
        }
        Asave_array[[s]] <- A_s
        Bsave_array[[s]] <- B_s
    }
    
    # prepare scalar parameters data frame
    pars_df <- data.frame(
        sigma2_proc = sigmasave,
        tau_A2 = tau_A_save,
        tau_B2 = tau_B_save,
        g2 = g2save
    )
    
    if (FAM$name == "gaussian") {
        pars_df$sigma2_obs <- sigma2_obs_save
    }
    
    # 
    draws <- list(
        theta = Thetasave,
        z = if (FAM$name %in% c("ordinal", "binary")) Zsave else NULL,
        pars = pars_df,
        misc = list(
            M = Msave,
            S = Ssave,
            A = Asave_array,
            B = Bsave_array,
            Pi = Pisave
        )
    )
    
    # 
    meta <- list(
        dims = dims,
        draws = n_keep,
        thin = odens,
        time_thin = time_thin,
        model = "hmm",
        R = R,
        pi0 = pi0
    )
    
    # 
    out <- list(
        # 
        draws = draws,
        meta = meta,
        family = FAM,
        
        # 
        model = "hmm",
        S = Ssave,
        A = Asave_array,
        B = Bsave_array,
        Pi = Pisave,
        pi0 = pi0,
        sigma2_proc = sigmasave,
        tau_A2 = tau_A_save,
        tau_B2 = tau_B_save,
        g2 = g2save,
        Y = Y,
        n_iter = n_iter,
        burn = burn,
        thin = odens,
        dims = dims,
        settings = list(
            R = R,
            nscan = nscan,
            burn = burn,
            odens = odens,
            delta = delta,
            time_thin = time_thin,
            family = family
        ),
        regime_counts = regime_counts
    )
    
    if (FAM$name == "gaussian") {
        out$sigma2_obs <- sigma2_obs_save
    }
    
    # if continuing from previous, append iteration counts
    if (!is.null(previous)) {
        prev_iter <- if (!is.null(previous$total_iter)) {
            previous$total_iter
        } else if (!is.null(previous$settings$nscan)) {
            previous$settings$burn + previous$settings$nscan
        } else {
            previous$n_iter
        }
        out$total_iter <- prev_iter + nscan
        out$continued_from <- prev_iter
    }
    
    class(out) <- "dbn"
    out
}

#' Step Log-likelihood
#'
#' @description Computes log-likelihood for one time step given regime
#' @param A Sender effects matrix
#' @param B Receiver effects matrix
#' @param Theta_prev Previous Theta
#' @param Theta_curr Current Theta
#' @param sigma2_proc Innovation variance
#' @return Log-likelihood value
#' @keywords internal
step_ll <- function(A, B, Theta_prev, Theta_curr, sigma2_proc) {
    # guard against sigma2_proc -> 0 blowup
    sigma2_proc <- max(sigma2_proc, 1e-6)
    
    # bilinear form: a theta_{t-1} b^t
    resid <- Theta_curr - bilinear_step(A, Theta_prev, B)
    
    # hopefully a more numerically stable log-likelihood computation
    -0.5 * sum(resid^2) / sigma2_proc - 0.5 * length(resid) * log(2 * pi * sigma2_proc)
}

#' Sample Dirichlet
#'
#' @description Samples from Dirichlet distribution
#' @param alpha Parameter vector
#' @return Sample vector
#' @keywords internal
rdirichlet <- function(alpha) {
    x <- rgamma(length(alpha), alpha, 1)
    x / sum(x)
}