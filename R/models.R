#' Dynamic Bilinear Network Analysis
#'
#' @description
#' Main wrapper function for Dynamic Bilinear Network (DBN) analysis. This is the primary
#' interface for fitting DBN models to network data. DBN models capture complex dependencies
#' in network data through bilinear interactions between latent sender and receiver effects.
#'
#' @param data Data array or path to .RData file containing Y array.
#'   Array should be 3-dimensional (actors x actors x time) for single relation 
#'   or 4-dimensional (actors x actors x relations x time) for multiple relations
#' @param family Character string specifying the data family/distribution:
#'   - "ordinal": For ordinal/ranked data (e.g., ratings 1-5). Data should be positive integers.
#'   - "gaussian": For continuous data. Data can be any real numbers.
#'   - "binary": For binary data. Data should be 0/1 or logical values.
#' @param model Character string specifying model type:
#'   - "static": Fixed sender/receiver effects across time
#'   - "dynamic": Time-varying sender/receiver effects
#'   - "lowrank": Low-rank factorization of sender effects
#'   - "hmm": Regime-switching model with hidden Markov states
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain (save every odens-th iteration)
#' @param verbose Logical or numeric. If TRUE, show progress. If numeric, print detailed info every n iterations (default: TRUE)
#' @param ... Additional model-specific parameters passed to model-specific functions
#' @return A list of class "dbn" containing:
#'   \item{B}{List of posterior samples for B matrices (static model)}
#'   \item{A}{List of posterior samples for time-varying A matrices (dynamic model)}
#'   \item{params}{Matrix of parameter traces (static model)}
#'   \item{sigma2}{Vector of sigma^2 samples (dynamic model)}
#'   \item{model}{Character string indicating which model was run}
#'   \item{dims}{List containing data dimensions}
#'   \item{settings}{List of model settings used}
#' @export
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_data)
#'
#' # Run static model with default settings
#' results <- dbn(example_data, model = "static")
#'
#' # Run dynamic model with custom MCMC settings
#' results <- dbn(example_data,
#'     model = "dynamic",
#'     nscan = 5000, burn = 1000, odens = 10
#' )
#'
#' # Run HMM model with 3 regimes
#' results <- dbn(example_data, model = "hmm", R = 3)
#'
#' # Run low-rank model with rank 2
#' results <- dbn(example_data, model = "lowrank", r = 2)
#'
#' # Run quietly without progress output
#' results <- dbn(example_data, model = "static", verbose = FALSE)
#'
#' # Run with detailed output every 100 iterations
#' results <- dbn(example_data, model = "dynamic", verbose = 100)
#' }
dbn <- function(data,
                family = c("ordinal", "gaussian", "binary"),
                model = c("static", "dynamic", "lowrank", "hmm"),
                nscan = 10000,
                burn = 1000,
                odens = 1,
                verbose = TRUE,
                ...) {
    family <- match.arg(family)
    model <- match.arg(model)

    # load data from file path if needed
    if (is.character(data)) {
        cli::cli_inform("Loading data from: {.path {data}}")
        env <- new.env()
        load(data, envir = env)
        Y <- env$Y
        if (is.null(Y)) cli::cli_abort("Data file must contain object {.var Y}")
    } else {
        Y <- data
    }

    # check that data has the right shape and convert 3D to 4D if needed
    if (length(dim(Y)) == 3) {
        # convert 3D array (actors x actors x time) to 4D with single relation
        dim_orig <- dim(Y)
        Y <- array(Y, dim = c(dim_orig[1], dim_orig[2], 1, dim_orig[3]))
        cli::cli_inform("Converting 3D array to 4D array with single relation")
    } else if (length(dim(Y)) != 4) {
        cli::cli_abort("Data must be a 3D array [actors x actors x time] or 4D array [actors x actors x relations x time]")
    }

    # make sure mcmc params are sensible
    if (nscan <= 0) {
        cli::cli_abort("nscan must be positive")
    }
    if (burn < 0) {
        cli::cli_abort("burn must be non-negative")
    }
    if (odens < 1) {
        cli::cli_abort("odens must be at least 1")
    }

    # show data dimensions if requested
    if (verbose) {
        cli::cli_h3("Data dimensions")
        cli::cli_bullets(c(
            " " = "Nodes: {dim(Y)[1]}",
            " " = "Relations: {dim(Y)[3]}",
            " " = "Time points: {dim(Y)[4]}"
        ))
    }

    # dispatch to a sub function call
    results <- switch(model,
        static = dbn_static(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, ...),
        dynamic = dbn_dynamic(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, ...),
        lowrank = dbn_lowrank(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, ...),
        hmm = dbn_hmm(Y, family = family, nscan = nscan, burn = burn, odens = odens, verbose = verbose, ...)
    )

    # tag results with model info and class
    results$model <- model
    class(results) <- "dbn"

    return(results)
}

#' Static DBN MCMC
#'
#' @description Fits DBN model with fixed sender/receiver effects
#' @param Y Data array (nodes x nodes x relations x time)
#' @param family Character string specifying the data family/distribution:
#'   - "ordinal": Ordinal data (ordered categories). Data should be positive integers.
#'   - "gaussian": Continuous data with Gaussian errors. Data can be any real numbers.
#'   - "binary": Binary (0/1) data. Data should be 0/1 or logical values.
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param seed Random seed for reproducibility
#' @param verbose Logical indicating whether to show progress
#' @param previous Previous dbn_static results to continue from (optional)
#' @param init List of initial values for parameters: B, s2, t2, g2, M, Z (optional)
#' @return List containing MCMC results
#' @export
dbn_static <- function(Y, family = c("ordinal", "gaussian", "binary"), 
                       nscan = 5000, burn = 1000, odens = 1,
                       seed = 6886, verbose = 100,
                       previous = NULL, init = NULL) {
    # setup
    family <- match.arg(family)
    FAM <- switch(family,
        ordinal = family_ordinal(),
        gaussian = family_gaussian(),
        binary = family_binary()
    )
    
    set.seed(seed)
    
    # use shared preprocessing
    pre <- shared_preprocess(Y, family = family)
    Z <- pre$Z
    R <- pre$R  
    M <- pre$M
    dims <- pre$dims
    m <- dims$m
    p <- dims$p
    n <- Tt <- dims$Tt
    
    # validate input
    if (m == 1) {
        stop("Single node network has zero variance - need at least 2 nodes")
    }
    
    # determine if we're dealing with a large network
    # Lower thresholds to use code paths more often
    is_large_network <- (m > 15) || (p > 1) || (Tt > 20) || (m * m * p * Tt > 10000)
    
    # for ordinal data, compute IR using the C++ function that matches rz_fc_batch expectations
    if (family == "ordinal") {
        R_flat <- array(R, c(m, m, p * Tt))
        IR <- precompute_rank_structure(R_flat, m, p, Tt)
    } else {
        IR <- pre$IR
    }
    K <- 3
    d <- c(m, m, p)
    
    # initialize from previous or defaults
    if (!is.null(previous)) {
        # warm start
        s2 <- tail(previous$draws$pars$s2, 1)
        t2 <- tail(previous$draws$pars$t2, 1)
        g2 <- tail(previous$draws$pars$g2, 1)
        M <- previous$draws$misc$M[[length(previous$draws$misc$M)]]
        B <- lapply(previous$draws$misc$B, function(b) b[, , dim(b)[3]])
    } else {
        s2 <- 1
        t2 <- 1
        g2 <- max(0.1, mean(M^2, na.rm = TRUE))
        B <- list(diag(m), diag(m), diag(p))
    }
    
    # total iterations
    n_iter <- burn + nscan
    keep_idx <- seq(burn + 1, n_iter, by = odens)
    n_keep <- length(keep_idx)
    
    # storage
    B_samples <- list()
    for (k in 1:K) B_samples[[k]] <- array(NA, c(d[k], d[k], n_keep))
    param_samples <- matrix(NA, n_keep, 3)
    colnames(param_samples) <- c("s2", "t2", "g2")
    Msave <- vector("list", n_keep)
    if (FAM$name %in% c("ordinal", "binary")) {
        Zsave <- vector("list", n_keep)
    }
    
    # flatten arrays for c++
    Z_flat <- array(Z, c(m, m, p * Tt))
    
    # progress
    if (verbose) {
        cli::cli_progress_step("Running static DBN MCMC")
        cli::cli_progress_bar("MCMC iterations", total = n_iter)
    }
    
    # prepare data for rcpp
    # use parallel version for large networks
    if (is_large_network) {
        Z_cube <- reshape_Z_to_cube_parallel(Z, m, p, Tt)
    } else {
        Z_cube <- reshape_Z_to_cube(Z, m, p, Tt)
    }
    
    # ensure M is always a proper 3D numeric array for C++
    if (!is.array(M) || length(dim(M)) != 3) {
        M <- array(as.numeric(M), dim = c(m, m, p))
    }
    # ensure both are numeric for C++
    storage.mode(Z_cube) <- "double"
    storage.mode(M) <- "double"
    
    # preallocate workspace for large networks
    if (is_large_network && family == "ordinal") {
        # precompute constants 
        Z_flat_mat <- matrix(0, nrow = m * m, ncol = p * Tt)
        EZ_flat_mat <- matrix(0, nrow = m * m, ncol = p * Tt)
    }
    
    for (iter in 1:n_iter) {
        # choose tiled version for large networks
        if (is_large_network) {
            B[[1]] <- update_B_static_tiled(Z_cube, M, s2, t2, m, p, Tt)
        } else {
            B[[1]] <- update_B_static(Z_cube, M, s2, t2, m, p, Tt)
        }
        
        # for static model, B[[2]] and B[[3]] remain as identity matrices
        # (they are not updated in the static model)
        
        # update t2
        sse <- compute_diagonal_sse(B, K)
        t2 <- safe_rinv_gamma((sum(d) + 1) / 2, (sse + 1) / 2)
        
        # z update for ordinal
        if (FAM$name == "ordinal") {
            # check if we should use gaussian approximation for speed
            use_approx <- should_use_gaussian_approximation(R) || 
                         (m * m * p * Tt > 5000)  # also use for large problems
            
            Z_flat <- array(Z, c(m, m, p * Tt))
            R_flat <- array(R, c(m, m, p * Tt))
            EZ_cube <- broadcast_M_and_compute_EZ(M, s2, m, p, Tt)
            
            if (use_approx) {
                # use fast gaussian approximation (C++ version if available)
                if (exists("rz_gaussian_approx_cpp", mode = "function")) {
                    Z_flat <- rz_gaussian_approx_cpp(R_flat, Z_flat, EZ_cube)
                } else {
                    Z_flat <- rz_gaussian_approx(R_flat, Z_flat, EZ_cube)
                }
            } else {
                # use exact truncated normal sampling
                Z_flat <- rz_fc_batch(R_flat, Z_flat, EZ_cube, IR, m, p, Tt)
            }
            
            # reshape back
            Z <- array(Z_flat, c(m, m, p, Tt))
            
            # update Z_cube for next B update
            if (is_large_network) {
                Z_cube <- reshape_Z_to_cube_parallel(Z, m, p, Tt)
            } else {
                Z_cube <- reshape_Z_to_cube(Z, m, p, Tt)
            }
        }
        
        # M update 
        if (FAM$name == "ordinal" || iter %% 5 == 0) {
            # already have Z_flat from ordinal update or update periodically
            if (FAM$name != "ordinal") {
                Z_flat <- array(Z, c(m, m, p * Tt))
            }
            # reshape Z_flat to matrix for C++ function
            Z_flat_mat <- matrix(Z_flat, nrow = m * m, ncol = p * Tt)
            
            # 
            if (is_large_network) {
                M <- compute_M_static_blocked(Z_flat_mat, m, p, Tt)
            } else {
                M <- compute_M_static(Z_flat_mat, m, p, Tt)
            }
            
            # update g2
            M_sum_sq <- sum(M^2)
            g2 <- (1 + M_sum_sq) / (2 * rgamma(1, shape = (1 + m*m*p)/2, rate = 1))
        }
        
        # s2 update
        if (FAM$name == "gaussian") {
            # compute residual sum of squares 
            if (is_large_network) {
                rss <- compute_rss_static_parallel(Z_cube, M, m, p, Tt)
            } else {
                rss <- compute_rss_static(Z, M, m, p, Tt)
            }
            s2 <- safe_rinv_gamma(1 + m * m * p * Tt / 2, 1 + rss / 2)
        } else {
            # for ordinal/binary, observation variance is fixed
            s2 <- 1
        }
        
        # save samples
        if (iter %in% keep_idx) {
            idx <- which(keep_idx == iter)
            # for static model, we only have one B matrix
            B_samples[[1]][, , idx] <- B[[1]]
            # B[[2]] and B[[3]] are identity matrices for static model
            if (K > 1) {
                if (length(B) >= 2 && !is.null(B[[2]])) {
                    B_samples[[2]][, , idx] <- B[[2]]
                } else {
                    B_samples[[2]][, , idx] <- diag(d[2])
                }
            }
            if (K > 2) {
                if (length(B) >= 3 && !is.null(B[[3]])) {
                    B_samples[[3]][, , idx] <- B[[3]]
                } else {
                    # ensure proper dimensions for small p
                    id_mat <- diag(d[3])
                    if (d[3] == 1) {
                        id_mat <- matrix(1, 1, 1)
                    }
                    B_samples[[3]][, , idx] <- id_mat
                }
            }
            param_samples[idx, ] <- c(s2, t2, g2)
            Msave[[idx]] <- M
            if (FAM$name %in% c("ordinal", "binary")) {
                Zsave[[idx]] <- Z
            }
        }
        
        if (verbose && iter %% verbose == 0) {
            cli::cli_progress_update()
        }
    }
    
    if (verbose) cli::cli_progress_done()
    
    # prepare output in standard format
    draws <- list(
        theta = NULL,
        z = if (FAM$name %in% c("ordinal", "binary")) Zsave else NULL,
        pars = data.frame(
            s2 = param_samples[, "s2"],
            t2 = param_samples[, "t2"],
            g2 = param_samples[, "g2"]
        ),
        misc = list(
            M = Msave,
            B = B_samples
        )
    )
    
    # calculate dic
    deviance <- -2 * sum(dnorm(as.vector(Y), 0, sqrt(param_samples[, "s2"]), log = TRUE))
    pd <- var(param_samples[, "s2"]) / 2
    dic <- deviance + 2 * pd
    
    out <- list(
        model = "static",
        family = family,
        Y = Y,  # Store Y for posterior predictions
        R = R,
        dims = list(m = m, p = p, n = n),
        settings = list(
            nscan = nscan,
            burn = burn,
            odens = odens,
            draws = n_keep
        ),
        meta = list(
            family = family,
            dims = list(m = m, p = p, Tt = Tt),
            draws = n_keep,
            settings = list(nscan = nscan, burn = burn, odens = odens)
        ),
        params = param_samples,
        M = M,
        B = B_samples,
        draws = draws,
        diagnostics = list(
            deviance = deviance,
            pd = pd,
            dic = dic
        ),
        # final states for warm start
        Z_final = if (FAM$name %in% c("ordinal", "binary")) Z else NULL,
        M_final = M
    )
    
    if (!is.null(previous)) {
        prev_total <- previous$total_iter %||% (previous$settings$burn + previous$settings$nscan)
        out$total_iter <- prev_total + nscan
        out$continued_from <- prev_total
    }
    
    class(out) <- "dbn"
    return(out)
}


#' Dynamic DBN MCMC
#'
#' @description Fits DBN model with time-varying sender/receiver effects
#' @param Y Data array (nodes x nodes x relations x time)
#' @param family Character string specifying the data family/distribution:
#'   - "ordinal": Ordinal data (ordered categories). Data should be positive integers.
#'   - "gaussian": Continuous data with Gaussian errors. Data can be any real numbers.
#'   - "binary": Binary (0/1) data. Data should be 0/1 or logical values.
#' @param nscan Number of iterations of the Markov chain (beyond burn-in)
#' @param burn Burn-in for the Markov chain
#' @param odens Output density for the Markov chain
#' @param ar1 Use AR(1) dynamics instead of random walk (default: FALSE)
#' @param update_rho Update AR coefficients (default: FALSE)
#' @param seed Random seed
#' @param verbose Print progress every n iterations (default: 100)
#' @param time_thin Save every nth time point to reduce memory (default: 1 = save all)
#' @param previous Previous dbn_dynamic results to continue from (optional)
#' @param init List of initial values: A, B, sigma2, tau_A2, tau_B2, g2, rho_A, rho_B, Theta, M, Z (optional)
#' @return List containing MCMC results
#' @export
dbn_dynamic <- function(Y, 
                        family = c("ordinal", "gaussian", "binary"),
                        nscan = 5000, 
                        burn = 1000, 
                        odens = 1,
                        ar1 = FALSE, 
                        update_rho = FALSE,
                        seed = 6886, 
                        verbose = 100, 
                        time_thin = 1,
                        previous = NULL,
                        init = NULL) {
    
    set.seed(seed)
    family <- match.arg(family)
    FAM <- switch(family,
        ordinal = family_ordinal(),
        gaussian = family_gaussian(),
        binary = family_binary()
    )
    
    # validate input for dynamic model
    if (dim(Y)[4] < 2) {
        stop("Dynamic model requires at least 2 time points")
    }
    
    # preprocess data
    pre <- shared_preprocess(Y, family = family)
    Z <- pre$Z
    R <- pre$R
    IR <- pre$IR
    M <- pre$M
    dims <- pre$dims
    m <- dims$m
    p <- dims$p
    Tt <- dims$Tt
    d <- m * m
    
    # determine if we're dealing with a large network
    # Lower thresholds to use code paths more often
    is_large_network <- (m > 15) || (p > 1) || (Tt > 20) || (m * m * p * Tt > 10000)
    
    # precompute time indices for ordinal data 
    if (family == "ordinal" && !is.null(IR)) {
        IR_time_indices <- precompute_time_indices(IR, m, p, Tt)
    }
    
    # flatten arrays 
    # maintain consistent memory layout throughout
    Z_4d <- matrix(Z, nrow = m * m, ncol = p * Tt)
    R_4d <- matrix(R, nrow = m * m, ncol = p * Tt)
    
    # initialize or continue from previous
    if (!is.null(previous)) {
        # validate previous results
        if (is.null(previous$A) || is.null(previous$B) || is.null(previous$sigma2)) {
            stop("previous must be results from dbn_dynamic")
        }
        # extract from previous run
        last_idx <- length(previous$sigma2)
        
        # Handle case where Theta might not be stored (removed to save memory)
        if (!is.null(previous$Theta)) {
            # Legacy case: extract last Theta from 5D array (m, m, p, time, iter)
            n_prev_iter <- dim(previous$Theta)[5]
            last_theta_idx <- n_prev_iter
            
            # handle Theta with time thinning
            prev_time_thin <- previous$time_thin %||% 1
            n_time_stored <- dim(previous$Theta)[4]
            
            if (n_time_stored < Tt) {
            # expand Theta back to full time
            time_indices <- seq(1, Tt, by = prev_time_thin)
            Theta_all <- array(0, dim = c(m, m, p, Tt))
            
            # fill in available time points from last iteration
            for (i in 1:min(length(time_indices), n_time_stored)) {
                if (time_indices[i] <= Tt) {
                    Theta_all[,,,time_indices[i]] <- previous$Theta[,,,i,last_theta_idx]
                }
            }
            
            # interpolate missing time points
            for (t in 1:Tt) {
                if (all(Theta_all[,,,t] == 0)) {
                    available_times <- time_indices[time_indices <= Tt]
                    nearest_idx <- which.min(abs(available_times - t))
                    nearest_stored <- nearest_idx
                    if (nearest_stored <= n_time_stored) {
                        Theta_all[,,,t] <- previous$Theta[,,,nearest_stored,last_theta_idx]
                    }
                }
            }
            } else {
                # extract last iteration from 5D array
                Theta_all <- previous$Theta[,,,1:Tt,last_theta_idx]
            }
        } else {
            # Theta not stored - reconstruct from previous A, B, and M
            # We'll initialize Theta_all from scratch (will be recalculated during MCMC)
            Theta_all <- pre$Theta
        }
        
        sigma2 <- previous$sigma2[last_idx]
        sigma2_obs <- if (!is.null(previous$sigma2_obs)) previous$sigma2_obs[last_idx] else 1
        tauA2 <- previous$tau_A2[last_idx]
        tauB2 <- previous$tau_B2[last_idx]
        g2 <- previous$g2[last_idx]
        rhoA <- if (!is.null(previous$rhoA)) previous$rhoA[last_idx] else 0
        rhoB <- if (!is.null(previous$rhoB)) previous$rhoB[last_idx] else 0
        
        # extract last a and b arrays
        last_A <- previous$A[[last_idx]]
        last_B <- previous$B[[last_idx]]
        
        # handle time thinning - expand back to full time if needed
        if (dim(last_A)[3] < Tt) {
            # previous results used time thinning
            prev_time_thin <- previous$time_thin %||% 1
            time_indices <- seq(1, Tt, by = prev_time_thin)
            
            # initialize full arrays
            Aarray <- array(0, dim = c(m, m, Tt))
            Barray <- array(0, dim = c(m, m, Tt))
            
            # fill in available time points
            for (i in 1:length(time_indices)) {
                if (time_indices[i] <= Tt && i <= dim(last_A)[3]) {
                    Aarray[,,time_indices[i]] <- last_A[,,i]
                    Barray[,,time_indices[i]] <- last_B[,,i]
                }
            }
            
            # interpolate missing time points using nearest available
            for (t in 1:Tt) {
                if (all(Aarray[,,t] == 0)) {
                    # find nearest available time point
                    available_times <- time_indices[time_indices <= Tt]
                    nearest_idx <- which.min(abs(available_times - t))
                    nearest_time <- available_times[nearest_idx]
                    nearest_stored <- which(time_indices == nearest_time)
                    
                    Aarray[,,t] <- last_A[,,nearest_stored]
                    Barray[,,t] <- last_B[,,nearest_stored]
                }
            }
        } else {
            Aarray <- last_A
            Barray <- last_B
        }
    } else {
        # initialize parameters
        Theta_all <- pre$Theta
        Aarray <- Barray <- array(0, dim = c(m, m, Tt))
        for (t in 1:Tt) {
            Aarray[, , t] <- Barray[, , t] <- diag(m)
        }
        
        sigma2 <- 1
        sigma2_obs <- FAM$init_pars$sigma2_obs %||% 1
        tauA2 <- tauB2 <- 1
        g2 <- 1
        rhoA <- rhoB <- 0
    }
    
    # flatten theta 
    Theta_4d <- matrix(Theta_all, nrow = m * m, ncol = p * Tt)
    
    # storage setup
    n_iter <- burn + nscan
    keep <- seq(burn + 1, n_iter, by = odens)
    n_keep <- length(keep)
    
    # preallocate storage arrays based on time_thin
    time_keep <- seq(1, Tt, by = time_thin)
    n_time_keep <- length(time_keep)
    
    # 
    Theta_store <- array(NA, dim = c(m, m, p, n_time_keep, n_keep))
    if (FAM$name %in% c("ordinal", "binary")) {
        Z_store <- array(NA, dim = c(m, m, p, n_time_keep, n_keep))
    }
    
    A_store <- B_store <- vector("list", n_keep)
    M_store <- array(NA, dim = c(m, m, p, n_keep))
    
    sigma2_store <- numeric(n_keep)
    sigma2_obs_store <- if (FAM$name == "gaussian") numeric(n_keep) else NULL
    tau_A2_store <- tau_B2_store <- g2_store <- numeric(n_keep)
    rhoA_store <- rhoB_store <- if (ar1) numeric(n_keep) else NULL
    
    keep_id <- 0
    a_sig <- b_sig <- 2
    a_g <- b_g <- 2
    
    # preallocate reusable workspace arrays for large networks
    if (is_large_network) {
        # precompute constants used in every iteration
        eye_m <- diag(m)
        shape_tauA <- (1 + m*m*(Tt-1))/2
        shape_tauB <- (1 + m*m*(Tt-1))/2
        shape_sigma_proc <- (a_sig + m * m * (Tt - 1) * p) / 2.0
        shape_sigma_obs <- (1.0 + m * m * Tt * p) / 2.0
        
        # preallocate workspace matrices to avoid repeated allocation
        workspace_mat <- matrix(0, m, m)
        workspace_vec <- numeric(m)
    }
    
    # progress tracking
    if (verbose) {
        cli::cli_progress_step("Running dynamic DBN MCMC")
        cli::cli_progress_bar("MCMC iterations", total = n_iter)
    }
    
    # main mcmc loop
    for (g in 1:n_iter) {
        # 1. batch update z for all relations (ordinal only)
        if (FAM$name == "ordinal") {
            # check if we should use gaussian approximation for speed
            use_approx <- should_use_gaussian_approximation(R) || 
                         (m * m * p * Tt > 5000)  # also use for large problems
            
            if (use_approx) {
                # use fast gaussian approximation
                EZ_cube <- array(Theta_4d + rep(as.vector(M), each = Tt), c(m, m, p * Tt))
                Z_cube <- array(Z_4d, c(m, m, p * Tt))
                R_cube <- array(R_4d, c(m, m, p * Tt))
                
                if (exists("rz_gaussian_approx_cpp", mode = "function")) {
                    Z_cube <- rz_gaussian_approx_cpp(R_cube, Z_cube, EZ_cube)
                } else {
                    Z_cube <- rz_gaussian_approx(R_cube, Z_cube, EZ_cube)
                }
                Z_4d <- matrix(Z_cube, nrow = m * m, ncol = p * Tt)
            } else {
                # use exact truncated normal sampling
                if (exists("IR_time_indices")) {
                    Z_4d <- batch_update_Z_ordinal_fast(R_4d, Z_4d, Theta_4d, M, IR, IR_time_indices, m, p, Tt)
                } else {
                    Z_4d <- batch_update_Z_ordinal(R_4d, Z_4d, Theta_4d, M, IR, m, p, Tt)
                }
            }
            # only reshape if needed for saving later
            if (g %in% keep) {
                Z <- array(Z_4d, dim = c(m, m, p, Tt))
            }
        } else if (FAM$name == "binary") {
            # use binary update
            Z <- update_Z_optimized(R, Z, Theta_all, M, IR = NULL, family = "binary")
            Z_4d <- matrix(Z, nrow = m * m, ncol = p * Tt)
        }
        
        # 2. mu update 
        mu_result <- update_mu_dynamic(Z_4d, Theta_4d, g2, a_g, b_g, m, p, Tt)
        M <- mu_result$M
        g2 <- mu_result$g2
        
        # 3. batch ffbs for all relations - use blocked version for large networks
        if (is_large_network && m > 100) {
            Theta_cube <- batch_ffbs_all_relations_blocked(Z_4d, M, Aarray, Barray, sigma2, m, p, Tt)
        } else {
            Theta_cube <- batch_ffbs_all_relations(Z_4d, M, Aarray, Barray, sigma2, m, p, Tt)
        }
        # only reshape to 4D array if needed for saving
        if (g %in% keep) {
            Theta_all <- array(Theta_cube, dim = c(m, m, p, Tt))
        }
        # convert cube back to matrix for update_AB_batch_extended
        Theta_4d <- matrix(Theta_cube, nrow = m * m, ncol = p * Tt)
        
        # 4. extended batch ab update
        if (is_large_network && m > 100) {
            AB_result <- update_AB_batch_large(
                Theta_4d, Aarray, Barray,
                sigma2, tauA2, tauB2,
                ar1, rhoA, rhoB,
                m, p, Tt
            )
        } else {
            AB_result <- update_AB_batch_extended(
                Theta_4d, Aarray, Barray,
                sigma2, tauA2, tauB2,
                ar1, rhoA, rhoB,
                m, p, Tt
            )
        }
        Aarray <- AB_result$Aarray
        Barray <- AB_result$Barray
        
        # 5. update hyperparameters
        # tau_a2 and tau_b2
        if (ar1) {
            # ar(1) case with vectorized operations
            if (is_large_network) {
                # vectorized computation for large networks
                innovA_ss <- 0
                innovB_ss <- 0
                for (t in 2:Tt) {
                    A_innov <- Aarray[,,t] - rhoA * Aarray[,,t-1] - (1-rhoA) * eye_m
                    B_innov <- Barray[,,t] - rhoB * Barray[,,t-1] - (1-rhoB) * eye_m
                    innovA_ss <- innovA_ss + sum(A_innov^2)
                    innovB_ss <- innovB_ss + sum(B_innov^2)
                }
                tauA2 <- safe_rinv_gamma(shape_tauA, (1 + innovA_ss)/2)
                tauB2 <- safe_rinv_gamma(shape_tauB, (1 + innovB_ss)/2)
            } else {
                # standard computation for smaller networks
                innovA_ss <- innovB_ss <- 0
                for (t in 2:Tt) {
                    innovA_ss <- innovA_ss + sum((Aarray[,,t] - rhoA * Aarray[,,t-1] - (1-rhoA) * diag(m))^2)
                    innovB_ss <- innovB_ss + sum((Barray[,,t] - rhoB * Barray[,,t-1] - (1-rhoB) * diag(m))^2)
                }
                tauA2 <- safe_rinv_gamma((1 + m*m*(Tt-1))/2, (1 + innovA_ss)/2)
                tauB2 <- safe_rinv_gamma((1 + m*m*(Tt-1))/2, (1 + innovB_ss)/2)
            }
        } else {
            # random walk case
            A_sum <- compute_deviation_sum(Aarray, m, Tt)
            B_sum <- compute_deviation_sum(Barray, m, Tt)
            if (is_large_network) {
                tauA2 <- safe_rinv_gamma(shape_tauA, (1 + A_sum)/2)
                tauB2 <- safe_rinv_gamma(shape_tauB, (1 + B_sum)/2)
            } else {
                tauA2 <- safe_rinv_gamma((1 + m*m*(Tt-1))/2, (1 + A_sum)/2)
                tauB2 <- safe_rinv_gamma((1 + m*m*(Tt-1))/2, (1 + B_sum)/2)
            }
        }
        
        # 6. update variances
        if (is_large_network && m > 100) {
            # use blocked variance computation for very large networks
            proc_rss <- compute_process_variance_blocked(Theta_4d, Aarray, Barray, m, p, Tt)
            if (exists("shape_sigma_proc")) {
                sigma2 <- (b_sig + proc_rss / 2.0) / rgamma(1, shape = shape_sigma_proc, rate = 1)
            } else {
                sigma2 <- (b_sig + proc_rss / 2.0) / rgamma(1, shape = (a_sig + m * m * (Tt - 1) * p) / 2.0, rate = 1)
            }
            
            if (FAM$name == "gaussian") {
                # compute observation variance for Gaussian
                obs_rss <- 0
                for (j in 1:p) {
                    M_j <- M[,,j]
                    for (t in 1:Tt) {
                        idx <- (j-1) * Tt + t
                        Z_jt <- matrix(Z_4d[,idx], m, m)
                        Theta_jt <- matrix(Theta_4d[,idx], m, m)
                        obs_rss <- obs_rss + sum((Z_jt - (Theta_jt + M_j))^2)
                    }
                }
                sigma2_obs <- (1 + obs_rss / 2) / rgamma(1, shape = (1 + m * m * Tt * p) / 2, rate = 1)
            }
        } else {
            # standard variance update for smaller networks
            var_result <- update_variances_dynamic(
                Theta_4d, Z_4d, M, Aarray, Barray,
                a_sig, b_sig, m, p, Tt,
                is_gaussian = (FAM$name == "gaussian")
            )
            sigma2 <- var_result$sigma2
            if (FAM$name == "gaussian") {
                sigma2_obs <- var_result$sigma2_obs
            }
        }
        
        # 7. update ar coefficients if requested
        if (ar1 && update_rho) {
            # ar coefficient updates with precomputed identity
            identity_mat <- if (is_large_network) eye_m else diag(m)
            
            # update rhoa - vectorized operations
            num <- denom <- 0
            for (t in 2:Tt) {
                diff_t <- Aarray[,,t] - identity_mat
                diff_tm1 <- Aarray[,,t-1] - identity_mat
                num <- num + sum(diff_t * diff_tm1)
                denom <- denom + sum(diff_tm1^2)
            }
            rho_mean <- num / (denom + 1e-10)
            rho_var <- tauA2 / (denom + 1e-10)
            rhoA <- truncnorm::rtruncnorm(1, a = -0.99, b = 0.99, mean = rho_mean, sd = sqrt(rho_var))
            
            # update rhob - vectorized operations
            num <- denom <- 0
            for (t in 2:Tt) {
                diff_t <- Barray[,,t] - identity_mat
                diff_tm1 <- Barray[,,t-1] - identity_mat
                num <- num + sum(diff_t * diff_tm1)
                denom <- denom + sum(diff_tm1^2)
            }
            rho_mean <- num / (denom + 1e-10)
            rho_var <- tauB2 / (denom + 1e-10)
            rhoB <- truncnorm::rtruncnorm(1, a = -0.99, b = 0.99, mean = rho_mean, sd = sqrt(rho_var))
        }
        
        # 8. save samples
        if (g %in% keep) {
            keep_id <- keep_id + 1
            
            # storage with time thinning
            # only reshape arrays if not already done above
            if (!exists("Theta_all") || !is.array(Theta_all) || length(dim(Theta_all)) != 4) {
                Theta_all <- array(Theta_cube, dim = c(m, m, p, Tt))
            }
            Theta_store[,,,,keep_id] <- Theta_all[,,,time_keep, drop = FALSE]
            
            if (FAM$name %in% c("ordinal", "binary")) {
                if (!is.array(Z) || length(dim(Z)) != 4) {
                    Z <- array(Z_4d, dim = c(m, m, p, Tt))
                }
                Z_store[,,,,keep_id] <- Z[,,,time_keep, drop = FALSE]
            }
            
            # store time-thinned a and b
            A_store[[keep_id]] <- Aarray[,,time_keep, drop = FALSE]
            B_store[[keep_id]] <- Barray[,,time_keep, drop = FALSE]
            
            M_store[,,,keep_id] <- M
            sigma2_store[keep_id] <- sigma2
            if (FAM$name == "gaussian") {
                sigma2_obs_store[keep_id] <- sigma2_obs
            }
            
            tau_A2_store[keep_id] <- tauA2
            tau_B2_store[keep_id] <- tauB2
            g2_store[keep_id] <- g2
            
            if (ar1) {
                rhoA_store[keep_id] <- rhoA
                rhoB_store[keep_id] <- rhoB
            }
        }
        
        if (verbose && (g %% verbose == 0)) cli::cli_progress_update()
    }
    
    if (verbose) cli::cli_progress_done()
    
    # prepare output
    out <- list(
        model = "dynamic",
        Y = Y,
        # Theta = Theta_store,  # Removed to save memory - use posterior_extract functions if needed
        A = A_store,
        B = B_store,
        M = M_store,  # keep 4D array structure
        sigma2 = sigma2_store,
        tau_A2 = tau_A2_store,
        tau_B2 = tau_B2_store,
        g2 = g2_store,
        n_iter = n_iter,
        burn = burn,
        thin = odens,
        time_thin = time_thin,
        dims = dims,
        family = family
    )
    
    if (FAM$name == "gaussian") {
        out$sigma2_obs <- sigma2_obs_store
    }
    
    if (FAM$name %in% c("ordinal", "binary")) {
        out$Z <- Z_store  # Keep 5D array structure
    }
    
    if (ar1) {
        out$rhoA <- rhoA_store
        out$rhoB <- rhoB_store
        out$ar1 <- TRUE
    }
    
    # add warm start info if continuing from previous
    if (!is.null(previous)) {
        prev_total <- previous$n_iter %||% (previous$burn + length(previous$sigma2))
        out$total_iter <- prev_total + nscan
        out$continued_from <- prev_total
    }
    
    class(out) <- "dbn"
    out
}