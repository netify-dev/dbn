#' Posterior Predictive Functions
#'
#' @description Functions for generating posterior predictive samples and checks
#' @name posterior_predict
#' @keywords internal
NULL

#' Generate posterior predictive samples
#'
#' @description Generate new observations from the posterior predictive distribution
#' @param fit A dbn model fit object
#' @param ndraws Number of posterior draws to use (default: 100)
#' @param seed Random seed for reproducibility
#' @param draws Specific draw indices to use (overrides ndraws)
#' @return List of predicted observations with class "dbn_ppd"
#' @export
posterior_predict_dbn <- function(fit, ndraws = 100, seed = NULL, draws = NULL) {
    # set seed if provided
    if (!is.null(seed)) set.seed(seed)

    # get family object
    fam <- get_family(fit)
    if (is.null(fam)) {
        stop("Model fit does not contain family information")
    }

    # determine which draws to use
    if (is.null(draws)) {
        # Determine total draws based on model type
        if (fit$model == "dynamic") {
            n_total_draws <- length(fit$sigma2)
        } else if (fit$model == "static") {
            n_total_draws <- nrow(fit$params)
        } else if (fit$model == "hmm") {
            n_total_draws <- length(fit$draws$pars$sigma2_proc)
        } else if (fit$model == "lowrank" || fit$model == "lowrank_accurate") {
            n_total_draws <- length(fit$sigma2)
        } else {
            n_total_draws <- fit$meta$draws %||% length(fit$draws$theta)
        }
        
        if (ndraws >= n_total_draws) {
            draws <- seq_len(n_total_draws)
        } else {
            draws <- sample(seq_len(n_total_draws), ndraws, replace = TRUE)
        }
    } else {
        ndraws <- length(draws)
    }

    # generate predictions for each draw
    Yrep <- vector("list", ndraws)

    for (d in seq_along(draws)) {
        draw_idx <- draws[d]

        # extract theta for this draw based on model type
        if (fit$model == "dynamic") {

            if (is.null(fit$A) || is.null(fit$B) || is.null(fit$M)) {
                cli::cli_abort(c(
                    "x" = "Cannot compute theta: A, B, or M matrices not found in fit object.",
                    "i" = "Ensure the model was fitted properly."
                ))
            }
            
            # get A and B arrays for this draw
            A_draw <- fit$A[[draw_idx]]
            B_draw <- fit$B[[draw_idx]]
            
            # handle M structure - could be 4D array or list
            if (is.array(fit$M) && length(dim(fit$M)) == 4) {
                # m stored as 4D array (m × m × p × draws)
                M_draw <- fit$M[, , , draw_idx, drop = FALSE]
                # ensure M_draw is 3D array
                dim(M_draw) <- dim(fit$M)[1:3]
            } else if (is.list(fit$M)) {
                # m stored as list
                M_draw <- fit$M[[draw_idx]]
            } else {
                cli::cli_abort("Unexpected M structure in dynamic model fit")
            }
            
            # get dims
            dims <- fit$dims
            m <- dims$m
            p <- dims$p
            Tt <- dims$Tt
            
            # init theta array
            th <- array(0, dim = c(m, m, p, Tt))
            
            # compute theta for each relation and time point
            # for dynamic model: theta_{j,t} = A_t * theta_{j,t-1} * B_t' + mu_j
            for (j in 1:p) {
                # first time point: theta_{j,1} = mu_j
                th[, , j, 1] <- M_draw[, , j]
                
                # subsequent time points: use bilinear dynamics
                for (t in 2:Tt) {
                    # Check if A and B have full time dimension or are thinned
                    t_idx <- if (dim(A_draw)[3] < Tt) {
                        # Time thinned - find appropriate index
                        time_thin <- fit$time_thin %||% 1
                        ceiling(t / time_thin)
                    } else {
                        t
                    }
                    
                    # apply bilinear transformation
                    th[, , j, t] <- A_draw[, , t_idx] %*% th[, , j, t-1] %*% t(B_draw[, , t_idx]) + M_draw[, , j]
                }
            }
        } else if (fit$model == "static") {
            # compute theta on-demand for static model
            if (is.null(fit$draws$misc$B) || is.null(fit$draws$misc$M)) {
                cli::cli_abort(c(
                    "x" = "Cannot compute theta: B matrices or M not found in fit object.",
                    "i" = "Ensure the model was fitted properly."
                ))
            }
            
            # get B matrix for this draw (only the first one is used in static model)
            B <- fit$draws$misc$B[[1]][, , draw_idx]
            
            # get M for this draw
            M <- fit$draws$misc$M[[draw_idx]]
            
            # get Z for this draw
            if (!is.null(fit$draws$z) && length(fit$draws$z) >= draw_idx) {
                # for ordinal/binary models, Z is stored
                Z <- fit$draws$z[[draw_idx]]
            } else if (fam$name == "gaussian") {
                # for Gaussian model, Z = Y (observed data)
                Z <- fit$Y
            } else {
                cli::cli_abort(c(
                    "x" = "Cannot determine Z values for theta computation.",
                    "i" = "Model family not supported for on-demand theta computation."
                ))
            }
            
            # compute theta
            th <- compute_theta_static(B, Z, M)
        } else if (fit$model == "hmm" && !is.null(fit$draws$theta)) {
            # HMM model may store theta (unless large network or time thinned)
            th <- fit$draws$theta[[draw_idx]]
        } else if (fit$model == "lowrank" || fit$model == "lowrank_accurate") {
            # lowrank models need to reconstruct theta from U, alpha, B
            if (is.null(fit$U) || is.null(fit$alpha) || is.null(fit$B)) {
                cli::cli_abort(c(
                    "x" = "Cannot reconstruct theta: U, alpha, or B not found in fit object.",
                    "i" = "Ensure the lowrank model was fitted properly."
                ))
            }
            
            # 
            dims <- fit$meta$dims %||% fit$dims
            m <- dims$m
            p <- dims$p
            T_keep <- dim(fit$alpha[[1]])[2]  # number of saved time points
            
            # 
            th <- array(0, dim = c(m, m, p, T_keep))
            
            # get components for this draw
            U_s <- fit$U[[draw_idx]]
            alpha_s <- fit$alpha[[draw_idx]]
            B_s <- fit$B[[draw_idx]]
            
            # reconstruct theta for each time point
            # Note: we need the initial theta or assume stationary start
            # straighfroward approach is to say that we
            # will compute forward from a zero initial state
            # this matches what's done in dyad_path for lowrank models
            
            for (t in 1:T_keep) {
                # reconstruct A_t from lowrank decomposition
                if (length(alpha_s[, t]) == 1) {
                    # rank-1 case
                    A_t <- alpha_s[, t] * U_s %*% t(U_s)
                } else {
                    # general rank-r case
                    A_t <- U_s %*% diag(alpha_s[, t]) %*% t(U_s)
                }
                
                # 
                B_t <- B_s[, , t]
                
                if (t == 1) {
                    # for first time point, we don't have previous theta
                    # use stationary approximation or zero
                    for (rel in 1:p) {
                        th[, , rel, t] <- matrix(0, m, m)  # or could use stationary solution
                    }
                } else {
                    # update theta using dynamic equation
                    for (rel in 1:p) {
                        th[, , rel, t] <- A_t %*% th[, , rel, t-1] %*% t(B_t)
                    }
                }
            }
        } else {
            # 
            if (!is.null(fit$draws$theta)) {
                th <- fit$draws$theta[[draw_idx]]
            } else {
                cli::cli_abort(c(
                    "x" = "Cannot extract theta for model type '{fit$model}'.",
                    "i" = "Model may not support posterior predictive checks yet."
                ))
            }
        }

        # extract miscellaneous parameters (eg m arrays, variance)
        misc <- list()

        # extract m (baseline mean)
        if (fit$model == "dynamic" && !is.null(fit$M)) {
            # dynamic model stores M as 4D array: m x m x p x iter
            misc$M <- fit$M[, , , draw_idx, drop = FALSE]
            # remove the iteration dimension to get 3D array
            dim(misc$M) <- dim(fit$M)[1:3]
        } else if ((fit$model == "lowrank" || fit$model == "lowrank_accurate") && !is.null(fit$M)) {
            # lowrank model stores M as a list
            if (is.list(fit$M)) {
                misc$M <- fit$M[[draw_idx]]
            } else {
                # 
                misc$M <- fit$M[, , , draw_idx, drop = FALSE]
                dim(misc$M) <- dim(fit$M)[1:3]
            }
        } else if (!is.null(fit$draws$misc$M) && length(fit$draws$misc$M) >= draw_idx) {
            misc$M <- fit$draws$misc$M[[draw_idx]]
        } else if (!is.null(fit$M)) {
            misc$M <- fit$M
        }

        if (fam$name == "gaussian") {
            if ((fit$model == "dynamic" || fit$model == "lowrank" || fit$model == "lowrank_accurate") 
                && !is.null(fit$sigma2_obs)) {
                # dynamic and lowrank models store sigma2_obs as a vector
                misc$sigma2_obs <- fit$sigma2_obs[draw_idx]
            } else if (!is.null(fit$draws$pars) && "sigma2_obs" %in% names(fit$draws$pars)) {
                # HMM and other models store in draws$pars
                misc$sigma2_obs <- fit$draws$pars$sigma2_obs[draw_idx]
            } else if (!is.null(fit$sigma2_obs) && length(fit$sigma2_obs) >= draw_idx) {
                misc$sigma2_obs <- fit$sigma2_obs[draw_idx]
            }
        }

        # generate observations
        Yrep[[d]] <- fam$rgen_obs(th, misc)
    }

    structure(Yrep,
        class = "dbn_ppd",
        family = fam$name,
        dims = fit$meta$dims %||% fit$dims
    )
}

#' Print method for posterior predictive distribution
#' @param x An object of class "dbn_ppd"
#' @param ... Additional arguments (currently unused)
#' @export
print.dbn_ppd <- function(x, ...) {
    cat("Posterior predictive distribution\n")
    cat("Family:", attr(x, "family"), "\n")
    cat("Number of draws:", length(x), "\n")

    dims <- attr(x, "dims")
    if (!is.null(dims)) {
        cat(
            "Data dimensions:",
            dims$m, "nodes x", dims$p, "relations x",
            dims$Tt %||% dims$T, "time points\n"
        )
    }

    invisible(x)
}

#' Get family object from model fit
#'
#' @description Extract or reconstruct family object from model fit
#' @param fit A dbn model fit object
#' @return A dbn_family object or NULL
#' @keywords internal
get_family <- function(fit) {
    # check if family is stored directly
    if (!is.null(fit$family) && inherits(fit$family, "dbn_family")) {
        return(fit$family)
    }

    # try to determine family from settings
    family_name <- fit$settings$family %||% fit$family

    if (is.null(family_name)) {
        # try to infer from model structure
        if ("IR" %in% names(fit)) {
            family_name <- "ordinal"
        } else if (!is.null(fit$sigma2_obs)) {
            family_name <- "gaussian"
        } else {
            warning("Cannot determine family type from model fit")
            return(NULL)
        }
    }

    # construct family object
    switch(family_name,
        ordinal = family_ordinal(),
        gaussian = family_gaussian(),
        binary = family_binary(),
        NULL
    )
}

#' Extract theta from legacy format
#'
#' @description Attempt to reconstruct theta array from legacy storage format
#' @param fit Model fit object
#' @param draw_idx Draw index
#' @return Theta array or NULL
#' @keywords internal
extract_theta_legacy <- function(fit, draw_idx) {
    # for lowrank models
    if (!is.null(fit$U) && !is.null(fit$alpha)) {
        U <- fit$U[[draw_idx]]
        alpha <- fit$alpha[[draw_idx]]

        # reconstruct a matrices and combine with b
        # this is model-specific and may need adjustment
        warning("Legacy theta extraction for lowrank models is approximate")
        return(NULL)
    }

    # for standard models with a and b
    if (!is.null(fit$A) && !is.null(fit$B)) {
        # this would need the actual theta samples, not just a/b
        warning("Cannot extract Theta from A/B matrices alone")
        return(NULL)
    }

    NULL
}

#' Posterior predictive ECDF plot
#'
#' @description Plot empirical CDFs of observed vs replicated data
#' @param fit A dbn model fit object
#' @param ppd Posterior predictive samples (from posterior_predict_dbn)
#' @param ndraws_plot Number of replicates to show (default: 20)
#' @param alpha Transparency for replicate lines (default: 0.3)
#' @param rel Relation index to plot (default: 1)
#' @param time Time index to plot (default: all)
#' @param Y_obs Observed data array (required)
#' @return A ggplot object
#' @export
plot_ppc_ecdf <- function(fit, ppd = NULL, ndraws_plot = 20,
                          alpha = 0.3, rel = 1, time = NULL, Y_obs = NULL) {
    # generate ppd if not provided
    if (is.null(ppd)) {
        ppd <- posterior_predict_dbn(fit, ndraws = ndraws_plot)
    }

    # extract observed data
    if (is.null(Y_obs)) {
        stop("Observed data (Y_obs) must be provided")
    }

    # subset to specific relation and time
    dims <- fit$meta$dims %||% fit$dims
    if (is.null(time)) {
        time <- seq_len(dims$Tt %||% dims$T)
    }

    # extract values for plotting
    obs_vals <- c(Y_obs[, , rel, time])
    obs_vals <- obs_vals[!is.na(obs_vals)]

    # subsample if too large for efficiency
    if (length(obs_vals) > 5e5) {
        obs_vals <- sample(obs_vals, 5e5)
    }

    # get unique values for ecdf
    vals <- sort(unique(obs_vals))

    # build data frame for plotting
    df <- data.frame(
        value = vals,
        ecdf = ecdf(obs_vals)(vals),
        type = "Observed",
        set = "Observed"
    )

    # add replicates
    draws_show <- sample(seq_along(ppd), min(ndraws_plot, length(ppd)))

    for (d in draws_show) {
        rep_vals <- c(ppd[[d]][, , rel, time])
        rep_vals <- rep_vals[!is.na(rep_vals)]

        # subsample if too large for efficiency
        if (length(rep_vals) > 5e5) {
            rep_vals <- sample(rep_vals, 5e5)
        }

        df_rep <- data.frame(
            value = vals,
            ecdf = ecdf(rep_vals)(vals),
            type = "Replicated",
            set = paste0("Rep", d)
        )

        df <- rbind(df, df_rep)
    }

    # create plot
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        p <- ggplot2::ggplot(df, ggplot2::aes(
            x = value, y = ecdf,
            group = set,
            color = type,
            size = type,
            alpha = type
        )) +
            ggplot2::geom_line() +
            ggplot2::scale_color_manual(values = c(
                "Observed" = "black",
                "Replicated" = "grey60"
            )) +
            ggplot2::scale_size_manual(values = c(
                "Observed" = 1.2,
                "Replicated" = 0.8
            )) +
            ggplot2::scale_alpha_manual(values = c(
                "Observed" = 1,
                "Replicated" = alpha
            )) +
            ggplot2::labs(
                title = paste0("Posterior Predictive Check: Relation ", rel),
                x = "Y",
                y = "ECDF",
                color = NULL, size = NULL, alpha = NULL
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "bottom")

        return(p)
    } else {
        # base r fallback
        plot(ecdf(obs_vals),
            main = paste0("PPC: Relation ", rel),
            xlab = "Y", ylab = "ECDF", lwd = 2
        )

        for (d in draws_show) {
            rep_vals <- c(ppd[[d]][, , rel, time])
            rep_vals <- rep_vals[!is.na(rep_vals)]
            lines(ecdf(rep_vals), col = adjustcolor("grey60", alpha.f = alpha))
        }

        legend("bottomright", c("Observed", "Replicated"),
            lwd = c(2, 1), col = c("black", "grey60")
        )

        invisible(NULL)
    }
}

#' Posterior predictive density plot
#'
#' @description Plot density of observed vs replicated data
#' @param fit A dbn model fit object
#' @param ppd Posterior predictive samples
#' @param rel Relation index
#' @param time Time indices
#' @param Y_obs Observed data array (required)
#' @return A ggplot object or NULL
#' @export
plot_ppc_density <- function(fit, ppd = NULL, rel = 1, time = NULL, Y_obs = NULL) {
    # generate ppd if not provided
    if (is.null(ppd)) {
        ppd <- posterior_predict_dbn(fit, ndraws = 50)
    }

    # extract observed data
    if (is.null(Y_obs)) {
        stop("Observed data (Y_obs) must be provided")
    }
    dims <- fit$meta$dims %||% fit$dims

    if (is.null(time)) {
        time <- seq_len(dims$Tt %||% dims$T)
    }

    # get observed values
    obs_vals <- c(Y_obs[, , rel, time])
    obs_vals <- obs_vals[!is.na(obs_vals)]

    # for discrete data, use bar plot instead
    if (length(unique(obs_vals)) <= 20) {
        return(plot_ppc_bars(fit, ppd, rel, time, Y_obs))
    }

    # continuous data - density plot
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        # build data frame
        df <- data.frame(value = obs_vals, type = "Observed")

        # add replicated data
        for (d in seq_along(ppd)) {
            rep_vals <- c(ppd[[d]][, , rel, time])
            rep_vals <- rep_vals[!is.na(rep_vals)]
            df <- rbind(df, data.frame(
                value = rep_vals,
                type = paste0("Rep", d)
            ))
        }

        # aggregate replicated densities
        p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = type == "Observed")) +
            ggplot2::geom_density(alpha = 0.5) +
            ggplot2::scale_fill_manual(
                values = c("TRUE" = "black", "FALSE" = "grey60"),
                labels = c("TRUE" = "Observed", "FALSE" = "Replicated")
            ) +
            ggplot2::labs(
                title = paste0("Posterior Predictive Density: Relation ", rel),
                x = "Y", y = "Density", fill = NULL
            ) +
            ggplot2::theme_minimal()

        return(p)
    }

    # base r fallback
    plot(density(obs_vals),
        main = paste0("PPC Density: Relation ", rel),
        xlab = "Y", ylab = "Density", lwd = 2
    )

    for (d in seq_along(ppd)) {
        rep_vals <- c(ppd[[d]][, , rel, time])
        rep_vals <- rep_vals[!is.na(rep_vals)]
        lines(density(rep_vals), col = adjustcolor("grey60", alpha.f = 0.3))
    }

    invisible(NULL)
}

#' Posterior predictive bar plot for discrete data
#'
#' @description Bar plot comparison for discrete outcomes
#' @param fit A dbn model fit object
#' @param ppd Posterior predictive samples
#' @param rel Relation index
#' @param time Time indices
#' @param Y_obs Observed data array (required)
#' @return A ggplot object or NULL
#' @keywords internal
plot_ppc_bars <- function(fit, ppd, rel = 1, time = NULL, Y_obs = NULL) {
    # extract observed data
    if (is.null(Y_obs)) {
        stop("Observed data (Y_obs) must be provided")
    }
    dims <- fit$meta$dims %||% fit$dims

    if (is.null(time)) {
        time <- seq_len(dims$Tt %||% dims$T)
    }

    # get observed frequencies
    obs_vals <- c(Y_obs[, , rel, time])
    obs_vals <- obs_vals[!is.na(obs_vals)]
    obs_freq <- table(obs_vals) / length(obs_vals)

    # get replicated frequencies
    rep_freq_list <- lapply(ppd, function(yrep) {
        rep_vals <- c(yrep[, , rel, time])
        rep_vals <- rep_vals[!is.na(rep_vals)]
        table(rep_vals) / length(rep_vals)
    })

    # aggregate across replicates
    all_vals <- sort(unique(c(
        names(obs_freq),
        unlist(lapply(rep_freq_list, names))
    )))

    rep_freq_mean <- numeric(length(all_vals))
    rep_freq_lower <- numeric(length(all_vals))
    rep_freq_upper <- numeric(length(all_vals))

    for (i in seq_along(all_vals)) {
        val <- all_vals[i]
        freqs <- sapply(rep_freq_list, function(f) f[val] %||% 0)
        rep_freq_mean[i] <- mean(freqs, na.rm = TRUE)
        rep_freq_lower[i] <- quantile(freqs, 0.05, na.rm = TRUE)
        rep_freq_upper[i] <- quantile(freqs, 0.95, na.rm = TRUE)
    }

    # create plot
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        df <- data.frame(
            value = rep(all_vals, 2),
            freq = c(as.numeric(obs_freq[all_vals]), rep_freq_mean),
            type = rep(c("Observed", "Replicated"), each = length(all_vals)),
            lower = c(rep(NA, length(all_vals)), rep_freq_lower),
            upper = c(rep(NA, length(all_vals)), rep_freq_upper)
        )
        df$freq[is.na(df$freq)] <- 0

        p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = freq, fill = type)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                position = ggplot2::position_dodge(0.9),
                width = 0.25, na.rm = TRUE
            ) +
            ggplot2::scale_fill_manual(values = c(
                "Observed" = "black",
                "Replicated" = "grey60"
            )) +
            ggplot2::labs(
                title = paste0("Posterior Predictive Check: Relation ", rel),
                x = "Y", y = "Frequency", fill = NULL
            ) +
            ggplot2::theme_minimal()

        return(p)
    }

    # base r fallback
    barplot(rbind(as.numeric(obs_freq[all_vals]), rep_freq_mean),
        beside = TRUE, names.arg = all_vals,
        col = c("black", "grey60"),
        main = paste0("PPC: Relation ", rel),
        xlab = "Y", ylab = "Frequency"
    )
    legend("topright", c("Observed", "Replicated"),
        fill = c("black", "grey60")
    )

    invisible(NULL)
}
