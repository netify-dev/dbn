#' Build shock matrix for IRF analysis
#'
#' @description Creates shock matrices for different types of network interventions
#' @param m Number of nodes in the network
#' @param type Type of shock: "unit_edge", "node_out", "node_in", or "density"
#' @param i Source node index (for unit_edge and node shocks)
#' @param j Target node index (for unit_edge shock)
#' @param magnitude Shock magnitude (default: 1)
#' @return m x m shock matrix
#' @export
build_shock <- function(m, type = c("unit_edge", "node_out", "node_in", "density"), 
                       i = 1, j = 2, magnitude = 1) {
    type <- match.arg(type)
    S <- matrix(0, m, m)
    
    switch(type,
        unit_edge = {
            if (i < 1 || i > m || j < 1 || j > m) {
                stop("Node indices must be between 1 and m")
            }
            S[i, j] <- magnitude
        },
        node_out = {
            if (i < 1 || i > m) {
                stop("Node index must be between 1 and m")
            }
            S[i, ] <- magnitude
        },
        node_in = {
            if (i < 1 || i > m) {
                stop("Node index must be between 1 and m")
            }
            S[, i] <- magnitude
        },
        density = {
            S[,] <- magnitude / (m * (m - 1))
            diag(S) <- 0
        }
    )
    
    S
}

#' Network statistic: density
#'
#' @description Compute network density (mean of off-diagonal elements)
#' @param X Network matrix
#' @return Scalar density value
#' @export
stat_density <- function(X) {
    mean(X[row(X) != col(X)])
}

#' Network statistic: in-degree
#'
#' @description Compute in-degree for all nodes (column sums)
#' @param X Network matrix
#' @return Vector of in-degrees
#' @export
stat_in_degree <- function(X) {
    colSums(X)
}

#' Network statistic: out-degree
#'
#' @description Compute out-degree for all nodes (row sums)
#' @param X Network matrix
#' @return Vector of out-degrees
#' @export
stat_out_degree <- function(X) {
    rowSums(X)
}

#' Network statistic: reciprocity
#'
#' @description Compute network reciprocity (correlation between X[i,j] and X[j,i])
#' @param X Network matrix
#' @return Scalar reciprocity value
#' @export
stat_reciprocity <- function(X) {
    upper_idx <- upper.tri(X)
    upper_vals <- X[upper_idx]
    lower_vals <- t(X)[upper_idx]
    
    if (length(unique(c(upper_vals, lower_vals))) == 1) {
        # All values are the same, correlation is undefined
        return(0)
    }
    
    cor(upper_vals, lower_vals)
}

#' Network statistic: transitivity
#'
#' @description Compute transitivity (clustering coefficient)
#' @param X Network matrix
#' @return Scalar transitivity value
#' @export
stat_transitivity <- function(X) {
    X_binary <- (X != 0) * 1
    diag(X_binary) <- 0
    
    # Count triangles
    triangles <- sum(diag(X_binary %*% X_binary %*% X_binary)) / 6
    
    # Count connected triples
    degrees <- rowSums(X_binary)
    triples <- sum(degrees * (degrees - 1)) / 2
    
    if (triples == 0) return(0)
    triangles / triples
}

#' Compute IRF for a single posterior draw
#'
#' @param fit dbn model fit object
#' @param draw_idx Index of posterior draw to use
#' @param shock Shock matrix
#' @param H Number of horizons
#' @param t0 Shock time (for dynamic models, 1-based)
#' @param stat_fun Network statistic function
#' @return Vector of IRF values at each horizon
#' @keywords internal
compute_irf_single <- function(fit, draw_idx, shock, H, t0 = 1, stat_fun = stat_density) {
    # Get dimensions safely
    dims <- fit$dims
    m <- as.integer(dims$m)
    p <- as.integer(dims$p)
    
    # Get A and B matrices for this draw
    if (fit$model == "dynamic") {
        # For dynamic model, we need time-varying A and B
        # Dynamic models store A and B as lists in fit$A and fit$B
        if (is.null(fit$A) || is.null(fit$B)) {
            stop("Dynamic model requires A and B matrices")
        }
        
        # Extract A and B arrays for this draw
        # Check time thinning
        time_thin <- fit$time_thin %||% 1
        T_stored <- dim(fit$A[[draw_idx]])[3]
        
        # Map t0 to stored time index if time thinning is used
        if (time_thin > 1) {
            t0_stored <- ceiling(t0 / time_thin)
        } else {
            t0_stored <- t0
        }
        
        # Check bounds
        if (t0_stored + H > T_stored) {
            stop("t0 + H exceeds stored time points. With time_thin = ", time_thin, 
                 ", maximum t0 + H is ", T_stored * time_thin)
        }
        
        # Extract the relevant time slices
        A_array <- fit$A[[draw_idx]]
        B_array <- fit$B[[draw_idx]]
        
        # Compute impulse response using stored indices
        Delta <- impulse_response_dynamic(A_array, B_array, shock, t0_stored - 1, H)
        
    } else if (fit$model == "static") {
        # For static model, A and B are constant
        if (is.null(fit$B)) {
            stop("Static model requires B matrix")
        }
        
        # Extract B for this draw
        B <- matrix(fit$B[draw_idx,,], m, m)
        
        # For static model, A is identity
        A <- diag(m)
        
        # Compute impulse response
        Delta <- impulse_response_const(A, B, shock, H)
        
    } else {
        stop("IRF computation not yet implemented for model type: ", fit$model)
    }
    
    # Compute network statistic at each horizon
    irf_vals <- numeric(H + 1)
    
    for (h in 0:H) {
        # For baseline, we assume steady state (M)
        baseline <- NULL
        
        if (fit$model == "dynamic" && !is.null(fit$M)) {
            # Dynamic models store M as 4D array [m, m, p, n_draws]
            if (length(dim(fit$M)) == 4) {
                # Average over relations for network-level statistics
                baseline <- matrix(0, m, m)
                for (rel in 1:p) {
                    baseline <- baseline + fit$M[,,rel,draw_idx]
                }
                baseline <- baseline / p
            } else if (length(dim(fit$M)) == 3 && p == 1) {
                baseline <- fit$M[,,draw_idx]
            }
        } else if (fit$model == "static" && !is.null(fit$M)) {
            # Static models store M differently
            if (is.matrix(fit$M)) {
                baseline <- fit$M
            } else if (length(dim(fit$M)) == 3) {
                baseline <- fit$M[,,1]  # First relation
            }
        }
        
        if (is.null(baseline)) {
            baseline <- matrix(0, m, m)
        }
        
        # Compute statistic on shocked and baseline networks
        shocked_net <- baseline + Delta[,,h+1]
        baseline_val <- stat_fun(baseline)
        shocked_val <- stat_fun(shocked_net)
        
        irf_vals[h+1] <- shocked_val - baseline_val
    }
    
    irf_vals
}

#' Compute network-level impulse response functions
#'
#' @description Computes IRFs for network-level statistics given a shock
#' @param fit A dbn model fit object
#' @param shock Shock matrix or shock type (see build_shock)
#' @param H Number of horizons to compute (default: 20)
#' @param t0 Shock time for dynamic models (default: 1)
#' @param stat_fun Network statistic function (default: stat_density)
#' @param n_draws Number of posterior draws to use (default: all)
#' @param shock_pars List of parameters for build_shock if shock is character
#' @param ... Additional arguments passed to stat_fun
#' @return Data frame with IRF results including posterior summaries
#' @export
#' @examples
#' \dontrun{
#' # Unit edge shock to density
#' irf_density <- compute_irf(fit, shock = "unit_edge", 
#'                           shock_pars = list(i = 1, j = 2))
#' 
#' # Node shock to out-degrees
#' irf_outdeg <- compute_irf(fit, shock = "node_out",
#'                          shock_pars = list(i = 1),
#'                          stat_fun = stat_out_degree)
#' }
compute_irf <- function(fit, shock, H = 20, t0 = 1, 
                       stat_fun = stat_density,
                       n_draws = NULL,
                       shock_pars = list(),
                       ...) {
    
    # Validate inputs
    if (!inherits(fit, "dbn")) {
        stop("fit must be a dbn object")
    }
    
    if (!fit$model %in% c("static", "dynamic")) {
        stop("IRF currently only implemented for static and dynamic models")
    }
    
    dims <- fit$dims
    m <- dims$m
    
    # Build shock matrix if needed
    if (is.character(shock)) {
        shock_args <- c(list(m = m, type = shock), shock_pars)
        shock <- do.call(build_shock, shock_args)
    } else if (!is.matrix(shock) || nrow(shock) != m || ncol(shock) != m) {
        stop("shock must be an m x m matrix or a valid shock type")
    }
    
    # Determine number of draws
    if (is.null(n_draws)) {
        if (fit$model == "dynamic") {
            # Dynamic models store draws differently
            if (!is.null(fit$A) && is.list(fit$A)) {
                n_draws <- length(fit$A)
            } else if (!is.null(fit$sigma2)) {
                n_draws <- length(fit$sigma2)
            } else {
                stop("Cannot determine number of posterior draws for dynamic model")
            }
        } else if (fit$model == "static") {
            if (!is.null(fit$B)) {
                n_draws <- dim(fit$B)[1]
            } else if (!is.null(fit$params)) {
                n_draws <- nrow(fit$params)
            } else {
                stop("Cannot determine number of posterior draws for static model")
            }
        }
    }
    
    # Ensure n_draws doesn't exceed available draws
    available_draws <- if (fit$model == "dynamic" && !is.null(fit$A)) {
        length(fit$A)
    } else if (fit$model == "static" && !is.null(fit$B)) {
        dim(fit$B)[1]
    } else {
        n_draws
    }
    n_draws <- min(n_draws, available_draws)
    
    # For dynamic models, check time bounds
    if (fit$model == "dynamic") {
        T_total <- dims$T %||% dims$Tt %||% dims$TT
        if (is.null(T_total)) {
            stop("Cannot determine time dimension T from model")
        }
        if (t0 < 1 || t0 + H > T_total) {
            stop("For dynamic models, t0 + H must not exceed T (", T_total, ")")
        }
    }
    
    # Storage for IRF values
    irf_array <- matrix(NA, n_draws, H + 1)
    
    # Progress message
    cli::cli_progress_bar("Computing IRFs", total = n_draws)
    
    # Compute IRF for each posterior draw
    for (s in 1:n_draws) {
        cli::cli_progress_update()
        
        tryCatch({
            irf_array[s, ] <- compute_irf_single(fit, s, shock, H, t0, stat_fun)
        }, error = function(e) {
            # More detailed error message
            if (grepl("aes", e$message, ignore.case = TRUE)) {
                # This is likely a namespace collision issue
                warning("Error in draw ", s, ": ", e$message, 
                       "\nThis might be a namespace issue. Ensure ggplot2 is not masking functions.")
            } else {
                warning("Error in draw ", s, ": ", e$message)
            }
            # Print traceback for debugging if verbose
            if (getOption("dbn.debug", FALSE)) {
                cat("Traceback for draw", s, ":\n")
                traceback()
            }
        })
    }
    
    cli::cli_progress_done()
    
    # Remove failed draws
    valid_draws <- complete.cases(irf_array)
    if (sum(valid_draws) < n_draws) {
        warning(n_draws - sum(valid_draws), " draws failed and were removed")
        irf_array <- irf_array[valid_draws, , drop = FALSE]
    }
    
    # Check if we have any valid draws
    if (nrow(irf_array) == 0) {
        stop("All IRF computations failed. Run debug_irf() for diagnostics.")
    }
    
    # Compute posterior summaries
    result <- data.frame(
        horizon = 0:H,
        mean = colMeans(irf_array),
        median = apply(irf_array, 2, median),
        sd = apply(irf_array, 2, sd),
        q025 = apply(irf_array, 2, quantile, 0.025),
        q975 = apply(irf_array, 2, quantile, 0.975),
        q10 = apply(irf_array, 2, quantile, 0.10),
        q90 = apply(irf_array, 2, quantile, 0.90)
    )
    
    # Add attributes
    attr(result, "irf_draws") <- irf_array
    attr(result, "shock") <- shock
    # Use class name instead of deparse which might cause issues
    if (is.function(stat_fun)) {
        attr(result, "stat_fun") <- "custom_function"
    } else {
        attr(result, "stat_fun") <- deparse(substitute(stat_fun))
    }
    attr(result, "model") <- fit$model
    attr(result, "t0") <- t0
    
    class(result) <- c("dbn_irf", "data.frame")
    
    result
}

#' Plot impulse response functions
#'
#' @description Creates a plot of IRF with credible intervals
#' @param x A dbn_irf object from compute_irf
#' @param ci_level Credible interval level (default: 0.95)
#' @param title Plot title (default: auto-generated)
#' @param ... Additional arguments (ignored)
#' @return A ggplot2 object
#' @export
plot.dbn_irf <- function(x, ci_level = 0.95, title = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required for plotting")
    }
    
    # Determine which quantiles to use
    if (ci_level == 0.95) {
        q_low <- "q025"
        q_high <- "q975"
    } else if (ci_level == 0.90) {
        q_low <- "q10" 
        q_high <- "q90"
    } else {
        stop("ci_level must be 0.90 or 0.95")
    }
    
    # Default title
    if (is.null(title)) {
        stat_name <- attr(x, "stat_fun")
        title <- paste0("Impulse Response Function: ", stat_name)
    }
    
    # Create plot
    p <- ggplot2::ggplot(x, ggplot2::aes(x = horizon)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = .data[[q_low]], ymax = .data[[q_high]]),
            alpha = 0.3, fill = "gray"
        ) +
        ggplot2::geom_line(ggplot2::aes(y = mean), linewidth = 1) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        ggplot2::labs(
            title = title,
            x = "Horizon",
            y = "Response"
        ) +
        ggplot2::theme_minimal()
    
    p
}

#' Print method for dbn_irf objects
#'
#' @param x A dbn_irf object
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @export
print.dbn_irf <- function(x, digits = 3, ...) {
    cat("Network Impulse Response Function\n")
    cat("Model:", attr(x, "model"), "\n")
    cat("Statistic:", attr(x, "stat_fun"), "\n")
    if (attr(x, "model") == "dynamic") {
        cat("Shock time:", attr(x, "t0"), "\n")
    }
    cat("\nSummary:\n")
    print(round(as.data.frame(x), digits))
    invisible(x)
}

#' Debug method for irf function
#'
#' @param fit A dbn model fit object
#' @param draw_idx Index of posterior draw to debug (default: 1)
#' @param shock_i Source node index for shock (default: 8)
#' @param shock_j Target node index for shock (default: 12)
#' @export
debug_irf <- function(fit, draw_idx = 1, shock_i = 8, shock_j = 12) {
    cat("=== IRF Debug Information ===\n")
    
    # Check model structure
    cat("\n1. Model Structure:\n")
    cat("   Model type:", fit$model, "\n")
    cat("   Dimensions:", names(fit$dims), "\n")
    cat("   m =", fit$dims$m, "\n")
    cat("   p =", fit$dims$p, "\n")
    cat("   T =", fit$dims$T %||% fit$dims$Tt, "\n")
    
    # Check A and B storage
    cat("\n2. A and B matrices:\n")
    cat("   A is list:", is.list(fit$A), "\n")
    cat("   B is list:", is.list(fit$B), "\n")
    if (is.list(fit$A)) {
        cat("   Length of A:", length(fit$A), "\n")
        cat("   Dim of A[[1]]:", dim(fit$A[[1]]), "\n")
    }
    if (is.list(fit$B)) {
        cat("   Length of B:", length(fit$B), "\n")
        cat("   Dim of B[[1]]:", dim(fit$B[[1]]), "\n")
    }
    
    # Check M storage
    cat("\n3. M matrix:\n")
    cat("   M dimensions:", dim(fit$M), "\n")
    
    # Try to extract matrices for first draw
    cat("\n4. Attempting to extract matrices for draw", draw_idx, ":\n")
    
    tryCatch({
        m <- fit$dims$m
        p <- fit$dims$p
        
        # Extract A and B
        A_array <- fit$A[[draw_idx]]
        B_array <- fit$B[[draw_idx]]
        cat("   Successfully extracted A and B arrays\n")
        cat("   A_array dim:", dim(A_array), "\n")
        cat("   B_array dim:", dim(B_array), "\n")
        
        # Extract M
        if (length(dim(fit$M)) == 4) {
            M_draw <- fit$M[,,1,draw_idx]
            cat("   M for draw", draw_idx, "extracted, dim:", dim(M_draw), "\n")
        }
        
        # Create simple shock
        S <- matrix(0, m, m)
        S[shock_i, shock_j] <- 1
        cat("   Created shock matrix\n")
        
        # Test C++ function directly
        cat("\n5. Testing C++ impulse_response_dynamic:\n")
        Delta <- impulse_response_dynamic(A_array, B_array, S, 46, 5)
        cat("   C++ function succeeded, Delta dim:", dim(Delta), "\n")
        cat("   Delta[,,2] non-zero elements:", sum(Delta[,,2] != 0), "\n")
        
        # Test statistic function
        cat("\n6. Testing stat_density:\n")
        test_mat <- matrix(runif(m*m), m, m)
        dens <- stat_density(test_mat)
        cat("   stat_density on random matrix:", dens, "\n")
        
        # Try M + Delta
        baseline <- M_draw
        shocked <- baseline + Delta[,,1]
        dens_baseline <- stat_density(baseline)
        dens_shocked <- stat_density(shocked)
        cat("   Baseline density:", dens_baseline, "\n")
        cat("   Shocked density:", dens_shocked, "\n")
        cat("   Difference:", dens_shocked - dens_baseline, "\n")
        
    }, error = function(e) {
        cat("   ERROR:", e$message, "\n")
        cat("   Call stack:\n")
        print(sys.calls())
    })
    
    cat("\n=== End Debug ===\n")
}

# Run debug
# debug_irf(drc_conf_dyn)