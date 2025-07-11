#' Posterior Visualization Functions
#'
#' @description Basic plotting utilities for posterior analysis
#' @name posterior_plots
#' @keywords internal
NULL

#' Plot parameter trace plots
#'
#' @description Create trace plots for scalar parameters
#' @param fit A dbn model fit object
#' @param pars Character vector of parameter names to plot
#' @param ncol Number of columns for multi-panel plot
#' @return A ggplot object or NULL
#' @export
plot_trace <- function(fit, pars = NULL, ncol = 2) {
    # Get parameter summary
    param_df <- param_summary(fit, probs = c(0.025, 0.5, 0.975))

    if (is.null(param_df)) {
        warning("No scalar parameters found to plot")
        return(invisible(NULL))
    }

    # Select parameters to plot
    if (is.null(pars)) {
        pars <- param_df$parameter
    } else {
        pars <- intersect(pars, param_df$parameter)
        if (length(pars) == 0) {
            warning("None of the requested parameters found")
            return(invisible(NULL))
        }
    }

    # Extract trace data
    if (!is.null(fit$draws$pars)) {
        trace_data <- fit$draws$pars[, pars, drop = FALSE]
    } else {
        # Fallback to legacy format
        trace_data <- data.frame(row.names = seq_len(length(fit[[pars[1]]])))
        for (p in pars) {
            if (!is.null(fit[[p]])) {
                trace_data[[p]] <- fit[[p]]
            }
        }
    }

    # Add iteration column
    trace_data$iteration <- seq_len(nrow(trace_data))

    # Create long format data for ggplot2
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        # Reshape to long format
        df_long <- data.frame()
        for (p in pars) {
            if (p %in% names(trace_data)) {
                # Calculate running mean
                run_mean <- cumsum(trace_data[[p]]) / seq_along(trace_data[[p]])
                post_mean <- param_df[param_df$parameter == p, "mean"]

                df_param <- data.frame(
                    iteration = trace_data$iteration,
                    value = trace_data[[p]],
                    running_mean = run_mean,
                    posterior_mean = post_mean,
                    parameter = p,
                    type = "trace"
                )

                df_long <- rbind(df_long, df_param)
            }
        }

        # Create the plot
        p <- ggplot2::ggplot(df_long, ggplot2::aes(x = iteration)) +
            ggplot2::geom_line(ggplot2::aes(y = value), color = "gray40", alpha = 0.7) +
            ggplot2::geom_line(ggplot2::aes(y = running_mean), color = "red", linewidth = 1) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = posterior_mean),
                color = "blue", linetype = "dashed"
            ) +
            ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = ncol) +
            ggplot2::labs(
                title = "Parameter Trace Plots",
                x = "Iteration", y = "Value"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(strip.text = ggplot2::element_text(size = 10, face = "bold"))

        return(p)
    } else {
        # Fallback to base R
        n_pars <- length(pars)
        nrow <- ceiling(n_pars / ncol)

        if (n_pars > 1) {
            oldpar <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))
            on.exit(par(oldpar))
        }

        # Create trace plots
        for (p in pars) {
            if (p %in% names(trace_data)) {
                plot(trace_data[[p]],
                    type = "l",
                    main = p, xlab = "Iteration", ylab = "Value",
                    col = "gray40"
                )

                # Add running mean
                run_mean <- cumsum(trace_data[[p]]) / seq_along(trace_data[[p]])
                lines(run_mean, col = "red", lwd = 2)

                # Add horizontal line at posterior mean
                abline(
                    h = param_df[param_df$parameter == p, "mean"],
                    col = "blue", lty = 2
                )
            }
        }

        invisible(NULL)
    }
}

# #' Plot Theta heatmap
# #'
# #' @description Visualize Theta matrix as heatmap
# #' @param fit A dbn model fit object
# #' @param time Time point to plot
# #' @param rel Relation to plot
# #' @param fun Summary function (default: mean)
# #' @param ... Additional arguments to theta_summary
# #' @return A ggplot object or NULL
# #' @export
# plot_theta <- function(fit, time = 1, rel = 1, fun = mean, ...) {
#     # Get theta summary for specific time and relation
#     theta_sum <- theta_summary(fit, fun = fun, rel = rel, time = time, ...)

#     if (is.null(theta_sum) || nrow(theta_sum) == 0) {
#         warning("No Theta values found for specified indices")
#         return(invisible(NULL))
#     }

#     # Convert to matrix format
#     dims <- fit$meta$dims %||% fit$dims
#     m <- dims$m

#     theta_mat <- matrix(NA, m, m)
#     for (row in seq_len(nrow(theta_sum))) {
#         i <- theta_sum$i[row]
#         j <- theta_sum$j[row]
#         theta_mat[i, j] <- theta_sum$value[row]
#     }

#     # Check for valid values
#     if (all(is.na(theta_mat))) {
#         return(invisible(NULL))
#     }

#     # Replace infinite values with NA
#     theta_mat[!is.finite(theta_mat)] <- NA

#     # Create heatmap
#     if (requireNamespace("ggplot2", quietly = TRUE)) {
#         # Convert matrix to long format for ggplot2
#         df <- data.frame()
#         for (i in 1:m) {
#             for (j in 1:m) {
#                 if (!is.na(theta_mat[i, j])) {
#                     df <- rbind(df, data.frame(
#                         sender = i,
#                         receiver = j,
#                         value = theta_mat[i, j]
#                     ))
#                 }
#             }
#         }

#         # Create the heatmap
#         p <- ggplot2::ggplot(df, ggplot2::aes(x = receiver, y = sender, fill = value)) +
#             ggplot2::geom_tile() +
#             ggplot2::scale_fill_gradient2(
#                 low = "blue", mid = "white", high = "red",
#                 midpoint = 0, na.value = "grey90"
#             ) +
#             ggplot2::scale_y_reverse() + # Flip y-axis to match matrix convention
#             ggplot2::coord_equal() + # Square tiles
#             ggplot2::labs(
#                 title = paste0("Theta: Relation ", rel, ", Time ", time),
#                 x = "Receiver (j)",
#                 y = "Sender (i)",
#                 fill = expression(theta[ij])
#             ) +
#             ggplot2::theme_minimal() +
#             ggplot2::theme(panel.grid = ggplot2::element_blank())

#         return(p)
#     } else if (requireNamespace("lattice", quietly = TRUE)) {
#         # Fallback to lattice
#         print(lattice::levelplot(theta_mat,
#             main = paste0("Theta: Relation ", rel, ", Time ", time),
#             xlab = "Receiver (j)",
#             ylab = "Sender (i)",
#             col.regions = heat.colors(100),
#             aspect = "iso"
#         ))
#         invisible(NULL)
#     } else {
#         # Base R heatmap
#         image(theta_mat,
#             main = paste0("Theta: Relation ", rel, ", Time ", time),
#             xlab = "Receiver (j)", ylab = "Sender (i)",
#             col = heat.colors(100)
#         )
#         invisible(NULL)
#     }
# }


#' Plot regime probabilities (HMM only)
#'
#' @description Plot posterior regime probabilities over time
#' @param fit A dbn_hmm model fit object
#' @return A ggplot object or NULL
#' @export
plot_regime_probs <- function(fit) {
    # Get regime probabilities
    probs <- regime_probs(fit)

    if (is.null(probs)) {
        warning("Not an HMM model or no regime information found")
        return(invisible(NULL))
    }

    # Convert to long format
    times <- seq_len(nrow(probs))
    R <- ncol(probs)

    if (requireNamespace("ggplot2", quietly = TRUE)) {
        df <- data.frame()
        for (r in 1:R) {
            df <- rbind(df, data.frame(
                time = times,
                prob = probs[, r],
                regime = paste("Regime", r)
            ))
        }

        p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = prob, fill = regime)) +
            ggplot2::geom_area(alpha = 0.7) +
            ggplot2::scale_fill_brewer(palette = "Set2") +
            ggplot2::labs(
                title = "Posterior Regime Probabilities",
                x = "Time", y = "Probability", fill = NULL
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "bottom")

        return(p)
    } else {
        # Base R stacked area plot
        plot(times, probs[, 1],
            type = "n", ylim = c(0, 1),
            main = "Posterior Regime Probabilities",
            xlab = "Time", ylab = "Probability"
        )

        # Build cumulative probabilities
        y_bottom <- rep(0, length(times))
        cols <- rainbow(R)

        for (r in 1:R) {
            y_top <- y_bottom + probs[, r]
            polygon(c(times, rev(times)),
                c(y_bottom, rev(y_top)),
                col = cols[r], border = NA
            )
            y_bottom <- y_top
        }

        legend("topright", paste("Regime", 1:R), fill = cols)

        invisible(NULL)
    }
}
