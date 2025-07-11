# note: generic plot.dbn method has been moved to zzz-dispatch.r
# which provides routing to all model types (static, dynamic, lowrank, hmm)

#' Plot Static DBN Results
#'
#' @description Creates diagnostic plots for static model results using ggplot2
#' @param results Output from dbn_static()
#' @param alpha Significance level for edge detection (default 0.01)
#' @return A ggplot2 object with multiple panels
#' @keywords internal
#' @import ggplot2 gridExtra
plot_static <- function(results, alpha = 0.01) {
    # define variables to avoid r cmd check notes
    iteration <- value <- parameter <- s2 <- t2 <- g2 <- NULL
    from_x <- from_y <- to_x <- to_y <- color <- x <- y <- label <- NULL

    # check for null params
    if (is.null(results$params)) {
        cli::cli_abort("results$params is NULL -- nothing to plot")
    }

    # prepare data for parameter traces
    params_df <- as.data.frame(results$params)
    params_df$iteration <- seq_len(nrow(params_df))

    # reshape for plotting
    params_long <- data.frame(
        iteration = rep(params_df$iteration, 3),
        value = c(params_df$s2, params_df$t2, params_df$g2),
        parameter = factor(rep(c("s2", "t2", "g2"), each = nrow(params_df)))
    )

    # create parameter trace plot
    p_traces <- ggplot2::ggplot(params_long, ggplot2::aes(x = iteration, y = value, color = parameter)) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 1) +
        ggplot2::labs(title = "Parameter Traces", x = "Iteration", y = "Value") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")

    # create histogram plots
    p_hist_s2 <- ggplot2::ggplot(params_df, ggplot2::aes(x = s2)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
        ggplot2::labs(title = "s2 Posterior", x = "s2", y = "Count") +
        ggplot2::theme_minimal()

    p_hist_t2 <- ggplot2::ggplot(params_df, ggplot2::aes(x = t2)) +
        ggplot2::geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
        ggplot2::labs(title = "t2 Posterior", x = "t2", y = "Count") +
        ggplot2::theme_minimal()

    p_hist_g2 <- ggplot2::ggplot(params_df, ggplot2::aes(x = g2)) +
        ggplot2::geom_histogram(bins = 30, fill = "darkred", alpha = 0.7) +
        ggplot2::labs(title = "g2 Posterior", x = "g2", y = "Count") +
        ggplot2::theme_minimal()

    # network plots for first b matrix
    if (length(results$B) >= 1 && !is.null(results$B[[1]])) {
        B_post <- results$B[[1]]
        EB <- apply(B_post, c(1, 2), mean)
        diag(EB) <- 0
        SB <- apply(B_post, c(1, 2), sd)

        # identify significant edges
        BIG <- abs(EB) > quantile(abs(EB), 1 - alpha, na.rm = TRUE)
        SIG <- abs(EB) > qnorm(1 - alpha / 2) * SB
        # handle case where all sb==0
        if (all(SB == 0, na.rm = TRUE)) {
            SIG[] <- FALSE
        }
        BSG <- BIG & SIG

        # filter nodes with connections
        csig <- (rowSums(BSG) + colSums(BSG)) > 0

        if (sum(csig) > 1) {
            # create network data
            BSG_sub <- BSG[csig, csig]
            EB_sub <- EB[csig, csig]
            nodes_sub <- which(csig)

            # get layout - simple circular layout ... this looks sooooo dumb
            # should just use netify
            n_nodes <- nrow(BSG_sub)
            angles <- seq(0, 2 * pi, length.out = n_nodes + 1)[1:n_nodes]
            xy <- cbind(x = cos(angles), y = sin(angles))

            # create edge data
            edges <- data.frame()
            for (i in 1:nrow(BSG_sub)) {
                for (j in 1:ncol(BSG_sub)) {
                    if (BSG_sub[i, j]) {
                        edges <- rbind(edges, data.frame(
                            from_x = xy[i, 1], from_y = xy[i, 2],
                            to_x = xy[j, 1], to_y = xy[j, 2],
                            weight = EB_sub[i, j],
                            color = ifelse(EB_sub[i, j] > 0, "Positive", "Negative")
                        ))
                    }
                }
            }

            # create node data
            nodes <- data.frame(
                x = xy[, 1], y = xy[, 2],
                label = nodes_sub
            )

            p_network <- ggplot2::ggplot() +
                ggplot2::geom_segment(
                    data = edges,
                    ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y, color = color),
                    arrow = grid::arrow(length = grid::unit(0.2, "cm"))
                ) +
                ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y), size = 4) +
                ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y, label = label), size = 3) +
                ggplot2::scale_color_manual(values = c("Positive" = "green3", "Negative" = "red3")) +
                ggplot2::labs(title = "Network B[[1]]", color = "Edge Type") +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    axis.text = ggplot2::element_blank(),
                    axis.title = ggplot2::element_blank(),
                    panel.grid = ggplot2::element_blank()
                ) +
                ggplot2::coord_fixed()
        } else {
            p_network <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No significant edges for B[[1]]") +
                ggplot2::theme_void()
        }
    } else {
        cli::cli_warn("No B samples available - skipping network panel")
        p_network <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No B samples available") +
            ggplot2::theme_void()
    }

    # combine plots using gridextra if available, otherwise return list
    if (requireNamespace("gridExtra", quietly = TRUE)) {
        tryCatch(
            {
                gridExtra::grid.arrange(p_network, p_traces, p_hist_s2, p_hist_t2, p_hist_g2,
                    layout_matrix = rbind(c(1, 1, 2), c(3, 4, 5))
                )
            },
            error = function(e) {
                # if grid.arrange fails, return the plots as a list
                cli::cli_warn("Could not arrange plots with gridExtra: {e$message}")
                list(
                    network = p_network, traces = p_traces,
                    hist_s2 = p_hist_s2, hist_t2 = p_hist_t2, hist_g2 = p_hist_g2
                )
            }
        )
    } else {
        list(
            network = p_network, traces = p_traces,
            hist_s2 = p_hist_s2, hist_t2 = p_hist_t2, hist_g2 = p_hist_g2
        )
    }
}

#' Plot Dynamic DBN Results
#'
#' @description Creates diagnostic plots for dynamic model results using ggplot2
#' @param results Output from dbn_dynamic()
#' @param time_points Which time points to display (default: start, middle, end)
#' @return A list of ggplot2 objects
#' @keywords internal
#' @import ggplot2 gridExtra
plot_dynamic <- function(results, time_points = NULL) {
    # define variables to avoid r cmd check notes
    iteration <- value <- parameter <- g2 <- NULL

    # get dims
    dims <- results$meta$dims %||% results$dims
    Tt <- dims$Tt %||% dims$T

    if (is.null(time_points)) {
        time_points <- c(1, floor(Tt / 2), Tt)
    }

    # extract params
    if (!is.null(results$draws$pars)) {
        # new format
        pars <- results$draws$pars
        trace_df <- data.frame(
            iteration = seq_len(nrow(pars)),
            sigma2 = pars$sigma2,
            tauA2 = pars$tau_A2,
            tauB2 = pars$tau_B2
        )
        if ("g2" %in% names(pars)) trace_df$g2 <- pars$g2
        if ("rho_A" %in% names(pars)) trace_df$rhoA <- pars$rho_A
        if ("rho_B" %in% names(pars)) trace_df$rhoB <- pars$rho_B
    } else {
        # 
        trace_df <- data.frame(
            iteration = seq_along(results$sigma2),
            sigma2 = results$sigma2,
            tauA2 = results$tau_A2 %||% results$tauA2,
            tauB2 = results$tau_B2 %||% results$tauB2
        )
        if (!is.null(results$g2)) trace_df$g2 <- results$g2
        if (!is.null(results$rho_A)) trace_df$rhoA <- results$rho_A
        if (!is.null(results$rho_B)) trace_df$rhoB <- results$rho_B
    }

    # reshape for plotting
    plist <- c("sigma^2", "tauA^2", "tauB^2")
    trace_long <- data.frame(
        iteration = rep(trace_df$iteration, 3),
        value = c(trace_df$sigma2, trace_df$tauA2, trace_df$tauB2),
        parameter = factor(rep(plist, each = nrow(trace_df)))
    )

    # 
    if (!is.null(results$g2)) {
        trace_long <- rbind(
            trace_long,
            data.frame(
                iteration = trace_df$iteration,
                value = trace_df$g2,
                parameter = "g^2"
            )
        )
        plist <- c(plist, "g^2")
    }

    p_traces <- ggplot2::ggplot(trace_long, ggplot2::aes(x = iteration, y = value)) +
        ggplot2::geom_line(color = "steelblue") +
        ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 1) +
        ggplot2::labs(title = "Parameter Traces", x = "Iteration", y = "Value") +
        ggplot2::theme_minimal()

    # 
    if (!is.null(results$draws$misc$A)) {
        last_A <- results$draws$misc$A[[length(results$draws$misc$A)]]
    } else {
        last_A <- results$A[[length(results$A)]]
    }

    # 
    time_idx_avail <- seq_len(dim(last_A)[3])
    time_points <- intersect(time_points, time_idx_avail)
    if (length(time_points) == 0) {
        # 
        time_points <- time_idx_avail[c(1, floor(length(time_idx_avail) / 2), length(time_idx_avail))]
    }

    # 
    A_hist_data <- data.frame()
    for (i in seq_along(time_points)) {
        t_idx <- time_points[i]
        # 
        orig_t <- if (!is.null(results$settings$time_thin) && results$settings$time_thin > 1) {
            (t_idx - 1) * results$settings$time_thin + 1
        } else {
            t_idx
        }
        A_hist_data <- rbind(A_hist_data, data.frame(
            value = as.vector(last_A[, , t_idx]),
            time = paste("t =", orig_t)
        ))
    }

    p_A_hist <- ggplot2::ggplot(A_hist_data, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(bins = 20, fill = "darkgreen", alpha = 0.7) +
        ggplot2::facet_wrap(~time, ncol = length(time_points)) +
        ggplot2::labs(title = "Distribution of A Matrix Elements", x = "Value", y = "Count") +
        ggplot2::theme_minimal()

    # ar1 parameters if present
    if ("rhoA" %in% names(trace_df) && "rhoB" %in% names(trace_df)) {
        rho_df <- data.frame(
            iteration = trace_df$iteration,
            rhoA = trace_df$rhoA,
            rhoB = trace_df$rhoB
        )

        rho_long <- data.frame(
            iteration = rep(rho_df$iteration, 2),
            value = c(rho_df$rhoA, rho_df$rhoB),
            parameter = factor(rep(c("rhoA", "rhoB"), each = nrow(rho_df)))
        )

        p_rho <- ggplot2::ggplot(rho_long, ggplot2::aes(x = iteration, y = value, color = parameter)) +
            ggplot2::geom_line() +
            ggplot2::facet_wrap(~parameter, ncol = 2) +
            ggplot2::labs(title = "AR(1) Parameters", x = "Iteration", y = "Value") +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none")
    } else if ("rhoA" %in% names(trace_df)) {
        # only rhoa available
        p_rho <- ggplot2::ggplot(
            data.frame(iteration = trace_df$iteration, rhoA = trace_df$rhoA),
            ggplot2::aes(x = iteration, y = rhoA)
        ) +
            ggplot2::geom_line(color = "darkblue") +
            ggplot2::labs(title = "AR(1) Parameter - rhoA", x = "Iteration", y = "rhoA") +
            ggplot2::theme_minimal()
    } else if ("rhoB" %in% names(trace_df)) {
        # only rhob available
        p_rho <- ggplot2::ggplot(
            data.frame(iteration = trace_df$iteration, rhoB = trace_df$rhoB),
            ggplot2::aes(x = iteration, y = rhoB)
        ) +
            ggplot2::geom_line(color = "darkred") +
            ggplot2::labs(title = "AR(1) Parameter - rhoB", x = "Iteration", y = "rhoB") +
            ggplot2::theme_minimal()
    } else {
        p_rho <- NULL
    }

    # g2 trace if available
    if ("g2" %in% names(trace_df)) {
        p_g2 <- ggplot2::ggplot(
            data.frame(iteration = trace_df$iteration, g2 = trace_df$g2),
            ggplot2::aes(x = iteration, y = g2)
        ) +
            ggplot2::geom_line(color = "darkred") +
            ggplot2::labs(title = "g^2 (tau_mu^2) Trace", x = "Iteration", y = "g^2") +
            ggplot2::theme_minimal()
    } else {
        p_g2 <- NULL
    }

    #  list of plots
    plots <- list(traces = p_traces, A_hist = p_A_hist)
    if (!is.null(p_rho)) plots$rho <- p_rho
    if (!is.null(p_g2)) plots$g2 <- p_g2

    # if gridextra available, arrange them
    if (requireNamespace("gridExtra", quietly = TRUE)) {
        do.call(gridExtra::grid.arrange, c(plots, ncol = 1))
    } else {
        plots
    }
}

#' Summary for static DBN fits
#'
#' @description Prints summary statistics for static DBN model results
#' @param object Object of class "dbn" with model="static"
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisible object
#' @keywords internal
summary_static <- function(object, digits = 3, ...) {
    if (object$model != "static") stop("Not a static model fit")

    cli::cli_h1("Static Bilinear Network Model")
    cli::cli_h3("Dimensions")

    # 
    dims <- object$meta$dims %||% object$dims
    cli::cli_bullets(c(
        " " = "Nodes: {dims$m}",
        " " = "Relations: {dims$p}",
        " " = "Time points: {dims$n %||% dims$Tt}"
    ))

    cli::cli_h3("Parameter estimates (mean [95% CI])")

    # 
    if (exists("param_summary")) {
        param_df <- param_summary(object, probs = c(0.025, 0.975))
        if (!is.null(param_df)) {
            for (i in seq_len(nrow(param_df))) {
                par <- param_df$parameter[i]
                # 
                if ("q2.5" %in% names(param_df)) {
                    cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q2.5[i], digits)}, {round(param_df$q97.5[i], digits)}]")
                } else if ("q5" %in% names(param_df)) {
                    cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q5[i], digits)}, {round(param_df$q95[i], digits)}]")
                } else if ("q50" %in% names(param_df)) {
                    cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q5[i], digits)}, {round(param_df$q95[i], digits)}]")
                } else {
                    # fallback to just mean
                    cli::cli_inform("  {par}: {round(param_df$mean[i], digits)}")
                }
            }
        }
    } else {
        # 
        params <- object$params
        for (par in colnames(params)) {
            vals <- params[, par]
            cli::cli_inform("  {par}: {round(mean(vals, na.rm = TRUE), digits)} [{round(quantile(vals, 0.025, na.rm = TRUE), digits)}, {round(quantile(vals, 0.975, na.rm = TRUE), digits)}]")
        }
    }

    cli::cli_h3("Settings")
    settings <- object$meta %||% object$settings
    for (s in names(settings)) {
        if (s != "dims") { # Skip dims as we already showed them
            cli::cli_inform("  {s}: {settings[[s]]}")
        }
    }

    invisible(object)
}

#' Summary for dynamic DBN fits
#'
#' @description Prints summary statistics for dynamic DBN model results
#' @param object Object of class "dbn" with model="dynamic"
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisible object
#' @keywords internal
summary_dynamic <- function(object, digits = 3, ...) {
    if (object$model != "dynamic") stop("Not a dynamic model fit")

    cli::cli_h1("Dynamic Bilinear Network Model")
    cli::cli_h3("Dimensions")

    # get dims
    dims <- object$meta$dims %||% object$dims
    cli::cli_bullets(c(
        " " = "Nodes: {dims$m}",
        " " = "Relations: {dims$p}",
        " " = "Time points: {dims$Tt %||% dims$T}"
    ))

    cli::cli_h3("Parameter estimates (mean [95% CI])")

    # 
    if (exists("param_summary")) {
        param_df <- param_summary(object, probs = c(0.025, 0.975))
        if (!is.null(param_df)) {
            for (i in seq_len(nrow(param_df))) {
                par <- param_df$parameter[i]
                # 
                if ("q2.5" %in% names(param_df)) {
                    cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q2.5[i], digits)}, {round(param_df$q97.5[i], digits)}]")
                } else if ("q5" %in% names(param_df)) {
                    cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q5[i], digits)}, {round(param_df$q95[i], digits)}]")
                } else if ("q50" %in% names(param_df)) {
                    cli::cli_inform("  {par}: {round(param_df$mean[i], digits)} [{round(param_df$q5[i], digits)}, {round(param_df$q95[i], digits)}]")
                } else {
                    # 
                    cli::cli_inform("  {par}: {round(param_df$mean[i], digits)}")
                }
            }
        }
    } else {
        # 
        for (par in c("sigma2", "tau_A2", "tau_B2", "g2")) {
            if (!is.null(object[[par]])) {
                vals <- object[[par]]
                cli::cli_inform("  {par}: {round(mean(vals, na.rm = TRUE), digits)} [{round(quantile(vals, 0.025, na.rm = TRUE), digits)}, {round(quantile(vals, 0.975, na.rm = TRUE), digits)}]")
            }
        }

        if (!is.null(object$rho_A)) {
            cli::cli_inform("  rho_A: {round(mean(object$rho_A, na.rm = TRUE), digits)} [{round(quantile(object$rho_A, 0.025, na.rm = TRUE), digits)}, {round(quantile(object$rho_A, 0.975, na.rm = TRUE), digits)}]")
            cli::cli_inform("  rho_B: {round(mean(object$rho_B, na.rm = TRUE), digits)} [{round(quantile(object$rho_B, 0.025, na.rm = TRUE), digits)}, {round(quantile(object$rho_B, 0.975, na.rm = TRUE), digits)}]")
        }
    }

    cli::cli_h3("Settings")
    settings <- object$meta %||% object$settings
    for (s in names(settings)) {
        if (s != "dims") { # Skip dims as we already showed them
            cli::cli_inform("  {s}: {settings[[s]]}")
        }
    }

    invisible(object)
}

#' Check MCMC Convergence
#'
#' @description Provides convergence diagnostics for MCMC chains
#' @param results Output from dbn()
#' @return Invisible NULL (diagnostics are printed and plotted)
#' @export
check_convergence <- function(results) {
    if (!requireNamespace("coda", quietly = TRUE)) {
        cli::cli_inform(c(
            "!" = "Package 'coda' is required for full convergence diagnostics.",
            "i" = "Install it with: {.code install.packages('coda')}",
            " " = "Showing basic diagnostics instead..."
        ))

        # basic diagnostics without coda
        if (results$model == "static") {
            params <- results$params
        } else {
            params <- cbind(
                sigma2 = results$sigma2,
                tauA2 = results$tauA2,
                tauB2 = results$tauB2
            )
            if (!is.null(results$g2)) {
                params <- cbind(params, g2 = results$g2)
            }
        }

        cli::cli_h3("Parameter Summary")
        print(summary(params))

        cli::cli_h3("Autocorrelation at lag 1")
        for (j in 1:ncol(params)) {
            ac1 <- cor(params[-nrow(params), j], params[-1, j])
            cli::cli_inform("  {colnames(params)[j]}: {round(ac1, 3)}")
        }

        return(invisible(NULL))
    }

    if (results$model == "static") {
        params_mcmc <- coda::mcmc(results$params)
    } else {
        params_df <- cbind(
            sigma2 = results$sigma2,
            tauA2 = results$tauA2,
            tauB2 = results$tauB2
        )
        if (!is.null(results$g2)) {
            params_df <- cbind(params_df, g2 = results$g2)
        }
        params_mcmc <- coda::mcmc(params_df)
    }

    # effective sample size
    cli::cli_h3("Effective Sample Sizes")
    print(coda::effectiveSize(params_mcmc))

    # geweke diagnostic
    cli::cli_h3("Geweke Diagnostic")
    print(coda::geweke.diag(params_mcmc))

    # plot autocorrelations
    par(mfrow = c(2, 2))
    coda::autocorr.plot(params_mcmc, auto.layout = FALSE)

    invisible(NULL)
}

#' Compare Multiple DBN Models
#'
#' @description Creates comparative plots for multiple DBN results using ggplot2
#' @param ... Multiple dbn objects to compare
#' @return A ggplot2 object or list of plots
#' @export
#' @import ggplot2 gridExtra
compare_dbn <- function(...) {
    # define variables to avoid r cmd check notes
    iteration <- value <- model <- NULL

    results_list <- list(...)
    n_models <- length(results_list)

    if (n_models < 2) {
        cli::cli_abort("Need at least 2 models to compare")
    }

    # prep comparison data
    compare_df <- data.frame()

    for (i in 1:n_models) {
        res <- results_list[[i]]

        if (res$model == "static") {
            temp_df <- data.frame(
                iteration = seq_len(nrow(res$params)),
                value = res$params[, "s2"],
                parameter = "s2",
                model = paste("Model", i)
            )
        } else {
            temp_df <- data.frame(
                iteration = seq_along(res$sigma2),
                value = res$sigma2,
                parameter = "sigma^2",
                model = paste("Model", i)
            )
        }

        compare_df <- rbind(compare_df, temp_df)
    }

    # get comparison plot
    p_compare <- ggplot2::ggplot(compare_df, ggplot2::aes(x = iteration, y = value, color = model)) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~parameter, scales = "free") +
        ggplot2::labs(title = "Model Comparison", x = "Iteration", y = "Value", color = "Model") +
        ggplot2::theme_minimal()

    # create density plots
    p_density <- ggplot2::ggplot(compare_df, ggplot2::aes(x = value, fill = model)) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::facet_wrap(~parameter, scales = "free") +
        ggplot2::labs(title = "Parameter Distributions", x = "Value", y = "Density", fill = "Model") +
        ggplot2::theme_minimal()

    if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(p_compare, p_density, ncol = 1)
    } else {
        list(traces = p_compare, density = p_density)
    }
}

#' Posterior Predictive ECDF Overlay
#'
#' @description Creates posterior predictive checks using ECDF overlays
#' @param fit Object of class "dbn"
#' @param n_rep Number of replicated datasets to draw (default: 20)
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
ppc_ecdf <- function(fit, n_rep = 20) {
    if (!inherits(fit, "dbn")) {
        cli::cli_abort("{.arg fit} must be of class {.cls dbn}")
    }
    if (!fit$model %in% c("static", "dynamic")) {
        cli::cli_abort("Model must be either 'static' or 'dynamic'")
    }

    # check if we need to load tprod
    if (!exists("tprod")) {
        cli::cli_abort("tprod function not found. Please load the dbn package.")
    }

    # orig data
    y_orig <- fit$R
    ecdf_orig <- ecdf(c(y_orig))
    stat_orig <- ecdf_orig(sort(unique(c(y_orig))))

    # helper to simulate one replicate
    draw_rep <- function() {
        if (fit$model == "static") {
            idx <- sample(seq_len(dim(fit$B[[1]])[3]), 1)
            Bsim <- lapply(fit$B, function(x) x[, , idx])
            # simple mean structure
            Yhat <- tprod(fit$M, Bsim) + array(rnorm(length(fit$M)), dim = dim(fit$M))
        } else {
            idx <- sample(seq_along(fit$A), 1)
            A <- fit$A[[idx]]
            B <- fit$B[[idx]]
            Theta <- array(0, dim = c(fit$dims$m, fit$dims$m, fit$dims$p, dim(A)[3]))
            # fwd simulate one theta path
            Theta[, , , 1] <- 0
            for (t in 2:dim(Theta)[4]) {
                for (rel in 1:fit$dims$p) {
                    Theta[, , rel, t] <- A[, , t] %*% Theta[, , rel, t - 1] %*% t(B[, , t]) +
                        sqrt(fit$sigma2[idx]) * matrix(rnorm(fit$dims$m^2), fit$dims$m)
                }
            }
            Yhat <- sweep(Theta, 1:3, fit$M[[idx]], "+")
        }
        # guard against values outside observed support
        vals <- sort(unique(c(y_orig)))
        Yrep <- pmax(min(vals), pmin(max(vals), round(Yhat)))
        ecdf(c(Yrep))(vals)
    }

    reps <- replicate(n_rep, draw_rep())

    df <- data.frame(
        quant = rep(sort(unique(c(y_orig))), n_rep + 1),
        ecdf = c(stat_orig, as.vector(reps)),
        sel = rep(c("Observed", paste0("Rep", 1:n_rep)), each = length(stat_orig))
    )

    ggplot2::ggplot(df, ggplot2::aes(
        x = quant, y = ecdf, group = sel,
        colour = sel == "Observed"
    )) +
        ggplot2::geom_line(alpha = 0.5) +
        ggplot2::scale_color_manual(values = c("grey70", "black"), guide = "none") +
        ggplot2::labs(
            title = "Posterior predictive ECDF overlay",
            x = "Ordinal category", y = "ECDF"
        ) +
        ggplot2::theme_minimal()
}

#' Plot Dyad Trajectory
#'
#' @description Plots posterior mean and 95% bands for a single dyad through time
#' @param fit Dynamic dbn object
#' @param i Actor i index
#' @param j Actor j index
#' @param rel Relation indices (default: NULL = all relations)
#' @param facet Whether to facet by relation (default: TRUE)
#' @param cred Credible interval quantiles (default: c(0.025, 0.975))
#' @return A ggplot2 object
#' @export
#' @import ggplot2
dyad_path <- function(fit, i, j, rel = NULL, facet = TRUE, cred=c(0.025, 0.975)) {
    if (!fit$model %in% c("dynamic", "lowrank", "hmm")) {
        cli::cli_abort("This function requires a time-varying model (dynamic, lowrank, or hmm)")
    }
    p <- fit$dims$p

    # default to all relations ... risky if someone is crazy
    if (is.null(rel)) {
        rel <- seq_len(p)
    }

    # helper function to compute trajectory for one relation
    compute_trajectory <- function(r) {
        if (fit$model == "dynamic") {
            n_keep <- length(fit$A)
            Tt <- dim(fit$A[[1]])[3]
            thetas <- matrix(NA, n_keep, Tt)

            for (s in seq_len(n_keep)) {
                for (t in seq_len(Tt)) {
                    thetas[s, t] <- sum(fit$A[[s]][i, , t] * fit$B[[s]][j, , t])
                }
            }
        } else if (fit$model == "lowrank") {
            n_keep <- length(fit$U)
            Tt <- ncol(fit$alpha[[1]])
            thetas <- matrix(NA, n_keep, Tt)

            for (s in seq_len(n_keep)) {
                for (t in seq_len(Tt)) {
                    U_s <- fit$U[[s]]
                    alpha_s <- fit$alpha[[s]][, t]
                    A_s <- U_s %*% diag(alpha_s) %*% t(U_s)
                    B_s <- fit$B[[s]][, , t]
                    thetas[s, t] <- sum(A_s[i, ] * B_s[j, ])
                }
            }
        } else if (fit$model == "hmm") {
            n_keep <- length(fit$S)
            Tt <- length(fit$S[[1]])
            thetas <- matrix(NA, n_keep, Tt)

            for (s in seq_len(n_keep)) {
                for (t in seq_len(Tt)) {
                    regime <- fit$S[[s]][t]
                    A_s <- fit$A[[s]][, , regime]
                    B_s <- fit$B[[s]][, , regime]
                    thetas[s, t] <- sum(A_s[i, ] * B_s[j, ])
                }
            }
        }

        # handle time_thin
        time_vals <- 1:Tt
        if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
            time_vals <- (time_vals - 1) * fit$settings$time_thin + 1
        }

        data.frame(
            relation = paste0("Relation ", r),
            time = time_vals,
            mean = colMeans(thetas),
            lo = apply(thetas, 2, quantile, cred[1]),
            hi = apply(thetas, 2, quantile, cred[2])
        )
    }

    # compute for all requested relations
    df_all <- do.call(rbind, lapply(rel, compute_trajectory))

    # viz
    g <- ggplot2::ggplot(df_all, ggplot2::aes(time, mean)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80") +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::labs(
            x = "Time", y = expression(theta[list(i, j, t)]),
            title = sprintf("Dyad (%d,%d) trajectory", i, j)
        ) +
        ggplot2::theme_minimal()

    # add faceting if multiple relations and facet = true
    if (facet && length(rel) > 1) {
        g <- g + ggplot2::facet_wrap(~relation, scales = "free_y")
    }

    #
    g
}

#' Track Role Evolution
#'
#' @description Tracks first left/right singular vector of A_t or B_t
#' @param fit Dynamic dbn object
#' @param mat "A" or "B"
#' @param comp Component index (default: 1)
#' @return Base R plot (invisible NULL)
#' @export
role_trajectory <- function(fit, mat = c("A", "B"), comp = 1) {
    mat <- match.arg(mat)
    n_keep <- length(fit[[mat]])
    Tt <- dim(fit[[mat]][[1]])[3]
    m <- fit$dims$m
    scores <- matrix(NA, Tt, m)

    # Posterior mean at each time
    for (t in 1:Tt) {
        Mbar <- Reduce(`+`, lapply(fit[[mat]], function(M) M[, , t])) / n_keep
        sv <- svd(Mbar)
        scores[t, ] <- if (mat == "A") sv$u[, comp] else sv$v[, comp]
    }

    matplot(1:Tt, scores,
        type = "l", lty = 1, col = 1:m,
        main = sprintf("%s_t singular vector %d", mat, comp),
        xlab = "Time", ylab = "Coordinate"
    )
    legend("topright", legend = 1:m, col = 1:m, lty = 1, cex = 0.6, ncol = 2)

    invisible(NULL)
}

#' Network Snapshot
#'
#' @description Heat map of Theta at given time (posterior mean)
#' @param fit Dynamic dbn object
#' @param t Time point
#' @param rel Relation index (default: 1)
#' @param sparse Auto-switch to sparse visualization for large networks
#' @param eps Threshold for sparse plotting
#' @param show_significant Logical, whether to show only significant effects (default: FALSE)
#' @param cred_level Credible level for significance (default corresponds to 95% CI)
#' @return ggplot2 object or base R plot
#' @export
#' @import ggplot2
net_snapshot <- function(fit, t, rel = 1, sparse = NULL, eps = 1e-4, 
                        show_significant = FALSE, cred_level = 0.025) {
    if (!fit$model %in% c("dynamic", "lowrank", "hmm")) {
        cli::cli_abort("This function requires a time-varying model (dynamic, lowrank, or hmm)")
    }
    m <- fit$dims$m
    Th <- matrix(0, m, m)
    
    # Initialize matrix to store all posterior draws if significance testing is requested
    if (show_significant) {
        Th_all <- array(NA, dim = c(m, m, length(fit$A)))
    }

    # auto-detect sparse mode
    if (is.null(sparse)) {
        sparse <- m > 150
    }

    # handle time_thin
    t_idx <- t
    if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
        t_idx <- ceiling(t / fit$settings$time_thin)
    }

    # compute Theta based on model type
    if (fit$model == "dynamic") {
        n_keep <- length(fit$A)
        if (t_idx > dim(fit$A[[1]])[3]) {
            cli::cli_abort("Time {t} not available (time_thin={fit$settings$time_thin})")
        }
        for (s in seq_len(n_keep)) {
            Th_s <- fit$A[[s]][, , t_idx] %*% t(fit$B[[s]][, , t_idx])
            if (show_significant) {
                Th_all[, , s] <- Th_s
            }
            Th <- Th + Th_s
        }
        Th <- Th / n_keep
    } else if (fit$model == "lowrank") {
        # for lowrank: A_t = U * diag(alpha_t) * U^T
        n_keep <- length(fit$U)
        for (s in seq_len(n_keep)) {
            U_s <- fit$U[[s]]
            alpha_s <- fit$alpha[[s]][, t_idx]
            A_s <- U_s %*% diag(alpha_s) %*% t(U_s)
            B_s <- fit$B[[s]][, , t_idx]
            Th_s <- A_s %*% t(B_s)
            if (show_significant) {
                Th_all[, , s] <- Th_s
            }
            Th <- Th + Th_s
        }
        Th <- Th / n_keep
    } else if (fit$model == "hmm") {
        # for HMM: use regime at time t
        n_keep <- length(fit$S)
        for (s in seq_len(n_keep)) {
            regime <- fit$S[[s]][t_idx]
            A_s <- fit$A[[s]][, , regime]
            B_s <- fit$B[[s]][, , regime]
            Th_s <- A_s %*% t(B_s)
            if (show_significant) {
                Th_all[, , s] <- Th_s
            }
            Th <- Th + Th_s
        }
        Th <- Th / n_keep
    }
    
    # Compute significance if requested
    if (show_significant) {
        # Compute credible intervals for each cell
        alpha <- 1 - cred_level
        lower_q <- alpha / 2
        upper_q <- 1 - alpha / 2
        
        # For each cell, check if credible interval includes zero
        is_significant <- matrix(TRUE, m, m)
        for (i in 1:m) {
            for (j in 1:m) {
                cell_vals <- Th_all[i, j, ]
                ci <- quantile(cell_vals, probs = c(lower_q, upper_q), na.rm = TRUE)
                # If CI includes zero, not significant
                if (ci[1] <= 0 && ci[2] >= 0) {
                    is_significant[i, j] <- FALSE
                    Th[i, j] <- NA  # Set to NA so it will be greyed out
                }
            }
        }
    }

    if (sparse) {
        # sparse visualization using ggplot2
        # threshold small values (but preserve NAs for significance)
        if (!show_significant) {
            Th[abs(Th) < eps] <- 0
        } else {
            # For significant mode, only threshold non-NA values
            non_na_mask <- !is.na(Th)
            Th[non_na_mask & abs(Th) < eps] <- 0
        }

        # find non-zero entries (including NAs if show_significant)
        if (show_significant) {
            idx <- which(Th != 0 | is.na(Th), arr.ind = TRUE)
        } else {
            idx <- which(Th != 0, arr.ind = TRUE)
        }
        if (nrow(idx) == 0) {
            cli::cli_warn("No edges above threshold eps = {eps}")
            df <- data.frame(i = 1, j = 1, val = 0)
        } else {
            df <- data.frame(
                i = idx[, 1],
                j = idx[, 2],
                val = Th[idx]
            )
        }

        ggplot2::ggplot(df, ggplot2::aes(i, j, color = val)) +
            ggplot2::geom_point(size = 1) +
            ggplot2::scale_color_gradient2(
                low = "blue", mid = "white", high = "red",
                midpoint = 0,
                name = expression(theta[ij]),
                na.value = "grey80"
            ) +
            ggplot2::scale_y_reverse() +
            ggplot2::coord_equal() +
            ggplot2::labs(
                title = paste("Network snapshot at t =", t),
                subtitle = paste("Relation", rel, "- Sparse view"),
                x = "Actor i", y = "Actor j"
            ) +
            ggplot2::theme_minimal()
    } else {
        # regular heatmap for smaller networks
        # convert to data frame for ggplot
        df <- expand.grid(i = 1:m, j = 1:m)
        df$val <- c(t(Th))

        ggplot2::ggplot(df, ggplot2::aes(i, j, fill = val)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient2(
                low = "navy", mid = "white", high = "firebrick",
                midpoint = 0,
                name = expression(theta[ij]),
                na.value = "grey80"
            ) +
            ggplot2::scale_y_reverse() +
            ggplot2::coord_equal() +
            ggplot2::labs(
                title = paste("Network snapshot at t =", t),
                subtitle = paste("Relation", rel),
                x = "Actor i", y = "Actor j"
            ) +
            ggplot2::theme_minimal()
    }
}

#' Tidy DBN Summary
#'
#' @description Extract posterior means in tidy format
#' @param fit DBN object
#' @param what Components to extract
#' @param time_subset Time points to include (dynamic model)
#' @return List of posterior mean arrays
#' @keywords internal
tidy_dbn <- function(fit, what = c("A", "B", "Theta"), time_subset = NULL) {
    what <- match.arg(what, several.ok = TRUE)
    n_keep <- ifelse(fit$model == "static", dim(fit$B[[1]])[3], length(fit$A))

    if (is.null(time_subset)) {
        time_subset <- if (fit$model == "static") 1 else seq_len(dim(fit$A[[1]])[3])
    }

    out <- list()

    if ("A" %in% what && fit$model == "dynamic") {
        # handle time_thin
        time_idx <- time_subset
        if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
            time_idx <- unique(ceiling(time_subset / fit$settings$time_thin))
            time_idx <- time_idx[time_idx <= dim(fit$A[[1]])[3]]
        }

        Amean <- Reduce(`+`, lapply(fit$A, function(a) a[, , time_idx, drop = FALSE])) / n_keep
        out$A <- Amean
    }

    if ("B" %in% what && fit$model == "dynamic") {
        # handle time_thin
        time_idx <- time_subset
        if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
            time_idx <- unique(ceiling(time_subset / fit$settings$time_thin))
            time_idx <- time_idx[time_idx <= dim(fit$B[[1]])[3]]
        }

        Bmean <- Reduce(`+`, lapply(fit$B, function(b) b[, , time_idx, drop = FALSE])) / n_keep
        out$B <- Bmean
    }

    if ("Theta" %in% what) {
        if (fit$model == "static") {
            out$Theta <- fit$M
        } else {
            # handle time_thin
            time_idx <- time_subset
            if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
                time_idx <- unique(ceiling(time_subset / fit$settings$time_thin))
                time_idx <- time_idx[time_idx <= dim(fit$A[[1]])[3]]
            }

            Th <- array(0, dim = c(fit$dims$m, fit$dims$m, fit$dims$p, length(time_idx)))
            for (s in seq_len(n_keep)) {
                for (i in seq_along(time_idx)) {
                    for (rel in 1:fit$dims$p) {
                        Th[, , rel, i] <- Th[, , rel, i] + fit$A[[s]][, , time_idx[i]] %*% t(fit$B[[s]][, , time_idx[i]])
                    }
                }
            }
            out$Theta <- Th / n_keep
        }
    }

    #
    out
}

#' Plot Group Influence Profile
#'
#' @description Plots posterior group influence over time for dynamic models
#' @param fit A "dbn" object from dbn_dynamic()
#' @param group Integer vector of actor indices
#' @param type "sender" (rows of A_t) or "target" (columns of B_t)
#' @param fun Aggregation across actors: "mean" or "sum"
#' @param measure Per-actor metric: "rowsum" (default), "rowmean", "l2"
#' @param cred Credible band level (0.95 gives 95% bands)
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' # Plot sender influence for actors 1, 3, 5
#' plot_group_influence(fit, group = c(1, 3, 5), type = "sender")
#'
#' # Plot target influence using L2 norm
#' plot_group_influence(fit,
#'     group = c(1, 3, 5), type = "target",
#'     fun = "sum", measure = "l2", cred = 0.8
#' )
#' }
plot_group_influence <- function(fit,
                                 group,
                                 type = c("sender", "target"),
                                 fun = c("mean", "sum"),
                                 measure = c("rowsum", "rowmean", "l2"),
                                 cred = 0.95) {
    if (fit$model != "dynamic") {
        cli::cli_abort("Group influence is defined for dynamic fit objects only.")
    }

    type <- match.arg(type)
    fun <- match.arg(fun)
    measure <- match.arg(measure)

    # validate group indices
    m <- fit$dims$m
    if (any(group < 1) || any(group > m)) {
        cli::cli_abort("Group indices must be between 1 and {m}")
    }

    S <- length(fit$A) # saved draws
    Tt <- dim(fit$A[[1]])[3] # time points
    influence <- matrix(NA, S, Tt) # s x t

    # define row/column functions based on measure
    row_fun <- switch(measure,
        rowsum = rowSums,
        rowmean = rowMeans,
        l2 = function(M) sqrt(rowSums(M^2))
    )

    col_fun <- switch(measure,
        rowsum = colSums,
        rowmean = colMeans,
        l2 = function(M) sqrt(colSums(M^2))
    )

    agg_fun <- if (fun == "mean") base::mean else base::sum

    # compute influence for each MCMC draw
    for (s in seq_len(S)) {
        if (type == "sender") {
            for (t in seq_len(Tt)) {
                a <- fit$A[[s]][, , t]
                influence[s, t] <- agg_fun(row_fun(a[group, , drop = FALSE]))
            }
        } else { # target
            for (t in seq_len(Tt)) {
                b <- fit$B[[s]][, , t]
                influence[s, t] <- agg_fun(col_fun(b[, group, drop = FALSE]))
            }
        }
    }

    # compute posterior summaries
    band <- apply(influence, 2, quantile, probs = c((1 - cred) / 2, 0.5, 1 - (1 - cred) / 2))

    # handle time_thin if present
    time_vals <- 1:Tt
    if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
        # map back to original time scale
        time_vals <- (time_vals - 1) * fit$settings$time_thin + 1
    }

    df <- data.frame(
        time = time_vals,
        lo = band[1, ],
        med = band[2, ],
        hi = band[3, ]
    )

    gtitle <- sprintf(
        "%s group influence (%s %s)",
        if (type == "sender") "Sender" else "Target",
        fun, measure
    )

    ggplot2::ggplot(df, ggplot2::aes(time, med)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80") +
        ggplot2::geom_line(linewidth = 0.8, colour = "steelblue") +
        ggplot2::labs(
            title = gtitle,
            subtitle = paste("Actors:", paste(group, collapse = ", ")),
            y = sprintf("Posterior median %s %d%% CI", "\u00B1", round(cred * 100)),
            x = "Time"
        ) +
        ggplot2::theme_minimal()
}

#' Extract Group Influence Trajectories
#'
#' @description Computes group influence trajectories with posterior quantiles
#' @param fit A "dbn" object from dbn_dynamic()
#' @param group Integer vector of actor indices
#' @param type "sender" or "target"
#' @param measure Per-actor metric: "rowsum", "rowmean", "l2"
#' @param fun Aggregation across actors: "mean" or "sum"
#' @param probs Quantile probabilities to compute
#' @return Data frame with time, posterior quantiles, and mean
#' @export
#' @examples
#' \dontrun{
#' # Get influence trajectory data
#' inf_data <- get_group_influence(fit, group = c(1, 3, 5), type = "sender")
#'
#' # Custom quantiles
#' inf_data <- get_group_influence(fit,
#'     group = c(1, 3, 5),
#'     probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
#' )
#' }
get_group_influence <- function(fit,
                                group,
                                type = c("sender", "target"),
                                measure = c("rowsum", "rowmean", "l2"),
                                fun = c("mean", "sum"),
                                probs = c(0.025, 0.5, 0.975)) {
    if (fit$model != "dynamic") {
        cli::cli_abort("Group influence is defined for dynamic fit objects only.")
    }

    type <- match.arg(type)
    measure <- match.arg(measure)
    fun <- match.arg(fun)

    # Validate inputs
    m <- fit$dims$m
    if (any(group < 1) || any(group > m)) {
        cli::cli_abort("Group indices must be between 1 and {m}")
    }

    S <- length(fit$A)
    Tt <- dim(fit$A[[1]])[3]
    influence <- matrix(NA, S, Tt)

    # define measurement functions
    row_fun <- switch(measure,
        rowsum = rowSums,
        rowmean = rowMeans,
        l2 = function(M) sqrt(rowSums(M^2))
    )

    col_fun <- switch(measure,
        rowsum = colSums,
        rowmean = colMeans,
        l2 = function(M) sqrt(colSums(M^2))
    )

    agg_fun <- if (fun == "mean") base::mean else base::sum

    # compute influence
    for (s in seq_len(S)) {
        if (type == "sender") {
            for (t in seq_len(Tt)) {
                a <- fit$A[[s]][, , t]
                influence[s, t] <- agg_fun(row_fun(a[group, , drop = FALSE]))
            }
        } else {
            for (t in seq_len(Tt)) {
                b <- fit$B[[s]][, , t]
                influence[s, t] <- agg_fun(col_fun(b[, group, drop = FALSE]))
            }
        }
    }

    # compute quantiles and mean
    quants <- apply(influence, 2, quantile, probs = probs)
    means <- colMeans(influence)

    # handle time_thin
    time_vals <- 1:Tt
    if (!is.null(fit$settings$time_thin) && fit$settings$time_thin > 1) {
        time_vals <- (time_vals - 1) * fit$settings$time_thin + 1
    }

    # create output data frame
    df <- data.frame(time = time_vals, mean = means)
    for (i in seq_along(probs)) {
        df[[paste0("q", probs[i])]] <- quants[i, ]
    }

    # metadata as attributes
    attr(df, "group") <- group
    attr(df, "type") <- type
    attr(df, "measure") <- measure
    attr(df, "fun") <- fun

    df
}

#' Compare Group Influences
#'
#' @description Compares influence trajectories of multiple groups
#' @param fit A "dbn" object from dbn_dynamic()
#' @param groups List of integer vectors, each defining a group
#' @param group_names Optional character vector of group names
#' @param type "sender" or "target"
#' @param measure Per-actor metric: "rowsum", "rowmean", "l2"
#' @param fun Aggregation: "mean" or "sum"
#' @param cred Credible band level
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' # Compare two groups
#' compare_group_influence(fit,
#'     groups = list(c(1, 3, 5), c(2, 4, 6)),
#'     group_names = c("Group A", "Group B")
#' )
#' }
compare_group_influence <- function(fit,
                                    groups,
                                    group_names = NULL,
                                    type = c("sender", "target"),
                                    measure = c("rowsum", "rowmean", "l2"),
                                    fun = c("mean", "sum"),
                                    cred = 0.95) {
    if (!is.list(groups)) {
        cli::cli_abort("groups must be a list of integer vectors")
    }

    n_groups <- length(groups)
    if (is.null(group_names)) {
        group_names <- paste("Group", 1:n_groups)
    } else if (length(group_names) != n_groups) {
        cli::cli_abort("group_names must have same length as groups")
    }

    # influence for each group
    all_data <- data.frame()

    for (i in seq_len(n_groups)) {
        inf_data <- get_group_influence(fit, groups[[i]], type, measure, fun,
            probs = c((1 - cred) / 2, 0.5, 1 - (1 - cred) / 2)
        )

        # get cols
        df_group <- data.frame(
            time = inf_data$time,
            median = inf_data[[paste0("q", 0.5)]],
            lo = inf_data[[paste0("q", (1 - cred) / 2)]],
            hi = inf_data[[paste0("q", 1 - (1 - cred) / 2)]],
            group = group_names[i]
        )

        all_data <- rbind(all_data, df_group)
    }

    # viz
    ggplot2::ggplot(all_data, ggplot2::aes(
        x = time, y = median,
        color = group, fill = group
    )) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.3) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::labs(
            title = sprintf(
                "%s group influence comparison (%s %s)",
                if (type == "sender") "Sender" else "Target", fun, measure
            ),
            x = "Time",
            y = sprintf("Posterior median %s %d%% CI", "\u00B1", round(cred * 100)),
            color = "Group",
            fill = "Group"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
}

#' Simulate data from static DBN model
#'
#' @description Generate posterior predictive samples from a fitted static DBN model
#' @param fit A fitted dbn object from dbn_static()
#' @param S Number of posterior samples to generate
#' @param summary Character string specifying summary type:
#'   - "none": Return full array of simulations (default)
#'   - "mean": Return posterior mean across simulations
#' @return If summary = "none", returns 5D array with dimensions nodes by nodes by relations by time by samples.
#'   If summary = "mean", returns 4D array with dimensions nodes by nodes by relations by time containing posterior means.
#' @keywords internal
simulate_static <- function(fit, S, summary = "none") {
    m <- fit$dims$m
    p <- fit$dims$p
    n <- fit$dims$n

    # check if we need to load tprod
    if (!exists("tprod")) {
        cli::cli_abort("tprod function not found. Please load the dbn package.")
    }

    out <- array(0, dim = c(m, m, p, n, S))

    # sample indices from posterior
    n_saved <- dim(fit$B[[1]])[3]
    idx <- sample(n_saved, S, replace = TRUE)

    for (s in seq_len(S)) {
        # extract B matrices for this draw
        Bdraw <- lapply(fit$B, function(b) b[, , idx[s]])

        # generate replicated data
        # Y = tprod(M, B) + noise
        Yrep <- tprod(fit$M, Bdraw)

        # add noise
        s2 <- fit$params[idx[s], "s2"]
        Yrep <- Yrep + array(rnorm(prod(dim(fit$M)), sd = sqrt(s2)),
            dim = dim(fit$M)
        )

        out[, , , , s] <- Yrep
    }

    if (summary == "mean") {
        apply(out, 1:4, mean)
    } else {
        out
    }
}

#' Forecast future network states from dynamic DBN model
#'
#' @description Generate H-step ahead forecasts from a fitted dynamic DBN model
#' @param fit A fitted dbn object from dbn_dynamic()
#' @param H Number of time steps to forecast ahead
#' @param S Number of posterior samples to generate
#' @param summary Character string specifying summary type:
#'   - "none": Return full array of forecasts (default)
#'   - "mean": Return posterior mean forecasts
#' @return If \code{summary = "none"}, returns 5D array with dimensions nodes by nodes by relations by horizon by samples.
#'   If \code{summary = "mean"}, returns 4D array with dimensions nodes by nodes by relations by horizon containing posterior means.
#' @keywords internal
simulate_dynamic <- function(fit, H, S, summary = "none") {
    m <- fit$dims$m
    p <- fit$dims$p
    Tt <- fit$dims$Tt

    # output array
    Theta_pred <- array(0, c(m, m, p, H, S))

    # sample from posterior
    n_saved <- length(fit$A)
    idx <- sample(n_saved, S, replace = TRUE)

    for (s in seq_len(S)) {
        # get last time point matrices
        A_last <- fit$A[[idx[s]]][, , Tt]
        B_last <- fit$B[[idx[s]]][, , Tt]
        sigma2 <- fit$sigma2[idx[s]]

        # initialize: could use last observed Theta or start from equilibrium
        # start from zero (equilibrium mean after centering)
        Theta_curr <- array(0, c(m, m, p))

        # forecast H steps
        for (h in seq_len(H)) {
            Theta_new <- array(0, c(m, m, p))

            for (rel in seq_len(p)) {
                # evolution equation: Theta_t = A_t * Theta_{t-1} * B_t' + noise
                Theta_new[, , rel] <- A_last %*% Theta_curr[, , rel] %*% t(B_last) +
                    sqrt(sigma2) * matrix(rnorm(m * m), m, m)
            }

            Theta_pred[, , , h, s] <- Theta_new
            Theta_curr <- Theta_new
        }
    }

    if (summary == "mean") {
        apply(Theta_pred, 1:4, mean)
    } else {
        Theta_pred
    }
}

#' Posterior Predictive Ordinal Data
#'
#' @description Generate ordinal data from posterior predictive distribution
#' @param fit A "dbn" object
#' @param draws Number of draws
#' @param H Forecast horizon (dynamic models)
#' @return Array of ordinal predictions
#' @keywords internal
predict_ordinal <- function(fit, draws = 100, H = NULL) {
    # get latent predictions
    if (fit$model == "static") {
        Z_pred <- predict(fit, S = draws, summary = "none")

        # add baseline mean
        for (s in seq_len(draws)) {
            Z_pred[, , , , s] <- sweep(Z_pred[, , , , s], 1:3, fit$M, "+")
        }

        # convert to ordinal using empirical CDF of original data
        vals <- sort(unique(c(fit$R)))
        Z_vec <- c(Z_pred)

        # map to ordinal
        R_pred <- vals[findInterval(
            Z_vec,
            quantile(c(fit$R),
                probs = seq(0, 1, length = length(vals) + 1)[-c(1, length(vals) + 1)],
                na.rm = TRUE
            )
        )]

        array(R_pred, dim = dim(Z_pred))
    } else {
        if (is.null(H)) H <- 1

        # get Theta predictions
        Theta_pred <- predict(fit, H = H, S = draws, summary = "none")

        # add baseline and convert to ordinal
        m <- fit$dims$m
        p <- fit$dims$p
        R_pred <- array(NA, dim = c(m, m, p, H, draws))

        vals <- sort(unique(c(fit$R)))

        for (s in seq_len(draws)) {
            # add mean structure
            M_sample <- fit$M[[sample(length(fit$M), 1)]]

            for (h in seq_len(H)) {
                Z_h <- sweep(Theta_pred[, , , h, s], 1:3, M_sample, "+")

                # convert to ordinal
                Z_vec <- c(Z_h)
                R_vec <- vals[findInterval(
                    Z_vec,
                    quantile(c(fit$R),
                        probs = seq(0, 1, length = length(vals) + 1)[-c(1, length(vals) + 1)],
                        na.rm = TRUE
                    )
                )]

                R_pred[, , , h, s] <- array(R_vec, dim = c(m, m, p))
            }
        }

        R_pred
    }
}
