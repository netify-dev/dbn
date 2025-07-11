#' HMM-DBN Posterior Analysis Methods
#'
#' @description Posterior analysis utilities for HMM-DBN models
#' @name hmm-methods
#' @keywords internal
NULL

#' Plot HMM-DBN results
#'
#' @description Creates diagnostic plots for HMM-DBN model fits
#' @param x A dbn object with model="hmm"
#' @return Either a list of ggplot objects or a combined plot (if gridExtra available)
#' @keywords internal
#' @import ggplot2 gridExtra
plot_hmm <- function(x) {
    if (x$model != "hmm") stop("Not an HMM fit")

    ## regime sequence heat-map
    S_mat <- do.call(cbind, lapply(x$S, identity)) # t* x saved
    R <- x$settings$R
    Tt_ <- nrow(S_mat)
    probs <- sapply(1:R, function(r) rowMeans(S_mat == r))

    # create melted data frame manually
    ### could bring in from netify if we add it in as a dependency
    df_reg <- data.frame()
    for (t in 1:Tt_) {
        for (r in 1:R) {
            df_reg <- rbind(df_reg, data.frame(
                time = t,
                regime = r,
                prob = probs[t, r]
            ))
        }
    }
    df_reg$regime <- factor(df_reg$regime)

    p_reg <- ggplot(df_reg, aes(time, regime, fill = prob)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(title = "Posterior regime probabilities", x = "Time", y = "Regime") +
        theme_minimal()

    ## transition matrix pi
    Pi_bar <- Reduce(`+`, x$Pi) / length(x$Pi)

    # create melted data frame manually
    df_Pi <- data.frame()
    for (i in 1:nrow(Pi_bar)) {
        for (j in 1:ncol(Pi_bar)) {
            df_Pi <- rbind(df_Pi, data.frame(
                from = i,
                to = j,
                prob = Pi_bar[i, j]
            ))
        }
    }
    df_Pi$from <- factor(df_Pi$from)
    df_Pi$to <- factor(df_Pi$to)

    p_Pi <- ggplot(df_Pi, aes(from, to, fill = prob)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "firebrick") +
        geom_text(aes(label = sprintf("%.2f", prob)), size = 3) +
        labs(title = "Posterior mean transition matrix Pi", x = "From", y = "To") +
        theme_minimal()

    ## parameter traces
    trace <- data.frame(
        iter = seq_along(x$sigma2),
        sigma2 = x$sigma2,
        tau_A2 = x$tau_A2,
        tau_B2 = x$tau_B2,
        g2 = x$g2
    )

    # convert to long format manually
    tl <- data.frame()
    for (var in c("sigma2", "tau_A2", "tau_B2", "g2")) {
        tl <- rbind(tl, data.frame(
            iter = trace$iter,
            variable = var,
            value = trace[[var]]
        ))
    }

    p_trace <- ggplot(tl, aes(iter, value)) +
        geom_line(colour = "darkgreen") +
        facet_wrap(~variable, scales = "free_y", ncol = 1) +
        theme_minimal() +
        labs(title = "MCMC traces", x = "Iteration", y = NULL)

    if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(p_reg, p_Pi, p_trace, ncol = 1)
    } else {
        list(regime = p_reg, Pi = p_Pi, trace = p_trace)
    }
}

#' Summary for HMM-DBN fits
#'
#' @description Prints posterior summaries for HMM-DBN models
#' @param object A dbn object with model="hmm"
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @keywords internal
summary_hmm <- function(object, digits = 3, ...) {
    if (object$model != "hmm") stop("Not an HMM fit")

    cat("Regime-switching (HMM) DBN model\n")
    cat(
        "  nodes     :", object$dims$m,
        "\n  relations :", object$dims$p,
        "\n  time pts  :", object$dims$Tt,
        "\n  regimes   :", object$settings$R, "\n\n"
    )

    ss <- function(x) {
        formatC(
            c(
                mean(x, na.rm = TRUE),
                quantile(x, c(.025, .975), na.rm = TRUE)
            ),
            digits = digits, format = "f"
        )
    }
    sm <- rbind(
        sigma2 = ss(object$sigma2),
        tau_A2 = ss(object$tau_A2),
        tau_B2 = ss(object$tau_B2),
        g2 = ss(object$g2)
    )
    colnames(sm) <- c("mean", "2.5%", "97.5%")
    print(sm, quote = FALSE)

    cat("\nPosterior mean transition matrix Pi:\n")
    print(round(Reduce(`+`, object$Pi) / length(object$Pi), 3))
    invisible(object)
}

#' Forecast / replicate from an HMM fit
#'
#' @description Generate predictions from an HMM-DBN model
#' @param object A dbn object with model="hmm"
#' @param H Number of time steps ahead to forecast
#' @param draws Number of posterior draws to use
#' @param summary Character string: "mean" for posterior mean, "none" for all draws
#' @return Array of predictions with dimensions appropriate to summary choice
#' @keywords internal
predict_hmm <- function(object, H = 1, draws = 100,
                        summary = c("mean", "none")) {
    if (object$model != "hmm") stop("Not an HMM fit")
    summary <- match.arg(summary)

    m <- object$dims$m
    p <- object$dims$p
    R <- object$settings$R
    S_saved <- length(object$A) # list of arrays per draw
    pick <- sample(S_saved, draws, TRUE)

    Theta_pred <- array(0, c(m, m, p, H, draws))

    for (d in seq_len(draws)) {
        s <- pick[d]
        A_list <- object$A[[s]] # m x m x R
        B_list <- object$B[[s]]
        Pi <- object$Pi[[s]]
        sigma2 <- object$sigma2[s]

        # start from stationary draw = 0 for simplicity
        regime <- sample(R, 1, prob = colMeans(Pi)) # initial
        Theta_now <- array(0, c(m, m, p))

        for (h in seq_len(H)) {
            # draw next regime
            regime <- sample(R, 1, prob = Pi[regime, ])
            A_t <- A_list[, , regime]
            B_t <- B_list[, , regime]

            for (rel in seq_len(p)) {
                Theta_now[, , rel] <- A_t %*% Theta_now[, , rel] %*% t(B_t) +
                    matrix(rnorm(m * m, 0, sqrt(sigma2)), m)
            }

            Theta_pred[, , , h, d] <- Theta_now
        }
    }
    if (summary == "mean") apply(Theta_pred, 1:4, mean) else Theta_pred
}
