#' Low-rank DBN Posterior Analysis Methods
#'
#' @description Posterior analysis utilities for low-rank DBN models
#' @name lowrank-methods
#' @keywords internal
NULL

#' Plot Low-rank DBN results
#'
#' @description Creates diagnostic plots for low-rank DBN model fits
#' @param x A dbn object with model="lowrank"
#' @param factors_show Number of factors to display (default: min(3, rank))
#' @param time_points Optional subset of time points to plot
#' @return Either a list of ggplot objects or a combined plot (if gridExtra available)
#' @keywords internal
#' @import ggplot2
plot_lowrank <- function(x,
                         factors_show = min(3, x$settings$r),
                         time_points = NULL) {
    if (x$model != "lowrank") stop("Input is not a low-rank fit")

    ##  parameter traces 
    df_trace <- data.frame(
        iter = seq_along(x$sigma2),
        sigma2 = x$sigma2,
        tau_alpha2 = x$tau_alpha2,
        tau_B2 = x$tau_B2,
        g2 = x$g2
    )

    # convert to long format for ggplot
    trace_long <- data.frame()
    for (par in c("sigma2", "tau_alpha2", "tau_B2", "g2")) {
        trace_long <- rbind(trace_long, data.frame(
            iter = df_trace$iter,
            par = par,
            val = df_trace[[par]]
        ))
    }

    p_trace <- ggplot(trace_long, aes(iter, val)) +
        geom_line(colour = "steelblue") +
        facet_wrap(~par, scales = "free_y", ncol = 1) +
        labs(title = "MCMC traces", x = "Iteration", y = NULL) +
        theme_minimal()

    ##  factor paths α_{k,t} 
    # dimensions: s (draw) x r x t*
    r <- x$settings$r
    S <- length(x$alpha)
    Tt_ <- ncol(x$alpha[[1]])

    # map thin back to full time for axis labels
    time_vals <- if (isTRUE(x$settings$time_thin > 1)) {
        (seq_len(Tt_) - 1) * x$settings$time_thin + 1
    } else {
        seq_len(Tt_)
    }

    # collect draws (k,t,s)
    arr <- array(NA_real_, c(r, Tt_, S))
    for (s in seq_len(S)) arr[, , s] <- x$alpha[[s]]

    # posterior summaries
    qfun <- function(a) apply(a, 1:2, quantile, probs = c(.025, .5, .975))
    qs <- qfun(arr)

    alpha_df <- data.frame()
    for (k in seq_len(factors_show)) {
        alpha_df <- rbind(
            alpha_df,
            data.frame(
                time = time_vals,
                lo = qs[1, k, ],
                med = qs[2, k, ],
                hi = qs[3, k, ],
                k = factor(k)
            )
        )
    }

    p_alpha <- ggplot(alpha_df, aes(time, med)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
        geom_line(colour = "firebrick") +
        facet_wrap(~k,
            scales = "free_y", ncol = 1,
            labeller = ggplot2::label_bquote(alpha[.(k)])
        ) +
        labs(title = "Latent factor trajectories", x = "Time", y = expression(alpha)) +
        theme_minimal()

    ##  node loadings u (posterior mean) 
    U_bar <- Reduce(`+`, x$U) / S

    # create melted data frame manually
    df_U <- data.frame()
    for (i in 1:nrow(U_bar)) {
        for (j in 1:ncol(U_bar)) {
            df_U <- rbind(df_U, data.frame(
                actor = i,
                factor = j,
                loading = U_bar[i, j]
            ))
        }
    }

    p_U <- ggplot(df_U, aes(factor, actor, fill = loading)) +
        geom_tile() +
        scale_fill_gradient2(
            low = "navy", mid = "white", high = "darkred",
            midpoint = 0
        ) +
        coord_equal() +
        labs(title = "Posterior mean of U", x = "Factor k", y = "Actor i") +
        theme_minimal()

    ##  combine 
    if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(p_trace, p_alpha, p_U, ncol = 1)
    } else {
        list(trace = p_trace, alpha = p_alpha, U = p_U)
    }
}

#' Summary for low-rank fits
#'
#' @description Prints posterior summaries for low-rank DBN models
#' @param object A dbn object with model="lowrank"
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @keywords internal
summary_lowrank <- function(object, digits = 3, ...) {
    if (object$model != "lowrank") stop("Not a low-rank fit")

    cat("Low-rank Dynamic Bilinear Network model\n")
    cat(
        "  nodes     :", object$dims$m,
        "\n  relations :", object$dims$p,
        "\n  time pts  :", object$dims$Tt,
        "\n  rank      :", object$settings$r, "\n\n"
    )

    sumstat <- function(x) {
        formatC(
            c(
                mean(x, na.rm = TRUE),
                quantile(x, c(.025, .975), na.rm = TRUE)
            ),
            digits = digits, format = "f"
        )
    }
    sm <- rbind(
        sigma2 = sumstat(object$sigma2),
        tau_alpha2 = sumstat(object$tau_alpha2),
        tau_B2 = sumstat(object$tau_B2),
        g2 = sumstat(object$g2)
    )
    colnames(sm) <- c("mean", "2.5%", "97.5%")
    print(sm, quote = FALSE)
    invisible(object)
}

#' Forecast / replicate from a low-rank fit
#'
#' @description Generate predictions from a low-rank DBN model
#' @param object A dbn object with model="lowrank"
#' @param H Number of time steps ahead to forecast
#' @param draws Number of posterior draws to use
#' @param summary Character string: "mean" for posterior mean, "none" for all draws
#' @return Array of predictions with dimensions appropriate to summary choice
#' @keywords internal
predict_lowrank <- function(object, H = 1, draws = 100,
                            summary = c("mean", "none")) {
    if (object$model != "lowrank") stop("Not a low-rank fit")
    summary <- match.arg(summary)

    m <- object$dims$m
    p <- object$dims$p
    r <- object$settings$r
    S_saved <- length(object$U)
    pick <- sample(S_saved, draws, TRUE)

    Theta_pred <- array(0, c(m, m, p, H, draws))

    for (d in seq_len(draws)) {
        s <- pick[d]
        sigma2 <- object$sigma2[s]
        U <- object$U[[s]]
        B_now <- object$B[[s]][, , ncol(object$B[[s]])] # last saved B_t
        alpha_path <- object$alpha[[s]]

        # extend alpha with rw step (assume random walk)
        for (h in seq_len(H)) {
            a_t <- alpha_path[, ncol(alpha_path)] +
                rnorm(r, 0, sqrt(object$tau_alpha2[s]))
            A_t <- U %*% diag(a_t) %*% t(U)

            for (rel in seq_len(p)) {
                Theta_pred[, , rel, h, d] <- A_t %*% Theta_pred[, , rel, max(h - 1, 1), d] %*% t(B_now) +
                    matrix(rnorm(m * m, 0, sqrt(sigma2)), m)
            }

            alpha_path <- cbind(alpha_path, a_t)
        }
    }
    if (summary == "mean") apply(Theta_pred, 1:4, mean) else Theta_pred
}

#' Tidy extractor for low-rank factor paths
#'
#' @description Extract factor trajectories in tidy format
#' @param fit A dbn object with model="lowrank"
#' @param factors Which factors to extract (default: all)
#' @return Data frame with columns: time, mean, lo, hi, factor
#' @keywords internal
tidy_dbn_lowrank <- function(fit, factors = NULL) {
    if (is.null(factors)) factors <- 1:fit$settings$r
    if (fit$model != "lowrank") stop("Not a low-rank fit")
    S <- length(fit$alpha)
    Tt_ <- ncol(fit$alpha[[1]])

    time_vals <- if (isTRUE(fit$settings$time_thin > 1)) {
        (seq_len(Tt_) - 1) * fit$settings$time_thin + 1
    } else {
        seq_len(Tt_)
    }

    out <- lapply(factors, function(k) {
        mat <- sapply(fit$alpha, `[`, k, ) # draws x time
        data.frame(
            time = time_vals,
            mean = colMeans(mat),
            lo = apply(mat, 2, quantile, .025),
            hi = apply(mat, 2, quantile, .975),
            factor = k
        )
    })
    do.call(rbind, out)
}
