#' Posterior Extraction and Summary Functions
#'
#' @description Core functions for extracting and summarizing posterior samples
#' @name posterior_extract
#' @keywords internal
NULL

#' Extract Theta slices from posterior draws
#'
#' @description Extract specific slices of Theta arrays from posterior draws
#' @param fit A dbn model fit object
#' @param draws Integer vector of draw indices to extract
#' @param i Row indices (sender nodes)
#' @param j Column indices (receiver nodes)
#' @param rel Relation indices
#' @param time Time indices
#' @return List of Theta slices
#' @keywords internal
theta_slice <- function(fit, draws = NULL, i = NULL, j = NULL, rel = NULL, time = NULL) {

    # validate draws exist
    if (is.null(fit$draws) || is.null(fit$draws$theta)) {
        warning("Model fit does not contain theta draws")
        return(NULL)
    }

    # get dims
    dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims

    # fefault to all indices if not specified
    if (is.null(i)) i <- seq_len(dims$m)
    if (is.null(j)) j <- seq_len(dims$m)
    if (is.null(rel)) rel <- seq_len(dims$p)
    if (is.null(time)) {
        time_dim <- dims$Tt %||% dims$T %||% dims$n %||% 1
        time <- seq_len(time_dim)
    }


    # default to all draws if not specified
    if (is.null(draws)) {
        n_draws <- length(fit$draws$theta)
        draws <- seq_len(n_draws)
    }

    # extract slices
    out <- lapply(fit$draws$theta[draws], function(th) {
        th[i, j, rel, time, drop = FALSE]
    })

    # if single value requested, return as vector
    if (length(i) == 1 && length(j) == 1 && length(rel) == 1 && length(time) == 1) {
        return(unlist(out))
    }

    # if multiple time points but single dyad/rel, return as matrix
    if (length(i) == 1 && length(j) == 1 && length(rel) == 1 && length(time) > 1) {
        return(do.call(rbind, lapply(out, as.vector)))
    }

    out
}

#' Summarize Theta over posterior draws
#'
#' @description Compute summary statistics for Theta parameters
#' @param fit A dbn model fit object
#' @param fun Summary function (default: mean)
#' @param draws Draw indices (default: all)
#' @param i Row indices
#' @param j Column indices
#' @param rel Relation indices
#' @param time Time indices
#' @param chunk Chunk size for memory-efficient processing
#' @return Data frame with summarized values
#' @keywords internal
theta_summary <- function(fit, fun = mean,
                          draws = NULL, i = NULL, j = NULL,
                          rel = NULL, time = NULL, chunk = 20) {
    # dfault to all draws
    if (is.null(draws)) {
        n_draws <- fit$meta$draws %||% length(fit$draws$theta)
        draws <- seq_len(n_draws)
    }

    # auto-disable chunking for non-linear functions
    if (!identical(fun, mean) && chunk < length(draws)) {
        chunk <- length(draws) # process all draws at once
    }

    # get dimensions
    dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims

    # process in chunks
    out <- NULL
    for (block in split(draws, ceiling(seq_along(draws) / chunk))) {
        # extract and combine slices
        slices <- theta_slice(fit, block, i, j, rel, time)

        # check if theta_slice returned null
        if (is.null(slices)) {
            return(NULL)
        }

        # if slices is a vector (single value requested), handle specially
        if (is.vector(slices) && !is.list(slices)) {
            # single value case
            res <- fun(slices)
            if (length(res) == 1) {
                df <- data.frame(
                    i = i,
                    j = j,
                    rel = rel,
                    time = time,
                    value = res,
                    .w_chunk = length(block)
                )
            } else {
                # function returned multiple values
                df <- data.frame(
                    i = rep(i, length(res)),
                    j = rep(j, length(res)),
                    rel = rep(rel, length(res)),
                    time = rep(time, length(res)),
                    value = res,
                    .w_chunk = length(block)
                )
            }
            out <- rbind(out, df)
            next
        }

        # if slices is a matrix (from single dyad/rel with multiple times), handle specially
        if (is.matrix(slices)) {
            # create a data frame directly
            vals <- apply(slices, 2, fun)
            df <- data.frame(
                i = rep(i, length(vals)),
                j = rep(j, length(vals)),
                rel = rep(rel, length(vals)),
                time = time,
                value = vals,
                .w_chunk = length(block)
            )
            out <- rbind(out, df)
            next
        }

        # stack along new dimension
        if (requireNamespace("abind", quietly = TRUE)) {
            arr <- do.call(abind::abind, c(slices, list(along = 5)))
        } else {
            # fallback without abind
            d <- dim(slices[[1]])
            arr <- array(unlist(slices), c(d, length(slices)))
        }

        # apply summary function
        res <- apply(arr, 1:4, fun)

        # handle case where fun returns multiple values
        if (is.matrix(res) || (is.array(res) && length(dim(res)) > 4)) {
            # if fun returns multiple values, melt appropriately
            df_list <- list()
            n_vals <- if (is.matrix(res)) ncol(res) else dim(res)[5]

            for (v in 1:n_vals) {
                if (is.matrix(res)) {
                    res_v <- res[, v]
                    dim(res_v) <- dim(arr)[1:4]
                } else {
                    res_v <- res[, , , , v]
                }
                df_v <- suppressWarnings(as.data.frame.table(res_v, responseName = "value"))
                # set names based on actual dimensions
                dim_names <- c("i", "j", "rel", "time")[1:length(dim(res_v))]
                names(df_v) <- c(dim_names, "value")
                # convert factors to numeric
                for (col in dim_names) {
                    if (is.factor(df_v[[col]])) {
                        df_v[[col]] <- as.numeric(df_v[[col]])
                    }
                }
                # add missing dimensions as constants
                if (!"i" %in% names(df_v)) df_v$i <- i
                if (!"j" %in% names(df_v)) df_v$j <- j
                if (!"rel" %in% names(df_v)) df_v$rel <- rel
                if (!"time" %in% names(df_v)) df_v$time <- time
                df_list[[v]] <- df_v
            }
            df <- do.call(rbind, df_list)
        } else {
            # convert to data frame
            df <- suppressWarnings(as.data.frame.table(res, responseName = "value"))
            names(df) <- c("i", "j", "rel", "time", "value")
            # convert factors to numeric
            for (col in c("i", "j", "rel", "time")) {
                if (col %in% names(df) && is.factor(df[[col]])) {
                    df[[col]] <- as.numeric(df[[col]])
                }
            }
        }

        df$.w_chunk <- length(block) # block is the current split(draws, …)

        out <- rbind(out, df)
    }

    # aggregate if processing multiple chunks
    if (!is.null(out) && nrow(out) > 0) {
        # only aggregate if we have duplicate rows AND we're processing in chunks
        # first check if aggregation is needed
        key_cols <- intersect(c("i", "j", "rel", "time"), names(out))

        if (length(key_cols) > 0 && "value" %in% names(out) && chunk < length(draws)) {
            # check for duplicates
            if (length(key_cols) == 1) {
                keys <- out[[key_cols[1]]]
            } else {
                keys <- do.call(paste, c(out[key_cols], sep = "_"))
            }

            if (anyDuplicated(keys)) {
                if (identical(fun, mean)) {
                    # weighted mean with exact draw counts
                    out$num <- out$value * out$.w_chunk
                    out$den <- out$.w_chunk
                    agg <- aggregate(cbind(num, den) ~ .,
                        data = out[c(key_cols, "num", "den")],
                        FUN = sum, na.action = na.pass
                    )
                    agg$value <- agg$num / agg$den
                    out <- agg[c(key_cols, "value")]
                } else {
                    cli::cli_abort("Chunked processing with non-linear `fun` not supported; set `chunk = length(draws)`.")
                }
            }
        }
    }

    # clean up chunk weights if present
    if (!is.null(out) && ".w_chunk" %in% names(out)) {
        out$.w_chunk <- NULL
    }

    out
}

#' Summarize scalar parameters
#'
#' @description Compute quantiles for scalar parameter traces
#' @param fit A dbn model fit object
#' @param probs Probability levels for quantiles
#' @return Data frame with parameter summaries
#' @keywords internal
param_summary <- function(fit, probs = c(0.05, 0.5, 0.95)) {
    # check for draws format
    if (!is.null(fit$draws$pars)) {
        pars <- fit$draws$pars
    } else {
        # fallback to legacy format
        pars <- NULL

        # try to collect scalar parameters
        scalar_pars <- c(
            "sigma2", "sigma2_proc", "sigma2_obs",
            "tau_A2", "tau_B2", "g2", "rho_A", "rho_B"
        )

        for (par in scalar_pars) {
            if (!is.null(fit[[par]])) {
                if (is.null(pars)) {
                    pars <- data.frame(fit[[par]])
                    names(pars) <- par
                } else {
                    pars[[par]] <- fit[[par]]
                }
            }
        }
    }

    if (is.null(pars)) {
        warning("No scalar parameters found in model fit")
        return(NULL)
    }

    # compute quantiles
    quants <- as.data.frame(t(apply(pars, 2, quantile, probs = probs, na.rm = TRUE)))

    # fix column names for quantiles immediately
    # use round to preserve decimal precision (e.g., q2.5, q50, q97.5)
    quant_names <- paste0("q", round(100 * probs, 1))
    names(quants) <- quant_names

    # add parameter name
    quants$parameter <- rownames(quants)
    rownames(quants) <- NULL

    # add mean and sd
    quants$mean <- apply(pars, 2, mean, na.rm = TRUE)
    quants$sd <- apply(pars, 2, sd, na.rm = TRUE)

    # reorder columns
    quants <- quants[, c("parameter", "mean", "sd", quant_names)]

    quants
}

#' Summarize latent means (M arrays)
#'
#' @description Compute summaries for latent mean arrays M
#' @param fit A dbn model fit object
#' @param fun Summary function
#' @param draws Draw indices
#' @param rel Relation indices (optional)
#' @param chunk Chunk size for processing
#' @return Data frame with M summaries
#' @keywords internal
latent_summary <- function(fit, fun = mean, draws = NULL, rel = NULL, chunk = 20) {
    # check if m exists in draws
    if (is.null(fit$draws$misc$M) && is.null(fit$M)) {
        warning("No M arrays found in model fit")
        return(NULL)
    }

    # default to all draws
    if (is.null(draws)) {
        n_draws <- fit$meta$draws %||% length(fit$draws$misc$M %||% 1)
        draws <- seq_len(n_draws)
    }

    # auto-disable chunking for non-linear functions
    if (!identical(fun, mean) && chunk < length(draws)) {
        chunk <- length(draws) # process all draws at once
    }

    # get m arrays
    if (!is.null(fit$draws$misc$M)) {
        M_list <- fit$draws$misc$M
    } else {
        # legacy format - replicate single m
        M_list <- replicate(length(draws), fit$M, simplify = FALSE)
    }

    # get dimensions
    dims <- if (!is.null(fit$meta$dims)) fit$meta$dims else fit$dims

    out <- NULL
    processed_any <- FALSE

    # debug: check split operation
    blocks <- split(draws, ceiling(seq_along(draws) / chunk))

    for (idx in seq_along(blocks)) {
        block <- blocks[[idx]]
        # extract m draws
        m_block <- M_list[block]

        # check if m_block is empty
        if (length(m_block) == 0) {
            next
        }

        # stack along new dimension
        if (requireNamespace("abind", quietly = TRUE)) {
            arr <- do.call(abind::abind, c(m_block, list(along = 4)))
        } else {
            d <- dim(m_block[[1]])
            arr <- array(unlist(m_block), c(d, length(m_block)))
        }

        # filter by relation if specified
        if (!is.null(rel)) {
            arr <- arr[, , rel, , drop = FALSE]
            # apply summary function to i,j dimensions only
            res <- apply(arr, 1:2, fun)
        } else {
            # apply summary function to i,j,rel dimensions
            res <- apply(arr, 1:3, fun)
        }

        # check if res is empty or has no data
        if (is.null(res) || length(res) == 0) {
            next
        }

        # convert to data frame
        df <- NULL
        tryCatch(
            {
                df <- as.data.frame.table(res, responseName = "value")
            },
            error = function(e) {
                warning("Failed to convert array to data frame: ", e$message)
                df <- NULL
            }
        )

        # check that we got the expected number of columns
        if (is.null(df) || ncol(df) == 0) {
            next
        }

        # set names based on dimensions
        n_vars <- ncol(df) - 1 # subtract 1 for the value column

        if (!is.null(rel)) {
            # when rel is specified, we filtered to a specific relation
            if (n_vars == 2) {
                names(df) <- c("i", "j", "value")
                df$rel <- rel # add rel as constant column
            } else {
                # unexpected structure, but try to handle it
                var_names <- paste0("Var", 1:n_vars)
                names(df) <- c(var_names, "value")
            }
        } else {
            # when rel is null, we have all relations
            if (n_vars == 3) {
                names(df) <- c("i", "j", "rel", "value")
            } else if (n_vars == 2) {
                # maybe p=1, so no rel dimension
                names(df) <- c("i", "j", "value")
                df$rel <- 1
            } else {
                # unexpected structure, but try to handle it
                var_names <- paste0("Var", 1:n_vars)
                names(df) <- c(var_names, "value")
            }
        }

        # convert factors to numeric
        for (col in c("i", "j", "rel")) {
            if (col %in% names(df) && is.factor(df[[col]])) {
                df[[col]] <- as.numeric(df[[col]])
            }
        }

        df$.w_chunk <- length(block) # block is the current split(draws, …)

        out <- rbind(out, df)
    }

    # aggregate if needed
    if (!is.null(out) && nrow(out) > 0) {
        # only aggregate if we have duplicate rows
        key_cols <- intersect(c("i", "j", "rel"), names(out))

        if (length(key_cols) > 0 && "value" %in% names(out) && chunk < length(draws)) {
            # check for duplicates
            if (length(key_cols) == 1) {
                keys <- out[[key_cols[1]]]
            } else {
                keys <- do.call(paste, c(out[key_cols], sep = "_"))
            }

            if (anyDuplicated(keys)) {
                if (identical(fun, mean)) {
                    # weighted mean with exact draw counts
                    out$num <- out$value * out$.w_chunk
                    out$den <- out$.w_chunk
                    agg <- aggregate(cbind(num, den) ~ .,
                        data = out[c(key_cols, "num", "den")],
                        FUN = sum, na.action = na.pass
                    )
                    agg$value <- agg$num / agg$den
                    out <- agg[c(key_cols, "value")]
                } else {
                    cli::cli_abort("Chunked processing with non-linear `fun` not supported; set `chunk = length(draws)`.")
                }
            }
        }
    }

    # clean up chunk weights if present
    if (!is.null(out) && ".w_chunk" %in% names(out)) {
        out$.w_chunk <- NULL
    }

    out
}

#' Extract regime probabilities for HMM models
#'
#' @description Compute posterior probabilities of regime assignments
#' @param fit A dbn_hmm model fit object
#' @return Matrix of regime probabilities (T x R) or NULL
#' @export
regime_probs <- function(fit) {
    # check if this is an hmm model
    if (fit$model != "hmm") {
        return(NULL)
    }

    # try new format first
    if (!is.null(fit$draws$misc$S)) {
        S_list <- fit$draws$misc$S
    } else if (!is.null(fit$S)) {
        # legacy format
        S_list <- fit$S
    } else {
        return(NULL)
    }

    if (length(S_list) == 0) {
        return(NULL)
    }

    # convert to matrix
    S_mat <- do.call(cbind, S_list)

    # get dimensions
    Tt <- nrow(S_mat)
    R <- fit$meta$R %||% fit$settings$R

    if (is.null(R)) {
        # infer r from unique states
        R <- max(unlist(S_list))
    }

    # compute probability of each regime at each time
    probs <- matrix(0, nrow = Tt, ncol = R)
    colnames(probs) <- paste0("Regime", 1:R)
    rownames(probs) <- paste0("Time", 1:Tt)

    for (r in 1:R) {
        probs[, r] <- rowMeans(S_mat == r)
    }

    probs
}

#' Derive new quantities from posterior draws
#'
#' @description Apply a transformation function to each posterior draw
#' @param fit A dbn model fit object
#' @param fun Function to apply to each Theta draw
#' @param draws Draw indices to process
#' @param chunk Chunk size for memory efficiency
#' @param name Name for the derived quantity
#' @return List of derived quantities with class "dbn_derived"
#' @keywords internal
derive_draws <- function(fit, fun, draws = NULL, chunk = 20, name = "derived") {
    # validate draws exist
    if (is.null(fit$draws)) {
        warning("Model fit does not contain posterior draws in expected format")
        return(NULL)
    }

    # default to all draws
    if (is.null(draws)) {
        n_draws <- fit$meta$draws %||%
            length(fit$draws$theta) %||%
            length(fit$draws$misc$A) %||%
            length(fit$draws$misc$B)
        if (is.null(n_draws) || n_draws == 0) {
            warning("Could not determine number of draws")
            return(NULL)
        }
        draws <- seq_len(n_draws)
    }

    # apply function to each draw
    out <- vector("list", length(draws))

    # process in chunks for memory efficiency
    for (block in split(draws, ceiling(seq_along(draws) / chunk))) {
        for (k_idx in seq_along(block)) {
            draw_idx <- block[k_idx]
            k <- which(draws == draw_idx)

            # create a draw object with all available components
            draw <- list()
            if (!is.null(fit$draws$theta)) draw$theta <- fit$draws$theta[[draw_idx]]
            if (!is.null(fit$draws$z)) draw$z <- fit$draws$z[[draw_idx]]
            if (!is.null(fit$draws$pars)) draw$pars <- fit$draws$pars[draw_idx, ]
            if (!is.null(fit$draws$misc)) {
                for (nm in names(fit$draws$misc)) {
                    if (length(fit$draws$misc[[nm]]) >= draw_idx) {
                        draw[[nm]] <- fit$draws$misc[[nm]][[draw_idx]]
                    }
                }
            }

            # apply transformation
            out[[k]] <- fun(draw)
        }
    }

    # try to simplify output if all elements have same structure
    if (length(out) > 0 && is.numeric(out[[1]])) {
        # if all outputs are numeric vectors of same length, combine as matrix
        lens <- sapply(out, length)
        if (all(lens == lens[1])) {
            out <- do.call(rbind, out)
        }
    }

    structure(out, class = "dbn_derived", name = name)
}

#' Print method for derived quantities
#' @param x An object of class "dbn_derived"
#' @param ... Additional arguments (currently unused)
#' @export
print.dbn_derived <- function(x, ...) {
    name <- attr(x, "name")
    cat("Derived posterior quantity:", name, "\n")
    cat("Number of draws:", length(x), "\n")

    if (length(x) > 0) {
        cat("First draw dimensions:", paste(dim(x[[1]]), collapse = " x "), "\n")
    }

    invisible(x)
}

#' Null-coalescing operator
#' @keywords internal
#' 