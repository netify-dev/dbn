# Minimal test of IRF computation to isolate the issue

test_irf_minimal <- function(fit, draw_idx = 1) {
    cat("\n=== Minimal IRF Test ===\n")
    
    # 1. Basic info
    m <- fit$dims$m
    p <- fit$dims$p
    cat("m =", m, ", p =", p, "\n")
    
    # 2. Try to run compute_irf_single in isolation
    cat("\nTesting compute_irf_single directly:\n")
    
    # Create simple shock
    S <- matrix(0, m, m) 
    S[1, 2] <- 1
    
    # Try with minimal parameters
    result <- tryCatch({
        compute_irf_single(fit, draw_idx, S, H = 5, t0 = 10)
    }, error = function(e) {
        cat("Error:", e$message, "\n")
        cat("Traceback:\n")
        print(traceback())
        return(NULL)
    })
    
    if (!is.null(result)) {
        cat("Success! IRF values:", result, "\n")
    }
    
    # 3. Test components individually
    cat("\n\nTesting components:\n")
    
    # Test A/B extraction
    tryCatch({
        A_test <- fit$A[[draw_idx]]
        B_test <- fit$B[[draw_idx]]
        cat("A and B extraction: OK\n")
        
        # Test C++ function
        Delta_test <- impulse_response_dynamic(A_test, B_test, S, 9, 5)
        cat("C++ impulse_response_dynamic: OK\n")
        
        # Test stat function
        test_mat <- matrix(runif(m*m), m, m)
        dens <- stat_density(test_mat)
        cat("stat_density function: OK, value =", dens, "\n")
        
    }, error = function(e) {
        cat("Component error:", e$message, "\n")
    })
    
    invisible(result)
}

# Alternative compute_irf without any potentially conflicting code
compute_irf_safe <- function(fit, shock_type = "unit_edge", i = 1, j = 2, 
                            H = 10, t0 = 10, n_draws = 10) {
    
    # Get dimensions
    m <- fit$dims$m
    p <- fit$dims$p
    
    # Build shock
    S <- matrix(0, m, m)
    if (shock_type == "unit_edge") {
        S[i, j] <- 1
    }
    
    # Storage
    irf_values <- matrix(NA, n_draws, H + 1)
    
    # Compute for each draw
    for (draw in 1:n_draws) {
        tryCatch({
            # Get A and B
            A_array <- fit$A[[draw]]
            B_array <- fit$B[[draw]]
            
            # Get M
            if (length(dim(fit$M)) == 4) {
                M_baseline <- fit$M[,,1,draw]
            } else {
                M_baseline <- matrix(0, m, m)
            }
            
            # Compute IRF
            Delta <- impulse_response_dynamic(A_array, B_array, S, t0 - 1, H)
            
            # Compute statistics
            for (h in 0:H) {
                baseline_val <- mean(M_baseline[row(M_baseline) != col(M_baseline)])
                shocked_net <- M_baseline + Delta[,,h+1]
                shocked_val <- mean(shocked_net[row(shocked_net) != col(shocked_net)])
                irf_values[draw, h+1] <- shocked_val - baseline_val
            }
            
        }, error = function(e) {
            cat("Draw", draw, "failed:", e$message, "\n")
        })
    }
    
    # Return summary
    list(
        mean = colMeans(irf_values, na.rm = TRUE),
        draws = irf_values
    )
}