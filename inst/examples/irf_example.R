# Example: Impulse Response Functions for Dynamic Gaussian Networks

library(dbn)

# Simulate a dynamic network
set.seed(123)
m <- 10  # Number of nodes
T <- 50  # Time points
p <- 1   # Relations

# Generate some dynamic network data
Y <- array(rnorm(m * m * p * T), dim = c(m, m, p, T))

# Fit a dynamic Gaussian model
# (This is a placeholder - you would use actual data and fitting function)
# fit <- dbn_dynamic(Y, model = "dynamic", family = "gaussian", ...)

# For demonstration, create a mock fit object with necessary components
# In practice, this would come from your actual model fitting
fit <- list(
    model = "dynamic",
    dims = list(m = m, p = p, T = T),
    draws = list(
        misc = list(
            A = array(rnorm(100 * m * m * T, 0, 0.1), dim = c(100, m, m, T)),
            B = array(rnorm(100 * m * m * T, 0, 0.1), dim = c(100, m, m, T)),
            M = array(rnorm(100 * m * m, 0, 0.1), dim = c(100, m, m))
        ),
        pars = matrix(rnorm(100 * 3), 100, 3)
    )
)
class(fit) <- c("dbn_dynamic", "dbn")

# Make A and B have appropriate spectral radii for stability
for (s in 1:100) {
    for (t in 1:T) {
        # Create stable transition matrices
        A_temp <- matrix(rnorm(m * m, 0, 0.1), m, m)
        B_temp <- matrix(rnorm(m * m, 0, 0.1), m, m)
        
        # Scale to ensure stability
        A_temp <- A_temp / (1.5 * max(abs(eigen(A_temp)$values)))
        B_temp <- B_temp / (1.5 * max(abs(eigen(B_temp)$values)))
        
        fit$draws$misc$A[s,,,t] <- A_temp
        fit$draws$misc$B[s,,,t] <- B_temp
    }
}

# Example 1: Unit edge shock to network density
# Shock edge (1,2) and see how network density responds
irf_density <- compute_irf(
    fit, 
    shock = "unit_edge",
    shock_pars = list(i = 1, j = 2),
    H = 20,
    t0 = 25,  # Middle of time series
    stat_fun = stat_density,
    n_draws = 50
)

print(irf_density)
plot(irf_density, title = "IRF: Unit Edge Shock → Network Density")

# Example 2: Node shock to out-degrees
# Shock all outgoing edges from node 1
irf_outdeg <- compute_irf(
    fit,
    shock = "node_out", 
    shock_pars = list(i = 1),
    H = 20,
    t0 = 25,
    stat_fun = function(X) stat_out_degree(X)[1],  # Track node 1's out-degree
    n_draws = 50
)

plot(irf_outdeg, title = "IRF: Node 1 Out-Shock → Node 1 Out-Degree")

# Example 3: Global density shock
# Shock all edges proportionally
irf_recip <- compute_irf(
    fit,
    shock = "density",
    H = 20,
    t0 = 25,
    stat_fun = stat_reciprocity,
    n_draws = 50
)

plot(irf_recip, title = "IRF: Density Shock → Network Reciprocity")

# Example 4: Custom shock matrix
# Create asymmetric shock affecting upper triangle
S_custom <- matrix(0, m, m)
S_custom[upper.tri(S_custom)] <- 0.1
irf_custom <- compute_irf(
    fit,
    shock = S_custom,
    H = 20,
    t0 = 25,
    stat_fun = stat_transitivity,
    n_draws = 50
)

plot(irf_custom, title = "IRF: Upper Triangle Shock → Transitivity")

# Compare multiple IRFs
if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    
    # Combine results
    df_compare <- rbind(
        cbind(irf_density[,c("horizon", "mean", "q025", "q975")], type = "Density"),
        cbind(irf_recip[,c("horizon", "mean", "q025", "q975")], type = "Reciprocity")
    )
    
    ggplot(df_compare, aes(x = horizon, y = mean, color = type, fill = type)) +
        geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.3) +
        geom_line(linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(~type, scales = "free_y") +
        labs(title = "Comparison of Network IRFs",
             x = "Horizon", y = "Response") +
        theme_minimal()
}