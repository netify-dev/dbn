####
# impulse response functions for dynamic gaussian networks
####

library(dbn)

####
# simulate data
set.seed(123)
m = 10
T = 50
p = 1

Y = array(rnorm(m * m * p * T), dim = c(m, m, p, T))
####

####
# mock fit object (in practice, use dbn_dynamic())
fit = list(
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
class(fit) = c("dbn_dynamic", "dbn")

# scale A and B for stability
for (s in 1:100) {
	for (t in 1:T) {
		A_temp = matrix(rnorm(m * m, 0, 0.1), m, m)
		B_temp = matrix(rnorm(m * m, 0, 0.1), m, m)
		A_temp = A_temp / (1.5 * max(abs(eigen(A_temp)$values)))
		B_temp = B_temp / (1.5 * max(abs(eigen(B_temp)$values)))
		fit$draws$misc$A[s,,,t] = A_temp
		fit$draws$misc$B[s,,,t] = B_temp
	}
}
####

####
# unit edge shock -> network density
irf_density = compute_irf(
	fit,
	shock = "unit_edge",
	shock_pars = list(i = 1, j = 2),
	H = 20,
	t0 = 25,
	stat_fun = stat_density,
	n_draws = 50
)

print(irf_density)
plot(irf_density, title = "IRF: Unit Edge Shock -> Network Density")
####

####
# node shock -> out-degree
irf_outdeg = compute_irf(
	fit,
	shock = "node_out",
	shock_pars = list(i = 1),
	H = 20,
	t0 = 25,
	stat_fun = function(X) stat_out_degree(X)[1],
	n_draws = 50
)

plot(irf_outdeg, title = "IRF: Node 1 Out-Shock -> Node 1 Out-Degree")
####

####
# density shock -> reciprocity
irf_recip = compute_irf(
	fit,
	shock = "density",
	H = 20,
	t0 = 25,
	stat_fun = stat_reciprocity,
	n_draws = 50
)

plot(irf_recip, title = "IRF: Density Shock -> Network Reciprocity")
####

####
# custom shock -> transitivity
S_custom = matrix(0, m, m)
S_custom[upper.tri(S_custom)] = 0.1
irf_custom = compute_irf(
	fit,
	shock = S_custom,
	H = 20,
	t0 = 25,
	stat_fun = stat_transitivity,
	n_draws = 50
)

plot(irf_custom, title = "IRF: Upper Triangle Shock -> Transitivity")
####

####
# comparison plot
if (requireNamespace("ggplot2", quietly = TRUE)) {
	library(ggplot2)

	df_compare = rbind(
		cbind(irf_density[, c("horizon", "mean", "q025", "q975")], type = "Density"),
		cbind(irf_recip[, c("horizon", "mean", "q025", "q975")], type = "Reciprocity")
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
####
