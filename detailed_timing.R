library(dbn)
set.seed(6886)

n <- 25
Tt <- 88

Y <- array(rnorm(n*n*Tt, sd=0.5), c(n, n, Tt))
for (t in 1:Tt) diag(Y[,,t]) <- 0

cat("===== SYMMETRIC (verbose timing) =====\n")
system.time(
    fit_sym <- dbn(Y, family="gaussian", model="dynamic", symmetric=TRUE,
                   nscan=20, burn=5, odens=2, seed=6886, verbose=200)
) -> timing_sym

