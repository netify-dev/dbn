# Development loader - loads all dbn functions for interactive development
# Usage: source("load_dev.R") from within the dbn directory

library(dbn)

# Create a wrapper environment that can access dbn functions
env <- new.env(parent = getNamespace("dbn"))

# Source all R files into this wrapper environment
# This allows them to access C++ functions from dbn while staying separate
r_files <- sort(list.files("R", pattern = "\\.R$", full.names = TRUE))

# Load non-dispatch files first, then dispatch last (which has bootstrap functions)
r_files <- c(setdiff(r_files, "R/dispatch.R"), "R/dispatch.R")

for (f in r_files) {
	tryCatch({
		source(f, verbose = FALSE, local = env)
	}, warning = function(w) {
		# Suppress warnings
	}, error = function(e) {
		# Skip files with sourcing errors
	})
}

# Export bootstrap functions to global environment for easy interactive access
bootstrap_funcs <- c("dbn_explore", "dbn_explore_bootstrap", "sampler_describe",
		  "print.dbn_boot", "summary.dbn_boot", "als_refit_warmstart",
		  "sampler_describe.dbn", "sampler_describe.dbn_boot")

for (fname in bootstrap_funcs) {
	if (exists(fname, envir = env)) {
		assign(fname, get(fname, envir = env), envir = .GlobalEnv)
	}
}

# Fix S3 method dispatch by creating a wrapper that explicitly calls the right method
sampler_describe_wrapper <- function(object, verbose = TRUE) {
	if (inherits(object, "dbn_boot")) {
		return(sampler_describe.dbn_boot(object, verbose = verbose))
	} else if (inherits(object, "dbn")) {
		return(sampler_describe.dbn(object, verbose = verbose))
	} else {
		stop("sampler_describe: object must be of class 'dbn' or 'dbn_boot'")
	}
}

# Replace the sampler_describe generic with our wrapper
if (exists("sampler_describe", envir = .GlobalEnv)) {
	assign("sampler_describe", sampler_describe_wrapper, envir = .GlobalEnv)
}

cat("✓ Development environment loaded\n")
cat("  Bootstrap functions available: dbn_explore, dbn_explore_bootstrap, sampler_describe\n")
cat("  Other functions available in global environment\n")
