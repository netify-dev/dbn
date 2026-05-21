####
# package initialization
####

#' @keywords internal
.onLoad <- function(libname, pkgname) {
	requireNamespace("Rcpp", quietly = TRUE)
	requireNamespace("RcppArmadillo", quietly = TRUE)

	# load compiled code
	tryCatch({
		library.dynam("dbn", pkgname, libname)
	}, error = function(e) {
		# Try alternative paths
		so_paths <- c(
			file.path(libname, pkgname, "libs", "dbn.so"),
			file.path(libname, pkgname, "libs", paste0("dbn", .Platform$dynlib.ext)),
			"src/dbn.so",
			"libs/dbn.so"
		)
		for (so_path in so_paths) {
			if (file.exists(so_path)) {
				tryCatch({
					dyn.load(so_path)
				}, error = function(e2) {
					# silently continue
				})
				break
			}
		}
	})

	# default options
	op <- options()
	op.dbn <- list(
		dbn.use_cpp_stability = TRUE,
		dbn.use_ffbs_dlm_cpp = TRUE,
		dbn.use_ffbs_cpp = TRUE,
		dbn.use_cpp_ranklik = TRUE,
		dbn.use_batch_ffbs = FALSE,
		dbn.n_threads = 1L
	)
	toset <- !(names(op.dbn) %in% names(op))
	if (any(toset)) options(op.dbn[toset])

	invisible()
}
