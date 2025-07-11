#ifndef THREAD_CONTROL_H
#define THREAD_CONTROL_H

#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

inline void set_dbn_threads() {
#ifdef _OPENMP
    // get the thread count from R options
    Rcpp::Environment base = Rcpp::Environment::base_env();
    Rcpp::Function getOption = base["getOption"];
    
    // try to get dbn.n_threads option
    SEXP result = getOption("dbn.n_threads");
    
    if (!Rf_isNull(result)) {
        int n_threads = Rcpp::as<int>(result);
        if (n_threads > 0) {
            omp_set_num_threads(n_threads);
        }
    }
#endif
}

inline int get_dbn_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

#endif // THREAD_CONTROL_H