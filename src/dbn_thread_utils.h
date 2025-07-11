#ifndef DBN_THREAD_UTILS_H
#define DBN_THREAD_UTILS_H

#include <Rcpp.h>
#include <omp.h>

// function to set the number of threads from r options
inline void set_dbn_threads() {
    // get thread count from r options
    Rcpp::Environment base("package:base");
    Rcpp::Function getOption = base["getOption"];
    
    int n_threads = 1;  // default
    
    try {
        SEXP opt = getOption("dbn.n_threads");
        if (!Rf_isNull(opt)) {
            n_threads = Rcpp::as<int>(opt);
            if (n_threads < 1) n_threads = 1;
            if (n_threads > omp_get_max_threads()) {
                n_threads = omp_get_max_threads();
            }
        }
    } catch(...) {
        // fallback to 1 thread if any error
        n_threads = 1;
    }
    
    omp_set_num_threads(n_threads);
}

#endif // DBN_THREAD_UTILS_H