#ifndef DBN_THREAD_UTILS_H
#define DBN_THREAD_UTILS_H

#include <Rcpp.h>
#include <omp.h>

// read dbn.n_threads from R options and apply it
inline void set_dbn_threads() {
    // look up the option
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
        // fallback to 1 thread
        n_threads = 1;
    }
    
    omp_set_num_threads(n_threads);
}

#endif // DBN_THREAD_UTILS_H