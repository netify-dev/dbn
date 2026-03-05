# Structured FFBS for Theta

Memory-efficient FFBS that avoids creating n^2 x n^2 matrices. Supports
rectangular (bipartite) Theta via row-covariance approximation.

## Usage

``` r
ffbs_theta_struct_5arg_cpp(Z, mu, A_array, B_array, sigma2)
```

## Arguments

- Z:

  Observations (n_row x n_col x T array as cube)

- mu:

  Baseline mean (n_row x n_col matrix)

- A_array:

  Time-varying A matrices (n_row x n_row x T)

- B_array:

  Time-varying B matrices (n_col x n_col x T)

- sigma2:

  Innovation variance (scalar)

## Value

Sampled Theta array (n_row x n_col x T)
