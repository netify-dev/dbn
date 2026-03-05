# Regularize Matrix for Cholesky

Adds small ridge to diagonal for numerical stability

## Usage

``` r
regularize_matrix(mat, eps = NULL)
```

## Arguments

- mat:

  Matrix to regularize

- eps:

  Regularization strength

## Value

Regularized matrix
