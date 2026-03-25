# Update B Matrices Using Tucker Products

updates all K influence matrices using Tucker product formulation

## Usage

``` r
update_B_tucker(Y, X, B, s2, t2, K = length(B))
```

## Arguments

- Y:

  response array (n_row x n_col x p x Tt)

- X:

  predictor array (n_row x n_col x p x Tt)

- B:

  list of K influence matrices

- s2:

  observation variance

- t2:

  prior precision for B

- K:

  number of modes (default 3)

## Value

list with updated B matrices and sse
