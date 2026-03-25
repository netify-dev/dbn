# Update Theta Using Tucker Products

updates theta array using B matrices via Tucker products

## Usage

``` r
update_theta_tucker(Z, M, B, s2)
```

## Arguments

- Z:

  centered data array (n_row x n_col x p x Tt)

- M:

  baseline mean (n_row x n_col x p)

- B:

  list of K influence matrices

- s2:

  observation variance

## Value

list with updated Y (theta deviations) and X (lagged theta)
