# Static Model MCMC Using Tucker Products

alternative static model implementation using full Tucker products

## Usage

``` r
static_gibbs_step_tucker(Y_obs, Z, M, B, s2, t2, g2)
```

## Arguments

- Y_obs:

  observed data array (n_row x n_col x p x Tt)

- Z:

  latent continuous array (n_row x n_col x p x Tt)

- M:

  baseline mean (n_row x n_col x p)

- B:

  list of K influence matrices

- s2:

  observation variance

- t2:

  prior precision for B

- g2:

  prior variance for M

## Value

list with updated B, M, s2, t2, g2, Y_theta, X_lag
