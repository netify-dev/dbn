# Compute Theta for Static Model

On-demand theta computation for static DBN model

## Usage

``` r
compute_theta_static(B, Z, M)
```

## Arguments

- B:

  Sender/receiver effect matrix (m x m)

- Z:

  Latent positions array (m x m x p x Tt) or specific slice

- M:

  Baseline mean array (m x m x p) or specific slice

## Value

Theta array with same dimensions as Z
