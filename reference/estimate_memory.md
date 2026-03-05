# Estimate Memory Usage for Dynamic DBN

Estimates RAM (in GB) needed for dynamic model output storage

## Usage

``` r
estimate_memory(
  n_row,
  n_col = n_row,
  p = 1,
  Tt = 50,
  nscan = 5000,
  burn = 1000,
  odens = 1,
  time_thin = 1,
  family = "ordinal",
  quiet = FALSE
)
```

## Arguments

- n_row:

  Number of sender nodes

- n_col:

  Number of receiver nodes (default: n_row)

- p:

  Number of relations

- Tt:

  Number of time points

- nscan:

  Number of MCMC iterations after burn-in

- burn:

  Number of burn-in iterations

- odens:

  Thinning interval

- time_thin:

  Time-thinning interval

- family:

  Data family ("ordinal", "gaussian", "binary")

- quiet:

  Suppress printed output

## Value

Estimated memory in GB (invisibly)

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`dbn_dynamic`](https://netify-dev.github.io/dbn/reference/dbn_dynamic.md)

## Examples

``` r
# Estimate memory for a moderate network
estimate_memory(n_row = 20, Tt = 30, nscan = 5000)
#> Dynamic DBN memory estimate:
#>   Network:   20 x 20, 1 relation(s), 30 time points
#>   MCMC:      5000 draws (nscan=5000, odens=1)
#>   Time thin: 1 (keeping 30 of 30 time points)
#>   --------------------------------
#>   Theta:    0.45 GB
#>   Z:        0.45 GB
#>   A:        0.45 GB
#>   B:        0.45 GB
#>   M:        0.01 GB
#>   --------------------------------
#>   TOTAL:    1.80 GB

# Bipartite network with time thinning
estimate_memory(n_row = 15, n_col = 25, Tt = 50,
                nscan = 10000, time_thin = 5)
#> Dynamic DBN memory estimate:
#>   Network:   15 x 25, 1 relation(s), 50 time points
#>   MCMC:      10000 draws (nscan=10000, odens=1)
#>   Time thin: 5 (keeping 10 of 50 time points)
#>   --------------------------------
#>   Theta:    0.28 GB
#>   Z:        0.28 GB
#>   A:        0.17 GB
#>   B:        0.47 GB
#>   M:        0.03 GB
#>   --------------------------------
#>   TOTAL:    1.22 GB

# Suppress printed output, just get the value
gb <- estimate_memory(n_row = 10, Tt = 20, nscan = 2000, quiet = TRUE)
gb
#> [1] 0.120759
```
