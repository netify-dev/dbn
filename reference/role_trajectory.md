# Track Role Evolution

Tracks first left/right singular vector of A_t or B_t

## Usage

``` r
role_trajectory(fit, mat = c("A", "B"), comp = 1)
```

## Arguments

- fit:

  Dynamic dbn object

- mat:

  "A" or "B"

- comp:

  Component index (default: 1)

## Value

Base R plot (invisible NULL)

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md),
[`dyad_path`](https://netify-dev.github.io/dbn/reference/dyad_path.md),
[`net_snapshot`](https://netify-dev.github.io/dbn/reference/net_snapshot.md)

## Examples

``` r
# \donttest{
sim <- simulate_dynamic_dbn(n = 6, time = 10, seed = 1)
fit <- dbn(sim$Y, model = "dynamic", nscan = 200, burn = 100, verbose = FALSE)
role_trajectory(fit, mat = "A", comp = 1)

# }
```
