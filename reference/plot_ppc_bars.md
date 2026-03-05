# Posterior predictive bar plot for discrete data

Bar plot comparison for discrete outcomes

## Usage

``` r
plot_ppc_bars(fit, ppd, rel = 1, time = NULL, Y_obs = NULL)
```

## Arguments

- fit:

  A dbn model fit object

- ppd:

  Posterior predictive samples

- rel:

  Relation index

- time:

  Time indices

- Y_obs:

  Observed data array (required)

## Value

A ggplot object or NULL
