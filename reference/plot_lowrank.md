# Plot Low-Rank DBN Results

Diagnostic plots for low-rank DBN model fits

## Usage

``` r
plot_lowrank(x, factors_show = min(3, x$settings$r), time_points = NULL)
```

## Arguments

- x:

  A dbn object with model="lowrank"

- factors_show:

  Number of factors to display

- time_points:

  Optional subset of time points

## Value

List of ggplot objects or combined plot
