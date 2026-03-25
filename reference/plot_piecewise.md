# Plot Piecewise DBN Results

diagnostic plots for piecewise-static model

## Usage

``` r
plot_piecewise(
  x,
  type = c("trace", "blocks", "stability", "influence"),
  block = NULL,
  ...
)
```

## Arguments

- x:

  piecewise dbn fit object

- type:

  plot type: "trace", "blocks", "stability", "influence"

- block:

  which block to plot (for block-specific plots)

- ...:

  additional arguments

## Value

plot object (invisibly)
