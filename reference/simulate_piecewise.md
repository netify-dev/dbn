# Simulate from Piecewise Model

generates forecasts or simulations from piecewise model

## Usage

``` r
simulate_piecewise(object, H = 10, ndraws = 100, seed = NULL, ...)
```

## Arguments

- object:

  piecewise dbn fit

- H:

  forecast horizon

- ndraws:

  number of posterior draws to use

- seed:

  random seed

- ...:

  additional arguments

## Value

list with simulated values
