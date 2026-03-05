# Forecast future network states from dynamic DBN model

Generate H-step ahead forecasts from a fitted dynamic DBN model

## Usage

``` r
simulate_dynamic(fit, H, S, summary = "none")
```

## Arguments

- fit:

  A fitted dbn object from dbn_dynamic()

- H:

  Number of time steps to forecast ahead

- S:

  Number of posterior samples to generate

- summary:

  Character string specifying summary type:

  - "none": Return full array of forecasts (default)

  - "mean": Return posterior mean forecasts

## Value

If `summary = "none"`, returns 5D array with dimensions nodes by nodes
by relations by horizon by samples. If `summary = "mean"`, returns 4D
array with dimensions nodes by nodes by relations by horizon containing
posterior means.
