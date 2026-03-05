# Simulate data from static DBN model

Generate posterior predictive samples from a fitted static DBN model

## Usage

``` r
simulate_static(fit, S, summary = "none")
```

## Arguments

- fit:

  A fitted dbn object from dbn_static()

- S:

  Number of posterior samples to generate

- summary:

  Character string specifying summary type:

  - "none": Return full array of simulations (default)

  - "mean": Return posterior mean across simulations

## Value

If summary = "none", returns 5D array with dimensions nodes by nodes by
relations by time by samples. If summary = "mean", returns 4D array with
dimensions nodes by nodes by relations by time containing posterior
means.
