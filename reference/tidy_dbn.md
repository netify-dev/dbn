# Tidy DBN Summary

Extract posterior means in tidy format

## Usage

``` r
tidy_dbn(fit, what = c("A", "B", "Theta"), time_subset = NULL)
```

## Arguments

- fit:

  DBN object

- what:

  Components to extract

- time_subset:

  Time points to include (dynamic model)

## Value

List of posterior mean arrays
