# Simulate Simple Test Data

Generates simple ordinal network data for quick testing

## Usage

``` r
simulate_test_data(n = 10, n_col = n, p = 2, time = 20, seed = NULL)
```

## Arguments

- n:

  Number of row actors / senders

- n_col:

  Number of column actors / receivers (default: n)

- p:

  Number of relations

- time:

  Number of time points

- seed:

  Random seed

## Value

Array of ordinal network data

## See also

[`dbn`](https://netify-dev.github.io/dbn/reference/dbn.md) for model
fitting,
[`simulate_static_dbn`](https://netify-dev.github.io/dbn/reference/simulate_static_dbn.md)
for full simulation with true parameters

## Examples

``` r
Y <- simulate_test_data(n = 10, time = 5, seed = 42)
str(Y)
#>  int [1:10, 1:10, 1:2, 1:5] NA 5 4 5 3 3 3 4 3 3 ...
```
