# Get the current number of threads used by dbn

Returns the number of OpenMP threads currently configured for parallel
MCMC computation. Defaults to 1 (single-threaded).

## Usage

``` r
get_dbn_threads()
```

## Value

Integer number of threads currently set for parallel computation

## See also

[`set_dbn_threads`](https://netify-dev.github.io/dbn/reference/set_dbn_threads.md)

## Examples

``` r
get_dbn_threads()
#> [1] 1
```
