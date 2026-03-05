# Set the number of threads used by dbn

Controls the number of OpenMP threads used for parallel MCMC updates in
the dynamic model. Parallelization applies to the row-wise A and B FFBS
updates and variance computations. The setting persists for the current
R session until changed.

For best performance, set to the number of physical (not logical) cores.
Speedup scales with network size: expect 2–4x improvement with 4–8 cores
for networks with 15+ actors. For small networks the overhead of thread
management can exceed the benefit.

Requires OpenMP support at compile time (standard on Linux and Windows;
on macOS install via `brew install libomp`). Without OpenMP the package
runs single-threaded regardless of this setting.

## Usage

``` r
set_dbn_threads(n_threads = NULL)
```

## Arguments

- n_threads:

  Integer number of threads to use for parallel computation. Use NULL to
  reset to default (1 thread).

## Value

The previous value of n_threads (invisibly)

## See also

[`get_dbn_threads`](https://netify-dev.github.io/dbn/reference/get_dbn_threads.md)

## Examples

``` r
set_dbn_threads(4)
set_dbn_threads(parallel::detectCores(logical = FALSE))
set_dbn_threads(NULL)
```
