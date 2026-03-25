# Parse Blocks Parameter

converts user-specified blocks into standardized format

## Usage

``` r
parse_blocks(blocks, Tt, n)
```

## Arguments

- blocks:

  block specification (integer, vector, or "auto")

- Tt:

  total number of time points

- n:

  number of actors (for minimum block length)

## Value

list with K, boundaries, lengths, block_indices, names
