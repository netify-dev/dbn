# Tidy Extractor for Low-Rank Factor Trajectories

Extracts the time-varying factor strengths \\\alpha_t\\ from a fitted
low-rank model in tidy data-frame format, ready for plotting with
ggplot2. Each factor captures how strongly one latent dimension of the
influence structure drives the dynamics at each time point. The returned
data frame includes posterior means and 95\\ credible intervals.

## Usage

``` r
tidy_dbn_lowrank(fit, factors = NULL)
```

## Arguments

- fit:

  A fitted `dbn` object with `model = "lowrank"`.

- factors:

  Integer vector specifying which factors to extract (default: all `r`
  factors). For example, `factors = 1:2` extracts the first two factors.

## Value

Data frame with columns: `time`, `mean` (posterior mean of
\\\alpha_k(t)\\), `lo` (2.5th percentile), `hi` (97.5th percentile), and
`factor` (factor index).

## See also

[`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md) with
`model = "lowrank"`,
[`plot_dbn()`](https://netify-dev.github.io/dbn/reference/plot_dbn.md)
for built-in diagnostic plots
