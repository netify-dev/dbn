# Predict from HMM-DBN Fit

Generates H-step-ahead forecasts from an HMM-DBN model

## Usage

``` r
predict_hmm(object, H = 1, draws = 100, summary = c("mean", "none"))
```

## Arguments

- object:

  A dbn object with model="hmm"

- H:

  Number of forecast steps

- draws:

  Number of posterior draws

- summary:

  "mean" for posterior mean, "none" for all draws

## Value

Prediction array
