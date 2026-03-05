# Predict from Low-Rank Fit

H-step-ahead forecasts using random-walk extension of alpha

## Usage

``` r
predict_lowrank(object, H = 1, draws = 100, summary = c("mean", "none"))
```

## Arguments

- object:

  A dbn object with model="lowrank"

- H:

  Number of forecast steps

- draws:

  Number of posterior draws

- summary:

  "mean" for posterior mean, "none" for all draws

## Value

Prediction array
