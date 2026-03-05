# Register DBN Model

Register a new DBN model variant

## Usage

``` r
register_dbn_model(name, fun_list)
```

## Arguments

- name:

  Model name

- fun_list:

  List of model functions

## Examples

``` r
if (FALSE) { # \dontrun{
custom_model <- list(
    init = function(data) list(theta = rnorm(10)),
    update_Z = function(state, data) state,
    update_Theta = function(state, data) {
        state$theta <- state$theta + rnorm(10, 0, 0.1)
        state
    },
    collect = function(state) list(theta = state$theta)
)
register_dbn_model("custom", custom_model)
} # }
```
