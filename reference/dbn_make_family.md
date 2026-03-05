# Family Constructor

Builds a dbn_family structure holding all family-specific callbacks

## Usage

``` r
dbn_make_family(
  name,
  draw_latent,
  ffbs_wrapper,
  loglik,
  linkinv,
  rgen_obs,
  init_pars = list()
)
```

## Arguments

- name:

  Family name string

- draw_latent:

  Latent Z sampler

- ffbs_wrapper:

  FFBS dispatch wrapper

- loglik:

  Log-likelihood function

- linkinv:

  Inverse link function

- rgen_obs:

  Observation generator

- init_pars:

  Initial parameter list

## Value

A dbn_family object
