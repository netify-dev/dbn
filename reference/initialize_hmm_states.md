# Initialize HMM States

Generates initial hidden state sequence for HMM model

## Usage

``` r
initialize_hmm_states(Tt, K, init_probs = NULL, trans_probs = NULL)
```

## Arguments

- Tt:

  Number of time points

- K:

  Number of states

- init_probs:

  Initial state probabilities

- trans_probs:

  Transition probability matrix

## Value

Integer vector of state assignments
