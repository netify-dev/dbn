# Example Dynamic Bilinear Network Data

A small dataset containing simulated ordinal relational data with two
relation types over 10 time points. This dataset is suitable for quick
examples and testing.

## Format

A 4-dimensional array with dimensions 10 x 10 x 2 x 10:

- dim 1:

  Sender actors (10)

- dim 2:

  Receiver actors (10)

- dim 3:

  Relation types (2): cooperation, conflict

- dim 4:

  Time points (10)

## Source

Generated using simulate_static_dbn(n=10, p=2, time=10, K=3, seed=6886)
