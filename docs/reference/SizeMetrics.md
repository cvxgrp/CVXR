# Problem Size Metrics

Reports scalar-counts and data-dimension metrics for a
[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md).
Constructed by
[`size_metrics`](https://www.cvxgrp.org/CVXR/reference/size_metrics.md);
end users normally call `size_metrics(prob)` rather than this
constructor directly.

## Usage

``` r
SizeMetrics(
  num_scalar_variables = 0L,
  num_scalar_data = 0L,
  num_scalar_eq_constr = 0L,
  num_scalar_leq_constr = 0L,
  max_data_dimension = 0L,
  max_big_small_squared = 0
)
```

## Arguments

- num_scalar_variables:

  Total scalar entries across all variables in the problem.

- num_scalar_data:

  Total scalar entries across all constants and parameters.

- num_scalar_eq_constr:

  Total scalar entries in equality (`Equality`, `Zero`) constraints.

- num_scalar_leq_constr:

  Total scalar entries in inequality (`Inequality`, `NonNeg`, `NonPos`)
  constraints.

- max_data_dimension:

  Largest single dimension across any data block (constant or
  parameter).

- max_big_small_squared:

  Maximum of `big * small^2` over all data blocks, where `big`/`small`
  are the larger/ smaller dimension of the block.

## Value

A `SizeMetrics` object.
