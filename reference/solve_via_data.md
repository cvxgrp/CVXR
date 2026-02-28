# Solve via Raw Data

Calls the solver on pre-compiled problem data (step 2 of the decomposed
solve pipeline). Dispatches on `x`: when `x` is a `SolvingChain`,
delegates to the terminal solver with proper cache management.

## Usage

``` r
solve_via_data(
  x,
  data,
  warm_start = FALSE,
  verbose = FALSE,
  solver_opts = list(),
  ...
)
```

## Arguments

- x:

  A `SolvingChain` (preferred) or `Solver` object.

- data:

  Named list of solver data from
  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md).

- warm_start:

  Logical; use warm-start if supported.

- verbose:

  Logical; print solver output.

- solver_opts:

  Named list of solver-specific options.

- ...:

  Additional arguments (e.g., `problem` for `SolvingChain` dispatch,
  `solver_cache` for `Solver` dispatch).

## Value

Solver-specific result (a named list).

## See also

[`problem_data`](https://www.cvxgrp.org/CVXR/reference/problem_data.md),
[`problem_unpack_results`](https://www.cvxgrp.org/CVXR/reference/problem_unpack_results.md)
