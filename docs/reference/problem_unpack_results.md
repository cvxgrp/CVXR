# Unpack Solver Results into a Problem

Inverts the reduction chain and unpacks the raw solver solution into the
original problem's variables and constraints. This is step 3 of the
decomposed solve pipeline:

1.  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md)
    – compile the problem

2.  [`solve_via_data`](https://www.cvxgrp.org/CVXR/reference/solve_via_data.md)`(chain, data)`
    – call the solver

3.  `problem_unpack_results()` – invert and unpack

## Usage

``` r
problem_unpack_results(problem, solution, chain, inverse_data)
```

## Arguments

- problem:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

- solution:

  The raw solver result from
  [`solve_via_data()`](https://www.cvxgrp.org/CVXR/reference/solve_via_data.md).

- chain:

  The `SolvingChain` from
  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md).

- inverse_data:

  The inverse data list from
  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md).

## Value

The problem object (invisibly), with solution unpacked.

## Details

After calling this function, variable values are available via
[`value()`](https://www.cvxgrp.org/CVXR/reference/value.md) and
constraint duals via
[`dual_value()`](https://www.cvxgrp.org/CVXR/reference/dual_value.md).

## See also

[`problem_data`](https://www.cvxgrp.org/CVXR/reference/problem_data.md),
[`solve_via_data`](https://www.cvxgrp.org/CVXR/reference/solve_via_data.md)
