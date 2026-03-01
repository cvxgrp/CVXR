# Get Problem Data for a Solver (deprecated)

**\[deprecated\]**

## Usage

``` r
get_problem_data(x, solver = NULL, ...)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

- solver:

  Character string naming solver, or `NULL` for automatic selection.

- ...:

  Additional arguments.

## Value

A list with components `data`, `chain`, and `inverse_data`.

## Details

Use
[`problem_data`](https://www.cvxgrp.org/CVXR/reference/problem_data.md)
instead.

## See also

[`problem_data`](https://www.cvxgrp.org/CVXR/reference/problem_data.md)
