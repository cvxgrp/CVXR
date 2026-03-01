# Unpack Results (backward-compatible alias)

**\[deprecated\]**

## Usage

``` r
unpack_results(problem, solution, chain, inverse_data)
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

Use
[`problem_unpack_results()`](https://www.cvxgrp.org/CVXR/reference/problem_unpack_results.md)
instead. This alias exists for backward compatibility with older CVXR
examples.

## See also

[`problem_unpack_results()`](https://www.cvxgrp.org/CVXR/reference/problem_unpack_results.md)
