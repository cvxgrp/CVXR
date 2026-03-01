# Get the Solution Status of a Problem

Returns the status string from the most recent solve, such as
`"optimal"`, `"infeasible"`, or `"unbounded"`.

## Usage

``` r
status(x)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

## Value

Character string, or `NULL` if the problem has not been solved.

## See also

[`OPTIMAL`](https://www.cvxgrp.org/CVXR/reference/status-constants.md),
[`INFEASIBLE`](https://www.cvxgrp.org/CVXR/reference/status-constants.md),
[`UNBOUNDED`](https://www.cvxgrp.org/CVXR/reference/status-constants.md)
