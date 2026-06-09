# List available solvers

Returns the names of installed solvers that are not currently excluded.
Use `exclude_solvers()` to temporarily disable solvers.

## Usage

``` r
available_solvers()

exclude_solvers(solvers)

include_solvers(solvers)

set_excluded_solvers(solvers)
```

## Arguments

- solvers:

  A character vector of solver names.

## Value

A character vector of solver names.

The current exclusion list (character vector), invisibly.

## Functions

- `exclude_solvers()`: Add solvers to the exclusion list

- `include_solvers()`: Remove solvers from the exclusion list

- `set_excluded_solvers()`: Replace the entire exclusion list

## See also

[`installed_solvers()`](https://www.cvxgrp.org/CVXR/reference/installed_solvers.md),
`exclude_solvers()`, `include_solvers()`, `set_excluded_solvers()`
