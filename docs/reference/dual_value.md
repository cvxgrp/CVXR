# Get the Dual Value of a Constraint

Returns the dual variable value(s) associated with a constraint after
solving. Returns `NULL` before the problem is solved.

## Usage

``` r
dual_value(x)
```

## Arguments

- x:

  A `Constraint` object.

## Value

A numeric matrix (single dual variable) or a list of numeric matrices
(multiple dual variables), or `NULL`.
