# Get the Numeric Value of an Expression

Returns the numeric value of a CVXR expression, variable, or constant.
For variables, the value is set after solving a problem.

## Usage

``` r
value(x, ...)
```

## Arguments

- x:

  An expression object.

- ...:

  Not used.

## Value

A numeric matrix, or `NULL` if no value has been set.
