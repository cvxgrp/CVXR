# Get Expression Name

Returns a human-readable string representation of a CVXR expression,
variable, or constraint.

## Usage

``` r
name(x)
```

## Arguments

- x:

  A CVXR expression, variable, parameter, constant, or constraint.

## Value

A character string.

## Examples

``` r
x <- Variable(2, name = "x")
name(x)  # "x"
#> [1] "x"
name(x + 1)  # "x + 1"
#> [1] "x + 1"
```
