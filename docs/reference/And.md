# Logical AND

Returns 1 if and only if all arguments equal 1, and 0 otherwise. For two
operands, can also be written with the `&` operator: `x & y`.

## Usage

``` r
And(..., id = NULL)
```

## Arguments

- ...:

  Two or more boolean
  [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md)s or
  logic expressions.

- id:

  Optional integer ID (internal use).

## Value

An `And` expression.

## See also

[`Not()`](https://www.cvxgrp.org/CVXR/reference/Not.md),
[`Or()`](https://www.cvxgrp.org/CVXR/reference/Or.md),
[`Xor()`](https://www.cvxgrp.org/CVXR/reference/Xor.md),
[`implies()`](https://www.cvxgrp.org/CVXR/reference/implies.md),
[`iff()`](https://www.cvxgrp.org/CVXR/reference/iff.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- Variable(boolean = TRUE)
y <- Variable(boolean = TRUE)
both <- x & y            # operator syntax
both <- And(x, y)        # functional syntax
all3 <- And(x, y, z)     # n-ary
} # }
```
