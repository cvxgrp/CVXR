# Logical OR

Returns 1 if and only if at least one argument equals 1, and 0
otherwise. For two operands, can also be written with the `|` operator:
`x | y`.

## Usage

``` r
Or(..., id = NULL)
```

## Arguments

- ...:

  Two or more boolean
  [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md)s or
  logic expressions.

- id:

  Optional integer ID (internal use).

## Value

An `Or` expression.

## See also

[`Not()`](https://www.cvxgrp.org/CVXR/reference/Not.md),
[`And()`](https://www.cvxgrp.org/CVXR/reference/And.md),
[`Xor()`](https://www.cvxgrp.org/CVXR/reference/Xor.md),
[`implies()`](https://www.cvxgrp.org/CVXR/reference/implies.md),
[`iff()`](https://www.cvxgrp.org/CVXR/reference/iff.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- Variable(boolean = TRUE)
y <- Variable(boolean = TRUE)
either <- x | y          # operator syntax
either <- Or(x, y)       # functional syntax
any3 <- Or(x, y, z)      # n-ary
} # }
```
