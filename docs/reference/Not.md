# Logical NOT

Returns `1 - x`, flipping 0 to 1 and 1 to 0. Can also be written with
the `!` operator: `!x`.

## Usage

``` r
Not(x, id = NULL)
```

## Arguments

- x:

  A boolean
  [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md) or logic
  expression.

- id:

  Optional integer ID (internal use).

## Value

A `Not` expression.

## See also

[`And()`](https://www.cvxgrp.org/CVXR/reference/And.md),
[`Or()`](https://www.cvxgrp.org/CVXR/reference/Or.md),
[`Xor()`](https://www.cvxgrp.org/CVXR/reference/Xor.md),
[`implies()`](https://www.cvxgrp.org/CVXR/reference/implies.md),
[`iff()`](https://www.cvxgrp.org/CVXR/reference/iff.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- Variable(boolean = TRUE)
not_x <- !x          # operator syntax
not_x <- Not(x)      # functional syntax
} # }
```
