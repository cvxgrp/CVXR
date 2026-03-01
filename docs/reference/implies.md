# Logical Implication

Logical implication: x =\> y. Returns 1 unless x = 1 and y = 0.
Equivalent to `Or(Not(x), y)`.

## Usage

``` r
implies(x, y)
```

## Arguments

- x, y:

  Boolean [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md)s
  or logic expressions.

## Value

An [Or](https://www.cvxgrp.org/CVXR/reference/Or.md) expression
representing `!x | y`.

## See also

[`iff()`](https://www.cvxgrp.org/CVXR/reference/iff.md),
[`Not()`](https://www.cvxgrp.org/CVXR/reference/Not.md),
[`And()`](https://www.cvxgrp.org/CVXR/reference/And.md),
[`Or()`](https://www.cvxgrp.org/CVXR/reference/Or.md),
[`Xor()`](https://www.cvxgrp.org/CVXR/reference/Xor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- Variable(boolean = TRUE)
y <- Variable(boolean = TRUE)
expr <- implies(x, y)
} # }
```
