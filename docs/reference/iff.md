# Logical Biconditional

Logical biconditional: x \<=\> y. Returns 1 if and only if x and y have
the same value. Equivalent to `Not(Xor(x, y))`.

## Usage

``` r
iff(x, y)
```

## Arguments

- x, y:

  Boolean [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md)s
  or logic expressions.

## Value

A [Not](https://www.cvxgrp.org/CVXR/reference/Not.md) expression
wrapping [Xor](https://www.cvxgrp.org/CVXR/reference/Xor.md).

## See also

[`implies()`](https://www.cvxgrp.org/CVXR/reference/implies.md),
[`Not()`](https://www.cvxgrp.org/CVXR/reference/Not.md),
[`And()`](https://www.cvxgrp.org/CVXR/reference/And.md),
[`Or()`](https://www.cvxgrp.org/CVXR/reference/Or.md),
[`Xor()`](https://www.cvxgrp.org/CVXR/reference/Xor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- Variable(boolean = TRUE)
y <- Variable(boolean = TRUE)
expr <- iff(x, y)
} # }
```
