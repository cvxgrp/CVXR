# Logical XOR

For two arguments: result is 1 iff exactly one is 1. For n arguments:
result is 1 iff an odd number are 1 (parity).

## Usage

``` r
Xor(..., id = NULL)
```

## Arguments

- ...:

  Two or more boolean
  [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md)s or
  logic expressions.

- id:

  Optional integer ID (internal use).

## Value

A `Xor` expression.

## Details

Note: R's `^` operator is used for
[`power()`](https://www.cvxgrp.org/CVXR/reference/power.md), so `Xor` is
functional syntax only.

## See also

[`Not()`](https://www.cvxgrp.org/CVXR/reference/Not.md),
[`And()`](https://www.cvxgrp.org/CVXR/reference/And.md),
[`Or()`](https://www.cvxgrp.org/CVXR/reference/Or.md),
[`implies()`](https://www.cvxgrp.org/CVXR/reference/implies.md),
[`iff()`](https://www.cvxgrp.org/CVXR/reference/iff.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- Variable(boolean = TRUE)
y <- Variable(boolean = TRUE)
exclusive <- Xor(x, y)
} # }
```
