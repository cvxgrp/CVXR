# Gradient of a Constant Expression

Returns a Variable -\> Jacobian list for an expression that is constant
in all of its variables. Each Jacobian is the appropriate-shape zero (a
scalar `0` for scalar/scalar; a sparse zero matrix otherwise).

## Usage

``` r
.constant_grad(expr)
```

## Arguments

- expr:

  An expression.

## Value

A list keyed by variable id (as character), each entry the zero
Jacobian.

## Details

Mirrors `cvxpy/utilities/grad.py:constant_grad`.
