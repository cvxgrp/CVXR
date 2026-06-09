# Compute the gradient of a solution with respect to Parameters

Differentiates through the solution map of `problem`: populates the
`gradient` slot of each `Parameter` with the sensitivity of a
scalar-valued function of the variables (defaulting to the sum-of-x
loss; override per variable by setting `gradient(variable) <-` before
calling) with respect to that parameter. Mirrors
`cvxpy.Problem.backward()`.

## Usage

``` r
backward(problem)
```

## Arguments

- problem:

  A solved `Problem`.

## Value

The `problem` (for piping); side-effect sets `gradient(param)` on each
parameter.

## Details

Must be called after
[`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) with
`requires_grad = TRUE`.

## See also

[`derivative()`](https://www.cvxgrp.org/CVXR/reference/derivative.md),
[`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md),
[`gradient()`](https://www.cvxgrp.org/CVXR/reference/gradient.md)
