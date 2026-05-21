# Apply the derivative of the solution map to perturbations

Forward-mode counterpart of
[`backward()`](https://www.cvxgrp.org/CVXR/reference/backward.md): reads
`delta(param)` for each parameter, applies the cone-program derivative,
and writes the predicted change in each variable's optimum to
`delta(var)`. Mirrors `cvxpy.Problem.derivative()`.

## Usage

``` r
derivative(problem)
```

## Arguments

- problem:

  A solved `Problem`.

## Value

The `problem` (for piping); side-effect sets `delta(variable)` on each
variable.

## Details

Must be called after
[`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) with
`requires_grad = TRUE`.

## See also

[`backward()`](https://www.cvxgrp.org/CVXR/reference/backward.md),
[`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md),
[`delta()`](https://www.cvxgrp.org/CVXR/reference/delta.md)
