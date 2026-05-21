# Adjoint of the parameter -\> (c, d, A, b) tensor map

Given derivatives `(delc, delA, delb)` of a downstream objective with
respect to the conic problem data, returns derivatives with respect to
each `Parameter` (keyed by parameter id). Mirrors
`ParamConeProg.apply_param_jac` in `cone_matrix_stuffing.py:242-280`.

## Usage

``` r
apply_param_jac(param_prog, delc, delA, delb, active_params = NULL)
```

## Arguments

- param_prog:

  A `ParamConeProg`.

- delc:

  Numeric vector of length `x_length`.

- delA:

  Sparse matrix of shape `(m, x_length)` (same shape as the conic
  constraint matrix `A`).

- delb:

  Numeric vector of length `m`.

- active_params:

  Optional character vector of parameter ids to restrict the output to.
  Default: all parameters.

## Value

A named list mapping `as.character(param_id)` to a numeric array of the
parameter's shape.

## Details

Reusing the tensor identity that `apply_parameters` exploits in the
forward direction: `c = c_tensor[1:x_length, ] %*% p` (omitting the
offset row), `vec(A) // b = A_tensor %*% p` (column-major; b is last
block), the adjoint is just the transpose of the same tensors applied to
the stacked (delc, vec(delA), delb).
