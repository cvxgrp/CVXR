# Split a primal solution into per-variable arrays

Mirrors `ParamConeProg.split_solution` in
`cone_matrix_stuffing.py:282-302`.

## Usage

``` r
split_solution(param_prog, sltn, active_vars = NULL)
```

## Arguments

- param_prog:

  A `ParamConeProg`.

- sltn:

  Numeric vector of length `x_length` – the primal `x` from the conic
  forward solve.

- active_vars:

  Optional character vector of variable ids to restrict the output to.
  Default: all variables.

## Value

A named list mapping `as.character(var_id)` to a numeric array of the
variable's shape.
