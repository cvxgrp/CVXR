# Adjoint of `split_solution`

Given a named list mapping variable ids to their (shape-aligned) deltas,
packs them into a single length-`x_length` numeric vector in
`var_id_to_col` order. Mirrors `ParamConeProg.split_adjoint` in
`cone_matrix_stuffing.py:304-318`.

## Usage

``` r
split_adjoint(param_prog, del_vars)
```

## Arguments

- param_prog:

  A `ParamConeProg`.

- del_vars:

  A named list of `as.character(var_id)` -\> array.

## Value

Numeric vector of length `x_length`.
