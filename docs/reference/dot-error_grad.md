# Error Gradient (None for every variable)

Returns a list keyed by variable id with every entry NULL. Used by the
chain-rule walker to propagate "could not compute" up the DAG when any
argument value is missing.

## Usage

``` r
.error_grad(expr)
```

## Arguments

- expr:

  An expression.

## Value

A list keyed by variable id (as character), each entry NULL.

## Details

Mirrors `cvxpy/utilities/grad.py:error_grad`.
