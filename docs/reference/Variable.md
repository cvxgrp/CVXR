# Create an Optimization Variable

Constructs a variable to be used in a CVXR optimization problem.
Variables are decision variables that the solver optimizes over.

## Usage

``` r
Variable(shape = c(1L, 1L), name = NULL, var_id = NULL, latex_name = NULL, ...)
```

## Arguments

- shape:

  Integer vector of length 1 or 2 giving the variable dimensions. A
  scalar `n` is interpreted as `c(n, 1)`. Defaults to `c(1, 1)`
  (scalar).

- name:

  Optional character string name for the variable. If `NULL`, an
  automatic name `"var<id>"` is generated.

- var_id:

  Optional integer ID. If `NULL`, a unique ID is generated.

- latex_name:

  Optional character string giving a custom LaTeX name for use in
  visualizations. For example, `"\\mathbf{x}"`. If `NULL` (default),
  visualizations auto-generate a LaTeX name.

- ...:

  Additional attributes: `nonneg`, `nonpos`, `PSD`, `NSD`, `symmetric`,
  `boolean`, `integer`, etc.

## Value

A `Variable` object (inherits from `Leaf` and `Expression`).

## Examples

``` r
x <- Variable(3)        # 3x1 column vector
X <- Variable(c(2, 3))  # 2x3 matrix
y <- Variable(2, nonneg = TRUE)  # non-negative variable
z <- Variable(3, name = "z", latex_name = "\\mathbf{z}")  # custom LaTeX
```
