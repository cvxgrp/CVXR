# Create a Parameter

Constructs a parameter whose numeric value can be changed without
re-canonicalizing the problem. Parameters are treated as constants for
DCP purposes but their value can be updated between solves.

## Usage

``` r
Parameter(
  shape = c(1L, 1L),
  name = NULL,
  value = NULL,
  id = NULL,
  latex_name = NULL,
  ...
)
```

## Arguments

- shape:

  Integer vector of length 1 or 2 giving the parameter dimensions. A
  scalar `n` is interpreted as `c(n, 1)`. Defaults to `c(1, 1)`
  (scalar).

- name:

  Optional character string name. If `NULL`, an automatic name
  `"param<id>"` is generated.

- value:

  Optional initial numeric value.

- id:

  Optional integer ID.

- latex_name:

  Optional character string giving a custom LaTeX name for use in
  visualizations. For example, `"\\gamma"`. If `NULL` (default),
  visualizations auto-generate a LaTeX name.

- ...:

  Additional attributes: `nonneg`, `nonpos`, etc.

## Value

A `Parameter` object (inherits from `Leaf` and `Expression`).

## Examples

``` r
p <- Parameter()
value(p) <- 5
p_vec <- Parameter(3, nonneg = TRUE)
gamma <- Parameter(1, name = "gamma", latex_name = "\\gamma")
```
