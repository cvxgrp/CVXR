# Callback Parameter

A `Parameter` whose value is computed by a user-supplied callback
function rather than stored explicitly. Mirrors CVXPY's
`cp.CallbackParam` (`cvxpy/expressions/constants/callback_param.py`).

## Usage

``` r
CallbackParam(
  callback,
  shape = c(1L, 1L),
  name = NULL,
  id = NULL,
  latex_name = NULL,
  ...
)
```

## Arguments

- callback:

  A function (no arguments) returning the parameter's numeric value.
  Re-evaluated on every read of `value(param)`.

- shape:

  Integer vector of length 1 or 2 giving parameter dimensions. Defaults
  to `c(1, 1)` (scalar).

- name:

  Optional character string name.

- id:

  Optional integer ID. If `NULL`, a unique ID is generated.

- latex_name:

  Optional LaTeX name for visualisation.

- ...:

  Other Parameter attributes (e.g., `nonneg`, `nonpos`).

## Value

A `CallbackParam` object (subclass of `Parameter`).

## Details

Each call to `value(param)` re-evaluates the callback and validates the
returned numeric against the parameter's shape and attribute domain.
Setting via `value(param) <- v` is not allowed and signals an error.

## DPP use

If `p` and `q` are scalar `Parameter`s, the expression `p * q` is not
DPP. Wrapping it in a `CallbackParam`,

    pq <- CallbackParam(callback = function() value(p) * value(q))

yields a DPP-compliant Parameter whose value tracks `p` and `q`
automatically.

## Examples

``` r
p <- Parameter(); value(p) <- 2
q <- Parameter(); value(q) <- 3
pq <- CallbackParam(callback = function() value(p) * value(q))
value(pq)   # evaluates the callback => 6
#>      [,1]
#> [1,]    6
```
