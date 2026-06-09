# Sampling Bounds for NLP Random Restarts

Get or set a
[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)'s
`sample_bounds` – a `(low, high)` region used to draw random initial
points in `best_of` NLP solves
([`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(prob, nlp = TRUE, best_of = n)`).
When set, it overrides the variable's value during random
initialization; when `NULL` (the default) finite variable bounds are
used instead. Supply a pair `c(low, high)` (scalars broadcast to the
variable shape) or a `list(low, high)` of per-entry vectors; set `NULL`
to clear.

## Usage

``` r
sample_bounds(x, ...)

sample_bounds(x) <- value
```

## Arguments

- x:

  A [`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md).

- ...:

  Not used.

- value:

  A `(low, high)` pair, or `NULL` to clear.

## Value

`sample_bounds(x)` returns the stored `list(low, high)` or `NULL`; the
setter returns the modified variable.
