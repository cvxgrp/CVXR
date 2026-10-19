# 1D discrete convolution (numpy-style)

`convolve()` is CVXPY 1.9's preferred name for
[`conv`](https://www.cvxgrp.org/CVXR/reference/conv.md) (CVXPY's `conv`
class is deprecated in favor of `convolve`). For a CVXR `Expression`
argument it builds a `Convolve` atom, the numpy-style convolution
(`numpy.convolve(x, y) == stats::convolve(x, rev(y), type = "open")`).
If neither argument is an Expression the call is forwarded as written to
[`convolve`](https://rdrr.io/r/stats/convolve.html), so the result is
exactly stats' (circular by default).

## Usage

``` r
convolve(x, y, conj = TRUE, type = c("circular", "open", "filter"))
```

## Arguments

- x, y:

  Expressions or numeric vectors; at least one must be constant when
  building an atom.

- conj, type:

  As in [`convolve`](https://rdrr.io/r/stats/convolve.html); only for
  plain numeric input, an error with an Expression.

## Value

A `Convolve` atom (for expressions) or a numeric vector.
