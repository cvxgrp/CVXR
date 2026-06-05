# 1D discrete convolution (numpy-style)

`convolve()` is CVXPY 1.9's preferred name for
[`conv`](https://www.cvxgrp.org/CVXR/reference/conv.md) (CVXPY's `conv`
class is deprecated in favor of `convolve`). For a CVXR `Expression`
argument it builds a `Convolve` atom; for plain numeric input it
computes the same numpy-style convolution the atom does
(`numpy.convolve(a, b) == stats::convolve(a, rev(b), type = "open")`),
so the numeric and Expression paths agree.

## Usage

``` r
convolve(a, b)
```

## Arguments

- a, b:

  Expressions or numeric vectors; at least one must be constant when
  building an atom.

## Value

A `Convolve` atom (for expressions) or a numeric vector.

## Details

CVXR masks [`stats::convolve`](https://rdrr.io/r/stats/convolve.html)
when attached (R reports this on load). If you specifically want stats'
circular cross-correlation or its other options, call
[`convolve`](https://rdrr.io/r/stats/convolve.html) directly.
