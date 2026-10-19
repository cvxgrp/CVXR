# Get Size Metrics for a Problem

Mirrors CVXPY's `Problem.size_metrics` property
(`cvxpy/problems/problem.py:515-521`, class at lines 1592-1654): returns
a `SizeMetrics` object summarizing the problem's scale.

## Usage

``` r
size_metrics(x, ...)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

- ...:

  Not used.

## Value

A `SizeMetrics` object with six numeric fields.
