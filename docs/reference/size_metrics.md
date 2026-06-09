# Get Size Metrics for a Problem

Mirrors CVXPY's `Problem.size_metrics` property
(`cvxpy/problems/problem.py:486-490`, class at lines 1690-1752): returns
a `SizeMetrics` object summarising the problem's scale.

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

A `SizeMetrics` object with seven numeric fields.
