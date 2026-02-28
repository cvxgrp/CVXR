# Standard R Functions for CVXR Expressions

CVXR registers methods so that standard R functions create the
appropriate atoms when applied to `Expression` objects.

For CVXR expressions, computes the matrix/vector norm atom. For other
inputs, falls through to [`norm`](https://rdrr.io/r/base/norm.html).

For CVXR expressions, computes the standard deviation atom (ddof=0 by
default, matching CVXPY/numpy convention). For numeric inputs, falls
through to [`sd`](https://rdrr.io/r/stats/sd.html).

For CVXR expressions, computes the variance atom (ddof=0 by default).
For numeric inputs, falls through to
[`var`](https://rdrr.io/r/stats/cor.html).

For CVXR expressions, computes the outer product of two vectors. For
other inputs, falls through to
[`outer`](https://rdrr.io/r/base/outer.html).

## Usage

``` r
norm(x, type = "2", ...)

sd(x, ...)

var(x, ...)

outer(X, Y, ...)
```

## Arguments

- x:

  An Expression or numeric.

- type:

  Norm type: `"1"`, `"2"` (default), `"I"`/`"i"` (infinity), `"F"`/`"f"`
  (Frobenius).

- ...:

  For non-Expression inputs: passed to
  [`outer`](https://rdrr.io/r/base/outer.html).

- X:

  An Expression or numeric.

- Y:

  An Expression or numeric.

## Value

An Expression or numeric value.

An Expression or numeric value.

An Expression or numeric value.

An Expression or matrix.

## Math group (elementwise, via S3 group generic)

- `abs(x)`:

  Absolute value (convex, nonneg)

- `exp(x)`:

  Exponential (convex, positive)

- `log(x)`:

  Natural logarithm (concave, domain x \>= 0)

- `sqrt(x)`:

  Square root via
  [`power`](https://www.cvxgrp.org/CVXR/reference/power.md)`(x, 0.5)`
  (concave)

- `log1p(x)`:

  log(1+x) compound expression (concave)

- `log2(x)`, `log10(x)`:

  Base-2/10 logarithm

- `cumsum(x)`:

  Cumulative sum (affine)

- `cummax(x)`:

  Cumulative max (convex)

- `cumprod(x)`:

  Cumulative product

- `ceiling(x)`, `floor(x)`:

  Round up/down (MIP)

## Summary group (via S3 group generic)

- `sum(x)`:

  Sum all entries (affine)

- `max(x)`:

  Maximum entry (convex)

- `min(x)`:

  Minimum entry (concave)

## S3 generic methods

- `mean(x)`:

  Arithmetic mean; pass `axis`/`keepdims` via `...`

- `diff(x)`:

  First-order differences; also
  [`cvxr_diff`](https://www.cvxgrp.org/CVXR/reference/cvxr_diff.md)

## Masking wrappers

These mask the base/stats versions and dispatch on argument type:

- `norm(x)`:

  2-norm; use `type` for "1", "I" (infinity), "F" (Frobenius)

- `sd(x)`:

  Standard deviation (ddof=0 for expressions)

- `var(x)`:

  Variance (ddof=0 for expressions)

- `outer(X, Y)`:

  Outer product of two vector expressions

## Advanced usage

For axis-aware reductions, keepdims, or other options not available
through the standard interface, use the explicit functions:
[`cvxr_norm`](https://www.cvxgrp.org/CVXR/reference/cvxr_norm.md),
[`cvxr_mean`](https://www.cvxgrp.org/CVXR/reference/cvxr_mean.md),
[`cvxr_diff`](https://www.cvxgrp.org/CVXR/reference/cvxr_diff.md),
[`cvxr_std`](https://www.cvxgrp.org/CVXR/reference/cvxr_std.md),
[`cvxr_var`](https://www.cvxgrp.org/CVXR/reference/cvxr_var.md),
[`cvxr_outer`](https://www.cvxgrp.org/CVXR/reference/cvxr_outer.md).

## See also

[`power`](https://www.cvxgrp.org/CVXR/reference/power.md),
[`sum_entries`](https://www.cvxgrp.org/CVXR/reference/sum_entries.md),
[`max_entries`](https://www.cvxgrp.org/CVXR/reference/max_entries.md),
[`min_entries`](https://www.cvxgrp.org/CVXR/reference/min_entries.md)

[`cvxr_norm`](https://www.cvxgrp.org/CVXR/reference/cvxr_norm.md) for
the full-featured version with `axis` and `keepdims` arguments

[`cvxr_std`](https://www.cvxgrp.org/CVXR/reference/cvxr_std.md) for the
full-featured version

[`cvxr_var`](https://www.cvxgrp.org/CVXR/reference/cvxr_var.md) for the
full-featured version

[`cvxr_outer`](https://www.cvxgrp.org/CVXR/reference/cvxr_outer.md) for
the CVXR-specific version
