# Set a Label on a Constraint

Attaches a human-readable label to a constraint for use in
visualizations and pretty-printing. Labels are visualization-only and
never affect the solver pipeline.

## Usage

``` r
set_label(constraint, label)
```

## Arguments

- constraint:

  A `Constraint` object.

- label:

  A character string label.

## Value

The modified constraint (invisibly).

## Details

Because R uses copy-on-modify semantics, you must either assign the
result back or use `set_label` fluently when building constraint lists.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- Variable(3, name = "x")

# Assign back
con <- (x >= 0)
con <- set_label(con, "non-negativity")

# Fluent use in constraint lists
constraints <- list(
  set_label(x >= 1, "lower bound"),
  set_label(sum_entries(x) <= 10, "budget")
)
} # }
```
