# Scalarize multiple objectives into a single objective

Transforms for combining several `Minimize`/`Maximize` objectives into
one objective for multi-objective optimization. Mirrors CVXPY's
`cvxpy.transforms.scalarize` submodule; access members with `$`:

## Usage

``` r
scalarize
```

## Format

A named list of four functions.

## Details

- `scalarize$weighted_sum(objectives, weights)` – weighted sum of
  objectives.

- `scalarize$targets_and_priorities(objectives, priorities, targets, limits = NULL, off_target = 1e-5)`
  – penalize each objective within a `[target, limit]` range; a negative
  priority flips the objective sense.

- `scalarize$max(objectives, weights)` – minimize the largest weighted
  objective term.

- `scalarize$log_sum_exp(objectives, weights, gamma = 1.0)` – smooth
  maximum; `gamma -> 0` approaches `weighted_sum`, `gamma -> Inf`
  approaches `max`.

## Examples

``` r
x <- Variable()
objs <- list(Minimize(square(x)), Minimize(square(x - 1)))
obj <- scalarize$weighted_sum(objs, c(1, 1))
if (FALSE) psolve(Problem(obj)) # \dontrun{}
```
