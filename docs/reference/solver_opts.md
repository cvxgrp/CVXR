# Create Solver Options

Constructs a structured list of solver options for use with
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md) and
[`problem_data`](https://www.cvxgrp.org/CVXR/reference/problem_data.md).
Known parameters are sorted into named slots; solver-specific parameters
are collected in `$solver_specific`.

## Usage

``` r
solver_opts(
  use_quad_obj = TRUE,
  feastol = NULL,
  reltol = NULL,
  abstol = NULL,
  num_iter = NULL,
  ...
)
```

## Arguments

- use_quad_obj:

  Logical. If `TRUE` (default), quadratic objectives use the QP matrix
  path. If `FALSE`, forces conic decomposition via `quad_form_canon`.

- feastol:

  Feasibility tolerance (solver-agnostic). Translated to solver-native
  name by internal mapping. `NULL` uses solver default.

- reltol:

  Relative tolerance. `NULL` uses solver default.

- abstol:

  Absolute tolerance. `NULL` uses solver default.

- num_iter:

  Maximum iterations. `NULL` uses solver default.

- ...:

  Solver-specific parameters passed directly to the solver (e.g.,
  `eps_abs`, `scip_params`, `mosek_params`).

## Value

A named list with class `"solver_opts"`.

## Examples

``` r
solver_opts(feastol = 1e-6)
#> $use_quad_obj
#> [1] TRUE
#> 
#> $feastol
#> [1] 1e-06
#> 
#> $reltol
#> NULL
#> 
#> $abstol
#> NULL
#> 
#> $num_iter
#> NULL
#> 
#> $solver_specific
#> list()
#> 
#> attr(,"class")
#> [1] "solver_opts"
solver_opts(use_quad_obj = FALSE, eps_abs = 1e-7)
#> $use_quad_obj
#> [1] FALSE
#> 
#> $feastol
#> NULL
#> 
#> $reltol
#> NULL
#> 
#> $abstol
#> NULL
#> 
#> $num_iter
#> NULL
#> 
#> $solver_specific
#> $solver_specific$eps_abs
#> [1] 1e-07
#> 
#> 
#> attr(,"class")
#> [1] "solver_opts"
solver_opts(scip_params = list("limits/time" = 10))
#> $use_quad_obj
#> [1] TRUE
#> 
#> $feastol
#> NULL
#> 
#> $reltol
#> NULL
#> 
#> $abstol
#> NULL
#> 
#> $num_iter
#> NULL
#> 
#> $solver_specific
#> $solver_specific$scip_params
#> $solver_specific$scip_params$`limits/time`
#> [1] 10
#> 
#> 
#> 
#> attr(,"class")
#> [1] "solver_opts"
```
