# What's New in CVXR 1.8

## Overview

CVXR 1.8 is a ground-up rewrite using R’s
[S7](https://rconsortium.github.io/S7/) object system, designed to
mirror [CVXPY 1.8](https://www.cvxpy.org/) for long-term
maintainability. It is approximately 4–5x faster than the previous
S4-based release and adds many new features.

This vignette summarizes the key changes from CVXR 1.x that may affect
existing users. For the introductory tutorial, see
[`vignette("cvxr_intro")`](https://www.cvxgrp.org/CVXR/articles/cvxr_intro.md).
For worked examples, visit the [CVXR website](https://cvxr.rbind.io).

## New Solve Interface

The primary solve function is now
[`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md), which
returns the optimal value directly. Variable values are extracted with
[`value()`](https://www.cvxgrp.org/CVXR/reference/value.md) and problem
status with
[`status()`](https://www.cvxgrp.org/CVXR/reference/status.md):

``` r

library(CVXR)
#> 
#> Attaching package: 'CVXR'
#> The following objects are masked from 'package:stats':
#> 
#>     power, sd, var
#> The following objects are masked from 'package:base':
#> 
#>     norm, outer
x <- Variable(2, name = "x")
prob <- Problem(Minimize(sum_squares(x)), list(x >= 1))
opt_val <- psolve(prob, solver = "CLARABEL")
opt_val
#> [1] 2
value(x)
#>      [,1]
#> [1,]    1
#> [2,]    1
status(prob)
#> [1] "optimal"
```

The old [`solve()`](https://rdrr.io/r/base/solve.html) still works and
returns a backward-compatible list:

``` r

result <- solve(prob, solver = "CLARABEL")
#> ℹ In a future CVXR release, `solve()` will return the optimal value directly
#>   (like `psolve()`).
#> ℹ The `$getValue()`/`$getDualValue()` interface will be removed.
#> ℹ New API: `psolve(prob)`, then `value(x)`, `dual_value(constr)`,
#>   `status(prob)`.
#> This message is displayed once per session.
result$value
#> [1] 2
result$status
#> [1] "optimal"
```

## API Changes

| Old (CVXR 1.x) | New (CVXR 1.8) |
|----|----|
| `solve(problem)` | `psolve(problem)` |
| `result$getValue(x)` | `value(x)` |
| `result$value` | return value of [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) |
| `result$status` | `status(problem)` |
| `result$getDualValue(con)` | `dual_value(con)` |
| `get_problem_data(prob, solver)` | `problem_data(prob, solver)` |

Old function names still work but emit once-per-session deprecation
warnings.

## 13 Solvers

CVXR 1.8 ships with four built-in solvers and supports nine additional
solvers via optional packages:

``` r

installed_solvers()
#>  [1] "CLARABEL" "SCS"      "OSQP"     "HIGHS"    "MOSEK"    "GUROBI"  
#>  [7] "GLPK"     "GLPK_MI"  "ECOS"     "ECOS_BB"  "CPLEX"    "CVXOPT"  
#> [13] "PIQP"
```

| Solver         | R Package   | Problem Classes                     |
|----------------|-------------|-------------------------------------|
| CLARABEL       | `clarabel`  | LP, QP, SOCP, SDP, ExpCone, PowCone |
| SCS            | `scs`       | LP, QP, SOCP, SDP, ExpCone, PowCone |
| OSQP           | `osqp`      | LP, QP                              |
| HIGHS          | `highs`     | LP, QP, MILP                        |
| MOSEK          | `Rmosek`    | LP, QP, SOCP, SDP, ExpCone, PowCone |
| GUROBI         | `gurobi`    | LP, QP, SOCP, MIP                   |
| ECOS / ECOS_BB | `ECOSolveR` | LP, SOCP, ExpCone (+ MIP)           |
| GLPK / GLPK_MI | `Rglpk`     | LP (+ MILP)                         |
| CPLEX          | `Rcplex`    | LP, QP, MIP                         |
| CVXOPT         | `cccp`      | LP, SOCP, SDP                       |
| PIQP           | `piqp`      | LP, QP                              |

CLARABEL is the default solver and handles the widest range of problem
types among the built-in solvers. You can specify a solver explicitly:

``` r

psolve(prob, solver = "SCS")
```

### Solver exclusion

You can temporarily exclude solvers from automatic selection without
uninstalling them:

``` r

available_solvers()
#>  [1] "CLARABEL" "SCS"      "OSQP"     "HIGHS"    "MOSEK"    "GUROBI"  
#>  [7] "GLPK"     "GLPK_MI"  "ECOS"     "ECOS_BB"  "CPLEX"    "CVXOPT"  
#> [13] "PIQP"
exclude_solvers("SCS")
available_solvers()
#>  [1] "CLARABEL" "OSQP"     "HIGHS"    "MOSEK"    "GUROBI"   "GLPK"    
#>  [7] "GLPK_MI"  "ECOS"     "ECOS_BB"  "CPLEX"    "CVXOPT"   "PIQP"
include_solvers("SCS")
```

## Disciplined Programming Extensions

### DGP (Geometric Programming)

Geometric programs are solved by converting to convex form via log-log
transformation:

``` r

x <- Variable(pos = TRUE, name = "x")
y <- Variable(pos = TRUE, name = "y")
obj <- Minimize(x * y)
constr <- list(x * y >= 40, x <= 20, y >= 2)
prob <- Problem(obj, constr)
cat("Is DGP:", is_dgp(prob), "\n")
#> Is DGP: TRUE
opt_val <- psolve(prob, gp = TRUE, solver = "CLARABEL")
cat("Optimal:", opt_val, " x =", value(x), " y =", value(y), "\n")
#> Optimal: 40  x = 9.96834  y = 4.012704
```

### DQCP (Quasiconvex Programming)

Quasiconvex problems are solved via bisection:

``` r

x <- Variable(name = "x")
prob <- Problem(Minimize(ceil_expr(x)), list(x >= 0.7, x <= 1.5))
cat("Is DQCP:", is_dqcp(prob), "\n")
#> Is DQCP: TRUE
opt_val <- psolve(prob, qcp = TRUE, solver = "CLARABEL")
cat("Optimal:", opt_val, " x =", value(x), "\n")
#> Optimal: 1  x = 0.8812532
```

### DPP (Parameterized Programming)

Parameters allow efficient re-solving when only data changes:

``` r

n <- 5
x <- Variable(n, name = "x")
lam <- Parameter(nonneg = TRUE, name = "lambda", value = 1.0)

set.seed(42)
A <- matrix(rnorm(10 * n), 10, n)
b <- rnorm(10)

prob <- Problem(Minimize(sum_squares(A %*% x - b) + lam * p_norm(x, 1)))
cat("Is DPP:", is_dpp(prob), "\n")
#> Is DPP: TRUE
psolve(prob, solver = "CLARABEL")
#> [1] 9.741033
value(x)
#>           [,1]
#> [1,] 0.1311826
#> [2,] 0.1125202
#> [3,] 0.3561130
#> [4,] 0.5394796
#> [5,] 0.7324974
```

Changing the parameter and re-solving reuses the cached compilation:

``` r

value(lam) <- 10.0
psolve(prob, solver = "CLARABEL")
#> [1] 26.58717
value(x)
#>           [,1]
#> [1,] 0.1311826
#> [2,] 0.1125202
#> [3,] 0.3561130
#> [4,] 0.5394796
#> [5,] 0.7324974
```

## Mixed-Integer Programming

Variables can be declared boolean or integer:

``` r

x <- Variable(3, integer = TRUE, name = "x")
prob <- Problem(Minimize(sum(x)), list(x >= 0.5, x <= 3.5))
psolve(prob, solver = "HIGHS")
#> [1] 3
value(x)
#>      [,1]
#> [1,]    1
#> [2,]    1
#> [3,]    1
```

## Complex Variables

CVXR supports complex-valued optimization:

``` r

z <- Variable(2, complex = TRUE, name = "z")
prob <- Problem(Minimize(p_norm(z, 2)), list(z[1] == 1 + 2i))
psolve(prob, solver = "CLARABEL")
#> [1] 2.236068
value(z)
#>      [,1]
#> [1,] 1+2i
#> [2,] 0+0i
```

## Boolean Logic

For mixed-integer formulations, boolean logic atoms are available:

``` r

b1 <- Variable(boolean = TRUE, name = "b1")
b2 <- Variable(boolean = TRUE, name = "b2")
## implies() returns an Expression; compare with == 1 to make a constraint
prob <- Problem(Maximize(b1 + b2),
                list(implies(b1, b2) == 1, b1 + b2 <= 1))
psolve(prob, solver = "HIGHS")
#> [1] 1
cat("b1 =", value(b1), " b2 =", value(b2), "\n")
#> b1 = 0  b2 = 1
```

## Decomposed Solve API

For bootstrapping or simulation, you can compile once and solve many
times:

``` r

pd <- problem_data(prob, solver = "CLARABEL")
chain <- pd$chain

# In a loop:
raw <- solve_via_data(chain, warm_start = FALSE, verbose = FALSE)
result <- problem_unpack_results(prob, raw$solution, chain, raw$inverse_data)
```

## Visualization

[`visualize()`](https://www.cvxgrp.org/CVXR/reference/visualize.md)
provides problem introspection in text or interactive HTML:

``` r

x <- Variable(3, name = "x")
prob <- Problem(Minimize(p_norm(x, 1) + sum_squares(x)),
                list(x >= -1, sum(x) == 1))
visualize(prob)
#> ── Expression Tree (MINIMIZE) ──────────────────────────────────────────────────
#> t_3 = AddExpression(...)  [convex, \mathbb{R}_+, 1x1]
#> |-- t_1 = Norm1(...)  [convex, \mathbb{R}_+, 1x1]
#> |   \-- x  [affine, 3x1]
#> \-- t_2 = QuadOverLin(...)  [convex, \mathbb{R}_+, 1x1]
#>     |-- x  [affine, 3x1]
#>     \-- 1  [constant, 1x1]
#> ── Constraints ─────────────────────────────────────────────────────────────────
#> ✓ [1] Inequality 3x1
#>       -1  [constant, 1x1]
#>       x  [affine, 3x1]
#> ✓ [2] Equality 1x1
#>       t_4 = SumEntries(...)  [affine, ?, 1x1]
#>       \-- x  [affine, 3x1]
#>       1  [constant, 1x1]
#> ── SMITH FORM ──────────────────────────────────────────────────────────────────
#>   t_{1} = phi^{Norm1}({x})
#>   t_{2} = phi^{qol}({x}, {1})
#>   t_{3} = {t_{1}} + {t_{2}}
#>   t_{4} = phi^{\Sigma}({x})
#> ── RELAXED SMITH FORM ──────────────────────────────────────────────────────────
#>   t_{1} >= phi^{Norm1}({x})
#>   t_{2} >= phi^{qol}({x}, {1})
#>   t_{3} = {t_{1}} + {t_{2}}
#>   t_{4} = 1^' {x}
#> ── CONIC FORM ──────────────────────────────────────────────────────────────────
#>   t_{1} >= phi^{Norm1}({x})    [conic form: see canonicalizer]
#>   (\frac{{1}+t_{2}}{2},  \frac{{1}-t_{2}}{2},  {x}) in Q^{n+2}
#>   {1} in Q^1
#>   t_{3} = {t_{1}} + {t_{2}}
#>   t_{4} = 1^' {x}
```

For interactive HTML output with collapsible expression trees, DCP
analysis, and curvature coloring:

``` r

visualize(prob, html = "my_problem.html")
```

## New Atoms

### Convenience atoms

| Function           | Description                             |
|--------------------|-----------------------------------------|
| `ptp(x)`           | Peak-to-peak (range): `max(x) - min(x)` |
| `cvxr_mean(x)`     | Arithmetic mean along an axis           |
| `cvxr_std(x)`      | Standard deviation                      |
| `cvxr_var(x)`      | Variance                                |
| `vdot(x, y)`       | Vector dot product                      |
| `cvxr_outer(x, y)` | Outer product                           |
| `inv_prod(x)`      | Reciprocal of product                   |
| `loggamma(x)`      | Log of gamma function                   |
| `log_normcdf(x)`   | Log of standard normal CDF              |
| `cummax_expr(x)`   | Cumulative maximum                      |
| `dotsort(X, W)`    | Weighted sorted dot product             |

### Math function dispatch

Standard R math functions work directly on CVXR expressions:

``` r

x <- Variable(3, name = "x")
abs(x)     # elementwise absolute value
#> CVXR::Abs(x)
sqrt(x)    # elementwise square root
#> Power(x, 0.5)
sum(x)     # sum of entries
#> SumEntries(x, NULL, FALSE)
max(x)     # maximum entry
#> CVXR::MaxEntries(x, NULL, FALSE)
norm(x, "2")  # Euclidean norm
#> CVXR::PnormApprox(x, 2, NULL, FALSE, 1024)
```

### Boolean logic atoms

[`Not()`](https://www.cvxgrp.org/CVXR/reference/Not.md),
[`And()`](https://www.cvxgrp.org/CVXR/reference/And.md),
[`Or()`](https://www.cvxgrp.org/CVXR/reference/Or.md),
[`Xor()`](https://www.cvxgrp.org/CVXR/reference/Xor.md),
[`implies()`](https://www.cvxgrp.org/CVXR/reference/implies.md),
[`iff()`](https://www.cvxgrp.org/CVXR/reference/iff.md) — for
mixed-integer programming with boolean variables.

### Other new atoms

- `perspective(f, s)` for perspective functions
- `FiniteSet(expr, values)` constraint for discrete optimization
- [`ceil_expr()`](https://www.cvxgrp.org/CVXR/reference/ceil_expr.md),
  [`floor_expr()`](https://www.cvxgrp.org/CVXR/reference/floor_expr.md)
  for DQCP problems
- [`condition_number()`](https://www.cvxgrp.org/CVXR/reference/condition_number.md),
  [`gen_lambda_max()`](https://www.cvxgrp.org/CVXR/reference/gen_lambda_max.md),
  [`dist_ratio()`](https://www.cvxgrp.org/CVXR/reference/dist_ratio.md)
  for DQCP

## Migration Guide

To migrate code from CVXR 1.x to 1.8:

1.  Replace `result <- solve(problem)` with `opt_val <- psolve(problem)`
2.  Replace `result$getValue(x)` with `value(x)`
3.  Replace `result$status` with `status(problem)`
4.  Replace `result$getDualValue(con)` with `dual_value(con)`
5.  Replace `axis = NA` with `axis = NULL` (axis values 1 and 2 are
    unchanged)
6.  Update solver preferences: the default is now CLARABEL (was ECOS)

All old function names continue to work with deprecation warnings.

## Further Reading

- [CVXR website](https://cvxr.rbind.io) — worked examples
- [Package reference](https://www.cvxgrp.org/CVXR/) — full API
  documentation
- [CVXPY documentation](https://www.cvxpy.org/) — mathematical framework
- Fu, Narasimhan, and Boyd (2020). “CVXR: An R Package for Disciplined
  Convex Optimization.” *Journal of Statistical Software*, 94(14).
  <doi:10.18637/jss.v094.i14>
