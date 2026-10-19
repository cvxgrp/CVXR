# Visualize the Canonicalization Pipeline of a CVXR Problem

Displays the Smith form decomposition of a convex optimization problem,
showing each stage of the DCP canonicalization pipeline: expression
tree, Smith form, relaxed Smith form, conic form, and (optionally)
standard cone form and solver data.

## Usage

``` r
visualize(
  problem,
  output = c("text", "json", "html", "latex", "tikz"),
  solver = NULL,
  digits = 4L,
  file = NULL,
  open = interactive(),
  doc_base = "https://cvxr.rbind.io/reference/"
)
```

## Arguments

- problem:

  A [Problem](https://www.cvxgrp.org/CVXR/reference/Problem.md) object.

- output:

  Character: output format.

  `"text"`

  : Console display (default).

  `"json"`

  : JSON data model (for interop with HTML/Python).

  `"html"`

  : Interactive D3+KaTeX HTML (Phase 2).

  `"latex"`

  : LaTeX align\* environments (Phase 3).

  `"tikz"`

  : TikZ forest tree diagrams (Phase 3).

- solver:

  Solver specification for matrix stuffing stages (4-5). `NULL`
  (default) shows only Stages 0-3 with zero overhead. `TRUE` uses the
  default solver (same as
  [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md)). A
  character string (e.g., `"Clarabel"`) uses that specific solver.

- digits:

  Integer: significant digits for displaying scalar constants.
  Integer-valued constants (0, 1, -3) always display without decimals
  regardless of this setting. Defaults to 4.

- file:

  Character: path for HTML output file. If `NULL` (default), a temporary
  file is used.

- open:

  Logical: whether to open the HTML file in a browser. Defaults to
  `TRUE` in interactive sessions.

- doc_base:

  Character: base URL for atom documentation links. Defaults to the CVXR
  pkgdown site.

## Value

For `"text"`: invisible model list. For `"json"`: a JSON string (or list
if jsonlite not available). For `"html"`: the file path (invisibly). For
other formats: the rendered output (Phase 2+).

## Examples

``` r
x <- Variable(3, name = "x")
prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))
visualize(prob)                          # Stages 0-3 only
#> ── Expression Tree (MINIMIZE) ──────────────────────────────────────────────────
#> t_1 = PnormApprox(...)  [convex, \mathbb{R}_+, 1x1]
#> \-- x  [affine, 3x1]
#> ── Constraints ─────────────────────────────────────────────────────────────────
#> ✓ [1] Inequality 3x1
#>       1  [constant, 1x1]
#>       x  [affine, 3x1]
#> ── SMITH FORM ──────────────────────────────────────────────────────────────────
#>   t_{1} = phi^{||.||_2}({x})
#> ── RELAXED SMITH FORM ──────────────────────────────────────────────────────────
#>   t_{1} >= phi^{||.||_2}({x})
#> ── CONIC FORM ──────────────────────────────────────────────────────────────────
#>   (t_{1}, {x}) in Q^{n+1}
visualize(prob, solver = TRUE)           # Stages 0-5, default solver
#> ── Expression Tree (MINIMIZE) ──────────────────────────────────────────────────
#> t_1 = PnormApprox(...)  [convex, \mathbb{R}_+, 1x1]
#> \-- x  [affine, 3x1]
#> ── Constraints ─────────────────────────────────────────────────────────────────
#> ✓ [1] Inequality 3x1
#>       1  [constant, 1x1]
#>       x  [affine, 3x1]
#> ── SMITH FORM ──────────────────────────────────────────────────────────────────
#>   t_{1} = phi^{||.||_2}({x})
#> ── RELAXED SMITH FORM ──────────────────────────────────────────────────────────
#>   t_{1} >= phi^{||.||_2}({x})
#> ── CONIC FORM ──────────────────────────────────────────────────────────────────
#>   (t_{1}, {x}) in Q^{n+1}
#> ── STANDARD CONE FORM (MOSEK) ──────────────────────────────────────────────────
#>   minimize  c'x
#>   s.t.      Ax + s = b,  s in K
#> 
#>   Variables (n):     4
#>   Constraints (m):   7
#>   nnz(A):           7
#>   Density:          25.0%
#> 
#>   K = R+^3 x Q^4
#> 
#>   Variable Mapping:
#>     x[1:1]  aux_781       auxiliary (canonicalization)
#>     x[2:4]  x             user variable (3x1)
#> 
#>   Block Structure of A:
#>     rows 1-3  (3 rows)  inequalities
#>     rows 4-7  (4 rows)  SOC(4)
#> ── SOLVER DATA MATRICES ────────────────────────────────────────────────────────
#>   A: 7 x 4, nnz 7 (25.0%)
#>   b: vector, length 7, nnz 3, range [   -1,     0]
#>   c: vector, length 4, nnz 1, range [    0,     1]
visualize(prob, solver = "CLARABEL")     # Stages 0-5, specific solver
#> ── Expression Tree (MINIMIZE) ──────────────────────────────────────────────────
#> t_1 = PnormApprox(...)  [convex, \mathbb{R}_+, 1x1]
#> \-- x  [affine, 3x1]
#> ── Constraints ─────────────────────────────────────────────────────────────────
#> ✓ [1] Inequality 3x1
#>       1  [constant, 1x1]
#>       x  [affine, 3x1]
#> ── SMITH FORM ──────────────────────────────────────────────────────────────────
#>   t_{1} = phi^{||.||_2}({x})
#> ── RELAXED SMITH FORM ──────────────────────────────────────────────────────────
#>   t_{1} >= phi^{||.||_2}({x})
#> ── CONIC FORM ──────────────────────────────────────────────────────────────────
#>   (t_{1}, {x}) in Q^{n+1}
#> ── STANDARD CONE FORM (CLARABEL) ───────────────────────────────────────────────
#>   minimize  c'x
#>   s.t.      Ax + s = b,  s in K
#> 
#>   Variables (n):     4
#>   Constraints (m):   7
#>   nnz(A):           7
#>   Density:          25.0%
#> 
#>   K = R+^3 x Q^4
#> 
#>   Variable Mapping:
#>     x[1:1]  aux_798       auxiliary (canonicalization)
#>     x[2:4]  x             user variable (3x1)
#> 
#>   Block Structure of A:
#>     rows 1-3  (3 rows)  inequalities
#>     rows 4-7  (4 rows)  SOC(4)
#> ── SOLVER DATA MATRICES ────────────────────────────────────────────────────────
#>   A: 7 x 4, nnz 7 (25.0%)
#>   b: vector, length 7, nnz 3, range [   -1,     0]
#>   c: vector, length 4, nnz 1, range [    0,     1]
if (FALSE) { # \dontrun{
visualize(prob, output = "html", solver = TRUE)
} # }
```
