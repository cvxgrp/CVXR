# Visualize the Canonicalization Pipeline of a CVXR Problem

Displays the Smith form decomposition of a convex optimization problem,
showing each stage of the DCP canonicalization pipeline: expression
tree, Smith form, relaxed Smith form, conic form, and (optionally)
standard cone form.

## Usage

``` r
visualize(
  problem,
  output = c("text", "json", "html", "latex", "tikz"),
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

  :   Console display (default).

  `"json"`

  :   JSON data model (for interop with HTML/Python).

  `"html"`

  :   Interactive D3+KaTeX HTML (Phase 2).

  `"latex"`

  :   LaTeX align\* environments (Phase 3).

  `"tikz"`

  :   TikZ forest tree diagrams (Phase 3).

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
if (FALSE) { # \dontrun{
x <- Variable(3, name = "x")
prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))
visualize(prob)
visualize(prob, output = "json")
visualize(prob, output = "html", open = FALSE)
} # }
```
