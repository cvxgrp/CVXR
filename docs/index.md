# CVXR

CVXR provides an object-oriented modeling language for convex
optimization, similar to [CVXPY](https://www.cvxpy.org/),
[CVX](https://cvxr.com/cvx/), [YALMIP](https://yalmip.github.io/), and
[Convex.jl](https://jump.dev/Convex.jl/stable/). It allows you to
formulate convex optimization problems in natural mathematical syntax
rather than the restrictive standard form required by most solvers. You
specify an objective and a set of constraints by combining constants,
variables, and parameters using a library of functions with known
mathematical properties. CVXR then applies signed [disciplined convex
programming
(DCP)](https://web.stanford.edu/~boyd/papers/pdf/disc_cvx_prog.pdf) to
verify the problem’s convexity. Once verified, the problem is converted
into standard conic form and passed to an appropriate backend solver.

This version is a ground-up rewrite built on the
[S7](https://rconsortium.github.io/S7/) object system, designed to
mirror CVXPY 1.9.1 closely. It is **~4–5x faster** than the previous
S4-based release, ships with 15 solvers (4 built-in), and supports DCP,
DGP, DQCP, disciplined nonlinear programming (DNLP), complex variables,
mixed-integer programming, warm-starting, and a derivative /
sensitivity-analysis API.

For tutorials, worked examples, and the full story, visit the [CVXR
website](https://cvxr.rbind.io).

## Installation

Install the released version from CRAN:

[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"CVXR"``)`

Or install the development version from GitHub:

`# install.packages("pak")`` ``pak``::`[`pak`](https://pak.r-lib.org/reference/pak.html)`(``"cvxgrp/CVXR"``)`

## Quick example

[`library`](https://rdrr.io/r/base/library.html)`(`[`CVXR`](https://cvxr.rbind.io)`)`` `` ``# Data`` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`` ``n`` ``<-`` ``50``; ``p`` ``<-`` ``10`` ``X`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n`` ``*`` ``p``)``, ``n``, ``p``)`` ``beta_true`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``1``, ``5``)``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``0``, ``5``)``)`` ``y`` ``<-`` ``X`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``beta_true`` ``+`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n``, sd ``=`` ``0.5``)`` `` ``# Problem`` ``beta`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``p``)`` ``objective`` ``<-`` `[`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md)`(`[`sum_squares`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md)`(``y`` ``-`` ``X`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``beta``)`` ``+`` ``0.1`` ``*`` `[`p_norm`](https://www.cvxgrp.org/CVXR/reference/p_norm.md)`(``beta``, ``1``)``)`` ``prob`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(``objective``)`` `` ``# Solve (Clarabel is the default solver)`` ``result`` ``<-`` `[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``prob``)`` ``result`` ``# optimal value`` ``estimated`` ``<-`` `[`value`](https://www.cvxgrp.org/CVXR/reference/value.md)`(``beta``)`` ``# coefficient estimates`

## Documentation

- **Tutorials and examples**: <https://cvxr.rbind.io>
- **Package reference**: <https://www.cvxgrp.org/CVXR/>
- **Paper**: Fu, Narasimhan, and Boyd (2020). “CVXR: An R Package for
  Disciplined Convex Optimization.” *Journal of Statistical Software*,
  94(14), 1–34.
  [doi:10.18637/jss.v094.i14](https://doi.org/10.18637/jss.v094.i14)

If you use CVXR in your work, please cite the paper above
(`citation("CVXR")`).

## License

Apache License 2.0
