# Introduction to CVXR

## Overview

`CVXR` is an R package that provides an object-oriented modeling
language for convex optimization, similar to
[CVXPY](https://www.cvxpy.org/) for Python. It allows you to formulate
and solve convex optimization problems in a natural mathematical syntax.

## A Simple Example: Least Squares

Consider a simple linear regression problem where we want to estimate
parameters using a least squares criterion.

We generate synthetic data where we know the true model:

``` math
 Y = X\beta + \epsilon 
```

where $`Y`$ is a $`100 \times 1`$ vector, $`X`$ is a $`100 \times 10`$
matrix, $`\beta = (-4, -3, \ldots, 5)^\top`$ is a $`10 \times 1`$
vector, and $`\epsilon \sim N(0, 1)`$.

[`set.seed`](https://rdrr.io/r/base/Random.html)`(``123``)`` ``n`` ``<-`` ``100`` ``p`` ``<-`` ``10`` ``beta`` ``<-`` ``-``4``:``5`` `` ``X`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n`` ``*`` ``p``)``, nrow ``=`` ``n``)`` ``Y`` ``<-`` ``X`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``beta`` ``+`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n``)`

Using base R, we can estimate $`\beta`$ via `lm`:

`ls.model`` ``<-`` `[`lm`](https://rdrr.io/r/stats/lm.html)`(``Y`` ``~`` ``0`` ``+`` ``X``)`

### The CVXR formulation

The same problem can be expressed as:

``` math
\underset{\beta}{\text{minimize}} \quad \|Y - X\beta\|_2^2
```

In CVXR, this translates directly:

[`library`](https://rdrr.io/r/base/library.html)`(`[`CVXR`](https://cvxr.rbind.io)`)`` ``#> `` ``#> Attaching package: 'CVXR'`` ``#> The following objects are masked from 'package:stats':`` ``#> `` ``#> convolve, power, sd, var`` ``#> The following objects are masked from 'package:base':`` ``#> `` ``#> diag, norm, outer`` `` ``betaHat`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``p``)`` ``objective`` ``<-`` `[`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md)`(`[`sum`](https://rdrr.io/r/base/sum.html)`(``(``Y`` ``-`` ``X`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``betaHat``)``^``2``)``)`` ``problem`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(``objective``)`` ``result`` ``<-`` `[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``problem``)`` ``## use default solver`

The optimal value and estimated coefficients:

[`cat`](https://rdrr.io/r/base/cat.html)`(``"Optimal value:"``, ``result``, ``"\n"``)`` ``#> Optimal value: 97.84759`` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(``CVXR ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`value`](https://www.cvxgrp.org/CVXR/reference/value.md)`(``betaHat``)``, ``3``)``,`` `` lm ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`coef`](https://rdrr.io/r/stats/coef.html)`(``ls.model``)``, ``3``)``)`` ``#> lm`` ``#> X1 -3.920 -3.920`` ``#> X2 -3.012 -3.012`` ``#> X3 -2.125 -2.125`` ``#> X4 -0.867 -0.867`` ``#> X5 0.091 0.091`` ``#> X6 0.949 0.949`` ``#> X7 2.076 2.076`` ``#> X8 3.127 3.127`` ``#> X9 3.961 3.961`` ``#> X10 5.135 5.135`

## Adding Constraints

The real power of CVXR is the ability to add constraints easily.

### Nonnegative Least Squares

Suppose we know the $`\beta`$s should be nonnegative:

`problem`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(``objective``, constraints ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``betaHat`` ``>=`` ``0``)``)`` ``result`` ``<-`` `[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``problem``, solver ``=`` ``"CLARABEL"``)`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`value`](https://www.cvxgrp.org/CVXR/reference/value.md)`(``betaHat``)``, ``3``)`` ``#> [,1]`` ``#> [1,] 0.000`` ``#> [2,] 0.000`` ``#> [3,] 0.000`` ``#> [4,] 0.000`` ``#> [5,] 1.237`` ``#> [6,] 0.623`` ``#> [7,] 2.123`` ``#> [8,] 2.804`` ``#> [9,] 4.445`` ``#> [10,] 5.207`

### Custom Constraints

Now suppose $`\beta_2 + \beta_3 \le 0`$ and all other $`\beta`$s are
nonnegative:

`A`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``1``, ``1``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``0``, ``7``)``)``, nrow ``=`` ``1``)`` ``B`` ``<-`` `[`diag`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md)`(`[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``0``, ``0``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``1``, ``7``)``)``)`` `` ``constraint1`` ``<-`` ``A`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``betaHat`` ``<=`` ``0`` ``constraint2`` ``<-`` ``B`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``betaHat`` ``>=`` ``0`` `` ``problem`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(``objective``, constraints ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``constraint1``, ``constraint2``)``)`` ``result`` ``<-`` `[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``problem``, solver ``=`` ``"CLARABEL"``, verbose ``=`` ``TRUE``)`` ``## verbose = TRUE for details`` ``#> ────────────────────────────────── CVXR v1.9.2 ─────────────────────────────────`` ``#> ``ℹ`` Problem: 1 variable, 2 constraints (QP)`` ``#> ``ℹ`` Compilation: ``"CLARABEL"`` via CVXR::Dcp2Cone -> CVXR::CvxAttr2Constr -> CVXR::ConeMatrixStuffing -> CVXR::Clarabel_Solver`` ``#> ``ℹ`` Compile time: 0.012s`` ``#> ─────────────────────────────── Numerical solver ───────────────────────────────`` ``#> ──────────────────────────────────── Summary ───────────────────────────────────`` ``#> ``✔`` Status: optimal`` ``#> ``✔`` Optimal value: 1287.63`` ``#> ``ℹ`` Compile time: 0.012s`` ``#> ``ℹ`` Solver time: 0.012s`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`value`](https://www.cvxgrp.org/CVXR/reference/value.md)`(``betaHat``)``, ``3``)`` ``#> [,1]`` ``#> [1,] 0.000`` ``#> [2,] -2.845`` ``#> [3,] -1.711`` ``#> [4,] 0.000`` ``#> [5,] 0.664`` ``#> [6,] 1.178`` ``#> [7,] 2.329`` ``#> [8,] 2.414`` ``#> [9,] 4.212`` ``#> [10,] 4.948`

This demonstrates the chief advantage of CVXR: *flexibility*. Users can
quickly modify and re-solve a problem, making the package ideal for
prototyping new statistical methods. Its syntax is simple and
mathematically intuitive.

## Available Solvers

CVXR supports 15 solvers, both open source and commercial: Clarabel,
SCS, OSQP, HiGHS, MOSEK, Gurobi, GLPK, GLPK_MI, ECOS, ECOS_BB, CPLEX,
CVXOPT, PIQP, SCIP, and XPRESS. Smooth nonlinear programs additionally
use the IPOPT and UNO solvers.

[`installed_solvers`](https://www.cvxgrp.org/CVXR/reference/installed_solvers.md)`(``)`` ``#> [1] "CLARABEL" "SCS" "OSQP" "HIGHS" "MOSEK" "GUROBI" `` ``#> [7] "GLPK" "GLPK_MI" "ECOS" "ECOS_BB" "CPLEX" "CVXOPT" `` ``#> [13] "PIQP" "SCIP" "XPRESS"`

You can specify a solver explicitly:

[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``problem``, solver ``=`` ``"CLARABEL"``)`

## What’s New

Recent releases add disciplined nonlinear programming
(`psolve(prob, nlp = TRUE)`), bounds propagation through expressions
([`get_bounds()`](https://www.cvxgrp.org/CVXR/reference/get_bounds.md)),
a derivative / sensitivity-analysis API (`requires_grad = TRUE`), and
new atoms such as
[`convolve()`](https://www.cvxgrp.org/CVXR/reference/convolve.md). CVXR
also supports element-wise matrix indexing using R’s native idioms:

`ind`` ``<-`` `[`which`](https://rdrr.io/r/base/which.html)`(``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``Rmiss``)``, arr.ind ``=`` ``TRUE``)`` ``prob`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(`[`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md)`(``obj``)``, `[`list`](https://rdrr.io/r/base/list.html)`(``X``[``ind``]`` ``==`` ``Rmiss``[``ind``]``)``)`

For a release-by-release summary see
[`vignette("whats_new")`](https://www.cvxgrp.org/CVXR/articles/whats_new.md),
and `news(package = "CVXR")` for the full details.

## Further Reading

- The [CVXR website](https://cvxr.rbind.io) has many worked examples
- The [CVXPY documentation](https://www.cvxpy.org/) covers the
  underlying mathematical framework
- [The published paper](https://doi.org/10.18637/jss.v094.i14): Fu,
  Narasimhan, and Boyd (2020). “CVXR: An R Package for Disciplined
  Convex Optimization.” *Journal of Statistical Software*, 94(14),
  <DOI:10.18637/jss.v094.i14>.

## Session Info

[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`` ``#> R version 4.6.1 (2026-06-24)`` ``#> Platform: aarch64-apple-darwin23`` ``#> Running under: macOS Tahoe 26.6.2`` ``#> `` ``#> Matrix products: default`` ``#> BLAS: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib `` ``#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib; LAPACK version 3.12.1`` ``#> `` ``#> locale:`` ``#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8`` ``#> `` ``#> time zone: America/Los_Angeles`` ``#> tzcode source: internal`` ``#> `` ``#> attached base packages:`` ``#> [1] stats graphics grDevices utils datasets methods base `` ``#> `` ``#> other attached packages:`` ``#> [1] CVXR_1.9.2`` ``#> `` ``#> loaded via a namespace (and not attached):`` ``#> [1] Matrix_1.7-6 piqp_0.6.2 jsonlite_2.0.0 compiler_4.6.1 `` ``#> [5] highs_1.14.0-2 Rcpp_1.1.2 slam_0.1-56 cccp_0.3-3 `` ``#> [9] jquerylib_0.1.4 systemfonts_1.3.2 textshaping_1.0.5 yaml_2.3.12 `` ``#> [13] fastmap_1.2.0 clarabel_0.11.3 lattice_0.22-9 R6_2.6.1 `` ``#> [17] scip_1.10.0-3 knitr_1.51 htmlwidgets_1.6.4 backports_1.5.1 `` ``#> [21] Rcplex_0.3-8 checkmate_2.3.4 gurobi_13.0-1 desc_1.4.3 `` ``#> [25] osqp_1.0.0 bslib_0.12.0 rlang_1.3.0 cachem_1.1.0 `` ``#> [29] xfun_0.60 fs_2.1.0 sass_0.4.10 S7_0.2.2 `` ``#> [33] otel_0.2.0 cli_3.6.6 pkgdown_2.2.1 Rglpk_0.6-5.1 `` ``#> [37] digest_0.6.39 grid_4.6.1 xpress_9.8.2 gmp_0.7-5.1 `` ``#> [41] lifecycle_1.0.5 ECOSolveR_0.6.1 scs_3.2.7 evaluate_1.0.5 `` ``#> [45] codetools_0.2-20 Rmosek_11.1.2 ragg_1.5.2 rmarkdown_2.31 `` ``#> [49] tools_4.6.1 htmltools_0.5.9`
