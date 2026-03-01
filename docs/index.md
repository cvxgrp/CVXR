# CVXR

CVXR provides an object-oriented modeling language for convex
optimization, similar to [CVXPY](https://www.cvxpy.org/),
[CVX](http://cvxr.com/cvx/), [YALMIP](https://yalmip.github.io/), and
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
mirror CVXPY 1.8 closely. It is **~4–5x faster** than the previous
S4-based release, ships with 13 solvers (4 built-in), and supports DCP,
DGP, DQCP, complex variables, mixed-integer programming, and
warm-starting.

For tutorials, worked examples, and the full story, visit the [CVXR
website](https://cvxr.rbind.io).

## Installation

Install the released version from CRAN:

``` r

install.packages("CVXR")
```

Or install the development version from GitHub:

``` r

# install.packages("pak")
pak::pak("cvxgrp/CVXR")
```

## Quick example

``` r

library(CVXR)

# Data
set.seed(42)
n <- 50; p <- 10
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(rep(1, 5), rep(0, 5))
y <- X %*% beta_true + rnorm(n, sd = 0.5)

# Problem
beta <- Variable(p)
objective <- Minimize(sum_squares(y - X %*% beta) + 0.1 * p_norm(beta, 1))
prob <- Problem(objective)

# Solve (Clarabel is the default solver)
result <- psolve(prob)
result                          # optimal value
estimated <- value(beta)        # coefficient estimates
```

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
