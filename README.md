# CVXR

An R modeling language for convex optimization problems.

## Note

This package is still under active development. You can install it, but crucial glue to the Canon interface which bridges any language---in our case R---and C-level libraries still requires work to do anything useful.

## Installation

To install, make sure you have the development tools for R available, including the C compilers etc.

1. Install the packages `devtools` `Rcpp`, `RcppEigen`, `BH`, `uuid`, `bitops`.
2. In R, execute
```
library(devtools)
install_github("anqif/cvxr")
```

## Support

You may experiment with the DSL expressing problems using the current version.
We are unable to provide any support until the package is officially released.

A simple (but otherwise useless) example for testing:

```
library(cvxr)

## Problem data.
m <- 1; n <- 2
A <- matrix(c(17, 19), nrow = m, byrow = TRUE)

## Construct the problem.
x <- Variable(n)
objective <- Minimize(A %*% x)
constraint <- list(x >= 1)
prob <- Problem(objective, constraint)
result <- solve(prob)

cat("Solver Status: ", result$status, "\n")
cat("Optimal Objective: ", result$value, "\n")
cat("Primal Solution:\n")
print(result$getValue(x))
cat("Dual Solution:\n")
print(result$getDualValues(constraint[[1]]))
```
