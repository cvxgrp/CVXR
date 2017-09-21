## Problem data
##set.seed(123)
m <- 3
n <- 2
## A <- matrix(c(-0.560475646552213,0.070508391424576,
##               -0.23017748948328,0.129287735160946,
##               1.55870831414912,1.71506498688328),
##             byrow = TRUE, ncol = 2)
## b <- matrix(c(0.460916205989202,
##               -1.26506123460653,
##               -0.686852851893526), ncol = 1)
A <- matrix(c(1, 1,
              -1, 2,
              -2, -3), byrow = TRUE, ncol = 2)

b <- matrix(c(1,2,3), nrow=3)

## Construct the problem
x <- Variable(n)
##objective <- Minimize(sum((A %*% x - b)^2))
objective <- Minimize(SumSquares(A %*% x - b))
constraints <- list(x >= 1, x <= 2)
p <- Problem(objective, constraints)

## The optimal objective is returned by solve(p)
result <- solve(p, verbose = TRUE)
## The optimal value for x is stored in result$x
## print(result$x)
## The optimal Lagrange multiplier for a constraint is stored in constraint$dual_value
## print(constraints[1].dual_value)
