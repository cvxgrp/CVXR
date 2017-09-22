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

n <- 100
p <- 10

X <- matrix(rnorm(p * n), nrow=n)
beta <- 1:p + rnorm(p)
y <- X %*% beta + rnorm(n)

beta_computed <- Variable(p)
objective <- Minimize(SumSquares(X %*% beta_computed - y))
p <- Problem(objective)

## The optimal objective is returned by solve(p)
result <- solve(p, verbose = TRUE)
## The optimal value for x is stored in result$x
## print(result$x)
## The optimal Lagrange multiplier for a constraint is stored in constraint$dual_value
## print(constraints[1].dual_value)

m <- lm(y ~ 0 + X)
p <- Problem(objective, list(beta_computed > 5))
