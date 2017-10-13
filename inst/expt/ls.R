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

library(cvxr)
set.seed(13893)
n <- 100
p <- 10

X <- matrix(rnorm(p * n), nrow=n)
beta <- c(1:2, rep(0.0, 4), rep(2, 4))
y <- X %*% beta + rnorm(n)

beta_computed <- Variable(p)
objective <- Minimize(SumSquares(X %*% beta_computed - y))
p <- Problem(objective)

## The optimal objective is returned by solve(p)
result <- solve(p, verbose = TRUE)g
## The optimal value for x is stored in result$x
## print(result$x)
## The optimal Lagrange multiplier for a constraint is stored in constraint$dual_value
## print(constraints[1].dual_value)

statistic <- function(data, indices) {
    d <- data[indices, ]
    p <- ncol(d) - 1L
    ys <- d[, 1, drop = FALSE]
    Xs <- d[, seq_len(p) + 1L, drop = FALSE]
    betaHat <- Variable(p)
    objective <- Minimize(SumSquares(Xs %*% betaHat - ys))
    prob <- Problem(objective)
    result <- solve(prob, verbose = FALSE)
    result$getValue(betaHat)
}
w <- cbind(y, X)
library(boot)
ans <- boot(data = w, statistic = statistic, R = 999)







m <- lm(y ~ 0 + X)
cbind(coef(m), result$getValue(beta_computed))

p <- Problem(objective, list(beta_computed > 5))
