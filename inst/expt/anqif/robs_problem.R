library(cvxr)
set.seed(34)
n <- 20
p <- 3
x <- abs(matrix(rnorm(n * p), n, p))
y <- 5 * (log(x[, 2]) - log(x[, 1])) + rnorm(n)
y <- y - mean(y)
z <- cbind(x, x)
lambda <- .1
y <- as.matrix(y)
theta <- Variable(2 * p)
objective <- 0.5 * SumSquares(z %*% theta - y) + lambda * (Norm(theta[1:2])+ Norm(theta[3:4]) + Norm(theta[5:6]))

p <- Problem(objective = Minimize(objective))
sol <- solve(p)
