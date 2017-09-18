library(cvxr)
m <- 3
n <- 2
A <- matrix(rnorm(m*n), nrow = m, ncol = n)
b <- matrix(rnorm(m), nrow = m, ncol = 1)
x <- Variable(n)
objective <- Minimize(sum((A %*% x - b)^2))
constraints <- list(x >= 0, x <= 1)
p <- Problem(objective, constraints)
result <- solve(p)
