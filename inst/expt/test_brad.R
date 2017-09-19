library(cvxr)
m <- 30
n <- 20
A <- matrix(rnorm(m*n), nrow = m, ncol = n)
b <- matrix(rnorm(m), nrow = m, ncol = 1)
x <- Variable(n)
objective <- Minimize(sum(Norm(A %*% x - b, 1)))
p <- Problem(objective)
result <- solve(p)
