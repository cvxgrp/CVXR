library(cvxr)

# Problem data
n <- 10
m <- 5
A <- matrix(rnorm(n*m), nrow = n, ncol = m)
b <- matrix(rnorm(n), nrow = n, ncol = 1)
# gamma <- Parameter(sign = "POSITIVE")
gamma <- 1

# Construct the problem
x <- Variable(rows = m)
objective <- Minimize(SumSquares(A %*% x - b) + gamma*Pnorm(x,1))
# objective <- Minimize(sum((A %*% x - b)^2) + gamma * Norm1(x))
p <- Problem(objective)
sol <- solve(p, solver = ECOS())
