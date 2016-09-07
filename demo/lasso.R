library(cvxr)

# Problem data
n <- 10
m <- 5
A <- matrix(rnorm(n*m), nrow = n, ncol = m)
b <- matrix(rnorm(n), nrow = n, ncol = 1)
gamma <- Parameter(sign = "POSITIVE")

# Construct the problem
x <- Variable(rows = m)
objective <- Minimize(SumSquares(A*x - b) + gamma*Pnorm(x,1))
p <- Problem(objective)
solve(p, solver = ECOS())
