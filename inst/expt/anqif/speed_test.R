library(cvxr)
library(MASS)

m <- 500
n <- 1000
mu <- 10
sigma <- 10
set.seed(12334)

A <- matrix(rnorm(n*m, mu, sigma), nrow = n, ncol = m)
x_true <- seq(0, 100, length.out = m)
eps <- rnorm(n, 0, 1)
b <- A %*% x_true + eps

##
## Write out output file to file so that you can
## compare against cvxpy
## speed_test.py reads these two output files for A and b.
##
write.matrix(format(A, digits = 14, scientific = TRUE), file = "a.txt")
write.matrix(format(b, digits = 14, scientific = TRUE), file = "b.txt")

x <- Variable(m)
objective <- Minimize(SumSquares(A %*% x - b))
prob <- Problem(objective)
system.time(result <- solve(prob))

result$optimal_value
result$primal_values
