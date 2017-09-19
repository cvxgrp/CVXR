library(cvxr)
N <- 5
y <- matrix(rnorm(N), nrow = N, ncol = 1)
h <- matrix(rnorm(2), nrow = 2, ncol = 1)
x <- Variable(N)
v <- Conv(h, x)
obj <- Minimize(sum(y * v[1:N]))
result <- solve(Problem(obj, list()))
