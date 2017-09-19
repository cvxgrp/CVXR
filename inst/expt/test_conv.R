library(cvxr)

n <- 3
x <- Variable(n)
f <- c(1, 2, 3)
g <- c(0, 1, 0.5)
f_conv_g <- c(0, 1, 2.5, 4, 1.5)
expr <- Conv(f, x)

## Matrix stuffing
prob <- Problem(Minimize(Norm(expr, 1)), list(x == g))
result <- solve(prob)

result$getValue(x)

result$getValue(expr)
