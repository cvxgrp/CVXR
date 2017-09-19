library(cvxr)

X <- Semidef(2)
Y <- Variable(2, 2)
F <- matrix(c(1, 0, 0, -1), byrow=TRUE, nrow=2)
obj <- Minimize(sum((X - F)^2))
p <- Problem(obj)
result <- solve(p, verbose = TRUE)
result$getValue(obj@expr)


