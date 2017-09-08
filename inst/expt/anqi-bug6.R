library(cvxr)

x <- matrix(c(1,-1), nrow = 2, ncol = 1)
p <- Problem(Minimize(MaxElemwise(t(x), 2, 2 + t(x))[2]))
solve(p)
