library(cvxr)
a <- Variable()
p <- Problem(Maximize(a), list(a >= 2))
result <- solve(p, solver = "ECOS", verbose = TRUE)
result$value
