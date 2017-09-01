library(cvxr)
A <- Variable(2,2)
obj <- Maximize(LogDet(A))
p <- Problem(obj)

debug(get_problem_matrix)
result <- solve(p, solver="SCS", verbose = TRUE, max_iters=20000)
result <- solve(p)
