library(cvxr)
A <- Variable(2,2)
obj <- Maximize(LogDet(A))
p <- Problem(obj)

#debug(get_problem_matrix)
result <- solve(p, solver="SCS", verbose = TRUE, scale = 1, cg_rate=2.0, max_iters = 10000)

