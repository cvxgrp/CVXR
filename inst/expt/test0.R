library(cvxr)

## Problem data.
n <- 2
A <- diag(c(1, 1))
# Construct the problem.
x <- Variable(n)
objective <- Minimize(sum(A %*% x))
#objective <- Minimize(SumSquares(A %*% x))
constraint <- list(1 <= x)
prob <- Problem(objective, constraint)

##base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("ECOS"))

##base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("ECOS"))

##debug(solve)
##debug(build_lin_op_tree)
##debug(get_matrix_data)

#base::trace("", tracer = browser, exit = browser, signature = c("ECOS"))
##debug(cvxr:::.lin_matrix)
#debug(cvxr:::build_lin_op_tree)
#debug(get_problem_matrix)
#base::trace("format_results", tracer = browser, exit = browser, signature = c("SCS"))
solve(prob, solver = "ECOS")

solve(prob, solver = "SCS")
