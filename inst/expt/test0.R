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

# The optimal objective is returned by prob.solve().
result = prob.solve()
# The optimal value for x is stored in x.value.
print x.value
# The optimal Lagrange multiplier for a constraint
# is stored in constraint.dual_value.
print constraints[0].dual_value

tmp <- list()
##C_objective <- build_lin_op_tree(objective, tmp)

root_linR <- objective

m <- 3
n <- 2
A <- matrix(c(1, 2, 3, 4, 2, 1), nrow=m, byrow=TRUE)
b <- matrix(c(1, 2, 3), nrow=m)
x <- Variable(n)
objective <- Minimize(SumSquares(A %*% x - b))
constraint <- list(1 <= x)
prob <- Problem(objective, constraint)
solve(prob)
