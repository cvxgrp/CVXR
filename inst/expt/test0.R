library(cvxr)

# Problem data.
n <- 2
A <- diag(c(1, 1))
# Construct the problem.
x <- Variable(n)
objective <- Minimize(SumEntries(A %*% x))
constraint <- list( 1 <= x)
prob <- Problem(objective, constraint)

##base::trace("Solver.solve", tracer=browser, exit = browser, signature = c("ECOS"))

#base::trace("Solver.get_problem_data", tracer=browser, exit = browser, signature = c("ECOS"))

##debug(solve)
##debug(build_lin_op_tree)
##debug(get_matrix_data)

#base::trace("", tracer=browser, exit = browser, signature = c("ECOS"))
##debug(cvxr:::.lin_matrix)
#debug(cvxr:::build_lin_op_tree)
debug(get_problem_matrix)
solve(prob)

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
