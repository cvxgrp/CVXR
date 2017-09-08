library(cvxr)

x <- Variable(2)
P <- diag(2)
obj <- MatrixFrac(x, P)
prob <- Problem(Minimize(obj))

base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("Solver"))
# debug(cvxr:::.lin_matrix)
##debug(set_slice_data)
debug(get_problem_matrix)
data <- get_problem_data(prob, "SCS")

##result <- solve(prob)

