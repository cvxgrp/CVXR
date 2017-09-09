library(cvxr)
x <- Variable(3)
P <- diag(3)
obj <- MatrixFrac(x, P)
prob <- Problem(Minimize(obj), list(x == c(1,5,3)))

# base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("Solver"))
# debug(cvxr:::.lin_matrix)
data <- get_problem_data(prob, "SCS")
result <- solve(prob, solver = "SCS", verbose=TRUE)

