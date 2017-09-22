library(cvxr)


x_bool <- Bool()
obj <- Minimize(Square(x_bool - 0.2))
p <- Problem(obj, list())
data <- get_problem_data(p, "ECOS_BB")
base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("ECOS_BB"))
#base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("ECOS_BB"))
soln <- solve(p, solver = "ECOS_BB", verbose = TRUE)

