library(cvxr)

# A <- Variable(2, 2, name = "A")
# obj <- Minimize(0)
# dom <- domain(LogDet(A))
# prob <- Problem(obj, dom)
# 
# base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("ECOS"))
# 
# base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("SCS"))
## base::trace("format_constr", tracer = browser, exit = browser, signature = c("SDP"))
## debug(get_problem_matrix)
## debug("SymData.get_var_offsets")
# solve(prob, solver = "SCS")

# TEST: test_atoms.R
# Passing wrong dims to ECOS.
# SCS matrices A and b are too large.
y <- Variable(1)
p1 <- Problem(Minimize(exp(y)))
# base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("SCS"))
base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("SCS"))
result1 <- solve(p1, solver = "SCS")

# TEST: test_ls.R
# Crashes in get_problem_matrix from SymData.presolve due to constraint mishandling
# n <- 100
# 
# # Specify the true value of the variable
# true_coeffs <- matrix(c(2, -2, 0.5), nrow = 3, ncol = 1)
# 
# # Generate data
# x_data <- matrix(runif(n) * 5, nrow = n, ncol = 1)
# x_data_expanded <- cbind(x_data, x_data^2, x_data^3)
# y_data <- x_data_expanded %*% true_coeffs + 0.5 * matrix(runif(n, 1), nrow = n, ncol = 1)
# 
# slope <- Variable()
# offset <- Variable()
# line <- offset + x_data * slope
# residuals <- line - y_data
# fit_error <- SumSquares(residuals)
# 
# debug("SymData.presolve")
# debug("get_problem_matrix")
# solve(Problem(Minimize(fit_error)), solver = "ECOS")

# TEST: test_non_optimal.R
# Should I pass axis = 1 or 2 into SOCAxis in Sqrt.graph_implementation?
# SCS matrices A and b are incorrect. Freezing in MatrixData initialization.
# x <- Variable(1)
# prob <- Problem(Maximize(sqrt(x)))
# # base::trace("format_constr", tracer = browser, exit = browser, signature = c("SOCAxis"))
# # base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("SCS"))
# debug(cvxr:::.lin_matrix)
# result <- solve(prob, solver = "SCS")
