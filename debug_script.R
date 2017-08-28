library(cvxr)

# TEST
# A <- Variable(2, 2, name = "A")
# obj <- Minimize(0)
# dom <- domain(LogDet(A))
# prob <- Problem(obj, dom)
 
# base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("SCS"))
# base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("SCS"))
# base::trace("format_constr", tracer = browser, exit = browser, signature = c("SDP"))
# base::trace("get_objective", tracer = browser, exit = browser, signature = c("MatrixData"))
# debug(get_problem_matrix)
# debug("SymData.get_var_offsets")
# result <- solve(prob, solver = "SCS")

# TEST: test_ls.R
# Crashes in get_problem_matrix from SymData.presolve due to constraint mishandling
n <- 100
 
# Specify the true value of the variable
true_coeffs <- matrix(c(2, -2, 0.5), nrow = 3, ncol = 1)
  
# Generate data
x_data <- matrix(runif(n) * 5, nrow = n, ncol = 1)
x_data_expanded <- cbind(x_data, x_data^2, x_data^3)
y_data <- x_data_expanded %*% true_coeffs + 0.5 * matrix(runif(n, 1), nrow = n, ncol = 1)
  
slope <- Variable()
offset <- Variable()
line <- offset + x_data * slope
residuals <- line - y_data
fit_error <- SumSquares(residuals)
# fit_error <- sum(residuals^2)
  
debug("SymData.presolve")
debug("get_problem_matrix")
result <- solve(Problem(Minimize(fit_error)), solver = "SCS")

# TEST: test_expressions.R
# C <- Variable(3, 2, name = "C")
# C[2,1]

