library(cvxr)


set.seed(123)
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

##debug("SymData.presolve")
##debug("get_problem_matrix")
##base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("ECOS"))
##base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("SCS"))
result <- solve(Problem(Minimize(fit_error)), solver = "SCS", verbose=1L)
result <- solve(Problem(Minimize(fit_error)), verbose=1L)
