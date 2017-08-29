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
# # fit_error <- SumSquares(residuals)
# fit_error <- sum(residuals^2)
#   
# debug(SymData.presolve)
# debug(get_problem_matrix)
# result <- solve(Problem(Minimize(fit_error)), solver = "SCS")

# TEST: test_quad_form.R
# n <- 3
# 
# # Construct a random 1-d finite distribution
# v <- exp(rnorm(n))
# v <- v / sum(v)
# 
# # Construct a random positive definite matrix
# A <- matrix(rnorm(n^2), nrow = n, ncol = n)
# Q <- A %*% t(A)
# 
# # Project onto the orthogonal complement of v
# # This turns Q into a singular matrix with a known nullspace
# E <- diag(rep(1,n)) - v %*% t(v) / as.numeric(t(v) %*% v)
# Q <- E %*% (Q %*% t(E))
# 
# x <- Variable(n)
# objective <- Minimize(QuadForm(x, Q))
# prob <- Problem(objective)
# solve(prob)

# TEST: test_nonlinear_atoms.R
# LinOp data field contains Parameter object rather than its value
# v <- Variable(1)
# p <- Parameter(1, sign = "positive")
# value(p) <- 1
# # p <- 1
#  
# obj <- Minimize(KLDiv(v, p))
# prob <- Problem(obj)
# result <- solve(prob)

# TEST: test_grad.R
x <- Variable(2, name = "x")
value(x) <- c(-1,0)
expr <- Pnorm(x, 1)

# base::trace("grad", tracer = browser, exit = browser, signature = c("Atom"))
grad(expr)
