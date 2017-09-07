library(cvxr)
library(testthat)
setwd("~/Documents/software/cvxr/tests/testthat")

# Fully passing
# test_file("test_constraints.R")
# test_file("test_curvature.R")
# test_file("test_expressions.R")     # Add more code for testing negative indices in the future
# test_file("test_lin_ops.R")
# test_file("test_matrices.R")
# test_file("test_monotonicity.R")
# test_file("test_non_optimal.R")
# test_file("test_objectives.R")
# test_file("test_quadratic.R")
# test_file("test_shape.R")
# test_file("test_sign.R")

# Passing with some comments
# test_file("test_atoms.R")             # No partial_optimize
# test_file("test_convolution.R")       # Need result retrieval function
# test_file("test_domain.R")            # No partial_optimize, need result retrieval function
# test_file("test_grad.R")              # No partial_optimize or linearize
# test_file("test_ls.R")                # No LS solver, need result retrieval function
# test_file("test_nonlinear_atoms.R")   # Need result retrieval function. Parameters unimplemented.
# test_file("test_quad_form.R")         # Need result retrieval function
# test_file("test_scs.R")               # Need result retrieval function
# test_file("test_semidefinite_vars.R") # Need result retrieval function

# Failing
# test_file("test_constant_atoms.R")
# test_file("test_examples.R")
# test_file("test_mip_vars.R")          # No ECOS BB solver
# test_file("test_problem.R")
# test_file("test_vignette.R")          # Need result retrieval function

# TEST: test_problem.R
# a <- Variable(name = "a")
# b <- Variable(name = "b")
# c <- Variable(name = "c")
# 
# x <- Variable(2, name = "x")
# y <- Variable(3, name = "y")
# z <- Variable(2, name = "z")
# 
# A <- Variable(2, 2, name = "A")
# B <- Variable(2, 2, name = "B")
# C <- Variable(3, 2, name = "C")

# TEST: test_problem.R
# cvxCanon's build_matrix throws a std::bad_alloc error after hanging for a minute
# require(gmp)
# x <- Variable(2)
# cost <- GeoMean(x)
# prob <- Problem(Maximize(cost), list(x <= 1))
# result <- solve(prob)

# TEST: test_ls.R
# SCS A matrix is incorrect. The get_problem_data returns incorrect indices (V, I, J). Correct ones should be:
# V = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.4142135623730951, -1.4142135623730951, -1.0, -1.4142135623730951, -1.0]
# I = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 5.0, 6.0, 6.0, 7.0, 8.0, 7.0, 9.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
# J = [2.0, 3.0, 5.0, 6.0, 0.0, 1.0, 8.0, 9.0, 10.0, 11.0, 3.0, 4.0, 5.0, 7.0, 8.0, 9.0, 2.0, 3.0, 4.0, 6.0, 7.0, 10.0]
# x <- Variable(2)
# P <- diag(2)
# obj <- MatrixFrac(x, P)
# prob <- Problem(Minimize(obj))
# 
# # base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("Solver"))
# # debug(cvxr:::.lin_matrix)
# result <- solve(prob)

# TEST: test_examples.R
# Problem in one of the lin_ops passed to get_problem_matrix
# x <- t(data.frame(c(0.55, 0.25, -0.2, -0.25, -0.0, 0.4),
#                   c(0.0, 0.35, 0.2, -0.1, -0.3, -0.2)))
# n <- nrow(x)
# m <- ncol(x)
# 
# # Create and solve the model
# A <- Variable(n, n)
# # b <- Variable(n)
# # constraints <- lapply(1:m, function(i) { norm2(A %*% as.matrix(x[,i]) + b) <= 1 })
# obj <- Maximize(LogDet(A))
# p <- Problem(obj)
# # p <- Problem(obj, constraints)
# 
# # debug(cvxr:::.lin_matrix)
# # debug(get_problem_matrix)
# result <- solve(p)

# TEST: test_constant_atoms.R
# A <- Variable(2,2)
# obj <- Minimize(LogSumExp(A))
# constraints <- list(A == cbind(c(5,7), c(0,-3)))
# prob <- Problem(obj, constraints)
# solve(prob, solver = "ECOS")

# A <- Variable(3,1)
# P <- Variable(3,3)
# obj <- Minimize(MatrixFrac(A, P))
# constraints <- list(A == matrix(1:3), P == diag(3))
# prob <- Problem(obj, constraints)
# solve(prob, solver = "SCS")

# TEST: test_constant_atoms.R
# Should not be infeasible.
# A <- Variable(2,2)
# obj <- Minimize(SigmaMax(A))
# constraints <- list(A == cbind(c(2,0), c(0,1)))
# prob <- Problem(obj, constraints)
# result <- solve(prob, solver = "SCS")
# print(result)

# TEST: test_constant_atoms.R
# x <- Variable(3)
# obj <- Maximize(HarmonicMean(x))
# constraints <- list(x == matrix(1:3))
# prob <- Problem(obj, constraints)
# solve(prob)
