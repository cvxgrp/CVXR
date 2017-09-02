library(cvxr)
library(testthat)
setwd("~/Documents/software/cvxr/tests/testthat")
# test_file("test_problem.R")

# TEST: test_problem.R
a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

# TEST: test_problem.R
# Problem data is correct. ECOSolveR should handle trivial constant problems.
# c <- matrix(c(1,-1), nrow = 2, ncol = 1)
# p <- Problem(Minimize(MaxElemwise(t(c), 2, 2 + t(c))[2]))
# base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("ECOS"))
# result <- solve(p)


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

# TEST: test_problem.R
# Problem in .cache_to_matrix because (V, I, J) are too short for dimensions. Variable size not promoted.
# c <- cbind(c(1,-1), c(2,-2))
# expr <- NormInf(MulElemwise(c, a))
# p <- Problem(Minimize(expr))
 
# base::trace(cvxr:::Solver.get_problem_data, tracer = browser, exit = browser, signature = c("Solver"))
# base::trace(cvxr::get_objective, tracer = browser, exit = browser, signature = c("MatrixData"))
# result <- solve(p, solver = "ECOS")

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

