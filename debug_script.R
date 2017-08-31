library(cvxr)
library(testthat)
setwd("~/Documents/software/cvxr/tests/testthat")

# TEST: Problem isn't DCP, but still goes through.
# Warning: m less than n, problem likely degenerate
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
# Problem in Index.block_eq
# x <- t(data.frame(c(0.55, 0.25, -0.2, -0.25, -0.0, 0.4),
#                   c(0.0, 0.35, 0.2, -0.1, -0.3, -0.2)))
# n <- nrow(x)
# m <- ncol(x)
# 
# # Create and solve the model
# A <- Variable(n, n)
# b <- Variable(n)
# obj <- Maximize(LogDet(A))
# constraints <- lapply(1:m, function(i) { norm2(A %*% as.matrix(x[,i]) + b) <= 1 })
# p <- Problem(obj, constraints)
# 
# result <- solve(p)

# TEST: test_ls.R
# Problem in Index.block_eq
# m <- 100
# n <- 80
# r <- 70
# 
# set.seed(1)
# A <- matrix(rnorm(m*n), nrow = m, ncol = n)
# b <- matrix(rnorm(m), nrow = m, ncol = 1)
# G <- matrix(rnorm(r*n), nrow = r, ncol = n)
# h <- matrix(rnorm(r), nrow = r, ncol = 1)
# 
# # ||Ax-b||^2 = x^T (A^T A) x - 2(A^T b)^T x + ||b||^2
# P <- t(A) %*% A
# q <- -2 * t(A) %*% b
# r <- t(b) %*% b
# Pinv <- base::solve(P)
# 
# x <- Variable(n)
# obj <- MatrixFrac(x, Pinv) + t(q) %*% x + r
# cons <- list(G %*% x == h)
# 
# solve(Problem(Minimize(obj), cons))

# TEST: test_quad_form.R
# Should throw DCP error since P is symmetric but not definite
# P <- rbind(c(1, 0), c(0, -1))
# x <- Variable(2)
# 
# # Forming quadratic form is okay
# expect_warning(cost <- QuadForm(x, P))
# prob <- Problem(Minimize(cost), list(x == c(1, 2)))
# expect_error(solve(prob))
