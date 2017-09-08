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
# test_file("test_examples.R")          # WARN: A->p (column pointers) not strictly increasing, column 9999 empty
# test_file("test_grad.R")              # No partial_optimize or linearize
# test_file("test_ls.R")                # No LS solver, need result retrieval function
# test_file("test_nonlinear_atoms.R")   # Need result retrieval function. Parameters unimplemented.
# test_file("test_problem.R")           # Need result retrieval function
# test_file("test_quad_form.R")         # Need result retrieval function
# test_file("test_scs.R")               # Need result retrieval function
# test_file("test_semidefinite_vars.R") # Need result retrieval function

# Failing
# test_file("test_constant_atoms.R")
# test_file("test_mip_vars.R")          # No ECOS BB solver
# test_file("test_vignette.R")          # Need result retrieval function

# TEST: test_problem.R
# cvxCanon's build_matrix throws a std::bad_alloc error after hanging for a minute
# require(gmp)
# x <- Variable(2)
# cost <- GeoMean(x)
# prob <- Problem(Maximize(cost), list(x <= 1))
# result <- solve(prob)

# TEST: test_constant_atoms.R
# Solver returns infeasible. SCS A matrix is incorrect.
# X <- Variable(3)
# P <- Variable(3,3)
# obj <- Minimize(MatrixFrac(X, P))
# constraints <- list(X == matrix(1:3), P == diag(3))
# prob <- Problem(obj, constraints)
# 
# base::trace("format_constr", tracer = browser, exit = browser, signature = c("SDP"))
# debug(cvxr:::.lin_matrix)
# result <- solve(prob, solver = "SCS")
