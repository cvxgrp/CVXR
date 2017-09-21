library(cvxr)
library(testthat)
setwd("~/Documents/software/cvxr/tests/testthat")

# Fully passing
# test_file("test_constraints.R")
# test_file("test_convolution.R")
# test_file("test_curvature.R")
# test_file("test_expressions.R")       # Add more code for testing negative indices in the future
# test_file("test_lin_ops.R")
# test_file("test_matrices.R")
# test_file("test_monotonicity.R")
# test_file("test_non_optimal.R")
# test_file("test_nonlinear_atoms.R")   # Parameters are unimplemented
# test_file("test_objectives.R")
# test_file("test_quad_form.R")
# test_file("test_quadratic.R")
# test_file("test_semidefinite_vars.R")
# test_file("test_scs.R")
# test_file("test_shape.R")
# test_file("test_sign.R")

# Passing with some comments
# test_file("test_atoms.R")             # No partial_optimize
# test_file("test_domain.R")            # No partial_optimize
# test_file("test_examples.R")          # Parameters and dual value retrieval are unimplemented
# test_file("test_grad.R")              # No partial_optimize or linearize
# test_file("test_problem.R")           # Need result retrieval function

# Failing
# test_file("test_constant_atoms.R")    # Failing on GeoMean (see below)
# test_file("test_ls.R")                # No LS solver
# test_file("test_mip_vars.R")          # No ECOS BB solver
# test_file("test_vignette.R")          # Need to finish adding examples from paper

# TEST: test_problem.R
# cvxCanon's build_matrix throws a std::bad_alloc error after hanging for a minute. This is a bug in the gmp package.
# require(gmp)
# x <- Variable(2)
# cost <- GeoMean(x)
# prob <- Problem(Maximize(cost), list(x <= 1))
# result <- solve(prob)

# TEST: test_problem.R
# Subscript out of bounds in saveValuesById
# a <- Variable(1)
# obj <- Minimize(0*a)
# p <- Problem(obj)
# result <- solve(p)

# TEST: test_problem.R
# A <- Variable(2, 2, name = "A")
# B <- Variable(2, 2, name = "B")
# C <- Variable(3, 2, name = "C")
# p <- Problem(Minimize(norm2(C[1:2,] + A)[,1]), list(C[2:3,] <= 2, C[1,] == 1, (A + B)[,1] == 3, (A + B)[,2] == 2, B == 1))
# result <- solve(p)
# result$value    # Should be exactly equal to 3

# TEST: test_examples.R
# WARN: A->p (column pointers) not strictly increasing, column 9999 empty.
# x <- Variable(5,5)
# obj <- Minimize(TotalVariation(x))
# prob <- Problem(obj)
# result <- solve(prob, solver = "SCS")
