# Tests for the NLP diff_engine adapter (CVXPY 1.9.0, Step 5a).
#
# The adapter (.de_C_problem) tree-walks a CVXR Problem onto the sparsediff C
# differentiation engine and exposes the value/gradient/Jacobian/Hessian oracle.
# Cross-checked against the sparsediff README worked example:
#   f(x) = sum(exp(x)) s.t. sum(x) == 1.5
#   at u = (0, 0.5, 1): objective = 5.367003, gradient = exp(u),
#   constraint Jacobian = (1,1,1), Hessian = diag(exp(u)).
#
# CVXPY SOURCE: reductions/solvers/nlp_solvers/diff_engine/{c_problem,converters,
#   registry,helpers}.py

library(CVXR)

.de_make_prob <- function() {
  x <- Variable(3)
  Problem(Minimize(sum_entries(Exp(x))), list(sum_entries(x) == 1.5))
}

## @cvxpy NONE
test_that("C_problem: objective value and gradient for sum(exp(x))", {
  skip_if_not_installed("sparsediff")
  cp <- .de_C_problem(.de_make_prob(), verbose = FALSE)
  cp$init_jacobian_coo()
  cp$init_hessian_coo_lower_tri()
  u <- c(0, 0.5, 1)
  expect_equal(cp$objective_forward(u), sum(exp(u)), tolerance = 1e-6)
  expect_equal(as.numeric(cp$gradient()), exp(u), tolerance = 1e-6)
})

## @cvxpy NONE
test_that("C_problem: constraint Jacobian of sum(x) is all ones", {
  skip_if_not_installed("sparsediff")
  cp <- .de_C_problem(.de_make_prob())
  cp$init_jacobian_coo()
  cp$init_hessian_coo_lower_tri()
  u <- c(0, 0.5, 1)
  cp$objective_forward(u)
  cp$constraint_forward(u)
  expect_equal(as.numeric(cp$jacobian_values()), c(1, 1, 1), tolerance = 1e-9)
})

## @cvxpy NONE
test_that("C_problem: Lagrangian Hessian of sum(exp(x)) is diag(exp(u))", {
  skip_if_not_installed("sparsediff")
  cp <- .de_C_problem(.de_make_prob())
  cp$init_jacobian_coo()
  cp$init_hessian_coo_lower_tri()
  u <- c(0, 0.5, 1)
  cp$objective_forward(u)
  cp$constraint_forward(u)
  ## obj_factor = 1, multiplier 0 on the (linear) constraint -> only obj Hessian
  h <- as.numeric(cp$hessian_values(obj_factor = 1, lagrange = c(0)))
  expect_equal(sort(h), sort(exp(u)), tolerance = 1e-6)
})

# ===================================================================
# Step 5a follow-on: matmul + structural atoms
# ===================================================================

## @cvxpy NONE
test_that("matmul A %*% x: objective and gradient", {
  skip_if_not_installed("sparsediff")
  A <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)  # 2x3, dense
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(A %*% x)), list())
  cp <- .de_C_problem(prob)
  cp$init_jacobian_coo(); cp$init_hessian_coo_lower_tri()
  u <- c(1, 2, 3)
  expect_equal(cp$objective_forward(u), sum(A %*% u), tolerance = 1e-9)
  ## d/dx_j sum(A x) = sum_i A_ij = colSums(A)
  expect_equal(as.numeric(cp$gradient()), as.numeric(colSums(A)), tolerance = 1e-9)
})

## @cvxpy NONE
test_that("quad_form x'Px: objective and gradient", {
  skip_if_not_installed("sparsediff")
  P <- diag(c(2, 3, 4))
  x <- Variable(3)
  prob <- Problem(Minimize(quad_form(x, P)), list())
  cp <- .de_C_problem(prob)
  cp$init_jacobian_coo(); cp$init_hessian_coo_lower_tri()
  u <- c(1, 2, 3)
  expect_equal(cp$objective_forward(u), as.numeric(t(u) %*% P %*% u), tolerance = 1e-9)  # 50
  expect_equal(as.numeric(cp$gradient()), as.numeric(2 * P %*% u), tolerance = 1e-9)       # (4,12,24)
})

## @cvxpy NONE
test_that("index x[1:2]: gradient selects the first two entries", {
  skip_if_not_installed("sparsediff")
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(x[1:2])), list())
  cp <- .de_C_problem(prob)
  cp$init_jacobian_coo(); cp$init_hessian_coo_lower_tri()
  u <- c(5, 7, 9)
  expect_equal(cp$objective_forward(u), u[1] + u[2], tolerance = 1e-9)
  expect_equal(as.numeric(cp$gradient()), c(1, 1, 0), tolerance = 1e-9)
})

## @cvxpy NONE
test_that("structural atoms convert and evaluate (reshape/transpose/hstack)", {
  skip_if_not_installed("sparsediff")
  x <- Variable(3)
  for (obj in list(
        sum_entries(reshape_expr(x, c(1L, 3L))),  # reshape (3,1)->(1,3)
        sum_entries(t(x)),                          # transpose vector
        sum_entries(hstack(x, x))                   # horizontal stack
      )) {
    cp <- .de_C_problem(Problem(Minimize(obj), list()))
    expect_true(is.finite(cp$objective_forward(c(1, 2, 3))))
  }
})

## @cvxpy NONE
test_that("parametrized matmul converts and evaluates via the dense param matmul", {
  skip_if_not_installed("sparsediff")
  P <- Parameter(c(2L, 3L))
  value(P) <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)  # [[1,3,5],[2,4,6]]
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(P %*% x)), list())
  cp <- .de_C_problem(prob)
  u <- c(1, 2, 3)
  ## sum(P %*% x) = (x1 + 3 x2 + 5 x3) + (2 x1 + 4 x2 + 6 x3) = 3 x1 + 7 x2 + 11 x3
  expect_equal(cp$objective_forward(u), 3 * u[1] + 7 * u[2] + 11 * u[3], tolerance = 1e-9)
  expect_equal(as.numeric(cp$gradient()), c(3, 7, 11), tolerance = 1e-9)
})
