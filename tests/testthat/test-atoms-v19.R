## CVXPY 1.9 test_atoms.py parity (2D-portable subset).
##
## The N-D / tuple-axis cases (test_*_nd_axis, test_sum_*tuple_axis,
## test_lambda_sum_largest_nd*) are N/A for CVXR's 2D-only model and are recorded
## in notes/cvxpy_na_tests.txt. This file covers the genuinely 2D ports:
## sum_largest/sum_smallest axis support (PR #3172), and a handful of singletons.
##
## Axis convention: CVXPY axis=0 (reduce rows) == CVXR axis=2; CVXPY axis=1
## (reduce cols) == CVXR axis=1; axis=None == NULL.

## @cvxpy test_atoms.py::TestAtoms::test_sum_largest_axis
test_that("sum_largest/sum_smallest support a 2D axis", {
  set.seed(42)
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  cst <- Constant(X)

  ## --- CVXR axis=2 (CVXPY axis=0): reduce rows -> row vector (1, 4) ---
  expr0 <- sum_largest(cst, 2, axis = 2)
  expect_equal(expr0@shape, c(1L, 4L))
  expected0 <- apply(X, 2, function(col) sum(sort(col, decreasing = TRUE)[1:2]))
  expect_equal(as.numeric(value(expr0)), expected0)

  ## --- CVXR axis=1 (CVXPY axis=1): reduce cols -> column vector (3, 1) ---
  expr1 <- sum_largest(cst, 2, axis = 1)
  expect_equal(expr1@shape, c(3L, 1L))
  expected1 <- apply(X, 1, function(row) sum(sort(row, decreasing = TRUE)[1:2]))
  expect_equal(as.numeric(value(expr1)), expected1)

  ## --- axis = NULL: scalar (same as no axis) ---
  expr_none <- sum_largest(cst, 2, axis = NULL)
  expect_equal(expr_none@shape, c(1L, 1L))
  expect_equal(as.numeric(value(expr_none)), sum(sort(as.numeric(X), decreasing = TRUE)[1:2]))

  ## --- keepdims keeps the (1, 4) shape and value ---
  expr_kd <- sum_largest(cst, 2, axis = 2, keepdims = TRUE)
  expect_equal(expr_kd@shape, c(1L, 4L))
  expect_equal(as.numeric(value(expr_kd)), expected0)

  ## --- fractional k with axis: top1 + 0.5 * top2 per column ---
  expr_frac <- sum_largest(cst, 1.5, axis = 2)
  expect_equal(expr_frac@shape, c(1L, 4L))
  expected_frac <- apply(X, 2, function(col) {
    s <- sort(col, decreasing = TRUE); s[1] + 0.5 * s[2]
  })
  expect_equal(as.numeric(value(expr_frac)), expected_frac)

  ## --- DCP properties ---
  xv <- Variable(c(3, 4))
  atom <- sum_largest(xv, 2, axis = 2)
  expect_true(is_convex(atom))
  expect_false(is_concave(atom))
  expect_true(is_pwl(atom))

  ## --- copy preserves type and get_data ([k, axis, keepdims]) ---
  cp <- tree_copy(atom)
  expect_s7_class(cp, SumLargest)
  expect_identical(get_data(cp), get_data(atom))
  expect_identical(get_data(atom), list(2, 2L, FALSE))

  ## --- solve: recovered atom value matches the per-column top-2 sums ---
  prob <- Problem(Minimize(sum_entries(sum_largest(xv, 2, axis = 2))),
                  list(xv == X))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(sum_largest(xv, 2, axis = 2))), expected0,
               tolerance = 1e-4)

  ## --- sum_smallest with axis: per-row sum of 2 smallest -> (3, 1) ---
  expr_sm <- sum_smallest(cst, 2, axis = 1)
  expect_equal(expr_sm@shape, c(3L, 1L))
  expected_sm <- apply(X, 1, function(row) sum(sort(row)[1:2]))
  expect_equal(as.numeric(value(expr_sm)), expected_sm)
})

## @cvxpy test_atoms.py::TestAtoms::test_lambda_max_value_none_parameter
test_that("lambda_max value is NULL when the parameter has no value", {
  P <- Parameter(c(3, 3), symmetric = TRUE)
  atom <- lambda_max(P)
  expect_null(value(atom))
  value(P) <- diag(3)
  expect_equal(as.numeric(value(atom)), 1.0)
})

## @cvxpy test_atoms.py::TestAtoms::test_length_all_zeros
test_that("length of an all-zero vector is 0", {
  x <- Variable(2)
  atom <- length_expr(x)
  value(x) <- rep(0, 2)
  expect_equal(as.numeric(value(atom)), 0)
})

## @cvxpy test_atoms.py::TestAtoms::test_one_minus_pos_grad
test_that("one_minus_pos per-atom gradient is a list of sparse -I", {
  atom <- one_minus_pos(Variable(2))
  g <- .grad(atom, list(c(0.5, 0.3)))
  expect_type(g, "list")
  expect_length(g, 1L)
  expect_equal(as.matrix(g[[1L]]), -diag(2))
})

## @cvxpy test_atoms.py::TestAtoms::test_sign_shape
test_that("sign preserves the shape of its argument", {
  x <- Variable(2)
  expect_equal(sign(x)@shape, c(2L, 1L))
  A <- Variable(c(2, 2))
  expect_equal(sign(A)@shape, c(2L, 2L))
})
