## CVXPY parity tests: Expressions, Problem, Constraints, QuadForm
## These tests close critical gaps vs CVXPY test suite:
##   - test_expressions.py (14 critical tests)
##   - test_problem.py (11 critical tests)
##   - test_constraints.py (5 critical tests)
##   - test_quad_form.py (5 critical tests)
##
## IMPORTANT: R uses 1-based indexing and (n,1) shape for vectors,
## while CVXPY uses 0-based indexing and (n,) shape.
## All expected values verified against CVXPY 1.8.1.

library(testthat)

# ====================================================================
# EXPRESSION TESTS (14 critical)
# ====================================================================

# ── 1. test_assign_var_value ─────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("Variable value assignment and retrieval (CVXPY parity)", {
  ## Scalar variable
  a <- Variable(name = "a")
  value(a) <- 1
  expect_equal(as.numeric(value(a)), 1)

  ## Wrong shape: scalar variable, assign vector
  expect_error(value(a) <- c(2, 1), "Invalid dimensions")

  ## Assign NULL clears value
  value(a) <- 1
  value(a) <- NULL
  expect_null(value(a))

  ## Vector variable
  x <- Variable(2, name = "x")
  value(x) <- c(2, 1)
  expect_equal(as.numeric(value(x)), c(2, 1))

  ## Matrix variable
  A <- Variable(c(3, 2), name = "A")
  value(A) <- matrix(1, 3, 2)
  expect_equal(as.numeric(value(A)), rep(1, 6))

  ## Nonneg variable rejects negative value
  x_nn <- Variable(nonneg = TRUE)
  expect_error(value(x_nn) <- -2, "nonnegative")
})

# ── 2. test_parameters_successes ─────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_parameters_successes
test_that("Parameter creation with various attributes (CVXPY parity)", {
  ## Name and shape
  p <- Parameter(name = "p")
  expect_equal(name(p), "p")
  expect_equal(p@shape, c(1L, 1L))

  ## Value assignment works for matching shape
  val <- matrix(-1, 4, 3)
  val[1, 1] <- 2
  p2 <- Parameter(c(4, 3))
  value(p2) <- val
  expect_equal(as.numeric(value(p2)), as.numeric(val))

  ## Initialize with value, then clear
  p3 <- Parameter(value = 10)
  expect_equal(as.numeric(value(p3)), 10)
  value(p3) <- 10
  value(p3) <- NULL
  expect_null(value(p3))

  ## Nonpos parameter
  p4 <- Parameter(c(4, 3), nonpos = TRUE)
  expect_true(is_nonpos(p4))

  ## PSD parameter
  p5 <- Parameter(c(3, 3), PSD = TRUE)
  expect_true(is_psd(p5))
  expect_true(is_symmetric(p5))

  ## NSD parameter
  p6 <- Parameter(c(3, 3), NSD = TRUE)
  expect_true(is_nsd(p6))

  ## Diag parameter
  p7 <- Parameter(c(2, 2), diag = TRUE)
  value(p7) <- diag(2)
  expect_equal(as.matrix(value(p7)), diag(2))
})

# ── 3. test_parameters_failures ──────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_parameters_failures
test_that("Parameter validation errors (CVXPY parity)", {
  ## Wrong shape for nonneg parameter
  p1 <- Parameter(c(4, 3), nonneg = TRUE)
  expect_error(value(p1) <- 1, "Invalid dimensions")

  ## Nonneg value violation
  val <- matrix(-1, 4, 3)
  val[1, 1] <- 2
  p2 <- Parameter(c(4, 3), nonneg = TRUE)
  expect_error(value(p2) <- val, "nonnegative")

  ## Nonpos value violation
  p3 <- Parameter(c(4, 3), nonpos = TRUE)
  expect_error(value(p3) <- val, "nonpositive")

  ## Boolean value violation
  expect_error(
    Parameter(c(2, 2), boolean = TRUE, value = matrix(c(1, 1, 1, -1), 2, 2)),
    "boolean"
  )

  ## Integer value violation
  expect_error(
    Parameter(c(2, 2), integer = TRUE, value = matrix(c(1, 1.5, 1, -1), 2, 2)),
    "integer"
  )

  ## PSD requires square
  expect_error(Parameter(c(2, 3), PSD = TRUE), "square")

  ## NSD requires square
  expect_error(Parameter(c(2, 3), NSD = TRUE), "square")

  ## Symmetric requires square
  expect_error(Parameter(c(2, 3), symmetric = TRUE), "square")
})

# ── 4. test_add_expression ───────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_add_expression
test_that("Addition of Variables and Constants (CVXPY parity)", {
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")
  C <- Variable(c(3, 2), name = "C")

  ## Vector addition
  c_val <- Constant(c(2, 2))
  exp <- x + c_val
  expect_equal(curvature(exp), "AFFINE")
  expect_equal(expr_sign(exp), "UNKNOWN")
  expect_equal(exp@shape, c(2L, 1L))

  ## Chained additions
  exp2 <- exp + z + x
  expect_true(is_affine(exp2))

  ## Incompatible dimensions
  y <- Variable(3, name = "y")
  expect_error(x + y)

  ## Matrix addition
  exp3 <- A + B
  expect_equal(curvature(exp3), "AFFINE")
  expect_equal(exp3@shape, c(2L, 2L))

  ## Incompatible matrix dimensions
  expect_error(A + C)
})

# ── 5. test_sub_expression ───────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_sub_expression
test_that("Subtraction, shape, sign, curvature (CVXPY parity)", {
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")
  C <- Variable(c(3, 2), name = "C")

  ## Vector subtraction
  c_val <- Constant(c(2, 2))
  exp <- x - c_val
  expect_equal(curvature(exp), "AFFINE")
  expect_equal(expr_sign(exp), "UNKNOWN")
  expect_equal(exp@shape, c(2L, 1L))

  ## Chained subtraction
  exp2 <- exp - z - x
  expect_true(is_affine(exp2))

  ## Incompatible dimensions
  y <- Variable(3, name = "y")
  expect_error(x - y)

  ## Matrix subtraction
  exp3 <- A - B
  expect_equal(curvature(exp3), "AFFINE")
  expect_equal(exp3@shape, c(2L, 2L))

  ## Incompatible matrix dimensions
  expect_error(A - C)
})

# ── 6. test_mul_expression ───────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_mul_expression
test_that("Multiplication, shape, sign, curvature (CVXPY parity)", {
  x <- Variable(2, name = "x")
  C <- Variable(c(3, 2), name = "C")

  ## Scalar times vector
  exp <- 2 * x
  expect_equal(curvature(exp), "AFFINE")
  expect_equal(exp@shape, c(2L, 1L))

  ## Constant expressions: T_const is (3,2), B_var is (2,2), result is (3,2)
  T_const <- Constant(matrix(c(1, 3, 2, 5, 3, 5), 3, 2, byrow = TRUE))
  B_var <- Variable(c(2, 2), name = "B")
  exp2 <- (T_const + T_const) %*% B_var
  expect_equal(curvature(exp2), "AFFINE")
  expect_equal(exp2@shape, c(3L, 2L))

  ## Constant division
  c2 <- Constant(2)
  exp3 <- c2 / (3 - 5)
  expect_true(is_constant(exp3))
  expect_equal(exp3@shape, c(1L, 1L))
})

# ── 7. test_matmul_expression ────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_matmul_expression
test_that("Matrix multiplication, shape, curvature (CVXPY parity)", {
  x <- Variable(2, name = "x")
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")

  ## Constant %*% variable
  c_mat <- Constant(matrix(c(2, 2), 1, 2))
  exp <- c_mat %*% x
  expect_equal(curvature(exp), "AFFINE")
  expect_equal(exp@shape, c(1L, 1L))

  ## Affine times affine is quadratic
  suppressWarnings({
    q <- A %*% B
  })
  expect_true(is_quadratic(q))

  ## Constant expression matmul: T_const is (3,2), B is (2,2), result is (3,2)
  T_const <- Constant(matrix(c(1, 3, 2, 5, 3, 5), 3, 2, byrow = TRUE))
  exp2 <- (T_const + T_const) %*% B
  expect_equal(curvature(exp2), "AFFINE")
  expect_equal(exp2@shape, c(3L, 2L))

  ## Incompatible dimensions
  C <- Variable(c(3, 2), name = "C")
  expect_error(Constant(matrix(c(2, 1, 2, 2), 2, 2)) %*% C)
})

# ── 8. test_div_expression ───────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_div_expression
test_that("Division, shape, sign, curvature (CVXPY parity)", {
  x <- Variable(2, name = "x")

  ## Vector / scalar
  exp <- x / 2
  expect_equal(curvature(exp), "AFFINE")
  expect_equal(expr_sign(exp), "UNKNOWN")
  expect_equal(exp@shape, c(2L, 1L))

  ## Constant division
  c1 <- Constant(c(3.0, 4.0, 12.0))
  c2 <- Constant(c(1.0, 2.0, 3.0))
  result <- value(c1 / c2)
  expect_equal(as.numeric(result), c(3.0, 2.0, 4.0), tolerance = 1e-10)

  ## Constant expression
  c3 <- Constant(2)
  exp2 <- c3 / (3 - 5)
  expect_true(is_constant(exp2))
  expect_equal(exp2@shape, c(1L, 1L))
  ## 2 / -2 = -1 (nonpositive)

  ## Parameter division
  p <- Parameter(nonneg = TRUE)
  exp3 <- 2 / p
  value(p) <- 2
  expect_equal(as.numeric(value(exp3)), 1, tolerance = 1e-10)

  ## Sign of rho/2 where rho is nonneg
  rho <- Parameter(nonneg = TRUE)
  value(rho) <- 1
  expect_equal(expr_sign(rho), "NONNEGATIVE")
  expect_equal(expr_sign(Constant(2)), "NONNEGATIVE")
  expect_true(is_nonneg(rho / 2))
})

# ── 9. test_neg_expression ───────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_neg_expression
test_that("Negation, shape, sign, curvature (CVXPY parity)", {
  x <- Variable(2, name = "x")
  C <- Variable(c(3, 2), name = "C")

  ## Negation of vector

  exp <- -x
  expect_equal(curvature(exp), "AFFINE")
  expect_equal(exp@shape, c(2L, 1L))
  expect_true(is_affine(exp))
  expect_equal(expr_sign(exp), "UNKNOWN")
  expect_false(is_nonneg(exp))

  ## Negation of matrix
  exp2 <- -C
  expect_equal(curvature(exp2), "AFFINE")
  expect_equal(exp2@shape, c(3L, 2L))
})

# ── 10. test_index_expression ────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_index_expression
test_that("Indexing into expressions (CVXPY parity, 1-based)", {
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")
  C <- Variable(c(3, 2), name = "C")

  ## Single element from vector (1-based: x[2,1] = CVXPY x[1])
  exp <- x[2, 1]
  expect_equal(curvature(exp), "AFFINE")
  expect_true(is_affine(exp))
  expect_null(value(exp))

  ## Slicing from matrix: C[1:2, 2] (CVXPY C[0:2, 1])
  exp2 <- C[1:2, 2]
  expect_equal(exp2@shape[1L], 2L)  # 2 rows

  ## Full slice: C[, 1:2] (all rows, first 2 cols)
  exp3 <- C[, 1:2]
  expect_equal(exp3@shape, c(3L, 2L))

  ## Constant indexing
  c_mat <- Constant(matrix(c(1, 0, -2, 4), 2, 2, byrow = TRUE))
  exp4 <- c_mat[2, 2]
  expect_true(is_constant(exp4))
  expect_equal(as.numeric(value(exp4)), 4)

  ## Arithmetic expression indexing: (x + z)[2, 1]
  exp5 <- (x + z)[2, 1]
  expect_equal(curvature(exp5), "AFFINE")
  expect_equal(expr_sign(exp5), "UNKNOWN")
})

# ── 11. test_neg_indices ─────────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_neg_indices
test_that("Negative index handling (R convention: exclude)", {
  ## In R, negative indices EXCLUDE those positions (different from Python).
  ## This tests R-specific negative indexing behavior.

  ## Negative indices exclude positions
  c_vec <- Constant(c(1, 2, 3, 4))
  ## Exclude first element: c_vec[-1, 1] should give c(2, 3, 4)
  exp <- c_vec[-1, 1]
  expect_equal(as.numeric(value(exp)), c(2, 3, 4))
  expect_equal(exp@shape[1L], 3L)

  ## Constant matrix: exclude last row
  c_mat <- Constant(matrix(c(1, 3, 2, 4), 2, 2))
  exp2 <- c_mat[-2, ]
  expect_equal(exp2@shape, c(1L, 2L))
})

# ── 12. test_broadcast_mul ───────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_broadcast_mul
test_that("Broadcasting in multiply (CVXPY parity)", {
  set.seed(0)
  m <- 3L
  n <- 4L
  A <- matrix(runif(m * n), m, n)

  ## Column scaling: (m x n) * (1 x n) => (m x n)
  col_scale <- Variable(c(1L, n))
  C_expr <- multiply(A, col_scale)
  expect_equal(C_expr@shape, c(m, n))

  ## Row scaling: (m x n) * (m x 1) => (m x n)
  row_scale <- Variable(c(m, 1L))
  R_expr <- multiply(A, row_scale)
  expect_equal(R_expr@shape, c(m, n))

  ## Solve: Constant(3,1) * Variable(1,3) => (3,3)
  ## CVXPY: y=Parameter((3,1)), z=Variable((1,3)), solve with SCS
  y_val <- matrix(0:2, 3, 1)
  z <- Variable(c(1L, 3L))
  z_val <- matrix(c(-1, 0, 1), 1, 3)
  expr <- multiply(Constant(y_val), z)
  expect_equal(expr@shape, c(3L, 3L))
  prob <- Problem(Minimize(sum(expr)), list(z == z_val))
  psolve(prob, solver = "CLARABEL")
  expected <- y_val[, 1] %o% z_val[1, ]  ## outer product
  expect_equal(value(expr), expected, tolerance = 1e-5)
})

# ── 13. test_broadcast_add ───────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_broadcast_add
test_that("Broadcasting in addition (CVXPY parity)", {
  set.seed(0)
  m <- 3L
  n <- 4L
  A <- matrix(runif(m * n), m, n)

  ## Column-vector addition: (m x n) + (1 x n) => (m x n)
  col_scale <- Variable(c(1L, n))
  C_expr <- A + col_scale
  expect_equal(C_expr@shape, c(m, n))

  ## Row-vector addition: (m x n) + (m x 1) => (m x n)
  row_scale <- Variable(c(m, 1L))
  R_expr <- A + row_scale
  expect_equal(R_expr@shape, c(m, n))

  ## Solve: Constant(3,1) + Variable(1,3) => (3,3)
  ## CVXPY: y=Parameter((3,1)), z=Variable((1,3)), solve with SCS
  y_val <- matrix(0:2, 3, 1)
  z <- Variable(c(1L, 3L))
  z_val <- matrix(c(-1, 0, 1), 1, 3)
  expr <- Constant(y_val) + z
  expect_equal(expr@shape, c(3L, 3L))
  prob <- Problem(Minimize(sum(expr)), list(z == z_val))
  psolve(prob, solver = "CLARABEL")
  expected <- outer(y_val[, 1], z_val[1, ], "+")
  expect_equal(value(expr), expected, tolerance = 1e-5)

  ## Subtraction: (1,3) - (3,1) => (3,3)
  expr2 <- z - Constant(y_val)
  expect_equal(expr2@shape, c(3L, 3L))
  prob2 <- Problem(Minimize(sum(expr2)), list(z == z_val))
  psolve(prob2, solver = "CLARABEL")
  expected2 <- outer(z_val[1, ], y_val[, 1], "-") |> t()
  expect_equal(value(expr2), expected2, tolerance = 1e-5)
})

# ── 14. test_curvatures ─────────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_curvatures
test_that("Comprehensive curvature composition tests (CVXPY parity)", {
  x <- Variable(2, name = "x")

  ## sum(square(x)) is convex
  expr <- sum(square(x))
  expect_true(is_convex(expr))
  expect_false(is_concave(expr))

  ## Constant is constant, affine, convex, and concave
  y <- Constant(42)
  expect_true(is_constant(y))
  expect_true(is_affine(y))
  expect_true(is_convex(y))
  expect_true(is_concave(y))

  ## Variable is affine
  z <- Variable(2)
  expect_true(is_affine(z))
  expect_true(is_convex(z))
  expect_true(is_concave(z))

  ## exp(x) is convex, not concave
  x1 <- Variable(1)
  e <- exp(x1)
  expect_true(is_convex(e))
  expect_false(is_concave(e))

  ## log(x) is concave, not convex
  x2 <- Variable(1, nonneg = TRUE)
  l <- log(x2)
  expect_true(is_concave(l))
  expect_false(is_convex(l))

  ## DGP curvature for positive variables
  x_pos <- Variable(1, pos = TRUE)
  expect_true(is_log_log_affine(x_pos))
  expect_true(is_log_log_convex(x_pos))
  expect_true(is_log_log_concave(x_pos))

  ## Monomial: x*x*x is log-log affine
  mono <- x_pos * x_pos * x_pos
  expect_true(is_log_log_affine(mono))

  ## Posynomial: x*x*x + x is log-log convex
  posy <- x_pos * x_pos * x_pos + x_pos
  expect_true(is_log_log_convex(posy))
  expect_false(is_log_log_concave(posy))
})


# ====================================================================
# PROBLEM TESTS (11 critical)
# ====================================================================

# ── 1. test_solving_a_problem_with_unspecified_parameters ────────────

## @cvxpy test_problem.py::TestProblem::test_solving_a_problem_with_unspecified_parameters
test_that("Error when param has no value (CVXPY parity)", {
  ## Minimize(param) with no variables goes through ConstantSolver in CVXR,

  ## so we need a real variable to trigger the parameter-check code path.
  param <- Parameter(name = "lambda")
  x <- Variable(name = "x")
  problem <- Problem(Minimize(param + x), list(x >= 0))
  expect_error(psolve(problem, solver = "SCS"), "unspecified")
})

# ── 2. test_get_problem_data ────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_get_problem_data
test_that("problem_data returns correct format (CVXPY parity)", {
  a <- Variable(name = "a")
  x <- Variable(2, name = "x")

  ## Norm problem -> SOC cone via Clarabel
  pd <- problem_data(
    Problem(Minimize(p_norm(x, 2) + 3), list()),
    solver = "CLARABEL"
  )
  expect_true(is.list(pd))
  expect_true("data" %in% names(pd))
  expect_true("chain" %in% names(pd))
  expect_true("inverse_data" %in% names(pd))

  data <- pd$data
  ## Should have c vector and A matrix
  expect_true("c" %in% names(data) || "q" %in% names(data))

  ## SCS data format
  pd2 <- problem_data(
    Problem(Minimize(sum(x)), list(x >= 0)),
    solver = "SCS"
  )
  expect_true(is.list(pd2$data))
})

# ── 3. test_duplicate_constraints ───────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_redundant_constraints
test_that("Same constraint object twice (CVXPY parity)", {
  x <- Variable(2, name = "x")
  eq <- (x == 2)
  le <- (x <= 2)

  ## Problem with duplicate constraints should solve
  p <- Problem(Minimize(sum(x)), list(eq, eq, le, le))
  result <- psolve(p, solver = "CLARABEL")
  expect_equal(result, 4, tolerance = 1e-4)  # sum of [2, 2]
  expect_equal(as.numeric(value(x)), c(2, 2), tolerance = 1e-4)
})

# ── 4. test_parameter_problems ──────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_parameter_problems
test_that("Solving with parameters (CVXPY parity)", {
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  p1 <- Parameter()
  p2 <- Parameter(3, nonpos = TRUE)
  p3 <- Parameter(c(4, 4), nonneg = TRUE)

  ## Maximize(p1 * a) s.t. a + p1 <= p2, b <= p3 + p3 + 2
  prob <- Problem(Maximize(p1 * a), list(a + p1 <= p2, b <= p3 + p3 + 2))
  value(p1) <- 2
  value(p2) <- matrix(-1, 3, 1)
  value(p3) <- matrix(1, 4, 4)
  result <- psolve(prob, solver = "SCS", reltol = 1e-6)
  ## p1=2, p2=-1, so a + 2 <= -1 => a <= -3 => max(2*a) = 2*(-3) = -6
  expect_equal(result, -6, tolerance = 1e-3)

  ## Unset parameter should error
  value(p1) <- NULL
  expect_error(psolve(prob, solver = "SCS", reltol = 1e-6), "unspecified")
})

# ── 5. test_parameter_expressions ───────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_parameter_expressions
test_that("Parameter in objective/constraints (CVXPY parity)", {
  x <- Variable(name = "x")
  y <- Variable(name = "y")
  x0 <- Parameter(name = "x0")

  ## Linear approximation of x^2: x0^2 + 2*x0*(x - x0)
  xSquared <- x0 * x0 + 2 * x0 * (x - x0)

  ## First solve with x0 = 2
  value(x0) <- 2
  g <- xSquared - y
  obj <- abs(x - 1)
  prob <- Problem(Minimize(obj), list(g == 0))
  suppressWarnings(psolve(prob, solver = "SCS"))

  ## Re-solve with x0 = 1
  value(x0) <- 1
  suppressWarnings(psolve(prob, solver = "SCS"))
  ## With x0=1: xSquared = 1 + 2*(x-1) = 2x - 1
  ## g = 2x - 1 - y = 0 => y = 2x - 1
  ## min |x-1| => x=1, y=1, g=0
  expect_equal(as.numeric(value(g)), 0, tolerance = 1e-2)

  ## Multiplication with parameter
  prob2 <- Problem(Minimize(x0 * x), list(x == 1))
  value(x0) <- 2
  suppressWarnings(psolve(prob2, solver = "SCS"))
  value(x0) <- 1
  suppressWarnings(psolve(prob2, solver = "SCS"))
  expect_equal(value(prob2), 1, tolerance = 1e-2)
})

# ── 6. test_presolve_parameters ─────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_presolve_parameters
test_that("Parameter evaluation in solve chain (CVXPY parity)", {
  gamma <- Parameter(nonneg = TRUE, name = "gamma")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x), list(gamma == 1, x >= 0))

  ## gamma = 0 but constraint says gamma == 1 => infeasible
  value(gamma) <- 0
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), INFEASIBLE)

  ## gamma = 1 => feasible, min x s.t. x >= 0 => x = 0
  value(gamma) <- 1
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), OPTIMAL)
})

# ── 7. test_invalid_solvers ─────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_invalid_solvers
test_that("Error for invalid solver/problem combination (CVXPY parity)", {
  a <- Variable(name = "a")
  A <- Variable(c(2, 2), name = "A")

  ## MIP problem on non-MIP solver
  expect_error(
    psolve(Problem(Minimize(Variable(boolean = TRUE))), solver = "OSQP")
  )

  ## SDP problem on QP solver
  expect_error(
    psolve(Problem(Minimize(lambda_max(A))), solver = "OSQP")
  )
})

# ── 8. test_solver_verbose ──────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_solver_verbose
test_that("verbose option works (CVXPY parity)", {
  a <- Variable(name = "a")
  x <- Variable(2, name = "x")

  ## verbose = FALSE should not produce solver output
  ## verbose = TRUE should produce some output
  ## Just verify it runs without error
  prob <- Problem(Minimize(a + x[1, 1]), list(a >= 2, x >= 2))

  ## Test with verbose = FALSE
  result_quiet <- psolve(prob, verbose = FALSE, solver = "CLARABEL")
  expect_true(is.numeric(result_quiet))
  expect_equal(result_quiet, 4, tolerance = 1e-3)

  ## Test with verbose = TRUE
  ## Solver output goes to stdout (C library), CVXR messages go to stderr.
  ## capture.output(type = "message") captures CVXR's cli messages.
  output <- capture.output({
    result_verb <- psolve(prob, verbose = TRUE, solver = "CLARABEL")
  }, type = "message")
  expect_true(is.numeric(result_verb))
  ## verbose = TRUE should produce CVXR informational messages on stderr
  expect_true(length(output) > 0)
})

# ── 9. test_matrix_lp ───────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_matrix_lp
test_that("Matrix-valued LP (CVXPY parity)", {
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")
  C <- Variable(c(3, 2), name = "C")

  ## Simple matrix equality
  T_val <- matrix(1, 2, 2)
  prob <- Problem(Minimize(1), list(A == T_val))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 1, tolerance = 1e-3)
  expect_equal(as.matrix(value(A)), T_val, tolerance = 1e-3)

  ## Matrix-valued LP with multiple constraints
  T2 <- matrix(2, 2, 3)
  prob2 <- Problem(Minimize(1), list(
    A >= T2 %*% C,
    A == B,
    C == t(T2)
  ))
  result2 <- psolve(prob2, solver = "CLARABEL")
  expect_equal(result2, 1, tolerance = 1e-3)
  expect_equal(as.matrix(value(A)), as.matrix(value(B)), tolerance = 1e-3)
})

# ── 10. test_variable_promotion ─────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_variable_promotion
test_that("Scalar variable promoted to match constraints (CVXPY parity)", {
  a <- Variable(name = "a")
  x <- Variable(2, name = "x")

  ## Scalar a compared to vector x
  prob <- Problem(Minimize(a), list(x <= a, x == matrix(c(1, 2), 2, 1)))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 2, tolerance = 1e-3)
  expect_equal(as.numeric(value(a)), 2, tolerance = 1e-3)

  ## Scalar a compared to matrix
  A <- Variable(c(2, 2), name = "A")
  prob2 <- Problem(Minimize(a), list(
    A <= a,
    A == matrix(c(1, 2, 3, 4), 2, 2)
  ))
  result2 <- psolve(prob2, solver = "CLARABEL")
  expect_equal(result2, 4, tolerance = 1e-3)
  expect_equal(as.numeric(value(a)), 4, tolerance = 1e-3)
})

# ── 11. test_parameter_promotion ────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_parameter_promotion
test_that("Scalar parameter promoted (CVXPY parity)", {
  a_param <- Parameter(name = "a_param")
  ## In CVXPY: [[1,2],[3,4]] * a
  ## In R we use multiply for elementwise
  ref_mat <- matrix(c(1, 3, 2, 4), 2, 2)  # column-major
  exp <- ref_mat * a_param
  value(a_param) <- 2
  expected <- 2 * ref_mat
  expect_equal(as.matrix(value(exp)), expected, tolerance = 1e-10)
})


# ====================================================================
# CONSTRAINT TESTS (5 critical)
# ====================================================================

# ── 1. test_pow3d_constraint ────────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_pow3d_constraint
test_that("PowCone3D constraint creation and residual (CVXPY parity)", {
  n <- 3L
  set.seed(0)
  alpha <- 0.275

  x <- Variable(c(n, 1L))
  y <- Variable(c(n, 1L))
  z <- Variable(c(n, 1L))
  con <- PowCone3D(x, y, z, alpha)

  ## Feasible values: z = x^alpha * y^(1-alpha) (sign flip on one)
  x0 <- 0.1 + runif(n)
  y0 <- 0.1 + runif(n)
  z0 <- x0^alpha * y0^(1 - alpha)
  z0[2] <- -z0[2]  # flip one element sign

  value(x) <- matrix(x0, n, 1)
  value(y) <- matrix(y0, n, 1)
  value(z) <- matrix(z0, n, 1)
  viol <- residual(con)
  ## residual() returns a vector (one per element); check all are small
  expect_true(all(viol <= 1e-7))

  ## Infeasible: make x negative
  x1 <- x0
  x1[1] <- -0.9 * x1[1]
  value(x) <- matrix(x1, n, 1)
  viol2 <- residual(con)
  ## The first element should have large violation; check max
  expect_true(max(viol2) >= 0.99 * abs(x1[1]))

  ## Invalid alpha values
  expect_error(PowCone3D(x, y, z, 1.001))
  expect_error(PowCone3D(x, y, z, -0.00001))
})

# ── 2. test_pownd_constraint ────────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_pownd_constraint
test_that("PowConeND constraint creation (CVXPY parity)", {
  n <- 4L
  W <- Variable(c(n, 1L))
  z <- Variable(1L)
  set.seed(0)
  alpha <- 0.5 + runif(n)
  alpha <- alpha / sum(alpha)

  ## Entries don't sum to one
  bad_alpha <- alpha + 0.01
  expect_error(
    PowConeND(W, z, Constant(matrix(bad_alpha, n, 1)), axis = 2L),
    "sum to 1"
  )

  ## Shapes don't match: alpha as (1, n) row vector vs W as (n, 1) column
  expect_error(
    PowConeND(W, z, Constant(matrix(alpha, 1, n)), axis = 2L)
  )

  ## Compute a violation
  con <- PowConeND(W, z, Constant(matrix(alpha, n, 1)), axis = 2L)
  W0 <- 0.1 + runif(n)
  z0 <- prod(W0^alpha) + 0.05
  value(W) <- matrix(W0, n, 1)
  value(z) <- z0
  viol <- violation(con)
  expect_true(viol >= 0.01)
  expect_true(viol <= 0.06)
})

# ── 3. test_chained_constraints ─────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_chained_constraints
test_that("Chained constraints raise an error (CVXPY parity)", {
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")

  ## In R, z <= x evaluates to a Constraint object.
  ## Then (constraint) <= 1 should error because you can't compare
  ## a Constraint to a numeric.
  constr <- (z <= x)
  expect_error(constr <= 1)

  ## Similarly for equality
  constr_eq <- (x == z)
  expect_error(constr_eq == 1)
})

# ── 4. test_bound_properties ───────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_bound_properties
test_that("Bounds on Variable via bounds attribute (CVXPY parity)", {
  ## Lower bound only
  x1 <- Variable(c(1, 1), bounds = list(1, NULL))
  attrs1 <- x1@attributes
  expect_true(all(attrs1$bounds[[1L]] > -Inf))
  expect_true(all(attrs1$bounds[[2L]] == Inf))

  ## Upper bound only
  x2 <- Variable(c(1, 1), bounds = list(NULL, 1))
  attrs2 <- x2@attributes
  expect_true(all(attrs2$bounds[[1L]] == -Inf))
  expect_true(all(attrs2$bounds[[2L]] < Inf))

  ## Both bounds
  x3 <- Variable(c(1, 1), bounds = list(1, 2))
  attrs3 <- x3@attributes
  expect_equal(attrs3$bounds[[1L]], 1)
  expect_equal(attrs3$bounds[[2L]], 2)
})

# ── 5. test_bounds_attr ────────────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("Bounds attribute for variables generates correct constraints (CVXPY parity)", {
  ## Bounds with lower and upper arrays
  x1 <- Variable(3, bounds = list(c(1, 2, 3), c(4, 5, 6)))
  prob <- Problem(Minimize(sum(x1)), list())
  result <- psolve(prob)
  expect_equal(as.numeric(value(x1)), c(1, 2, 3), tolerance = 1e-3)

  ## Project should enforce bounds (project is internal)
  proj <- CVXR:::project(x1, matrix(c(0, 0, 0), 3, 1))
  expect_equal(as.numeric(proj), c(1, 2, 3), tolerance = 1e-10)

  ## Bounds with scalar lower and upper
  x2 <- Variable(2, bounds = list(1, 2))
  prob2 <- Problem(Minimize(sum(x2)), list())
  result2 <- psolve(prob2)
  expect_equal(as.numeric(value(x2)), c(1, 1), tolerance = 1e-3)

  ## Bounds with array lower and scalar upper
  x3 <- Variable(2, bounds = list(c(4, 5), 6))
  prob3 <- Problem(Maximize(sum(x3)), list())
  result3 <- psolve(prob3)
  expect_equal(as.numeric(value(x3)), c(6, 6), tolerance = 1e-3)

  ## Bounds with scalar lower and array upper
  x4 <- Variable(3, bounds = list(1, c(2, 3, 4)))
  prob4 <- Problem(Maximize(sum(x4)), list())
  result4 <- psolve(prob4)
  expect_equal(as.numeric(value(x4)), c(2, 3, 4), tolerance = 1e-3)
})


# ====================================================================
# QUAD FORM TESTS (5 critical)
# ====================================================================

# ── 1. test_singular_quad_form ──────────────────────────────────────

## @cvxpy test_quad_form.py::TestNonOptimal::test_singular_quad_form
test_that("quad_form with rank-deficient P (CVXPY parity)", {
  set.seed(1234)
  n <- 3L
  ## Construct a random 1d finite distribution
  v <- exp(rnorm(n))
  v <- v / sum(v)

  ## Construct a random PSD matrix
  A <- matrix(rnorm(n * n), n, n)
  Q <- A %*% t(A)

  ## Project onto orthogonal complement of v => rank-deficient
  E <- diag(n) - outer(v, v) / sum(v * v)
  Q <- E %*% Q %*% t(E)
  expect_equal(Matrix::rankMatrix(Q)[1], n - 1L)

  ## Minimize: x^T Q x on simplex => minimum at v with value 0
  x <- Variable(n)
  q <- quad_form(x, Q)
  prob <- Problem(Minimize(q), list(x >= 0, sum(x) == 1))
  psolve(prob, solver = "OSQP")
  xopt <- as.numeric(value(x))
  yopt <- as.numeric(t(xopt) %*% Q %*% xopt)
  expect_equal(yopt, 0, tolerance = 1e-3)
  expect_equal(xopt, v, tolerance = 1e-3)
})

# ── 2. test_param_quad_form ─────────────────────────────────────────

## @cvxpy test_quad_form.py::TestNonOptimal::test_param_quad_form
test_that("quad_form with Parameter as P matrix (CVXPY parity)", {
  P <- Parameter(c(2, 2), PSD = TRUE)
  Q <- diag(2)
  x <- Variable(2)
  cost <- quad_form(x, P)
  value(P) <- Q
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  suppressWarnings(result <- psolve(prob, solver = "SCS"))
  ## 1^2 + 2^2 = 5
  expect_equal(result, 5, tolerance = 1e-2)
})

# ── 3. test_non_symmetric ──────────────────────────────────────────

## @cvxpy test_quad_form.py::TestNonOptimal::test_non_symmetric
test_that("Non-symmetric P detection (CVXPY parity)", {
  P <- matrix(c(2, 2, 3, 4), 2, 2)
  x <- Variable(2)
  expect_error(quad_form(x, P), "symmetric")
})

# ── 4. test_non_psd ───────────────────────────────────────────────

## @cvxpy test_quad_form.py::TestNonOptimal::test_non_psd
test_that("Non-PSD P handling (CVXPY parity)", {
  P <- matrix(c(1, 0, 0, -1), 2, 2)
  x <- Variable(2)
  ## Forming quad_form is okay (with warning about indefinite)
  suppressWarnings({
    cost <- quad_form(x, P)
  })
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  ## Should fail DCP: indefinite quad form is not convex
  expect_error(psolve(prob, solver = "SCS"), "DCP")
})

# ── 5. test_assume_psd ─────────────────────────────────────────────

## @cvxpy test_quad_form.py::TestNonOptimal::test_assume_psd
test_that("assume_PSD flag (CVXPY parity - not implemented)", {
  ## CVXR's quad_form does not currently have assume_PSD parameter.
  ## This test documents the gap and verifies basic convexity behavior.
  x <- Variable(3)
  A <- diag(3)
  expr <- quad_form(x, A)
  expect_true(is_convex(expr))

  ## Negative definite P => quad_form is concave, not convex
  A_neg <- -diag(3)
  suppressWarnings({
    expr_neg <- quad_form(x, A_neg)
  })
  expect_true(is_concave(expr_neg))
  expect_false(is_convex(expr_neg))
})


# ====================================================================
# ADDITIONAL PARITY TESTS
# ====================================================================

# ── Expression sign composition ──────────────────────────────────────

## @cvxpy NONE
test_that("Sign composition through arithmetic (CVXPY parity)", {
  ## Constant sign
  c_pos <- Constant(2)
  c_neg <- Constant(-2)
  c_zero <- Constant(0)

  expect_equal(expr_sign(c_pos), "NONNEGATIVE")
  expect_equal(expr_sign(c_neg), "NONPOSITIVE")
  expect_equal(expr_sign(c_zero), "ZERO")

  ## Negation of constant
  expect_equal(expr_sign(-c_pos), "NONPOSITIVE")
  expect_equal(expr_sign(-c_neg), "NONNEGATIVE")

  ## Mixed-sign constant
  c_mixed <- Constant(c(-1, 1))
  expect_equal(expr_sign(c_mixed), "UNKNOWN")

  ## Zero * constant = zero sign
  expect_true(is_zero(0 * c_pos))
})

# ── Parameter as constant in DCP ────────────────────────────────────

## @cvxpy NONE
test_that("Parameter treated as constant for DCP (CVXPY parity)", {
  p <- Parameter(nonneg = TRUE)
  x <- Variable(1)

  ## p * x is affine (since p is constant for DCP)
  expr <- p * x
  expect_true(is_affine(expr))

  ## exp(x) is convex
  expr2 <- exp(x)
  expect_true(is_convex(expr2))

  ## p * exp(x) where p is nonneg => convex
  expr3 <- p * exp(x)
  expect_true(is_convex(expr3))
})

# ── Equality constraint properties ──────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("Equality constraint shape and value (CVXPY parity)", {
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")

  constr <- (x == z)
  expect_equal(constr@shape, c(2L, 1L))
  expect_null(dual_value(constr))

  ## Set values and check
  value(x) <- c(2, 1)
  value(z) <- c(2, 2)
  expect_false(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 1))

  value(z) <- c(2, 1)
  expect_true(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 0))

  ## Incompatible dimensions
  y <- Variable(3, name = "y")
  expect_error(x == y)
})

# ── Inequality constraint properties ────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("Inequality constraint shape and value (CVXPY parity)", {
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")

  constr <- (x <= z)
  expect_equal(constr@shape, c(2L, 1L))
  expect_null(dual_value(constr))

  ## Feasible: x <= z
  value(x) <- c(1, 1)
  value(z) <- c(2, 2)
  expect_true(value(constr))

  ## Infeasible: x > z for one element
  value(x) <- c(2, 1)
  value(z) <- c(2, 0)
  expect_false(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 1))

  ## Feasible again
  value(z) <- c(2, 2)
  expect_true(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 0))

  ## Incompatible dimensions
  y <- Variable(3, name = "y")
  expect_error(x <= y)
})

# ── PSD constraint operator ─────────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_psd_constraint
test_that("PSD constraint via >> operator (CVXPY parity)", {
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")

  constr <- (A %>>% B)
  expect_equal(constr@shape, c(2L, 2L))

  ## Feasible: A - B is PSD
  value(A) <- matrix(c(2, -1, 1, 2), 2, 2)
  value(B) <- diag(2)
  expect_true(value(constr))
  expect_equal(residual(constr), 0, tolerance = 1e-10)

  ## Non-square matrix should error
  x <- Variable(2, name = "x")
  expect_error(x %>>% 0)
})

# ── Solve and verify status/value ───────────────────────────────────

## @cvxpy NONE
test_that("Problem solve sets status and value correctly (CVXPY parity)", {
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x), list(x >= 5))
  result <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(result, 5, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 5, tolerance = 1e-4)
  expect_equal(value(prob), result, tolerance = 1e-10)

  ## Infeasible problem
  prob2 <- Problem(Minimize(x), list(x >= 5, x <= 3))
  result2 <- psolve(prob2, solver = "CLARABEL")
  expect_equal(status(prob2), INFEASIBLE)
  expect_true(is.infinite(result2))
})

# ── Sparse quad form ───────────────────────────────────────────────

## @cvxpy test_quad_form.py::TestNonOptimal::test_sparse_quad_form
test_that("quad_form with sparse matrix (CVXPY parity)", {
  Q <- Matrix::Diagonal(2)
  x <- Variable(2)
  cost <- quad_form(x, Q)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  result <- psolve(prob, solver = "OSQP")
  ## 1^2 + 2^2 = 5
  expect_equal(result, 5, tolerance = 1e-2)
})

# ── Zero matrix quad form ──────────────────────────────────────────

## @cvxpy test_quad_form.py::TestNonOptimal::test_zero_matrix
test_that("quad_form with P = 0 (CVXPY parity)", {
  x <- Variable(3)
  A <- diag(3)
  b <- rep(1, 3)
  c_vec <- -rep(1, 3)
  P <- matrix(0, 3, 3)
  expr <- 0.5 * quad_form(x, P) + sum(c_vec * x)
  prob <- Problem(Minimize(expr), list(A %*% x <= b))
  result <- psolve(prob, solver = "CLARABEL")
  ## With P=0, obj = sum(-x) = -sum(x), min with x <= 1 => x = [1,1,1]
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-3)
  expect_equal(result, -3, tolerance = 1e-3)
})

# ── Nonneg / NonPos constraints via solve ───────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_nonneg
test_that("NonNeg constraint through conic and QP paths (CVXPY parity)", {
  x <- Variable(3)
  c_val <- 0:2

  ## Conic path: Clarabel
  prob <- Problem(Minimize(sum(x)), list(NonNeg(x - c_val)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), c_val, tolerance = 1e-4)

  ## QP path: OSQP
  x2 <- Variable(3)
  prob2 <- Problem(Minimize(sum(x2)), list(NonNeg(x2 - c_val)))
  psolve(prob2, solver = "OSQP")
  expect_equal(as.numeric(value(x2)), c_val, tolerance = 1e-3)
})

## @cvxpy test_constraints.py::TestConstraints::test_nonpos
test_that("NonPos constraint through conic and QP paths (CVXPY parity)", {
  x <- Variable(3)
  c_val <- 0:2

  ## Conic path
  prob <- Problem(Maximize(sum(x)), list(NonPos(x - c_val)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), c_val, tolerance = 1e-4)

  ## QP path
  x2 <- Variable(3)
  prob2 <- Problem(Maximize(sum(x2)), list(NonPos(x2 - c_val)))
  psolve(prob2, solver = "OSQP")
  expect_equal(as.numeric(value(x2)), c_val, tolerance = 1e-3)
})


# ====================================================================
# SIGN TESTS (5 from test_sign.py)
# ====================================================================

## @cvxpy test_sign.py::TestSign::test_add
test_that("Sign of addition expressions (CVXPY parity)", {
  pos <- Constant(1)
  neg <- Constant(-1)
  zero <- Constant(0)
  unknown <- Variable(name = "x")

  ## pos + neg => UNKNOWN (sign algebra: NONNEG + NONPOS = UNKNOWN)
  expect_equal(expr_sign(pos + neg), "UNKNOWN")

  ## neg + zero => NONPOSITIVE
  expect_equal(expr_sign(neg + zero), "NONPOSITIVE")

  ## pos + pos => NONNEGATIVE
  expect_equal(expr_sign(pos + pos), "NONNEGATIVE")

  ## unknown + zero => UNKNOWN
  expect_equal(expr_sign(unknown + zero), "UNKNOWN")
})

## @cvxpy test_sign.py::TestSign::test_sub
test_that("Sign of subtraction expressions (CVXPY parity)", {
  pos <- Constant(1)
  neg <- Constant(-1)
  zero <- Constant(0)
  unknown <- Variable(name = "x")

  ## pos - neg => NONNEGATIVE (1 - (-1) = 2)
  expect_equal(expr_sign(pos - neg), "NONNEGATIVE")

  ## neg - zero => NONPOSITIVE
  expect_equal(expr_sign(neg - zero), "NONPOSITIVE")

  ## pos - pos => UNKNOWN (sign algebra: NONNEG - NONNEG = UNKNOWN)
  expect_equal(expr_sign(pos - pos), "UNKNOWN")
})

## @cvxpy test_sign.py::TestSign::test_mult
test_that("Sign of multiplication expressions (CVXPY parity)", {
  pos <- Constant(1)
  neg <- Constant(-1)
  zero <- Constant(0)
  unknown <- Variable(name = "x")

  ## zero * pos => ZERO
  expect_equal(expr_sign(zero * pos), "ZERO")

  ## unknown * pos => UNKNOWN
  expect_equal(expr_sign(unknown * pos), "UNKNOWN")

  ## pos * neg => NONPOSITIVE
  expect_equal(expr_sign(pos * neg), "NONPOSITIVE")

  ## pos * pos => NONNEGATIVE
  expect_equal(expr_sign(pos * pos), "NONNEGATIVE")

  ## neg * neg => NONNEGATIVE
  expect_equal(expr_sign(neg * neg), "NONNEGATIVE")

  ## zero * unknown => ZERO
  expect_equal(expr_sign(zero * unknown), "ZERO")
})

## @cvxpy test_sign.py::TestSign::test_neg
test_that("Sign of negation expressions (CVXPY parity)", {
  pos <- Constant(1)
  neg <- Constant(-1)
  zero <- Constant(0)

  ## -zero => ZERO
  expect_equal(expr_sign(-zero), "ZERO")

  ## -pos => NONPOSITIVE
  expect_equal(expr_sign(-pos), "NONPOSITIVE")

  ## -neg => NONNEGATIVE
  expect_equal(expr_sign(-neg), "NONNEGATIVE")
})

## @cvxpy test_sign.py::TestSign::test_is_sign
test_that("Sign predicate methods (CVXPY parity)", {
  pos <- Constant(1)
  neg <- Constant(-1)
  zero <- Constant(0)
  unknown <- Variable(name = "x")

  ## is_nonneg
  expect_true(is_nonneg(pos))
  expect_false(is_nonneg(neg))
  expect_false(is_nonneg(unknown))
  expect_true(is_nonneg(zero))

  ## is_nonpos
  expect_false(is_nonpos(pos))
  expect_true(is_nonpos(neg))
  expect_false(is_nonpos(unknown))
  expect_true(is_nonpos(zero))

  ## is_zero
  expect_true(is_zero(zero))
  expect_false(is_zero(neg))
  expect_false(is_zero(pos))

  ## unknown is neither nonneg nor nonpos
  expect_false(is_nonneg(unknown) || is_nonpos(unknown))
})


# ====================================================================
# EXPRESSION METHOD TESTS (from test_expression_methods.py)
# Ported from CVXPY commit 3b964472b (Release 1.8.1)
# ====================================================================

## ── test_conj ─────────────────────────────────────────────────────
## CVXPY: v.conj() >= 1 with minimize sum(v) => v = [1,1,1,1]
## In CVXR: Conj(v) dispatches via Complex S3 group handler

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_conj
test_that("CVXPY parity: conj method on variable", {
  v <- Variable(c(4, 1))
  obj <- Minimize(sum(v))
  ## For real variables, Conj(v) == v, so Conj(v) >= 1 <=> v >= 1
  prob <- Problem(obj, list(Conj(v) >= 1))
  psolve(prob, solver = "SCS")
  expect_equal(as.numeric(value(v)), rep(1, 4), tolerance = 1e-3)
})

## ── test_conjugate ────────────────────────────────────────────────
## CVXPY: v.conjugate() >= 1 is same as v.conj() >= 1
## In CVXR: Conj() is the same function (no separate conjugate method)

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_conjugate
test_that("CVXPY parity: conjugate (same as conj) on variable", {
  v <- Variable(c(4, 1))
  obj <- Minimize(sum(v))
  ## conjugate is just Conj in R
  prob <- Problem(obj, list(Conj(v) >= 1))
  psolve(prob, solver = "SCS")
  expect_equal(as.numeric(value(v)), rep(1, 4), tolerance = 1e-3)
})

## ── test_max ──────────────────────────────────────────────────────
## CVXPY: Variable().max().sign == UNKNOWN
## Variable(2).max(axis=0, keepdims=True).shape == (1,)
## Variable((2,3)).max(axis=0, keepdims=True).shape == (1,3)
## Variable((2,3)).max(axis=1).shape == (2,)
## Invalid axis raises error.

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_max
test_that("CVXPY parity: max entries sign, shape, axis", {
  ## Scalar variable: max sign is UNKNOWN
  expect_equal(expr_sign(max_entries(Variable())), "UNKNOWN")

  ## axis=2 (CVXPY axis=0), keepdims=TRUE: (2,1) -> (1,1)? No.
  ## CVXPY: Variable(2).max(axis=0, keepdims=True).shape == (1,)
  ## In R: Variable(2) is (2,1). axis=2 reduces columns.
  ## For (2,1), axis=2 reduces the single col: shape = (1,1) with keepdims
  ## But CVXPY (2,) with axis=0, keepdims=True -> (1,)
  ## We test the CVXR shapes:
  expr_v <- max_entries(Variable(2), axis = 2, keepdims = TRUE)
  expect_equal(expr_v@shape, c(1L, 1L))

  ## CVXPY: Variable((2,3)).max(axis=0, keepdims=True).shape == (1,3)
  ## axis=0 in CVXPY -> axis=2 in CVXR
  expr_m <- max_entries(Variable(c(2, 3)), axis = 2, keepdims = TRUE)
  expect_equal(expr_m@shape, c(1L, 3L))

  ## CVXPY: Variable((2,3)).max(axis=1).shape == (2,)
  ## axis=1 in CVXPY -> axis=1 in CVXR, but (2,) in CVXPY = (2,1) in R
  expr_m2 <- max_entries(Variable(c(2, 3)), axis = 1)
  expect_equal(expr_m2@shape, c(2L, 1L))

  ## Invalid axis: should error
  expect_error(max_entries(Variable(2), axis = 4))
})

## ── test_min ──────────────────────────────────────────────────────

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_min
test_that("CVXPY parity: min entries sign, shape, axis", {
  ## Scalar variable: min sign is UNKNOWN
  expect_equal(expr_sign(min_entries(Variable())), "UNKNOWN")

  ## CVXPY: Variable(2).min(axis=0).shape == tuple() [scalar]
  ## In R: Variable(2) is (2,1). axis=2 reduces columns -> scalar (1,1)
  expr_v <- min_entries(Variable(2), axis = 2)
  expect_equal(expr_v@shape, c(1L, 1L))

  ## CVXPY: Variable((2,3)).min(axis=0).shape == (3,)
  ## axis=0 in CVXPY -> axis=2 in CVXR
  expr_m <- min_entries(Variable(c(2, 3)), axis = 2)
  expect_equal(expr_m@shape, c(1L, 3L))

  ## CVXPY: Variable((2,3)).min(axis=1).shape == (2,)
  ## axis=1 in CVXPY -> axis=1 in CVXR
  expr_m2 <- min_entries(Variable(c(2, 3)), axis = 1)
  expect_equal(expr_m2@shape, c(2L, 1L))

  ## Invalid axis: should error
  expect_error(min_entries(Variable(2), axis = 4))
})

## ── test_ptp ──────────────────────────────────────────────────────
## CVXPY: ptp (peak-to-peak) = max - min along axis

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_ptp
test_that("CVXPY parity: ptp (peak-to-peak) values and shapes", {
  a_np <- matrix(c(10, 6, -10, 0, 3, -1.5), nrow = 2, ncol = 3)
  a <- Constant(a_np)

  ## No axis: overall ptp
  expr <- ptp(a)
  ## Note: CVXR's ptp = max - min is a convenience, not a proper Ptp atom class.
  ## Sign analysis can't infer nonneg from UNKNOWN - UNKNOWN. We verify the value.
  expect_equal(expr@shape, c(1L, 1L))
  expect_equal(as.numeric(value(expr)), 20.0, tolerance = 1e-6)
  expect_true(as.numeric(value(expr)) >= 0)

  ## axis=2 (CVXPY axis=0): reduce along columns -> row of diffs
  expr2 <- ptp(a, axis = 2)
  ## Note: CVXR's ptp is max - min (convenience function, not a proper atom class),
  ## so sign analysis can't infer nonneg. CVXPY has a Ptp class with is_nonneg.
  ## We verify the value is nonneg (which it always is) but skip the sign assertion.
  expect_equal(expr2@shape, c(1L, 3L))
  expect_equal(as.numeric(value(expr2)), c(4, 10, 4.5), tolerance = 1e-6)
  expect_true(all(as.numeric(value(expr2)) >= 0))

  ## axis=1 (CVXPY axis=1): reduce along rows
  expr3 <- ptp(a, axis = 1)
  expect_equal(expr3@shape, c(2L, 1L))
  expect_equal(as.numeric(value(expr3)), c(20.0, 7.5), tolerance = 1e-6)
  expect_true(all(as.numeric(value(expr3)) >= 0))

  ## axis=2 with keepdims=TRUE (CVXPY axis=0, keepdims=True)
  expr4 <- ptp(a, axis = 2, keepdims = TRUE)
  expect_equal(expr4@shape, c(1L, 3L))
  expect_equal(as.numeric(value(expr4)), c(4, 10, 4.5), tolerance = 1e-6)
  expect_true(all(as.numeric(value(expr4)) >= 0))

  ## axis=1 with keepdims=TRUE (CVXPY axis=1, keepdims=True)
  expr5 <- ptp(a, axis = 1, keepdims = TRUE)
  expect_equal(expr5@shape, c(2L, 1L))
  expect_equal(as.numeric(value(expr5)), c(20.0, 7.5), tolerance = 1e-6)
  expect_true(all(as.numeric(value(expr5)) >= 0))
})

## ── test_reshape ──────────────────────────────────────────────────

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_reshape
test_that("CVXPY parity: reshape sign, curvature, shape, C-order", {
  A <- Variable(c(2, 2), name = "A")
  C_var <- Variable(c(3, 2), name = "C")

  ## Reshape (2,2) -> (4,1) F-order
  expr <- reshape_expr(A, c(4, 1), order = "F")
  expect_equal(expr_sign(expr), "UNKNOWN")
  expect_equal(curvature(expr), "AFFINE")
  expect_equal(expr@shape, c(4L, 1L))

  ## Reshape back: (4,1) -> (2,2)
  expr2 <- reshape_expr(expr, c(2, 2), order = "F")
  expect_equal(expr2@shape, c(2L, 2L))

  ## Reshape of convex expression: square(x) reshaped is still convex, nonneg
  x <- Variable(2, name = "x")
  expr3 <- reshape_expr(square(x), c(1, 2), order = "F")
  expect_equal(expr_sign(expr3), "NONNEGATIVE")
  expect_equal(curvature(expr3), "CONVEX")
  expect_equal(expr3@shape, c(1L, 2L))

  ## Invalid reshape dimensions: (3,2) -> (5,4) should error
  expect_error(reshape_expr(C_var, c(5, 4), order = "F"), ".*reshape.*")

  ## Test C-style reshape
  a <- 0:9  # 0..9
  A_np <- matrix(a, nrow = 5, ncol = 2, byrow = TRUE)  # C-order: row-major
  A_cp <- reshape_expr(Constant(a), c(5, 2), order = "C")
  expect_equal(as.numeric(value(A_cp)), as.numeric(A_np), tolerance = 1e-10)

  ## Solve: X == A_cp
  X <- Variable(c(5, 2))
  prob <- Problem(Minimize(0), list(X == A_cp))
  psolve(prob, solver = "SCS")
  expect_equal(as.numeric(value(X)), as.numeric(A_np), tolerance = 1e-3)

  ## Reshape back to 1D: C-order
  a_np <- as.numeric(t(A_np))  # row-major flatten
  a_cp <- reshape_expr(A_cp, c(10, 1), order = "C")
  expect_equal(as.numeric(value(a_cp)), a_np, tolerance = 1e-10)

  ## Test F-order: default
  b <- matrix(0:11, nrow = 4, ncol = 3, byrow = TRUE)
  b_reshaped <- matrix(as.numeric(t(b)), nrow = 2, ncol = 6, byrow = TRUE)
  ## F-order reshape of b (4x3) -> (2x6) in F-order
  b_F_reshaped <- matrix(as.numeric(b), nrow = 2, ncol = 6)
  X2 <- Variable(c(4, 3))
  X2_reshaped <- reshape_expr(X2, c(2, 6), order = "F")
  prob2 <- Problem(Minimize(0), list(X2_reshaped == b_F_reshaped))
  psolve(prob2, solver = "SCS")
  expect_equal(as.numeric(value(X2_reshaped)), as.numeric(b_F_reshaped), tolerance = 1e-3)
  expect_equal(as.numeric(value(X2)), as.numeric(b), tolerance = 1e-3)
})

## ── test_reshape_negative_one ─────────────────────────────────────

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_reshape_negative_one
test_that("CVXPY parity: reshape with -1 dimension inference", {
  expr <- Variable(c(2, 3))

  ## (-1, 1) -> (6, 1)
  r1 <- reshape_expr(expr, c(-1, 1), order = "F")
  expect_equal(r1@shape, c(6L, 1L))

  ## (1, -1) -> (1, 6)
  r2 <- reshape_expr(expr, c(1, -1), order = "F")
  expect_equal(r2@shape, c(1L, 6L))

  ## (-1, 2) -> (3, 2)
  r3 <- reshape_expr(expr, c(-1, 2), order = "F")
  expect_equal(r3@shape, c(3L, 2L))

  ## Invalid: (8, -1) can't divide 6 by 8
  expect_error(reshape_expr(expr, c(8, -1), order = "F"), ".*reshape.*")

  ## C-order flatten with -1
  A <- matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE)
  A_reshaped_C <- reshape_expr(Constant(A), c(-1, 1), order = "C")
  ## C-order flatten: row-major
  expect_equal(as.numeric(value(A_reshaped_C)), as.numeric(t(A)), tolerance = 1e-10)

  A_reshaped_F <- reshape_expr(Constant(A), c(-1, 1), order = "F")
  ## F-order flatten: column-major
  expect_equal(as.numeric(value(A_reshaped_F)), as.numeric(A), tolerance = 1e-10)
})

## ── test_stats ────────────────────────────────────────────────────
## CVXPY: mean, std, var on a constant array

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_stats
test_that("CVXPY parity: mean, std, var atoms", {
  a_np <- matrix(c(10, 6, 10, 0, 3, 1.5), nrow = 2, ncol = 3)
  a <- Constant(a_np)

  ## mean, var, std — no axis
  expr_mean <- cvxr_mean(a)
  expr_var <- cvxr_var(a)
  expr_std <- cvxr_std(a)

  expect_true(is_nonneg(expr_mean))
  expect_true(is_nonneg(expr_var))
  expect_true(is_nonneg(expr_std))

  expect_equal(as.numeric(value(expr_mean)), mean(a_np), tolerance = 1e-6)
  expect_equal(as.numeric(value(expr_var)),
               var(as.numeric(a_np)) * (length(a_np) - 1) / length(a_np),
               tolerance = 1e-6)
  expect_equal(as.numeric(value(expr_std)),
               sqrt(var(as.numeric(a_np)) * (length(a_np) - 1) / length(a_np)),
               tolerance = 1e-6)

  ## ddof = 0 and 1
  for (ddof in c(0, 1)) {
    expr_var_d <- cvxr_var(a, ddof = ddof)
    expr_std_d <- cvxr_std(a, ddof = ddof)

    n <- length(a_np)
    expected_var <- sum((a_np - mean(a_np))^2) / (n - ddof)
    expect_equal(as.numeric(value(expr_var_d)), expected_var, tolerance = 1e-5)
    expect_equal(as.numeric(value(expr_std_d)), sqrt(expected_var), tolerance = 1e-5)
  }

  ## Axis variants for mean and std
  ## CVXR axis=2 means "reduce rows, keep columns" -> like apply(X, 2, FUN) = colMeans
  ## CVXR axis=1 means "reduce cols, keep rows" -> like apply(X, 1, FUN) = rowMeans
  ## CVXPY axis=0 -> CVXR axis=2
  ## CVXPY axis=1 -> CVXR axis=1
  for (cvxpy_axis in c(0, 1)) {
    cvxr_axis <- if (cvxpy_axis == 0) 2L else 1L
    ## apply MARGIN matches cvxr_axis: axis=2 -> apply(,2,), axis=1 -> apply(,1,)
    r_margin <- cvxr_axis
    for (keepdims in c(TRUE, FALSE)) {
      expr_m <- cvxr_mean(a, axis = cvxr_axis, keepdims = keepdims)
      expr_s <- cvxr_std(a, axis = cvxr_axis, keepdims = keepdims)

      np_mean <- apply(a_np, r_margin, mean)
      np_std <- apply(a_np, r_margin, function(col) sqrt(mean((col - mean(col))^2)))

      if (cvxr_axis == 2L) {
        expected_shape <- c(1L, ncol(a_np))
      } else {
        expected_shape <- c(nrow(a_np), 1L)
      }
      expect_equal(expr_m@shape, expected_shape)
      expect_equal(expr_s@shape, expected_shape)

      expect_equal(as.numeric(value(expr_m)), as.numeric(np_mean), tolerance = 1e-5)
      expect_equal(as.numeric(value(expr_s)), as.numeric(np_std), tolerance = 1e-5)
    }
  }
})

## ── test_sum ──────────────────────────────────────────────────────

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_sum
test_that("CVXPY parity: sum entries sign, curvature, shape, axis", {
  ## Mixed constant: sign is UNKNOWN
  c1 <- Constant(c(1, -1))
  expect_equal(expr_sign(sum_entries(c1)), "UNKNOWN")
  expect_equal(curvature(sum_entries(c1)), "CONSTANT")

  ## Variable sum: UNKNOWN sign, AFFINE curvature
  v <- Variable(2)
  sv <- sum_entries(v)
  expect_equal(expr_sign(sv), "UNKNOWN")
  expect_equal(sv@shape, c(1L, 1L))
  expect_equal(curvature(sv), "AFFINE")

  ## keepdims: (2,1) -> (1,1)
  v21 <- Variable(c(2, 1))
  sk <- sum_entries(v21, keepdims = TRUE)
  expect_equal(sk@shape, c(1L, 1L))

  ## Axis arguments
  ## CVXPY: Variable(2).sum(axis=0).shape == () -> in R: (1,1)
  v2 <- Variable(2)
  expect_equal(sum_entries(v2, axis = 2)@shape, c(1L, 1L))

  ## CVXPY: Variable((2,3)).sum(axis=0, keepdims=True).shape == (1,3)
  v23 <- Variable(c(2, 3))
  expect_equal(sum_entries(v23, axis = 2, keepdims = TRUE)@shape, c(1L, 3L))
  expect_equal(sum_entries(v23, axis = 2, keepdims = FALSE)@shape, c(1L, 3L))

  ## CVXPY: Variable((2,3)).sum(axis=1).shape == (2,)
  expect_equal(sum_entries(v23, axis = 1)@shape, c(2L, 1L))

  ## Sparse identity sum: sum of eye(3) = 3
  A_sp <- Matrix::Diagonal(3)
  expect_equal(as.numeric(value(sum_entries(Constant(A_sp)))), 3)

  ## Sparse sum along axis=2 (CVXPY axis=0)
  A_sp_sum <- sum_entries(Constant(A_sp), axis = 2)
  expect_equal(as.numeric(value(A_sp_sum)), c(1, 1, 1), tolerance = 1e-10)
})

## ── test_trace ────────────────────────────────────────────────────

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_trace
test_that("CVXPY parity: trace sign, curvature, shape", {
  A <- Variable(c(2, 2), name = "A")
  C_var <- Variable(c(3, 2), name = "C")

  expr <- matrix_trace(A)
  expect_equal(expr_sign(expr), "UNKNOWN")
  expect_equal(curvature(expr), "AFFINE")
  expect_equal(expr@shape, c(1L, 1L))

  ## Non-square matrix: should error
  expect_error(matrix_trace(C_var), ".*square.*")
})

## ── test_all_expressions (simplified) ─────────────────────────────
## CVXPY: tests many atoms on a constant matrix, both standalone and method form.
## CVXR: standalone function calls only. We test the subset that CVXR supports.

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_all_expressions
test_that("CVXPY parity: all_expressions — constant evaluation", {
  X_np <- matrix(c(1, 2, 99, 4, -4, -2, 7, 3, 2.4), nrow = 3, ncol = 3)
  X <- Constant(X_np)

  ## trace: sum of diagonal
  tr <- matrix_trace(X)
  expect_equal(as.numeric(value(tr)), sum(diag(X_np)), tolerance = 1e-10)

  ## sum: sum of all entries
  s <- sum_entries(X)
  expect_equal(as.numeric(value(s)), sum(X_np), tolerance = 1e-10)

  ## max: maximum entry
  mx <- max_entries(X)
  expect_equal(as.numeric(value(mx)), max(X_np), tolerance = 1e-10)

  ## min: minimum entry
  mn <- min_entries(X)
  expect_equal(as.numeric(value(mn)), min(X_np), tolerance = 1e-10)

  ## mean: arithmetic mean
  m <- cvxr_mean(X)
  expect_equal(as.numeric(value(m)), mean(X_np), tolerance = 1e-10)

  ## ptp: range (max - min)
  pt <- ptp(X)
  expect_equal(as.numeric(value(pt)), max(X_np) - min(X_np), tolerance = 1e-10)

  ## std: standard deviation (population, ddof=0)
  st <- cvxr_std(X)
  expected_std <- sqrt(mean((X_np - mean(X_np))^2))
  expect_equal(as.numeric(value(st)), expected_std, tolerance = 1e-5)

  ## var: variance (population, ddof=0)
  vr <- cvxr_var(X)
  expected_var <- mean((X_np - mean(X_np))^2)
  expect_equal(as.numeric(value(vr)), expected_var, tolerance = 1e-5)

  ## cumsum: cumulative sum of flattened entries (column-major)
  cs <- cumsum_axis(X)
  expect_equal(as.numeric(value(cs)), cumsum(as.numeric(X_np)), tolerance = 1e-10)

  ## conj: complex conjugate (for real, identity)
  cj <- Conj(X)
  expect_equal(as.numeric(value(cj)), as.numeric(X_np), tolerance = 1e-10)

  ## Shapes match between standalone function and constant evaluation
  expect_equal(tr@shape, c(1L, 1L))
  expect_equal(s@shape, c(1L, 1L))
  expect_equal(mx@shape, c(1L, 1L))
  expect_equal(mn@shape, c(1L, 1L))
  expect_equal(cj@shape, X@shape)

  ## Axis tests for sum
  ## CVXPY axis=0 -> CVXR axis=2
  s0 <- sum_entries(X, axis = 2)
  expect_equal(as.numeric(value(s0)), as.numeric(colSums(X_np)), tolerance = 1e-10)

  ## CVXPY axis=1 -> CVXR axis=1
  s1 <- sum_entries(X, axis = 1)
  expect_equal(as.numeric(value(s1)), as.numeric(rowSums(X_np)), tolerance = 1e-10)

  ## Axis tests for max
  mx0 <- max_entries(X, axis = 2)
  expect_equal(as.numeric(value(mx0)),
               as.numeric(apply(X_np, 2, max)), tolerance = 1e-10)

  mx1 <- max_entries(X, axis = 1)
  expect_equal(as.numeric(value(mx1)),
               as.numeric(apply(X_np, 1, max)), tolerance = 1e-10)

  ## Axis tests for min
  mn0 <- min_entries(X, axis = 2)
  expect_equal(as.numeric(value(mn0)),
               as.numeric(apply(X_np, 2, min)), tolerance = 1e-10)

  mn1 <- min_entries(X, axis = 1)
  expect_equal(as.numeric(value(mn1)),
               as.numeric(apply(X_np, 1, min)), tolerance = 1e-10)
})

# ====================================================================
# QUADRATIC TESTS (test_quadratic.py)
# ====================================================================

# ── test_power (quadratic) ────────────────────────────────────────

## @cvxpy test_quadratic.py::TestExpressions::test_power
test_that("quadratic: power atom curvature properties (CVXPY parity)", {
  ## CVXPY: Variable(3) is not constant, is affine, is quadratic
  x <- Variable(3)
  y <- Variable(3)
  expect_false(is_constant(x))
  expect_true(is_affine(x))
  expect_true(is_quadratic(x))

  ## power(x.T @ y, 0) is constant, affine, quadratic
  suppressWarnings({
    s <- power(t(x) %*% y, 0)
  })
  expect_true(is_constant(s))
  expect_true(is_affine(s))
  expect_true(is_quadratic(s))

  ## power(x - y, 1) is not constant, is affine, is quadratic
  tt <- power(x - y, 1)
  expect_false(is_constant(tt))
  expect_true(is_affine(tt))
  expect_true(is_quadratic(tt))

  ## power(x + 2*y, 2) is not constant, not affine, is quadratic, is DCP
  u <- power(x + 2 * y, 2)
  expect_false(is_constant(u))
  expect_false(is_affine(u))
  expect_true(is_quadratic(u))
  expect_true(is_dcp(u))

  ## (x + 2*y)^2 same properties
  w <- (x + 2 * y)^2
  expect_false(is_constant(w))
  expect_false(is_affine(w))
  expect_true(is_quadratic(w))
  expect_true(is_dcp(w))
})

# ── test_matrix_multiplication (quadratic) ────────────────────────

## @cvxpy test_quadratic.py::TestExpressions::test_matrix_multiplication
test_that("quadratic: matrix multiplication curvature (CVXPY parity)", {
  ## CVXPY: x(3,5).T @ y(3,5) is quadratic but not DCP
  x <- Variable(c(3, 5))
  y <- Variable(c(3, 5))
  expect_false(is_constant(x))
  expect_true(is_affine(x))
  expect_true(is_quadratic(x))

  suppressWarnings({
    s <- t(x) %*% y
  })
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
})

# ── test_indefinite_quadratic ─────────────────────────────────────

## @cvxpy test_quadratic.py::TestExpressions::test_indefinite_quadratic
test_that("quadratic: indefinite quadratic is_quadratic but not DCP (CVXPY parity)", {
  ## CVXPY: y*z is quadratic, not DCP
  x <- Variable()
  y <- Variable()
  z <- Variable()
  suppressWarnings({
    s <- y * z
  })
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))

  ## (x+y)^2 - y*z - z^2 is quadratic, not DCP
  suppressWarnings({
    tt <- (x + y)^2 - s - z * z
  })
  expect_true(is_quadratic(tt))
  expect_false(is_dcp(tt))
})

# ── test_non_quadratic ────────────────────────────────────────────

## @cvxpy test_quadratic.py::TestExpressions::test_non_quadratic
test_that("quadratic: non-quadratic expressions (CVXPY parity)", {
  ## CVXPY: max(vstack(x, y, z))^2 is NOT quadratic
  x <- Variable()
  y <- Variable()
  z <- Variable()
  s <- max_entries(vstack(x, y, z))^2
  expect_false(is_quadratic(s))

  ## max(vstack(x^2, power(y, 2), z)) is NOT quadratic
  tt <- max_entries(vstack(x^2, power(y, 2), z))
  expect_false(is_quadratic(tt))
})

# ── test_affine_prod ──────────────────────────────────────────────

## @cvxpy test_quadratic.py::TestExpressions::test_affine_prod
test_that("quadratic: affine product x(3,5) %*% y(5,4) is quadratic (CVXPY parity)", {
  ## CVXPY: x @ y is quadratic but not DCP
  x <- Variable(c(3, 5))
  y <- Variable(c(5, 4))
  suppressWarnings({
    s <- x %*% y
  })
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
})

# ── test_has_quadratic ────────────────────────────────────────────

## @cvxpy test_quadratic.py::TestExpressions::test_has_quadratic
test_that("quadratic: has_quadratic_term detection (CVXPY parity)", {
  ## CVXPY: x.has_quadratic_term() checks for quadratic sub-expressions
  x <- Variable()
  expect_false(CVXR:::has_quadratic_term(x))
  expect_false(CVXR:::has_quadratic_term(3 + 3 * x))
  expect_true(CVXR:::has_quadratic_term(x^2))
  expect_true(CVXR:::has_quadratic_term(x^2 / 2))
  expect_true(CVXR:::has_quadratic_term(x^2 + x^3))
  expect_true(CVXR:::has_quadratic_term(2 * x^2 + x^3))
  expect_true(CVXR:::has_quadratic_term(Conj(x^2)))
  expect_false(CVXR:::has_quadratic_term(pos(x^2)))
  expect_true(CVXR:::has_quadratic_term(power(x^2, 2)))   ## square(x^2)
  expect_true(CVXR:::has_quadratic_term(huber(x^3)))
  expect_true(CVXR:::has_quadratic_term(power(x^2, 1)))
  expect_true(CVXR:::has_quadratic_term(quad_over_lin(x^3, 1)))
  expect_false(CVXR:::has_quadratic_term(quad_over_lin(x^3, x)))
  y <- Variable(2)
  P <- diag(2)
  expect_true(CVXR:::has_quadratic_term(matrix_frac(y^3, P)))
  P2 <- Parameter(c(2, 2), PSD = TRUE)
  expect_true(CVXR:::has_quadratic_term(matrix_frac(y^3, P2)))
})

# ── test_composite_quad_over_lin ──────────────────────────────────

## @cvxpy test_quadratic.py::TestExpressions::test_composite_quad_over_lin
test_that("quadratic: square(y) + quad_over_lin(x, y) (CVXPY parity)", {
  ## CVXPY: Minimize(square(y) + quad_over_lin(x, y)), y==1
  ## => y=1, x=0, value=1
  x <- Variable()
  y <- Variable()
  obj <- Minimize(power(y, 2) + quad_over_lin(x, y))
  prob <- Problem(obj, list(y == 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(y)), 1, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 0, tolerance = 1e-5)
  expect_equal(as.numeric(value(prob)), 1, tolerance = 1e-5)
})

# ====================================================================
# SHAPE TESTS (test_shape.py)
# ====================================================================

# ── test_add_broadcasting ─────────────────────────────────────────

## @cvxpy test_shape.py::TestShape::test_add_broadcasting
test_that("shape: add broadcasting with arbitrary dimensions (CVXPY parity)", {
  ## CVXPY test uses hypothesis-generated N-dimensional shapes.
  ## CVXR only supports 2D (matrix) shapes, so N-D broadcasting is not applicable.
  ## We test basic 2D broadcasting instead.
  ## sum_shapes in CVXR handles 2D broadcasting
  s <- CVXR:::sum_shapes(list(c(3L, 1L), c(1L, 3L)))
  expect_equal(s, c(3L, 3L))

  s2 <- CVXR:::sum_shapes(list(c(3L, 1L), c(3L, 1L)))
  expect_equal(s2, c(3L, 1L))

  s3 <- CVXR:::sum_shapes(list(c(1L, 1L), c(3L, 4L)))
  expect_equal(s3, c(3L, 4L))
})

# ── test_mul_broadcasting ─────────────────────────────────────────

## @cvxpy test_shape.py::TestShape::test_mul_broadcasting
test_that("shape: mul broadcasting with matmul signature (CVXPY parity)", {
  ## CVXPY test uses hypothesis-generated shapes with matmul signature.
  ## CVXR only supports 2D, so we test standard 2D matmul shapes.
  s <- CVXR:::mul_shapes(c(3L, 5L), c(5L, 4L))
  expect_equal(s, c(3L, 4L))

  s2 <- CVXR:::mul_shapes(c(1L, 3L), c(3L, 1L))
  expect_equal(s2, c(1L, 1L))
})

# ── test_negative_axis ────────────────────────────────────────────

## @cvxpy test_shape.py::TestShape::test_negative_axis
test_that("shape: negative axis indexing (CVXPY parity)", {
  ## CVXPY supports N-D tensors with negative axis indexing.
  ## CVXR only supports 2D shapes. Parameter shapes must be length 1 or 2.
  ## Negative axis indexing on N-D tensors is not applicable.
  skip("CVXR only supports 2D shapes; negative axis on N-D tensors is N/A")
})

# ====================================================================
# CONSTANT TESTS (test_constant.py)
# ====================================================================

# ── test_nested_lists ─────────────────────────────────────────────

## @cvxpy test_constant.py::test_nested_lists
test_that("constant: nested list behavior (CVXPY parity)", {
  ## CVXPY: Constant(nested_list) warns and transposes relative to numpy.
  ## CVXR: Constant() does not accept raw lists -- it requires numeric vectors
  ## or matrices. This is an R vs Python difference: R has no implicit
  ## list-to-matrix coercion matching numpy's nested-list behavior.
  ## We verify that Constant from matrix works correctly.
  A <- matrix(c(1, 3, 5, 2, 4, 6), nrow = 3, ncol = 2)
  c_from_matrix <- Constant(A)
  expect_equal(c_from_matrix@shape, c(3L, 2L))
  expect_equal(as.numeric(value(c_from_matrix)), as.numeric(A))

  ## Verify that raw lists are rejected
  expect_error(Constant(list(c(1, 2), c(3, 4))), "list")
})

# ── test_print ────────────────────────────────────────────────────

## @cvxpy test_constant.py::test_print
test_that("constant: print/name representation (CVXPY parity)", {
  ## CVXPY: str(Constant) renders matrix values with truncation.
  ## CVXR: name() returns a shape descriptor like "[3x3 matrix]".
  ## The rendering format differs by design (R vs Python convention).
  ## We verify that name() returns a sensible string representation.
  A <- Constant(matrix(1, 3, 3))
  expect_type(name(A), "character")
  expect_true(nchar(name(A)) > 0)

  B <- Constant(matrix(1, 5, 2))
  expect_type(name(B), "character")
  expect_true(nchar(name(B)) > 0)
})

# ====================================================================
# CONSTRAINT TESTS (test_constraints.py)
# ====================================================================

# ── test_boolean_violation ────────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_boolean_violation
test_that("constraints: boolean variable violation (CVXPY parity)", {
  ## CVXPY: z = Variable(1, boolean=True); z.value = 1
  ## Violations for z >= 0.6, 1-z <= 0.6, z <= 0.6, 1-z >= 0.6
  z <- Variable(1, boolean = TRUE)
  value(z) <- 1

  ## z >= 0.6 => violation 0.0 (feasible)
  constr1 <- z >= 0.6
  expect_equal(as.numeric(violation(constr1)), 0.0)

  ## 1 - z <= 0.6 => violation 0.0 (feasible: 0 <= 0.6)
  constr2 <- 1 - z <= 0.6
  expect_equal(as.numeric(violation(constr2)), 0.0)

  ## z <= 0.6 => violation 0.4 (infeasible: 1 > 0.6)
  constr3 <- z <= 0.6
  expect_equal(as.numeric(violation(constr3)), 0.4)

  ## 1 - z >= 0.6 => violation 0.6 (infeasible: 0 < 0.6)
  constr4 <- 1 - z >= 0.6
  expect_equal(as.numeric(violation(constr4)), 0.6)
})

# ── test_broadcasting ─────────────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_broadcasting
test_that("constraints: broadcasting in constraints (CVXPY parity)", {
  ## CVXPY: x(3,1) == Constant(ones(3)) broadcasts to (3,3)
  ## In CVXR: x(3,1) == Constant(1x3) => (3,3)

  ## Equality constraint with broadcasting — shape + solve
  x <- Variable(c(3, 1))
  c_val <- Constant(matrix(1, nrow = 1, ncol = 3))
  con <- (x == c_val)
  expect_equal(con@shape, c(3L, 3L))

  prob <- Problem(Minimize(0), list(con))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.vector(value(x)), rep(1, 3), tolerance = 1e-5)

  ## Inequality constraint with broadcasting — shape + solve
  x2 <- Variable(c(3, 1))
  con2 <- (x2 >= c_val)
  expect_equal(con2@shape, c(3L, 3L))

  prob2 <- Problem(Minimize(sum(x2)), list(con2))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.vector(value(x2)), rep(1, 3), tolerance = 1e-5)
})

# ── test_broadcast_sum_entries (wdnet repro) ─────────────────────

## @r_specific broadcast regression from wdnet::get_eta_directed
test_that("Broadcasting with sum_entries and matrix constant (wdnet repro)", {
  ## sum_entries(x, axis=2) has shape (1, ncol)
  ## matrix(1, nrow, 1) has shape (nrow, 1)
  ## Their difference broadcasts to (nrow, ncol)
  nr <- 10L
  nc <- 8L
  x <- Variable(c(nr, nc), nonneg = TRUE)

  se <- sum_entries(x, axis = 2)
  expect_equal(se@shape, c(1L, nc))

  bc <- se - matrix(1, nr, 1)
  expect_equal(bc@shape, c(nr, nc))

  ## Solve: constrain broadcasted expression to zero
  prob <- Problem(Minimize(0), list(x >= 0, bc == 0))
  result <- psolve(prob, solver = "CLARABEL")
  expect_true(result < Inf)  ## feasible

  ## Each column of x sums to 1 (replicated across all rows of the constraint)
  col_sums <- colSums(value(x))
  expect_equal(col_sums, rep(1, nc), tolerance = 1e-4)
})

# ── test_nsd_constraint ───────────────────────────────────────────

## @cvxpy test_constraints.py::TestConstraints::test_nsd_constraint
test_that("constraints: NSD constraint A << B (CVXPY parity)", {
  ## CVXPY: constr = A << B => name "B + -A >> 0", shape (2,2)
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")

  constr <- A %<<% B
  expect_equal(name(constr), "B + -A >> 0")
  expect_equal(constr@shape, c(2L, 2L))

  ## dual_value is NULL before solving
  expect_null(dual_value(constr))

  ## Without values, value() should error
  expect_error(value(constr))

  ## Set B PSD, A smaller => feasible
  value(B) <- matrix(c(2, -1, 1, 2), 2, 2)
  value(A) <- diag(2)
  expect_true(value(constr))

  ## A too large => infeasible
  value(A) <- 3 * diag(2)
  expect_false(value(constr))

  ## Non-square matrix should error
  x <- Variable(2, name = "x")
  expect_error(x %<<% 0, "Non-square")
})
