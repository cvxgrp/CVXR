## CVXPY parity tests: MEDIUM-priority — Problem and Expression gaps
## These tests close MEDIUM gaps vs CVXPY test suite:
##   - test_problem.py (17 tests)
##   - test_expressions.py (10 tests)
##
## IMPORTANT: R uses 1-based indexing and (n,1) shape for vectors,
## while CVXPY uses 0-based indexing and (n,) shape.
## All expected values verified against CVXPY 1.8.1.

library(testthat)

# ====================================================================
# PROBLEM TESTS (17 medium)
# ====================================================================

# ── 1. CLARABEL_noineq ──────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_CLARABEL_noineq
test_that("Problem with no inequality constraints, Clarabel solver (CVXPY parity)", {
  ## CVXPY: test_CLARABEL_noineq
  ## T = ones(2,2), Problem(Minimize(1), [A == T]) => result=1, A=T
  A <- Variable(c(2, 2), name = "A")
  T_val <- matrix(1, 2, 2)
  prob <- Problem(Minimize(1), list(A == T_val))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 1, tolerance = 1e-4)
  expect_equal(as.matrix(value(A)), T_val, tolerance = 1e-4)
})

# ── 2. add_problems ─────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_add_problems
test_that("Adding two problems together (CVXPY parity)", {
  ## CVXPY: test_add_problems — problem arithmetic (prob1 + prob2)
  x <- Variable()
  y <- Variable()
  p1 <- Problem(Minimize(x), list(x >= 3))
  p2 <- Problem(Minimize(y), list(y >= 5))
  p3 <- p1 + p2
  psolve(p3)
  expect_equal(as.numeric(value(p3)), 8.0, tolerance = 1e-4)
})

# ── 3. mul_problems ─────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_mul_problems
test_that("Multiplying problem by scalar (CVXPY parity)", {
  ## CVXPY: test_mul_problems — f * prob and prob * f
  x <- Variable()
  p1 <- Problem(Minimize(x), list(x >= 3))
  p2 <- p1 * 2
  psolve(p2)
  expect_equal(as.numeric(value(p2)), 6.0, tolerance = 1e-4)
  ## Commutative
  p3 <- 2 * p1
  psolve(p3)
  expect_equal(as.numeric(value(p3)), 6.0, tolerance = 1e-4)
})

# ── 4. lin_combination_problems ──────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_lin_combination_problems
test_that("Linear combination of problems (CVXPY parity)", {
  ## CVXPY: test_lin_combination_problems — prob1 + 2*prob2
  x <- Variable()
  y <- Variable()
  p1 <- Problem(Minimize(x), list(x >= 3))
  p2 <- Problem(Minimize(y), list(y >= 5))
  p3 <- p1 + 2 * p2
  psolve(p3)
  expect_equal(as.numeric(value(p3)), 13.0, tolerance = 1e-4)
  ## Division
  p4 <- p1 / 2
  psolve(p4)
  expect_equal(as.numeric(value(p4)), 1.5, tolerance = 1e-4)
  ## Negation
  p5 <- -p1
  psolve(p5)
  expect_equal(as.numeric(value(p5)), -3.0, tolerance = 1e-4)
})

# ── 5. mult_constant_atoms ──────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_mult_constant_atoms
test_that("Multiple constant atoms in objective (CVXPY parity)", {
  ## CVXPY: test_mult_constant_atoms
  ## pnorm([3,4], p=2) * a s.t. a >= 2 => 5*2 = 10
  a <- Variable(name = "a")
  prob <- Problem(
    Minimize(p_norm(matrix(c(3, 4), 2, 1), p = 2) * a),
    list(a >= 2)
  )
  result <- psolve(prob, solver = "SCS", reltol = 1e-6)
  expect_equal(result, 10, tolerance = 1e-3)
  expect_equal(as.numeric(value(a)), 2, tolerance = 1e-3)
})

# ── 6. mult_by_zero ─────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_mult_by_zero
test_that("Multiplying expression by zero (CVXPY parity)", {
  ## CVXPY: test_mult_by_zero
  ## 0 * a => value is 0, problem value is 0
  a <- Variable(name = "a")
  value(a) <- 1
  exp_zero <- 0 * a
  expect_equal(as.numeric(value(exp_zero)), 0)
  prob <- Problem(Minimize(exp_zero))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 0, tolerance = 1e-6)
  ## Variable should still have a value after solving
  expect_false(is.null(value(a)))
})

# ── 7. multiplication_on_left ────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_multiplication_on_left
test_that("Left-multiply: constant * variable (CVXPY parity)", {
  ## CVXPY: test_multiplication_on_left
  A <- Variable(c(2, 2), name = "A")
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")

  ## c^T A c where c = [1, 2]^T and A >= 2
  ## min c^T A c = [1,2] %*% A %*% [1,2] with A_ij >= 2
  ## => all entries of A = 2, c^T A c = [1,2] %*% [[2,2],[2,2]] %*% [1,2]
  ## = [1,2] %*% [6, 6] = 6 + 12 = 18
  c_vec <- matrix(c(1, 2), 2, 1)
  prob <- Problem(Minimize(t(c_vec) %*% A %*% c_vec), list(A >= 2))
  result <- psolve(prob, solver = "SCS", reltol = 1e-6)
  expect_equal(result, 18, tolerance = 1e-2)

  ## a * 2 s.t. a >= 2 => 4
  a <- Variable(name = "a2")
  prob2 <- Problem(Minimize(a * 2), list(a >= 2))
  result2 <- psolve(prob2, solver = "SCS", reltol = 1e-6)
  expect_equal(result2, 4, tolerance = 1e-3)

  ## x^T c s.t. x >= 2 => [2,2]^T [1,2] = 2 + 4 = 6
  prob3 <- Problem(Minimize(t(x) %*% c_vec), list(x >= 2))
  result3 <- psolve(prob3, solver = "SCS", reltol = 1e-6)
  expect_equal(result3, 6, tolerance = 1e-3)

  ## (x^T + z^T) c s.t. x >= 2, z >= 1 => [3,3]^T [1,2] = 3 + 6 = 9
  prob4 <- Problem(Minimize((t(x) + t(z)) %*% c_vec), list(x >= 2, z >= 1))
  result4 <- psolve(prob4, solver = "SCS", reltol = 1e-6)
  expect_equal(result4, 9, tolerance = 1e-3)
})

# ── 8. diag_prob ─────────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_diag_prob
test_that("Problem using DiagVec/DiagMat atoms (CVXPY parity)", {
  ## CVXPY: test_diag_prob
  ## C is 3x3, diag(C) == 1, C[0,1]==0.6, C[1,2]==-0.3, C is PSD
  ## Maximize C[0,2] => approximately 0.583
  C <- Variable(c(3, 3), name = "C")
  C_psd <- Variable(c(3, 3), PSD = TRUE, name = "C_psd")
  prob <- Problem(
    Maximize(C[1, 3]),
    list(
      DiagMat(C) == matrix(1, 3, 1),
      C[1, 2] == 0.6,
      C[2, 3] == -0.3,
      C == C_psd
    )
  )
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 0.583, tolerance = 0.02)
})

## Helper: build a matrix with vector v on k-th diagonal of an n x n zero matrix
.make_diag_matrix <- function(v, k, n) {
  m <- matrix(0, n, n)
  len <- length(v)
  if (k >= 0L) {
    for (i in seq_len(len)) m[i, i + k] <- v[i]
  } else {
    for (i in seq_len(len)) m[i - k, i] <- v[i]
  }
  m
}

# ── 9. diag_offset_problem ──────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_diag_offset_problem
test_that("DiagVec/DiagMat with offset k != 0 (CVXPY parity)", {
  ## CVXPY: test_diag_offset_problem
  ## For different k offsets, test DiagVec and DiagMat atoms
  n <- 4L
  A <- matrix(0:(n^2 - 1), n, n, byrow = TRUE)  # np.arange(16).reshape(4,4)

  ## Helper to extract k-th diagonal from a matrix (like np.diag(A, k))
  .extract_diag <- function(M, k) {
    nr <- nrow(M)
    len <- nr - abs(k)
    v <- numeric(len)
    if (k >= 0L) {
      for (i in seq_len(len)) v[i] <- M[i, i + k]
    } else {
      for (i in seq_len(len)) v[i] <- M[i - k, i]
    }
    v
  }

  ## Test a subset of offsets to keep test brief
  for (k in c(-2L, -1L, 0L, 1L, 2L)) {
    ## DiagVec: x -> diag matrix with x on k-th diagonal
    diag_vals <- .extract_diag(A, k)
    ## CVXPY diag(x, k) produces a matrix. DiagVec in R does the same.
    x <- Variable(n - abs(k))
    prob <- Problem(
      Minimize(sum(x)),
      list(DiagVec(x, k) == .make_diag_matrix(diag_vals, k, n))
    )
    result <- psolve(prob, solver = "SCS", reltol = 1e-6)
    expect_equal(result, sum(diag_vals), tolerance = 1e-3,
                 label = paste0("DiagVec k=", k))
    expect_equal(as.numeric(value(x)), diag_vals, tolerance = 1e-3,
                 label = paste0("DiagVec x values k=", k))

    ## DiagMat: extract k-th diagonal from matrix variable
    X <- Variable(c(n, n), nonneg = TRUE)
    prob2 <- Problem(
      Minimize(sum(X)),
      list(DiagMat(X, k) == matrix(diag_vals, ncol = 1))
    )
    result2 <- psolve(prob2, solver = "SCS", reltol = 1e-6)
    expect_equal(result2, sum(diag_vals), tolerance = 1e-3,
                 label = paste0("DiagMat k=", k))
  }
})

# ── 10. neg_slice ────────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_neg_slice
test_that("Negative indexing on expressions (CVXPY parity)", {
  ## CVXPY: test_neg_slice
  ## x[-2:] >= 1 in Python means last 2 elements.
  ## In R, negative indexing EXCLUDES elements. The CVXPY intent is
  ## to constrain the last 2 elements. In R: x[1:2, 1] for a 2-vector.
  ## For a 2-vector, x[-2:] in Python means x[0:], i.e. all elements.
  ## Let's test: Minimize x[1] + x[2] s.t. x >= 1 for 2-element vector.
  x <- Variable(2)
  prob <- Problem(Minimize(x[1, 1] + x[2, 1]), list(x >= 1))
  psolve(prob, solver = "SCS", reltol = 1e-6)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-3)

  ## Test R-style negative indexing: exclude position 1, keep the rest
  ## For a 4-vector, x[-1, 1] gives elements 2,3,4
  y <- Variable(4)
  prob2 <- Problem(Minimize(sum(y)), list(y[-1, 1] >= 5, y[1, 1] >= 10))
  psolve(prob2, solver = "SCS", reltol = 1e-6)
  expect_equal(as.numeric(value(y))[1], 10, tolerance = 1e-3)
  expect_equal(as.numeric(value(y))[2:4], c(5, 5, 5), tolerance = 1e-3)
})

# ── 11. bool_constr ──────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_bool_constr
test_that("Boolean constraints in optimization (CVXPY parity)", {
  ## CVXPY: test_bool_constr — supports TRUE/FALSE in constraints list.
  x <- Variable()
  ## TRUE is trivially satisfied, problem is feasible
  p1 <- Problem(Minimize(x), list(TRUE, x >= 3))
  psolve(p1)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-4)
  ## FALSE makes problem infeasible
  p2 <- Problem(Minimize(x), list(FALSE, x >= 3))
  psolve(p2)
  expect_equal(status(p2), "infeasible")
})

# ── 12. indicator ────────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_indicator
test_that("Indicator constraints (CVXPY parity)", {
  ## CVXPY: test_indicator — indicator(constraints)
  ## Converts constraints into an indicator function for the objective.
  n <- 5L
  m <- 2L
  q <- seq(0, n - 1)  # 0:4 like np.arange(5)
  a <- matrix(1, m, n)
  b <- matrix(1, m, 1)
  x <- Variable(c(n, 1))
  constr <- list(a %*% x == b)
  objective <- Minimize(0.5 * square(t(q) %*% x) + indicator(constr))
  problem <- Problem(objective)
  solution1 <- psolve(problem, solver = "SCS", eps = 1e-5)

  ## Without indicators.
  objective2 <- Minimize(0.5 * square(t(q) %*% x))
  problem2 <- Problem(objective2, constr)
  solution2 <- psolve(problem2, solver = "SCS", eps = 1e-5)

  expect_equal(as.numeric(solution1), as.numeric(solution2), tolerance = 1e-2)
})

# ── 13. rmul_scalar_mats ────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_rmul_scalar_mats
test_that("Right-multiply scalar by matrix (CVXPY parity)", {
  ## CVXPY: test_rmul_scalar_mats
  ## Two formulations should give the same result:
  ## (1) quad_form(z, [[x]]) - 2 * z^T @ [[y]] with z as (1,1) Variable
  ## (2) x * z^2 - 2 * z * y with z as scalar Variable
  x_val <- 4144.30127531
  y_val <- 7202.52114311

  ## Formulation 1: 1x1 matrix
  z1 <- Variable(c(1, 1))
  obj1 <- Minimize(quad_form(z1, matrix(x_val, 1, 1)) - 2 * t(z1) %*% matrix(y_val, 1, 1))
  prob1 <- Problem(obj1)
  psolve(prob1, solver = "OSQP")
  result1 <- value(prob1)

  ## Formulation 2: scalar
  z2 <- Variable()
  obj2 <- Minimize(x_val * square(z2) - 2 * z2 * y_val)
  prob2 <- Problem(obj2)
  psolve(prob2, solver = "OSQP")

  expect_equal(as.numeric(value(prob2)), as.numeric(result1), tolerance = 1e-2)
})

# ── 14. min_with_axis ────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_min_with_axis
test_that("MinEntries/MaxEntries with axis parameter (CVXPY parity)", {
  ## CVXPY: test_min_with_axis
  ## In CVXR, axis=0 reduces rows, keeping ncol as a column vector: (ncol, 1)
  ## axis=1 reduces columns, keeping nrow as a column vector: (nrow, 1)
  x <- Variable(c(3, 2))

  ## min_entries with axis=2 (reduce rows) => shape (1, 2)
  ## CVXPY: shape (2,); CVXR: shape (1, 2) — row vector
  expr_min_col <- min_entries(x, axis = 2L)
  expect_equal(expr_min_col@shape[1L], 1L)
  expect_equal(expr_min_col@shape[2L], 2L)

  ## min_entries with axis=1 (reduce cols) => shape (3, 1)
  expr_min_row <- min_entries(x, axis = 1L)
  expect_equal(expr_min_row@shape[1L], 3L)
  expect_equal(expr_min_row@shape[2L], 1L)

  ## max_entries with axis=2 (reduce rows) => shape (1, 2)
  expr_max_col <- max_entries(x, axis = 2L)
  expect_equal(expr_max_col@shape[1L], 1L)
  expect_equal(expr_max_col@shape[2L], 2L)

  ## max_entries with axis=1 (reduce cols) => shape (3, 1)
  expr_max_row <- max_entries(x, axis = 1L)
  expect_equal(expr_max_row@shape[1L], 3L)
  expect_equal(expr_max_row@shape[2L], 1L)

  ## Solve a simple problem with elementwise minimum
  ## x is (5,2), y is (5,2), min_elemwise(x, y) with x==1, y==2
  ## => result is min(1, 2) = 1 for each entry, sum = 10
  x2 <- Variable(c(5, 2))
  y2 <- Variable(c(5, 2))
  elem_min <- min_elemwise(x2, y2)
  prob <- Problem(Maximize(sum(elem_min)), list(x2 == 1, y2 == 2))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 10, tolerance = 1e-3)
})

# ── 15. huber_scs ────────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_huber_scs
test_that("Huber atom solved with SCS (CVXPY parity)", {
  ## CVXPY: test_huber_scs — regression test for issue #1370
  ## Huber regression with outliers
  set.seed(1)
  m <- 5L
  n <- 2L

  x0 <- rnorm(n)
  A <- matrix(rnorm(m * n), m, n)
  b <- as.numeric(A %*% x0) + 0.01 * rnorm(m)
  ## Add outlier noise
  k <- max(1L, as.integer(0.02 * m))
  idx <- sample.int(m, k)
  b[idx] <- b[idx] + 10 * rnorm(k)

  x <- Variable(n)
  prob <- Problem(Minimize(sum(huber(A %*% x - matrix(b, m, 1)))))
  ## This should solve without error
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), OPTIMAL)
})

# ── 16. rmul_param ──────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_rmul_param
test_that("Right-multiply with Parameter (CVXPY parity)", {
  ## CVXPY: test_rmul_param — issue #1555
  ## (2 * b) @ param with b as (1,) Variable, param as Parameter(1)
  ## Unbounded problem => -Inf
  b <- Variable(c(1, 1))
  param <- Parameter(1)

  obj <- Minimize((2 * b) %*% param)
  prob <- Problem(obj)
  value(param) <- matrix(1, 1, 1)
  psolve(prob, solver = "CLARABEL")
  ## Should be unbounded below (Clarabel reports DUAL_INFEASIBLE for unbounded)
  expect_true(
    value(prob) == -Inf ||
    status(prob) %in% c(UNBOUNDED, INFEASIBLE_OR_UNBOUNDED, "infeasible_or_unbounded")
  )
})

# ── 17. multiply_by_scalar ──────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_multiply_by_scalar
test_that("Multiply expression by scalar (CVXPY parity)", {
  ## CVXPY: test_multiply_by_scalar
  ## In CVXPY: multiply(2, x) is elementwise multiply.
  ## In CVXR, multiply() is deprecated; use `*` operator for scalar * expression.
  ## 2 * x s.t. x == 2 => sum = 2*2 + 2*2 = 8
  x <- Variable(2, name = "x")
  obj <- Minimize(sum(2 * x))
  prob <- Problem(obj, list(x == 2))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 8, tolerance = 1e-3)

  ## Also test that multiply function still works (deprecated but functional)
  x2 <- Variable(2, name = "x2")
  suppressWarnings({
    obj2 <- Minimize(sum(multiply(2, x2)))
  })
  prob2 <- Problem(obj2, list(x2 == 2))
  result2 <- psolve(prob2, solver = "SCS")
  expect_equal(result2, 8, tolerance = 1e-3)
})


# ====================================================================
# EXPRESSION TESTS (10 medium)
# ====================================================================

# ── 1. transpose_variable ───────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_transpose_variable
test_that("Transpose of variable: shape and solving (CVXPY parity)", {
  ## CVXPY: test_transpose_variable
  ## Scalar transpose is a no-op
  a <- Variable(name = "a")
  var_t <- t(a)
  ## Scalar: shape remains (1, 1)
  expect_equal(var_t@shape, c(1L, 1L))

  ## Column vector (2,1) transposed => (1,2)
  x <- Variable(c(2, 1), name = "x")
  var_t2 <- t(x)
  expect_equal(var_t2@shape, c(1L, 2L))

  ## Set value and check transpose value
  value(x) <- matrix(c(1, 2), 2, 1)
  expect_equal(as.numeric(value(var_t2)), c(1, 2))
  expect_equal(dim(value(var_t2)), c(1L, 2L))

  ## Matrix (3,2) transposed => (2,3)
  C <- Variable(c(3, 2), name = "C")
  var_t3 <- t(C)
  expect_equal(var_t3@shape, c(2L, 3L))

  ## Index into transpose
  idx <- var_t3[2, 1]
  expect_equal(idx@shape, c(1L, 1L))

  ## Double transpose: back to original shape
  var_t4 <- t(t(x))
  expect_equal(var_t4@shape, c(2L, 1L))

  ## Use in a problem: min t(x) %*% c s.t. x >= 1
  c_vec <- matrix(c(3, 4), 2, 1)
  x2 <- Variable(c(2, 1))
  prob <- Problem(Minimize(t(x2) %*% c_vec), list(x2 >= 1))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 7, tolerance = 1e-4)
  expect_equal(as.numeric(value(x2)), c(1, 1), tolerance = 1e-4)
})

# ── 2. constant_psd_nsd ─────────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_constant_psd_nsd
test_that("Constant.is_psd() and Constant.is_nsd() queries (CVXPY parity)", {
  ## CVXPY: test_constant_psd_nsd
  ## Test indefinite matrices
  set.seed(0)
  n <- 5L
  U <- matrix(rnorm(n * n), n, n)
  U <- U %*% t(U)
  eig <- eigen(U, symmetric = TRUE)
  U <- eig$vectors  # orthogonal matrix

  ## Indefinite: has both positive and negative eigenvalues
  v1 <- c(3, 2, 1, 1e-8, -1)
  P1 <- Constant(U %*% diag(v1) %*% t(U))
  expect_false(is_psd(P1))
  expect_false(is_nsd(P1))

  v2 <- c(3, 2, 2, 1e-6, -1)
  P2 <- Constant(U %*% diag(v2) %*% t(U))
  expect_false(is_psd(P2))
  expect_false(is_nsd(P2))

  v3 <- c(3, 2, 2, 1e-4, -1e-6)
  P3 <- Constant(U %*% diag(v3) %*% t(U))
  expect_false(is_psd(P3))
  expect_false(is_nsd(P3))

  v4 <- c(-1, 3, 0, 0, 0)
  P4 <- Constant(U %*% diag(v4) %*% t(U))
  expect_false(is_psd(P4))
  expect_false(is_nsd(P4))

  ## GitHub issue 1451 equivalent: [[1,2],[2,1]] is indefinite
  P_indef <- Constant(matrix(c(1, 2, 2, 1), 2, 2))
  x <- Variable(2)
  expr <- quad_form(x, matrix(c(1, 2, 2, 1), 2, 2))
  expect_false(is_dcp(expr))
  suppressWarnings(expect_false(is_dcp(-expr)))

  ## Large eigenvalue spread: diag(9 x 1e-4, -1e4)
  P_large <- Constant(diag(c(rep(1e-4, 9), -1e4)))
  expect_false(is_psd(P_large))
  expect_false(is_nsd(P_large))

  ## Actually PSD: ones(5,5) has eigenvalues {5, 0, 0, 0, 0}
  P_psd <- Constant(matrix(1, 5, 5))
  expect_true(is_psd(P_psd))
  expect_false(is_nsd(P_psd))

  ## Sparse PSD: identity matrix
  P_eye <- Constant(diag(10))
  expect_true(is_psd(P_eye))
  neg_eye <- Constant(-diag(10))
  expect_true(is_nsd(neg_eye))
})

# ── 3. constant_skew_symmetric ──────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_constant_skew_symmetric
test_that("Constant.is_skew_symmetric() (CVXPY parity)", {
  ## CVXPY: test_constant_skew_symmetric
  ## is_skew_symmetric is an internal generic; access via CVXR:::
  is_skew <- CVXR:::is_skew_symmetric

  ## Identity: NOT skew-symmetric
  M1 <- diag(3)
  expect_false(is_skew(Constant(M1)))

  ## Zero matrix: IS skew-symmetric (A = -A^T trivially)
  M2 <- matrix(0, 3, 3)
  expect_true(is_skew(Constant(M2)))

  ## [[0, 1], [-1, 0]]: IS skew-symmetric
  ## R column-major: matrix(c(0, -1, 1, 0), 2, 2) gives row1=[0,1], row2=[-1,0]
  M3 <- matrix(c(0, -1, 1, 0), 2, 2)
  expect_true(is_skew(Constant(M3)))

  ## [[0, -1], [1, 0]]: IS skew-symmetric
  M4 <- matrix(c(0, 1, -1, 0), 2, 2)
  expect_true(is_skew(Constant(M4)))

  ## [[0, 1], [1, 0]]: NOT skew-symmetric (symmetric instead)
  M5 <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_false(is_skew(Constant(M5)))

  ## [[1, 1], [-1, 0]]: NOT skew-symmetric (nonzero diagonal)
  M6 <- matrix(c(1, -1, 1, 0), 2, 2)
  expect_false(is_skew(Constant(M6)))

  ## [[0, 1], [-1.1, 0]]: NOT skew-symmetric (off-diagonals don't cancel)
  M7 <- matrix(c(0, -1.1, 1, 0), 2, 2)
  expect_false(is_skew(Constant(M7)))
})

# ── 4. psd_nsd_parameters ──────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_psd_nsd_parameters
test_that("Parameter with PSD/NSD, use in problem (CVXPY parity)", {
  ## CVXPY: test_psd_nsd_parameters
  ## PSD parameter accepts a valid PSD matrix
  set.seed(42)
  n <- 10L
  a <- matrix(rnorm(n * (n - 1)), n, n - 1)
  a2 <- a %*% t(a)  # PSD matrix (n x n)
  p <- Parameter(c(n, n), PSD = TRUE)
  value(p) <- a2
  expect_equal(as.matrix(value(p)), a2, tolerance = 1e-10)

  ## Arithmetic preserves PSD/NSD
  p2 <- Parameter(c(2, 2), PSD = TRUE)
  expect_true(is_psd(2 * p2))
  expect_true(is_psd(p2 + p2))
  expect_true(is_nsd(-p2))

  ## Invalid PSD parameter: indefinite matrix should error
  n2 <- 5L
  P <- Parameter(c(n2, n2), PSD = TRUE)
  N <- Parameter(c(n2, n2), NSD = TRUE)
  set.seed(0)
  U <- matrix(rnorm(n2 * n2), n2, n2)
  U <- U %*% t(U)
  eig <- eigen(U, symmetric = TRUE)
  U <- eig$vectors

  ## Each of these has a negative eigenvalue, so should fail for PSD
  v_bad <- c(3, 2, 1, 1e-8, -1)
  expect_error(
    value(P) <- U %*% diag(v_bad) %*% t(U),
    "semidefinite"
  )
  ## And the negation should fail for NSD
  expect_error(
    value(N) <- -(U %*% diag(v_bad) %*% t(U)),
    "semidefinite"
  )
})

# ── 5. round_attr ───────────────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_round_attr
test_that("Rounding attributes via project (CVXPY parity)", {
  ## CVXPY: test_round_attr — tests Variable.project() for various attributes
  ## In CVXR, project() is an internal function but exported for tests.

  ## Nonpos: project positive values to 0
  v1 <- Variable(1, nonpos = TRUE)
  proj1 <- CVXR:::project(v1, matrix(1, 1, 1))
  expect_equal(as.numeric(proj1), 0, tolerance = 1e-10)

  v1b <- Variable(c(2, 1), nonpos = TRUE)
  proj1b <- CVXR:::project(v1b, matrix(c(1, -1), 2, 1))
  expect_equal(as.numeric(proj1b), c(0, -1), tolerance = 1e-10)

  ## Nonneg: project negative values to 0
  v2 <- Variable(1, nonneg = TRUE)
  proj2 <- CVXR:::project(v2, matrix(-1, 1, 1))
  expect_equal(as.numeric(proj2), 0, tolerance = 1e-10)

  v2b <- Variable(c(2, 1), nonneg = TRUE)
  proj2b <- CVXR:::project(v2b, matrix(c(1, -1), 2, 1))
  expect_equal(as.numeric(proj2b), c(1, 0), tolerance = 1e-10)

  ## Boolean: project to {0, 1}
  v3 <- Variable(c(2, 2), boolean = TRUE)
  proj3 <- CVXR:::project(v3, matrix(c(1, 1, -1, 0), 2, 2))
  expect_equal(as.numeric(proj3), c(1, 1, 0, 0), tolerance = 1e-10)

  ## Integer: project to nearest integer
  v4 <- Variable(c(2, 2), integer = TRUE)
  proj4 <- CVXR:::project(v4, matrix(c(1, 1, -1.6, 0), 2, 2))
  expect_equal(as.numeric(proj4), c(1, 1, -2, 0), tolerance = 1e-10)

  ## Symmetric: project (A + A^T)/2
  v5 <- Variable(c(2, 2), symmetric = TRUE)
  proj5 <- CVXR:::project(v5, matrix(c(1, 1, -1, 0), 2, 2))
  expect_equal(as.numeric(proj5), c(1, 0, 0, 0), tolerance = 1e-10)

  ## PSD: project to positive semidefinite cone
  v6 <- Variable(c(2, 2), PSD = TRUE)
  proj6 <- CVXR:::project(v6, matrix(c(1, 1, -1, -1), 2, 2))
  ## Eigenvalues of [[1,-1],[1,-1]] are 0 and 0, but the matrix is not symmetric.
  ## PSD projection first symmetrizes, then projects eigenvalues.
  ## Expected: [[1, 0], [0, 0]] (matching CVXPY)
  expect_equal(as.numeric(proj6), c(1, 0, 0, 0), tolerance = 1e-3)

  ## NSD: project to negative semidefinite cone
  v7 <- Variable(c(2, 2), NSD = TRUE)
  proj7 <- CVXR:::project(v7, matrix(c(1, 1, -1, -1), 2, 2))
  ## Expected: [[0, 0], [0, -1]] (matching CVXPY)
  expect_equal(as.numeric(proj7), c(0, 0, 0, -1), tolerance = 1e-3)

  ## Diag: project to diagonal matrices
  v8 <- Variable(c(2, 2), diag = TRUE)
  proj8 <- CVXR:::project(v8, matrix(c(1, 1, -1, 0), 2, 2))
  expect_equal(as.numeric(as.matrix(proj8)), c(1, 0, 0, 0), tolerance = 1e-10)
})

# ── 6. scalar_const_promotion ────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_scalar_const_promotion
test_that("Scalar constant auto-promoted in expressions (CVXPY parity)", {
  ## CVXPY: test_scalar_const_promotion
  x <- Variable(2, name = "x")
  A <- Variable(c(2, 2), name = "A")

  ## Vector + scalar => vector shape
  exp1 <- x + 2
  expect_equal(curvature(exp1), "AFFINE")
  expect_true(is_affine(exp1))
  expect_equal(exp1@shape, c(2L, 1L))

  ## Scalar - vector => vector shape
  exp2 <- 4 - x
  expect_equal(exp2@shape, c(2L, 1L))

  ## Scalar * vector => vector shape
  exp3 <- 4 * x
  expect_equal(exp3@shape, c(2L, 1L))

  ## Scalar <= vector => constraint shape
  constr1 <- (4 <= x)
  expect_equal(constr1@shape, c(2L, 1L))

  ## Scalar == vector => constraint shape
  constr2 <- (4 == x)
  expect_equal(constr2@shape, c(2L, 1L))

  ## Vector >= scalar => constraint shape
  constr3 <- (x >= 4)
  expect_equal(constr3@shape, c(2L, 1L))

  ## Matrix + scalar + scalar => matrix shape
  exp4 <- (A + 2) + 4
  expect_equal(curvature(exp4), "AFFINE")
  expect_equal(exp4@shape, c(2L, 2L))

  ## Scalar * matrix => matrix shape
  exp5 <- 3 * A
  expect_equal(exp5@shape, c(2L, 2L))
})

# ── 7. var_copy ──────────────────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_var_copy
test_that("Copying/cloning a Variable (CVXPY parity)", {
  ## CVXPY: test_var_copy — Leaf.copy() returns self
  x <- Variable(c(3, 4), name = "x")
  y <- expr_copy(x)
  expect_identical(y, x)
  expect_equal(dim(y), c(3L, 4L))
  expect_equal(expr_name(y), "x")
})

# ── 8. param_copy ────────────────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_param_copy
test_that("Copying/cloning a Parameter (CVXPY parity)", {
  ## CVXPY: test_param_copy — Leaf.copy() returns self
  p <- Parameter(nonneg = TRUE, name = "p")
  q <- expr_copy(p)
  expect_identical(q, p)
  expect_true(is_nonneg(q))
  expect_equal(expr_name(q), "p")
})

# ── 9. constant_copy ─────────────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_constant_copy
test_that("Copying/cloning a Constant (CVXPY parity)", {
  ## CVXPY: test_constant_copy — Leaf.copy() returns self
  c1 <- Constant(42)
  c2 <- expr_copy(c1)
  expect_identical(c2, c1)
  expect_equal(as.numeric(value(c2)), 42)
})

# ── 10. matmul_scalars ──────────────────────────────────────────────

## @cvxpy test_expressions.py::TestExpressions::test_matmul_scalars
test_that("%*% with scalar operands (CVXPY parity)", {
  ## CVXPY: test_matmul_scalars
  ## quad_form(x, I) has shape () — then multiplying by array([2]) should work
  x <- Variable(2)
  quad <- quad_form(x, diag(2))
  ## quad should be scalar-shaped
  expect_true(CVXR:::expr_is_scalar(quad))

  ## Set value and check computation
  value(x) <- c(1, 2)
  ## quad_form = x^T I x = 1 + 4 = 5
  expect_equal(as.numeric(value(quad)), 5, tolerance = 1e-10)

  ## Multiply quad_form * 2
  expr <- quad * 2
  expect_equal(as.numeric(value(expr)), 10, tolerance = 1e-10)

  ## Use in a problem: Minimize(quad_form(x, I) * 2) s.t. sum(x) == 3
  x2 <- Variable(2)
  prob <- Problem(Minimize(quad_form(x2, diag(2)) * 2), list(sum(x2) == 3))
  psolve(prob, solver = "OSQP")
  ## Minimum of 2*(x1^2 + x2^2) s.t. x1+x2=3 => x1=x2=1.5, obj=2*(2.25+2.25)=9
  expect_equal(value(prob), 9, tolerance = 1e-2)
  expect_equal(as.numeric(value(x2)), c(1.5, 1.5), tolerance = 1e-2)
})
