## Wave 2: Problem Solve-Level Tests — ports from CVXPY's test_problem.py
## Tests operator patterns, validation, SDP, and solve patterns in solve context.
##
## Expected values verified via `uv run python` against CVXPY.

# ══════════════════════════════════════════════════════════════════
# Operator patterns in solve context
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_slicing
test_that("slicing: basic row slicing with constraints", {
  C <- Variable(c(3, 2))
  prob <- Problem(Maximize(sum_entries(C)), list(C[2:3, ] <= 2, C[1, ] == 1))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 10, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_slicing
test_that("slicing: transpose of slice", {
  C <- Variable(c(3, 2))
  prob <- Problem(Maximize(sum_entries(C)), list(t(C[2:3, ]) <= 2, t(C[1, ]) == 1))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 10, tolerance = 1e-3)
})
## @cvxpy test_problem.py::TestProblem::test_transpose

test_that("transpose: A[i,j] == t(A)[j,i]", {
  A <- Variable(c(2, 2))
  ## A = [1 2; 3 4] column-major: c(1,3,2,4)
  prob <- Problem(Minimize(A[1, 2] - t(A)[2, 1]),
                  list(A == matrix(c(1, 3, 2, 4), 2, 2)))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 0, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_transpose
test_that("transpose: negate then transpose", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(t(-x) <= 1))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, -2, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_transpose
test_that("transpose: slice of transpose", {
  C <- Variable(c(3, 2))
  ## t(C) is (2,3); t(C)[, 2:3] selects cols 2,3
  prob <- Problem(Maximize(sum_entries(C)), list(t(C)[, 2:3] <= 2, t(C)[, 1] == 1))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 10, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_hstack
test_that("hstack: basic vector hstack", {
  x <- Variable(2)
  y <- Variable(3)
  cvec <- matrix(1, 1, 5)
  ## R shapes: x is (2,1), t(x) is (1,2); hstack of (1,2) and (1,3) → (1,5)
  prob <- Problem(Minimize(cvec %*% t(hstack(t(x), t(y)))),
                  list(x == matrix(c(1, 2), 2, 1),
                       y == matrix(c(3, 4, 5), 3, 1)))
  result <- psolve(prob, solver = "SCS", eps = 1e-8)
  expect_equal(result, 15, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_hstack
test_that("hstack: same variable twice", {
  x <- Variable(2)
  cvec <- matrix(1, 1, 4)
  prob <- Problem(Minimize(cvec %*% t(hstack(t(x), t(x)))),
                  list(x == matrix(c(1, 2), 2, 1)))
  result <- psolve(prob, solver = "SCS", eps = 1e-8)
  expect_equal(result, 6, tolerance = 1e-3)
})
## @cvxpy test_problem.py::TestProblem::test_div

test_that("div: division by scalar", {
  A <- Variable(c(2, 2))
  prob <- Problem(Minimize(norm_inf(A / 5)), list(A >= 5))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 1, tolerance = 1e-4)
})
## @cvxpy test_problem.py::TestProblem::test_div

test_that("div: scalar promotion", {
  a <- Variable(1)
  cvec <- matrix(c(1, 2, -1, -2), 2, 2)
  expr <- a * cvec  ## equivalent to a/(1/c) when c is constant
  prob <- Problem(Minimize(norm_inf(expr)), list(a == 5))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 10, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_multiplication_on_left
test_that("left multiply: c' A c quadratic form", {
  A <- Variable(c(2, 2))
  cvec <- matrix(c(1, 2), 2, 1)
  prob <- Problem(Minimize(t(cvec) %*% A %*% cvec), list(A >= 2))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 18, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_multiplication_on_left
test_that("left multiply: scalar * variable", {
  a <- Variable(1)
  prob <- Problem(Minimize(a * 2), list(a >= 2))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 4, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_multiplication_on_left
test_that("left multiply: t(x) %*% c", {
  x <- Variable(2)
  cvec <- matrix(c(1, 2), 2, 1)
  prob <- Problem(Minimize(t(x) %*% cvec), list(x >= 2))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 6, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_multiplication_on_left
test_that("left multiply: (t(x) + t(z)) %*% c", {
  x <- Variable(2)
  z <- Variable(2)
  cvec <- matrix(c(1, 2), 2, 1)
  prob <- Problem(Minimize((t(x) + t(z)) %*% cvec),
                  list(x >= 2, z >= 1))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 9, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_vector_lp
test_that("vector LP: basic", {
  x <- Variable(2)
  cvec <- matrix(c(1, 2), 2, 1)
  prob <- Problem(Minimize(t(cvec) %*% x), list(x >= cvec))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 5, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_vector_lp
test_that("vector LP: complex with matrix constraints", {
  x <- Variable(2)
  a <- Variable(1)
  z <- Variable(2)
  cvec <- matrix(c(1, 2), 2, 1)
  ## A.T in Python = t(A) in R. Python A = [[3,5],[1,2]], A.T = [[3,1],[5,2]]
  A <- matrix(c(3, 1, 5, 2), 2, 2)  ## R column-major: same as Python A.T
  Imat <- diag(2)
  prob <- Problem(Minimize(t(cvec) %*% x + a),
                  list(A %*% x >= c(-1, 1),
                       4 * Imat %*% z == x,
                       z >= c(2, 2),
                       a >= 2))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 26, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(8, 8), tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# Validation & error paths
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_bad_objective
test_that("bad objective: expression not Minimize/Maximize", {
  x <- Variable(2)
  expect_error(Problem(x + 2))
})

## @cvxpy test_problem.py::TestProblem::test_invalid_constr
test_that("invalid constraint: expression in constraint list", {
  x <- Variable(1)
  expect_error(Problem(Minimize(x), list(sum_entries(x))))
})

## @cvxpy test_problem.py::TestProblem::test_invalid_solvers
test_that("invalid solver for problem type: OSQP rejects SOCP", {
  skip_if_not_installed("osqp")
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))
  expect_error(psolve(prob, solver = "OSQP"))
})

## @cvxpy test_problem.py::TestProblem::test_invalid_solvers
test_that("invalid solver for problem type: HiGHS rejects SDP", {
  skip_if_not_installed("highs")
  A <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(A)), list(A >= 1))
  expect_error(psolve(prob, solver = "HIGHS"))
})

## @cvxpy test_problem.py::TestProblem::test_redundant_constraints
test_that("redundant constraints: duplicate constraints", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)),
                  list(x == 2, x == 2, t(x) == 2, x[1] == 2))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 4, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_redundant_constraints
test_that("redundant constraints: tautology x == x", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(square(x))), list(x == x))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 0, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_mult_by_zero
test_that("multiply by zero: 0 * a produces 0", {
  a <- Variable(1)
  value(a) <- 1
  expr <- 0 * a
  expect_equal(as.numeric(value(expr)), 0)
  prob <- Problem(Minimize(expr))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 0, tolerance = 1e-6)
})

# ══════════════════════════════════════════════════════════════════
# SDP-specific tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_sdp_symmetry
test_that("SDP symmetry: lambda_max forces symmetric solution", {
  A <- Variable(c(2, 2))
  prob <- Problem(Minimize(lambda_max(A)), list(A >= 2))
  psolve(prob, solver = "SCS")
  A_val <- value(A)
  expect_equal(A_val, t(A_val), tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_sdp_symmetry
test_that("SDP symmetry: asymmetric constraint makes infeasible", {
  A <- Variable(c(2, 2))
  prob <- Problem(Minimize(lambda_max(A)),
                  list(A == matrix(c(1, 3, 2, 4), 2, 2)))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "infeasible")
})

## @cvxpy test_problem.py::TestProblem::test_psd_duals
test_that("PSD duals: basic", {
  C <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Maximize(C[1, 1]),
                  list(PSD(matrix(c(2, 0, 0, 2), 2, 2) - C)))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 2, tolerance = 1e-2)
})

## @cvxpy test_problem.py::TestProblem::test_psd_duals
test_that("PSD duals: off-diagonal", {
  C <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Maximize(C[1, 2] + C[2, 1]),
                  list(PSD(matrix(c(2, 0, 0, 2), 2, 2) - C), C >= 0))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 4, tolerance = 1e-2)
})

## @cvxpy test_problem.py::TestProblem::test_psd_constraints
test_that("PSD constraints: correlation matrix", {
  C <- Variable(c(3, 3))
  prob <- Problem(Maximize(C[1, 3]),
                  list(DiagMat(C) == 1,
                       C[1, 2] == 0.6,
                       C[2, 3] == -0.3,
                       C == t(C),
                       PSD(C)))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 0.583, tolerance = 1e-2)
})

## @cvxpy test_problem.py::TestProblem::test_psd_constraints
test_that("PSD constraints: infeasible", {
  C <- Variable(c(2, 2))
  ## C == 1 (all ones) cannot be PSD >> [[2,0],[0,2]]
  prob <- Problem(Maximize(C[1, 2]),
                  list(C == 1,
                       PSD(C - matrix(c(2, 0, 0, 2), 2, 2))))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "infeasible")
})

## @cvxpy test_problem.py::TestProblem::test_psd_constraints
test_that("PSD constraints: unbounded NSD", {
  C <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(C[1, 1]),
                  list(PSD(matrix(c(2, 0, 0, 2), 2, 2) - C)))
  psolve(prob, solver = "SCS")
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
})

# ══════════════════════════════════════════════════════════════════
# Solve patterns
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_cumsum_axis
test_that("cumsum axis: axis=0 on row matrix", {
  n <- 5L
  x1 <- Variable(c(1L, n))
  expr1 <- cumsum(x1)  ## default axis=0
  prob1 <- Problem(Minimize(0), list(expr1 == 1))
  psolve(prob1)
  expect_equal(as.numeric(value(expr1)), rep(1, n), tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_cummax_axis
test_that("cummax axis: axis=0 on row matrix", {
  n <- 5L
  x1 <- Variable(c(1L, n))
  expr1 <- cummax_expr(x1, axis = 2L)
  prob1 <- Problem(Maximize(sum_entries(x1)), list(expr1 <= 1))
  psolve(prob1)
  expect_equal(as.numeric(value(expr1)), rep(1, n), tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_cummax_axis
test_that("cummax axis: axis=1 on column matrix", {
  n <- 5L
  x2 <- Variable(c(n, 1L))
  expr2 <- cummax_expr(x2, axis = 1L)
  prob2 <- Problem(Maximize(sum_entries(x2)), list(expr2 <= 1))
  psolve(prob2)
  expect_equal(as.numeric(value(expr2)), rep(1, n), tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_pnorm_axis
test_that("pnorm axis=1 with binding constraints", {
  X <- Variable(c(2, 10))
  b <- c(0, 1)
  expr <- p_norm(X, 2, axis = 1L) - b
  prob <- Problem(Maximize(sum_entries(X)), list(expr <= 0))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(expr)), c(0, 0), tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_pnorm_concave
test_that("pnorm concave: negative entry forces domain constraint", {
  x <- Variable(3)
  a <- c(-1, 2, 3)
  ## p=0.5 is concave, domain requires x >= 0
  prob <- Problem(Minimize(sum_entries(abs(x - a))), list(p_norm(x, 0.5) >= 0))
  result <- psolve(prob, solver = "CLARABEL")
  ## x can't go negative with p < 1 domain; best is x = [0, 2, 3], obj = 1
  expect_equal(result, 1, tolerance = 1e-3)
})
## @cvxpy test_problem.py::TestProblem::test_pnorm_concave

test_that("pnorm concave: all positive allows exact match", {
  x <- Variable(3)
  a <- c(1, 2, 3)
  prob <- Problem(Minimize(sum_entries(abs(x - a))), list(p_norm(x, 0.5) >= 0))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 0, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_mixed_atoms
test_that("mixed atoms: norm1 + norm_inf + pnorm in one problem", {
  x <- Variable(2)
  z <- Variable(2)
  prob <- Problem(
    Minimize(p_norm(5 + norm1(z) + norm1(x) + norm_inf(x - z), 2)),
    list(x >= c(2, 3), z <= c(-1, -4), p_norm(x + z, 2) <= 2))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, 22, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(2, 3), tolerance = 1e-3)
  expect_equal(as.numeric(value(z)), c(-1, -4), tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_diag_offset_problem
test_that("diag offset: DiagMat extracts main diagonal (k=0)", {
  ## C++ backend (get_diag_matrix_mat) only supports k=0
  ## Off-diagonal k != 0 causes segfault (wrong coefficient matrix dimensions)
  n <- 4L
  A <- matrix(0:(n^2 - 1), n, n, byrow = TRUE)
  d <- diag(A)  # main diagonal: 0, 5, 10, 15

  X <- Variable(c(n, n), nonneg = TRUE)
  prob <- Problem(Minimize(sum_entries(X)), list(DiagMat(X) == d))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, sum(d), tolerance = 1e-2)
  expect_equal(as.numeric(value(DiagMat(X))), d, tolerance = 1e-2)
})

## @cvxpy test_problem.py::TestProblem::test_diag_offset_problem
test_that("diag offset: DiagVec creates diagonal matrix (k=0)", {
  ## C++ backend (get_diag_vec_mat) only supports k=0
  n <- 4L
  d <- c(0, 5, 10, 15)

  x <- Variable(n)
  target <- diag(d)
  prob <- Problem(Minimize(sum_entries(x)), list(DiagVec(x) == target))
  result <- psolve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result, sum(d), tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), d, tolerance = 1e-2)
})

## @cvxpy test_problem.py::TestProblem::test_diag_offset_problem
test_that("diag offset: DiagVec k=1 (super-diagonal)", {
  ## CVXPY: x=[1,2,3], diag(x,k=1) -> 4x4 matrix, sum(x)=6
  x <- Variable(3)
  M <- matrix(0, 4, 4)
  M[1, 2] <- 1; M[2, 3] <- 2; M[3, 4] <- 3
  prob <- Problem(Minimize(sum_entries(x)), list(DiagVec(x, k = 1L) == M))
  psolve(prob, solver = "SCS", eps = 1e-8)
  expect_equal(value(prob), 6, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 2, 3), tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_diag_offset_problem
test_that("diag offset: DiagVec k=-1 (sub-diagonal)", {
  ## CVXPY: x=[10,20,30], diag(x,k=-1) -> 4x4 matrix, sum(x)=60
  x <- Variable(3)
  M <- matrix(0, 4, 4)
  M[2, 1] <- 10; M[3, 2] <- 20; M[4, 3] <- 30
  prob <- Problem(Minimize(sum_entries(x)), list(DiagVec(x, k = -1L) == M))
  psolve(prob, solver = "SCS", eps = 1e-8)
  expect_equal(value(prob), 60, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(10, 20, 30), tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_diag_offset_problem
test_that("diag offset: DiagMat k=1 (extract super-diagonal)", {
  ## CVXPY: target = matrix(1:16, 4, 4), diag(target, k=1) = [5, 10, 15], sum=30
  X <- Variable(c(4, 4))
  target <- matrix(1:16, 4, 4)
  prob <- Problem(Minimize(sum_entries(DiagMat(X, k = 1L))), list(X == target))
  psolve(prob, solver = "SCS", eps = 1e-8)
  expect_equal(value(prob), 30, tolerance = 1e-4)
  expect_equal(as.numeric(value(DiagMat(X, k = 1L))), c(5, 10, 15), tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_diag_offset_problem
test_that("diag offset: DiagMat k=-1 (extract sub-diagonal)", {
  ## CVXPY: target = matrix(1:16, 4, 4), diag(target, k=-1) = [2, 7, 12], sum=21
  X <- Variable(c(4, 4))
  target <- matrix(1:16, 4, 4)
  prob <- Problem(Minimize(sum_entries(DiagMat(X, k = -1L))), list(X == target))
  psolve(prob, solver = "SCS", eps = 1e-8)
  expect_equal(value(prob), 21, tolerance = 1e-4)
  expect_equal(as.numeric(value(DiagMat(X, k = -1L))), c(2, 7, 12), tolerance = 1e-4)
})
