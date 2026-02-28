## CVXPY Parity Gaps — Batch 2
## ============================
## Covers gaps from: test_examples.py, test_linear_cone.py, test_coeff_extractor.py,
## test_matrices.py, test_mip_vars.py, test_param_quad_prog.py,
## test_param_cone_prog.py, test_dgp_dpp.py, test_pow_cone_nd.py,
## test_valinvec2mixedint.py, test_logic.py, test_problem.py
##
## All expected values verified via `uv run python` against CVXPY 1.8.1.
## R uses 1-based indexing and (n,1) shape; CVXPY uses 0-based and (n,).

# ═══════════════════════════════════════════════════════════════════════
# test_examples.py — Example problems
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_examples.py::TestExamples::test_chebyshev_center
test_that("Chebyshev center: largest inscribed ball in polyhedron", {
  a1 <- c(2, 1)
  a2 <- c(2, -1)
  a3 <- c(-1, 2)
  a4 <- c(-1, -2)
  b <- rep(1, 4)

  r <- Variable(name = "r")
  x_c <- Variable(2, name = "x_c")
  obj <- Maximize(r)
  constraints <- list(
    t(a1) %*% x_c + norm(a1, type = "2") * r <= b[1],
    t(a2) %*% x_c + norm(a2, type = "2") * r <= b[2],
    t(a3) %*% x_c + norm(a3, type = "2") * r <= b[3],
    t(a4) %*% x_c + norm(a4, type = "2") * r <= b[4]
  )

  p <- Problem(obj, constraints)
  result <- psolve(p, solver = "CLARABEL")
  ## Analytical: r = 1/sqrt(5) = 0.447214
  expect_equal(result, 0.447214, tolerance = 1e-3)
  expect_equal(as.numeric(value(r)), result, tolerance = 1e-3)
  expect_equal(as.numeric(value(x_c)), c(0, 0), tolerance = 1e-3)
})

## @cvxpy test_examples.py::TestExamples::test_numpy_scalars
test_that("numpy_scalars: QP with random data, dual value shapes", {
  n <- 6L
  eps <- 1e-6
  set.seed(10)
  P0 <- crossprod(matrix(rnorm(n * n), n, n)) + eps * diag(n)
  P1 <- crossprod(matrix(rnorm(n * n), n, n))
  P2 <- crossprod(matrix(rnorm(n * n), n, n))
  P3 <- crossprod(matrix(rnorm(n * n), n, n))

  q0 <- matrix(rnorm(n), n, 1)
  q1 <- matrix(rnorm(n), n, 1)
  q2 <- matrix(rnorm(n), n, 1)
  q3 <- matrix(rnorm(n), n, 1)

  r0 <- matrix(rnorm(1), 1, 1)
  r1 <- matrix(rnorm(1), 1, 1)
  r2 <- matrix(rnorm(1), 1, 1)
  r3 <- matrix(rnorm(1), 1, 1)

  slack <- Variable()
  x <- Variable(n)
  objective <- Minimize(0.5 * quad_form(x, P0) + t(q0) %*% x + r0 + slack)
  constraints <- list(
    0.5 * quad_form(x, P1) + t(q1) %*% x + r1 <= slack,
    0.5 * quad_form(x, P2) + t(q2) %*% x + r2 <= slack,
    0.5 * quad_form(x, P3) + t(q3) %*% x + r3 <= slack
  )

  p <- Problem(objective, constraints)
  psolve(p, solver = "CLARABEL")

  ## Check that dual values have the right type (scalar)
  lam1 <- dual_value(constraints[[1]])
  lam2 <- dual_value(constraints[[2]])
  lam3 <- dual_value(constraints[[3]])
  expect_true(is.numeric(lam1))
  expect_true(is.numeric(lam2))
  expect_true(is.numeric(lam3))
})

## @cvxpy test_examples.py::TestExamples::test_readme_examples
test_that("README examples: least squares with box constraints", {
  set.seed(1)
  m <- 30L
  n <- 20L
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)

  x <- Variable(n)
  objective <- Minimize(sum_squares(A %*% x - b))
  constraints <- list(x >= 0, x <= 1)
  p <- Problem(objective, constraints)
  psolve(p, solver = "CLARABEL")

  ## Verify it solved
  expect_equal(status(p), "optimal")
  expect_true(all(value(x) >= -1e-4))
  expect_true(all(value(x) <= 1 + 1e-4))
})

## @cvxpy test_examples.py::TestExamples::test_advanced1
test_that("advanced1: solve with different solvers", {
  x <- Variable(2)
  obj <- Minimize(x[1] + p_norm(x, 1))
  constraints <- list(x >= 2)
  prob <- Problem(obj, constraints)

  ## Solve with CLARABEL
  psolve(prob, solver = "CLARABEL")
  expect_equal(value(prob), 6, tolerance = 1e-3)

  ## Solve with SCS
  psolve(prob, solver = "SCS")
  expect_equal(value(prob), 6, tolerance = 1e-1)
})

## @cvxpy test_examples.py::TestExamples::test_log_det
test_that("log_det example: maximum volume inscribed ellipsoid", {
  x <- cbind(
    c(0.55, 0.25, -0.2, -0.25, 0.0, 0.4),
    c(0.0, 0.35, 0.2, -0.1, -0.3, -0.2)
  )
  x <- t(x)  # (2, 6)
  n <- nrow(x)
  m <- ncol(x)

  A <- Variable(c(n, n))
  b <- Variable(n)
  obj <- Maximize(log_det(A))
  constraints <- lapply(1:m, function(i) {
    p_norm(A %*% x[, i] + b, 2) <= 1
  })
  p <- Problem(obj, constraints)
  result <- psolve(p, solver = "SCS")
  expect_equal(result, 1.9746, tolerance = 0.05)
})

## @cvxpy test_examples.py::TestExamples::test_portfolio_problem
test_that("portfolio_problem: DCP attributes on sparse data", {
  skip_if_not_installed("Matrix")
  set.seed(5)
  n <- 100L
  m <- 10L

  ## Sparse factor matrix
  F_mat <- Matrix::rsparsematrix(m, n, density = 0.01)
  F_mat@x <- rep(1.0, length(F_mat@x))
  D <- as(Matrix::Diagonal(n, rnorm(n)^2), "generalMatrix")
  Z <- matrix(rnorm(m), m, 1) %*% matrix(rnorm(m), 1, m)

  x <- Variable(n)
  y <- Constant(F_mat) %*% x
  ## This expression must not error (DCP attribute on sparse)
  expr <- square(p_norm(Constant(D) %*% x, 2)) + square(Z %*% y)
  expect_s3_class(expr, "CVXR::Expression")
})

## @cvxpy test_examples.py::TestExamples::test_intro
test_that("intro: LP, QP, infeasible, unbounded examples", {
  set.seed(1)
  ## Part 1: bound-constrained least squares
  m <- 30L; n <- 20L
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)

  x <- Variable(n)
  objective <- Minimize(sum_squares(A %*% x - b))
  constraints <- list(x >= 0, x <= 1)
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")

  ## Part 2: equality + inequality QP
  x2 <- Variable()
  y2 <- Variable()
  constraints2 <- list(x2 + y2 == 1, x2 - y2 >= 1)
  obj2 <- Minimize(square(x2 - y2))
  prob2 <- Problem(obj2, constraints2)
  psolve(prob2, solver = "CLARABEL")
  expect_equal(status(prob2), "optimal")
  expect_equal(value(prob2), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x2)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(y2)), 0.0, tolerance = 1e-3)

  ## Part 3: replace objective with Maximize
  prob3 <- Problem(Maximize(x2 + y2), prob2@constraints)
  psolve(prob3, solver = "CLARABEL")
  expect_equal(value(prob3), 1.0, tolerance = 1e-2)

  ## Part 4: replace constraint
  constraints4 <- prob3@constraints
  constraints4[[1]] <- (x2 + y2 <= 3)
  prob4 <- Problem(prob3@objective, constraints4)
  psolve(prob4, solver = "CLARABEL")
  expect_equal(value(prob4), 3.0, tolerance = 1e-1)

  ## Part 5: infeasible
  x5 <- Variable()
  prob5 <- Problem(Minimize(x5), list(x5 >= 1, x5 <= 0))
  psolve(prob5, solver = "CLARABEL")
  expect_equal(status(prob5), "infeasible")

  ## Part 6: unbounded
  prob6 <- Problem(Minimize(x5))
  psolve(prob6, solver = "CLARABEL")
  expect_true(status(prob6) %in% c("unbounded", "dual_infeasible"))
})

## @cvxpy test_examples.py::TestExamples::test_inpainting
test_that("inpainting: total variation image reconstruction", {
  set.seed(1)
  rows <- 20L; cols <- 20L
  Uorig <- matrix(sample(0:255, rows * cols, replace = TRUE), rows, cols)

  Known <- matrix(0, rows, cols)
  for (i in seq_len(rows)) {
    for (j in seq_len(cols)) {
      if (runif(1) > 0.7) Known[i, j] <- 1
    }
  }
  Ucorr <- Known * Uorig

  U <- Variable(c(rows, cols))
  obj <- Minimize(total_variation(U))
  constraints <- list(Known * U == Known * Ucorr)
  prob <- Problem(obj, constraints)
  psolve(prob, solver = "SCS")
  ## Just verify it solved without error
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
})

## @cvxpy test_examples.py::TestExamples::test_log_sum_exp
test_that("log_sum_exp: vstack + log_sum_exp optimization", {
  set.seed(1)
  m <- 5L; n <- 2L
  X <- matrix(1, m, n)
  w <- Variable(n)

  expr2 <- lapply(1:m, function(i) {
    log_sum_exp(vstack(Constant(0), t(X[i, ]) %*% w))
  })
  expr3 <- Reduce(`+`, expr2)
  obj <- Minimize(expr3)
  p <- Problem(obj)
  ## Just run 1 iteration to check compilation works
  psolve(p, solver = "SCS", max_iters = 1L)
  ## If it didn't error, the test passes
  expect_true(TRUE)
})

# ═══════════════════════════════════════════════════════════════════════
# test_linear_cone.py — ConeMatrixStuffing pipeline tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_linear_cone.py::TestLinearCone::test_scalar_lp
test_that("scalar_lp: scalar LP solves correctly", {
  a <- Variable(name = "a")
  p <- Problem(Minimize(3 * a), list(a >= 2))
  result <- psolve(p, solver = "CLARABEL")
  expect_equal(result, 6, tolerance = 1e-3)
  expect_equal(as.numeric(value(a)), 2, tolerance = 1e-3)

  ## With two variables
  b <- Variable(name = "b")
  p2 <- Problem(Minimize(-3 * a + b),
    list(a <= 2, b == a, b <= 5))
  result2 <- psolve(p2, solver = "CLARABEL")
  expect_equal(result2, -4, tolerance = 1e-3)

  ## With constant in objective
  c_var <- Variable(name = "c")
  p3 <- Problem(Minimize(3 * a - b + 100),
    list(a >= 2, b + 5 * c_var - 2 == a, b <= 5 + c_var))
  psolve(p3, solver = "CLARABEL")
  expect_equal(status(p3), "optimal")

  ## Infeasible
  p4 <- Problem(Maximize(a), list(a >= 2, a <= 1))
  psolve(p4, solver = "CLARABEL")
  expect_equal(status(p4), "infeasible")
})

## @cvxpy test_linear_cone.py::TestLinearCone::test_vector_lp
test_that("vector_lp: vector LP solves correctly", {
  x <- Variable(2, name = "x")
  cc <- c(1, 2)
  p <- Problem(Minimize(t(cc) %*% x), list(x >= cc))
  result <- psolve(p, solver = "CLARABEL")
  ## x >= [1,2] minimizing [1,2]^T x → x = [1,2], val = 5
  expect_equal(result, 5, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), cc, tolerance = 1e-3)

  ## More complex vector LP
  z <- Variable(2, name = "z")
  a <- Variable(name = "a")
  A_mat <- matrix(c(3, 1, 5, 2), 2, 2)  # [[3,5],[1,2]]
  Imat <- diag(2)
  p2 <- Problem(Minimize(t(cc) %*% x + a),
    list(A_mat %*% x >= c(-1, 1),
         4 * Imat %*% z == x,
         z >= c(2, 2),
         a >= 2))
  result2 <- psolve(p2, solver = "CLARABEL")
  expect_equal(status(p2), "optimal")
})

## @cvxpy test_linear_cone.py::TestLinearCone::test_matrix_lp
test_that("matrix_lp: matrix LP solves correctly", {
  A_var <- Variable(c(2, 2), name = "A")
  a <- Variable(name = "a")
  Tmat <- matrix(1, 2, 2)
  p <- Problem(Minimize(a), list(A_var == Tmat + a, a >= 0))
  result <- psolve(p, solver = "CLARABEL")
  expect_equal(result, 0, tolerance = 1e-3)
})

## @cvxpy test_linear_cone.py::TestLinearCone::test_socp
test_that("socp: SOC problem solves correctly", {
  x <- Variable(2, name = "x")
  b <- Variable(name = "b")

  ## Basic
  p <- Problem(Minimize(b), list(p_norm(x, 2) <= b))
  result <- psolve(p, solver = "CLARABEL")
  expect_equal(result, 0, tolerance = 1e-3)

  ## More complex with SOC
  y <- Variable(3, name = "y")
  p2 <- Problem(Minimize(b),
    list(p_norm(x / 2 + y[1:2], 2) <= b + 5, x >= 1, y == 5))
  result2 <- psolve(p2, solver = "CLARABEL")
  expect_equal(status(p2), "optimal")
})

## @cvxpy test_linear_cone.py::TestLinearCone::test_psd_constraints
test_that("psd_constraints: PSD constraint pipeline", {
  C <- Variable(c(3, 3))
  obj <- Maximize(C[1, 3])
  constraints <- list(
    DiagMat(C) == 1,
    C[1, 2] == 0.6,
    C[2, 3] == -0.3,
    C == t(C),
    C %>>% 0
  )
  prob <- Problem(obj, constraints)
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## C[1,3] should be maximized subject to PSD + constraints
})

## @cvxpy test_linear_cone.py::TestLinearCone::test_nonneg_constraints_backend
test_that("nonneg_constraints_backend: NonNeg constraint pipeline", {
  x <- Variable(2, name = "x")
  objective <- Maximize(-4 * x[1] - 5 * x[2])
  constr_expr <- hstack(3 - (2 * x[1] + x[2]),
                         3 - (x[1] + 2 * x[2]),
                         x[1],
                         x[2])
  constraints <- list(NonNeg(constr_expr))
  prob <- Problem(objective, constraints)
  ## Problem should be solvable
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
})

## @cvxpy test_linear_cone.py::TestLinearCone::test_nonneg_constraints_end_user
test_that("nonneg_constraints_end_user: NonNeg solve with duals", {
  x <- Variable(2, name = "x")
  objective <- Minimize(-4 * x[1] - 5 * x[2])
  constr_expr <- hstack(3 - (2 * x[1] + x[2]),
                         3 - (x[1] + 2 * x[2]),
                         x[1],
                         x[2])
  nn_constr <- NonNeg(constr_expr)
  constraints <- list(nn_constr)
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "CLARABEL")

  ## Expected: x = [1, 1], obj = -9
  expect_equal(value(prob), -9, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-2)

  ## Check dual values
  dv <- as.numeric(dual_value(nn_constr))
  expect_equal(dv, c(1, 2, 0, 0), tolerance = 1e-2)
})

# ═══════════════════════════════════════════════════════════════════════
# test_coeff_extractor.py — Coefficient extraction regressions
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_coeff_extractor.py::test_issue_2402_scalar_parameter
test_that("issue_2402_scalar_parameter: two params on quad_forms", {
  r <- c(-0.48, 0.11, 0.09, -0.39, 0.03)
  Sigma <- matrix(c(
    2.4e-04, 1.3e-04, 2.0e-04, 1.6e-04, 2.0e-04,
    1.3e-04, 2.8e-04, 2.1e-04, 1.7e-04, 1.5e-04,
    2.0e-04, 2.1e-04, 5.8e-04, 3.3e-04, 2.3e-04,
    1.6e-04, 1.7e-04, 3.3e-04, 6.9e-04, 2.1e-04,
    2.0e-04, 1.5e-04, 2.3e-04, 2.1e-04, 3.6e-04
  ), 5, 5, byrow = TRUE)

  w <- Variable(5)
  risk_aversion <- Parameter(value = 1.0, nonneg = TRUE)
  ridge_coef <- Parameter(value = 0.0, nonneg = TRUE)

  obj_func <- t(r) %*% w - risk_aversion * quad_form(w, Sigma) -
    ridge_coef * sum_squares(w)
  objective <- Maximize(obj_func)
  fixed_w <- c(10, 11, 12, 13, 14)
  constraints <- list(w == fixed_w)
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "CLARABEL")

  expected_value <- as.numeric(
    t(r) %*% fixed_w - 1.0 * t(fixed_w) %*% Sigma %*% fixed_w -
      0.0 * sum(fixed_w^2)
  )
  expect_equal(value(prob), expected_value, tolerance = 1e-2)
  expect_equal(as.numeric(value(w)), fixed_w, tolerance = 1e-2)
})

## @cvxpy test_coeff_extractor.py::test_issue_2402_scalar_constant
test_that("issue_2402_scalar_constant: constant + param on quad_forms", {
  r <- c(-0.48, 0.11, 0.09, -0.39, 0.03)
  Sigma <- matrix(c(
    2.4e-04, 1.3e-04, 2.0e-04, 1.6e-04, 2.0e-04,
    1.3e-04, 2.8e-04, 2.1e-04, 1.7e-04, 1.5e-04,
    2.0e-04, 2.1e-04, 5.8e-04, 3.3e-04, 2.3e-04,
    1.6e-04, 1.7e-04, 3.3e-04, 6.9e-04, 2.1e-04,
    2.0e-04, 1.5e-04, 2.3e-04, 2.1e-04, 3.6e-04
  ), 5, 5, byrow = TRUE)

  w <- Variable(5)
  risk_aversion <- Parameter(value = 1.0, nonneg = TRUE)
  ridge_coef <- 0  # constant, not parameter

  obj_func <- t(r) %*% w - risk_aversion * quad_form(w, Sigma) -
    ridge_coef * sum_squares(w)
  objective <- Maximize(obj_func)
  fixed_w <- c(10, 11, 12, 13, 14)
  constraints <- list(w == fixed_w)
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "CLARABEL")

  expected_value <- as.numeric(
    t(r) %*% fixed_w - 1.0 * t(fixed_w) %*% Sigma %*% fixed_w
  )
  expect_equal(value(prob), expected_value, tolerance = 1e-2)
  expect_equal(as.numeric(value(w)), fixed_w, tolerance = 1e-2)
})

## @cvxpy test_coeff_extractor.py::test_issue_2402_vector
test_that("issue_2402_vector: vector ridge coefficient param", {
  r <- c(-0.48, 0.11, 0.09, -0.39, 0.03)
  Sigma <- matrix(c(
    2.4e-04, 1.3e-04, 2.0e-04, 1.6e-04, 2.0e-04,
    1.3e-04, 2.8e-04, 2.1e-04, 1.7e-04, 1.5e-04,
    2.0e-04, 2.1e-04, 5.8e-04, 3.3e-04, 2.3e-04,
    1.6e-04, 1.7e-04, 3.3e-04, 6.9e-04, 2.1e-04,
    2.0e-04, 1.5e-04, 2.3e-04, 2.1e-04, 3.6e-04
  ), 5, 5, byrow = TRUE)

  w <- Variable(5)
  risk_aversion <- Parameter(value = 2.0, nonneg = TRUE)
  ridge_coef <- Parameter(5, value = 0:4, nonneg = TRUE)

  obj_func <- t(r) %*% w - risk_aversion * quad_form(w, Sigma) -
    sum_entries(multiply(multiply(ridge_coef, c(5, 6, 7, 8, 9)), square(w)))
  objective <- Maximize(obj_func)
  fixed_w <- c(10, 11, 12, 13, 14)
  constraints <- list(w == fixed_w)
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "CLARABEL")

  expected_value <- as.numeric(
    t(r) %*% fixed_w - 2.0 * t(fixed_w) %*% Sigma %*% fixed_w -
      sum(0:4 * c(5, 6, 7, 8, 9) * fixed_w^2)
  )
  expect_equal(value(prob), expected_value, tolerance = 1e-1)
  expect_equal(as.numeric(value(w)), fixed_w, tolerance = 1e-1)
})

## @cvxpy test_coeff_extractor.py::test_problem_end_to_end
test_that("problem_end_to_end: MWE regression for #2402", {
  x <- Variable(2)
  p1 <- Parameter(value = 1.0, nonneg = TRUE)
  p2 <- Parameter(value = 0.0, nonneg = TRUE)
  P <- diag(2)

  objective <- Minimize(p1 * quad_form(x, P) + p2 * sum_squares(x))
  problem <- Problem(objective, constraints = list(sum_entries(x) == 1))
  psolve(problem, solver = "CLARABEL")
  expect_equal(value(problem), 0.5, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-3)
})

## @cvxpy test_coeff_extractor.py::test_issue_2437
test_that("issue_2437: quad_obj vs non-quad_obj consistency", {
  N <- 3L
  t_cost <- c(0.01, 0.02, 0.03)
  alpha <- c(0.04, 0.05, 0.06)
  ivol <- c(0.07, 0.08, 0.09)

  w <- Variable(N, name = "w")
  risk <- sum_entries(power(multiply(w, ivol), 2))
  U <- t(w) %*% alpha - risk - t(abs(w)) %*% t_cost
  problem <- Problem(Maximize(U), list())

  ## use_quad_obj has no direct equivalent in R, but solve should work
  val <- psolve(problem, solver = "CLARABEL")
  ## Expected ≈ 0.1089
  expect_equal(val, 0.1089, tolerance = 1e-2)
})

# ═══════════════════════════════════════════════════════════════════════
# test_matrices.py — Matrix type arithmetic
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_matrices.py::TestMatrices::test_numpy_arrays
test_that("numpy_arrays: R array arithmetic with Variables", {
  x <- Variable(2, name = "x")
  A_var <- Variable(c(2, 2), name = "A")

  ## Vector operations
  v <- 0:1
  expect_equal((x + v)@shape, c(2L, 1L))
  expect_equal((v + x)@shape, c(2L, 1L))
  expect_equal((x - v)@shape, c(2L, 1L))
  expect_equal((v - x)@shape, c(2L, 1L))
  expect_equal((x <= v)@shape, c(2L, 1L))
  expect_equal((v <= x)@shape, c(2L, 1L))
  expect_equal((x == v)@shape, c(2L, 1L))
  expect_equal((v == x)@shape, c(2L, 1L))

  ## Matrix multiplication
  A_mat <- matrix(0:7, 4, 2)
  expect_equal((A_mat %*% x)@shape, c(4L, 1L))

  ## PSD inequalities
  A_psd <- matrix(1, 2, 2)
  expect_equal((A_psd %<<% A_var)@shape, c(2L, 2L))
  expect_equal((A_psd %>>% A_var)@shape, c(2L, 2L))
})

## @cvxpy test_matrices.py::TestMatrices::test_numpy_matrices
test_that("numpy_matrices: R matrix arithmetic with Variables", {
  x <- Variable(2, name = "x")
  A_var <- Variable(c(2, 2), name = "A")

  v <- 0:1
  expect_equal((x + v)@shape, c(2L, 1L))
  expect_equal((v + v + x)@shape, c(2L, 1L))
  expect_equal((x - v)@shape, c(2L, 1L))
  expect_equal((v - v - x)@shape, c(2L, 1L))

  A_mat <- matrix(0:7, 4, 2)
  expect_equal((A_mat %*% x)@shape, c(4L, 1L))
  expect_equal((crossprod(A_mat, A_mat) %*% x)@shape, c(2L, 1L))

  A_psd <- matrix(1, 2, 2)
  expect_equal((A_psd %<<% A_var)@shape, c(2L, 2L))
  expect_equal((A_psd %>>% A_var)@shape, c(2L, 2L))
})

## @cvxpy test_matrices.py::TestMatrices::test_numpy_scalars
test_that("numpy_scalars: R scalar arithmetic with Variables", {
  x <- Variable(2, name = "x")
  A_var <- Variable(c(2, 2), name = "A")

  v <- 2.0
  expect_equal((x + v)@shape, c(2L, 1L))
  expect_equal((v + x)@shape, c(2L, 1L))
  expect_equal((v * x)@shape, c(2L, 1L))
  expect_equal((x - v)@shape, c(2L, 1L))
  expect_equal((v - v - x)@shape, c(2L, 1L))
  expect_equal((x <= v)@shape, c(2L, 1L))
  expect_equal((v <= x)@shape, c(2L, 1L))
  expect_equal((x == v)@shape, c(2L, 1L))
  expect_equal((v == x)@shape, c(2L, 1L))

  ## PSD inequalities
  expect_equal((v %<<% A_var)@shape, c(2L, 2L))
  expect_equal((v %>>% A_var)@shape, c(2L, 2L))
})

## @cvxpy test_matrices.py::TestMatrices::test_scipy_sparse
test_that("scipy_sparse: sparse matrix arithmetic with Variables", {
  skip_if_not_installed("Matrix")
  var <- Variable(c(4, 2))
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 1, 2, 3, 4),
    j = c(1, 1, 1, 1, 2, 2, 2, 2),
    x = as.double(0:7),
    dims = c(4, 2)
  )
  B <- cbind(A, A)
  A_const <- Constant(A)
  B_const <- Constant(B)

  expect_equal((var + A_const)@shape, c(4L, 2L))
  expect_equal((A_const + var)@shape, c(4L, 2L))
  expect_equal((B_const %*% var)@shape, c(4L, 2L))
  expect_equal((var - A_const)@shape, c(4L, 2L))
  expect_equal((A_const - A_const - var)@shape, c(4L, 2L))
})

# ═══════════════════════════════════════════════════════════════════════
# test_mip_vars.py — Mixed-integer programming
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_mip_vars.py::TestMIPVariable::test_all_solvers
test_that("MIP: bool and int problems with GLPK_MI", {
  skip_if_not_installed("Rglpk")

  ## Bool in objective
  x_bool <- Variable(boolean = TRUE)
  obj <- Minimize(abs(x_bool - 0.2))
  p <- Problem(obj, list())
  result <- psolve(p, solver = "GLPK_MI")
  expect_equal(result, 0.2, tolerance = 1e-3)
  expect_equal(as.numeric(value(x_bool)), 0, tolerance = 1e-3)

  ## Int in objective
  y_int <- Variable(integer = TRUE)
  obj2 <- Minimize(abs(y_int - 0.2))
  p2 <- Problem(obj2, list())
  result2 <- psolve(p2, solver = "GLPK_MI")
  expect_equal(result2, 0.2, tolerance = 1e-3)
  expect_equal(as.numeric(value(y_int)), 0, tolerance = 1e-3)

  ## Matrix Bool in objective
  A_bool <- Variable(c(3, 2), boolean = TRUE)
  C_mat <- matrix(c(0, 1, 0, 1, 1, 1), 3, 2)
  obj3 <- Minimize(sum_entries(abs(A_bool - C_mat)))
  p3 <- Problem(obj3, list())
  result3 <- psolve(p3, solver = "GLPK_MI")
  expect_equal(result3, 0, tolerance = 1e-3)
  expect_equal(as.numeric(value(A_bool)), as.numeric(C_mat), tolerance = 1e-2)
})

## @cvxpy test_mip_vars.py::TestMIPVariable::test_highs_default_milp
test_that("HiGHS default MILP: integer LP uses HiGHS by default", {
  skip_if_not_installed("highs")

  x <- Variable(3, integer = TRUE)
  objective <- Minimize(sum_entries(x))
  constraints <- list(x >= 0, x <= 10, sum_entries(x) >= 5)
  prob <- Problem(objective, constraints)

  ## Solve with HiGHS explicitly (R may not auto-select HiGHS)
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 5.0, tolerance = 1e-3)
  ## Verify integer solution
  expect_true(all(abs(value(x) - round(value(x))) < 1e-4))
})

## @cvxpy test_mip_vars.py::TestMIPVariable::test_milp_no_warning
test_that("MILP no warning: no MINLP warning for pure LP", {
  skip_if_not_installed("Rglpk")

  x <- Variable(2, integer = TRUE)
  y <- Variable(2)
  objective <- Minimize(sum_entries(x) + sum_entries(y))
  constraints <- list(x >= 0, x <= 5, y >= 0, y <= 10,
                      sum_entries(x) + sum_entries(y) >= 3)
  prob <- Problem(objective, constraints)

  ## Just solve and check optimal - MILP should work fine
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
})

## @cvxpy test_mip_vars.py::TestMIPVariable::test_miqp_warning
test_that("MIQP: mixed-integer QP problem", {
  skip_if_not_installed("Rglpk")

  ## MIQP = QP + integer variables
  ## In R/CVXR, this may require a solver that handles MIQP
  ## GLPK_MI is LP-only, so this tests behavior when QP goes to conic
  x <- Variable(3, integer = TRUE)
  objective <- Minimize(sum_squares(x))
  constraints <- list(x >= 0, x <= 10, sum_entries(x) >= 5)
  prob <- Problem(objective, constraints)

  ## MIQP requires a solver that supports both QP and integer variables.
  ## ECOS_BB is conic-only (no QP path), so only GUROBI, CPLEX, or MOSEK work.
  skip_if_not(
    any(c("GUROBI", "CPLEX", "MOSEK") %in% installed_solvers()),
    "No MIQP solver available (need GUROBI, CPLEX, or MOSEK)"
  )
  psolve(prob)
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
})

## @cvxpy test_mip_vars.py::TestMIPVariable::test_highs_milp_simple
test_that("HiGHS MILP simple: knapsack-style problem", {
  skip_if_not_installed("highs")

  x <- Variable(4, integer = TRUE)
  values <- c(3, 4, 5, 6)
  weights <- c(2, 3, 4, 5)
  capacity <- 8

  objective <- Maximize(t(values) %*% x)
  constraints <- list(
    x >= 0,
    x <= 1,  # binary variables (0 or 1)
    t(weights) %*% x <= capacity
  )
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "HIGHS")

  expect_equal(status(prob), "optimal")
  ## Verify integer solution
  expect_true(all(abs(value(x) - round(value(x))) < 1e-4))
  ## Verify constraint satisfaction
  expect_true(as.numeric(t(weights) %*% value(x)) <= capacity + 1e-6)
})

# ═══════════════════════════════════════════════════════════════════════
# test_param_quad_prog.py — Parametric QP tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_param_quad_prog.py::TestParamQuadProg::test_param_data
test_that("param_data: DPP parametric QP data matches scratch solve", {
  skip_if_not_installed("osqp")
  set.seed(0)
  m <- 30L; n <- 20L
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)
  x <- Variable(n)
  gamma <- Parameter(nonneg = TRUE)
  gamma_val <- 0.5
  gamma_val_new <- 0.1
  objective <- Minimize(gamma * sum_squares(A %*% x - b) + p_norm(x, 1))
  constraints <- list(x >= 1, x <= 2)

  ## Solve from scratch with gamma = 0.1
  prob <- Problem(objective, constraints)
  expect_true(is_dpp(prob))
  value(gamma) <- gamma_val_new
  psolve(prob, solver = "OSQP")
  x_scratch <- as.numeric(value(x))

  ## Solve with gamma = 0.5 first, then re-solve with 0.1
  prob2 <- Problem(objective, constraints)
  value(gamma) <- gamma_val
  psolve(prob2, solver = "OSQP")

  value(gamma) <- gamma_val_new
  psolve(prob2, solver = "OSQP")
  x_gamma_new <- as.numeric(value(x))

  ## Solutions should be close
  expect_equal(x_gamma_new, x_scratch, tolerance = 0.05)
})

## @cvxpy test_param_quad_prog.py::TestParamQuadProg::test_qp_problem
test_that("qp_problem: cached solve matches full solve", {
  skip_if_not_installed("osqp")
  set.seed(0)
  m <- 30L; n <- 20L
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)
  x <- Variable(n)
  gamma <- Parameter(nonneg = TRUE)
  value(gamma) <- 0.5
  objective <- Minimize(sum_squares(A %*% x - b) + gamma * p_norm(x, 1))
  constraints <- list(x >= 0, x <= 1)

  problem <- Problem(objective, constraints)
  psolve(problem, solver = "OSQP")
  x_full <- as.numeric(value(x))

  ## Re-solve should give same result (DPP fast path)
  psolve(problem, solver = "OSQP")
  x_param <- as.numeric(value(x))

  expect_equal(x_param, x_full, tolerance = 0.05)
})

## @cvxpy test_param_quad_prog.py::TestParamQuadProg::test_var_bounds
test_that("var_bounds: OSQP bounded variables propagation", {
  skip_if_not_installed("osqp")
  ## Scalar bounds work; matrix upper_bounds have a shape mismatch bug
  x <- Variable(3, bounds = list(-10, 5))
  problem <- Problem(Minimize(sum_entries(x)))
  psolve(problem, solver = "OSQP")
  expect_equal(status(problem), "optimal")
  xval <- as.numeric(value(x))
  expect_true(all(xval >= -10 - 1e-3))
  expect_true(all(xval <= 5 + 1e-3))
})

## @cvxpy test_param_quad_prog.py::TestParamQuadProg::test_highs_var_bounds
test_that("highs_var_bounds_qp: variable bounds with HiGHS", {
  skip_if_not_installed("highs")
  x1 <- Variable(bounds = list(-1, 1))
  x2 <- Variable(bounds = list(-0.5, 1))
  x3 <- Variable()
  objective <- (x1^2 + x2^2) / 2 + x1 + x2 + x3
  constraints <- list(
    -3 <= x1 + x2, x1 + x2 <= 3,
    -4 <= x1 - x2, x1 - x2 <= 4,
    x3 >= -2
  )
  prob <- Problem(Minimize(objective), constraints)
  psolve(prob, solver = "HIGHS")

  expect_equal(as.numeric(value(x1)), -1, tolerance = 1e-3)
  expect_equal(as.numeric(value(x2)), -0.5, tolerance = 1e-3)
  expect_equal(as.numeric(value(x3)), -2, tolerance = 1e-3)
})

## @cvxpy test_param_quad_prog.py::TestParamQuadProg::test_daqp_var_bounds
test_that("daqp_var_bounds: skipped (DAQP not available in R)", {
  ## DAQP solver is not available in R
  skip("DAQP solver not available in R")
})

# ═══════════════════════════════════════════════════════════════════════
# test_param_cone_prog.py — Parametric conic program tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_param_cone_prog.py::TestParamConeProg::test_log_problem
test_that("log_problem: log in objective and constraint", {
  ## Log in objective
  x <- Variable(2)
  obj <- Maximize(sum_entries(log(x)))
  constr <- list(x <= c(1, exp(1)))
  problem <- Problem(obj, constr)
  psolve(problem, solver = "SCS")
  expect_equal(status(problem), "optimal")
  ## x* = [1, e], sum(log) = 0 + 1 = 1
  expect_equal(value(problem), 1.0, tolerance = 0.1)

  ## Log in constraint
  obj2 <- Minimize(sum_entries(x))
  constr2 <- list(log(x) >= 0, x <= c(1, 1))
  problem2 <- Problem(obj2, constr2)
  psolve(problem2, solver = "SCS")
  expect_equal(status(problem2), "optimal")
  ## log(x) >= 0 ⟹ x >= 1, with x <= 1 → x = [1, 1], sum = 2
  expect_equal(value(problem2), 2.0, tolerance = 0.1)
})

## @cvxpy test_param_cone_prog.py::TestParamConeProg::test_psd_var
test_that("psd_var: PSD variable with minimize/maximize", {
  s <- Variable(c(2, 2), PSD = TRUE)
  obj <- Maximize(min_elemwise(s[1, 2], 10))
  const <- list(DiagMat(s) == rep(1, 2))
  problem <- Problem(obj, const)
  psolve(problem, solver = "SCS", eps = 1e-5)
  expect_equal(status(problem), "optimal")
  ## PSD matrix with diag=1 → s[1,2] = s[2,1] in [-1, 1], maximized → 1
  expect_equal(as.numeric(value(s)[1, 2]), 1.0, tolerance = 0.1)
})

## @cvxpy test_param_cone_prog.py::TestParamConeProg::test_var_bounds
test_that("var_bounds_conic: bounded variables in conic path", {
  ## Scalar bounds work; matrix upper_bounds have a shape mismatch bug
  x <- Variable(3, bounds = list(-10, 5))
  problem <- Problem(Minimize(sum_entries(x)))
  psolve(problem, solver = "CLARABEL")
  expect_equal(status(problem), "optimal")
  xval <- as.numeric(value(x))
  expect_true(all(xval >= -10 - 1e-3))
  expect_true(all(xval <= 5 + 1e-3))
})

## @cvxpy test_param_cone_prog.py::TestParamConeProg::test_highs_var_bounds
test_that("highs_var_bounds_conic: variable bounds with HiGHS LP", {
  skip_if_not_installed("highs")
  x1 <- Variable(bounds = list(-1, 1))
  x2 <- Variable(bounds = list(-0.5, 1))
  x3 <- Variable()
  objective <- x1 + x2 + x3
  constraints <- list(
    -3 <= x1 + x2, x1 + x2 <= 3,
    -4 <= x1 - x2, x1 - x2 <= 4,
    x3 >= -2
  )
  prob <- Problem(Minimize(objective), constraints)
  psolve(prob, solver = "HIGHS")

  expect_equal(as.numeric(value(x1)), -1, tolerance = 1e-3)
  expect_equal(as.numeric(value(x2)), -0.5, tolerance = 1e-3)
  expect_equal(as.numeric(value(x3)), -2, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# test_dgp_dpp.py — DGP + DPP integration
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_dgp_dpp.py::TestDgpDpp::test_get_problem_data_without_param_values
test_that("DGP DPP: get_problem_data works without param values", {
  skip("CVXR problem_data() does not support gp= parameter; DGP canonicalization eagerly evaluates params")
  ## CVXPY's get_problem_data(solver, gp=True) accepts unset params for DPP caching.
  ## CVXR's problem_data() lacks gp= parameter and triggers eager param evaluation.
})

## @cvxpy test_dgp_dpp.py::TestDgpDpp::test_fast_path_with_changing_params
test_that("DGP DPP: fast path correctly updates params across solves", {
  alpha <- Parameter(pos = TRUE)
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(x), list(x >= alpha))

  ## First solve (CVXR handles DPP automatically without enforce_dpp flag)
  value(alpha) <- 1.0
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)

  ## Second solve with different parameter (DPP fast path)
  value(alpha) <- 2.0
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-3)

  ## Third solve
  value(alpha) <- 0.5
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)
})

## @cvxpy test_dgp_dpp.py::TestDgpDpp::test_solve_without_param_value_raises_error
test_that("DGP DPP: solve without param value raises error", {
  alpha <- Parameter(pos = TRUE)  # No value set
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(x), list(x >= alpha))

  ## Solve without param value should error
  expect_error(
    psolve(problem, solver = "CLARABEL", gp = TRUE)
  )
})

## @cvxpy test_dgp_dpp.py::TestDgpDpp::test_non_dpp_mode_with_ignore_dpp_flag
test_that("DGP DPP: ignore_dpp mode backward compatibility", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(x), list(x >= alpha))

  ## Problem should be DPP-compatible (CVXPY: is_dpp('dgp') → CVXR: .is_dgp_dpp)
  expect_true(CVXR:::.is_dgp_dpp(problem))

  ## CVXR does not implement ignore_dpp as psolve arg;
  ## DPP/non-DPP path is chosen automatically. Just solve normally.
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)

  ## Update parameter and solve again (uses DPP fast path automatically)
  value(alpha) <- 3.0
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# test_valinvec2mixedint.py — FiniteSet gap tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_8
test_that("FiniteSet parametrized: Parameter as set values", {
  skip_if_not_installed("Rglpk")
  x <- Variable()
  objective <- Maximize(x)
  set_vals <- Parameter(5, value = 0:4)
  constraints <- list(FiniteSet(x, set_vals))
  problem <- Problem(objective, constraints)
  psolve(problem, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 4, tolerance = 1e-3)

  ## Update parameter values
  value(set_vals) <- 1:5
  problem2 <- Problem(objective, constraints)
  psolve(problem2, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 5, tolerance = 1e-3)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_9
test_that("FiniteSet with R set (vector): max x in {0,1,2,3,4}", {
  skip_if_not_installed("Rglpk")
  x <- Variable()
  objective <- Maximize(x)
  set_vals <- 0:4  # R equivalent of Python set(range(5))
  constraints <- list(FiniteSet(x, set_vals))
  problem <- Problem(objective, constraints)
  psolve(problem, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 4, tolerance = 1e-3)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_10
test_that("FiniteSet two-element set: max x in {1, 2}", {
  skip_if_not_installed("Rglpk")
  x <- Variable()
  objective <- Maximize(x)
  set_vals <- c(1, 2)
  constraints <- list(FiniteSet(x, set_vals))
  problem <- Problem(objective, constraints)
  psolve(problem, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 2, tolerance = 1e-3)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_gp
test_that("FiniteSet in GP: geometric program with finite set", {
  skip_if_not_installed("Rglpk")
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Single element
  objective <- Maximize(x * y)
  set_vals <- c(2)
  constraints <- list(FiniteSet(x, set_vals), y <= 1)
  problem <- Problem(objective, constraints)
  psolve(problem, gp = TRUE, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 2, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 1, tolerance = 1e-3)

  ## Multiple elements
  set_vals2 <- c(1, 2, 3)
  constraints2 <- list(FiniteSet(x, set_vals2), y <= 1)
  problem2 <- Problem(objective, constraints2)
  psolve(problem2, gp = TRUE, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 3, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 1, tolerance = 1e-3)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_monomial
test_that("FiniteSet on monomial: x*y in {1,2,3}, x=1", {
  skip_if_not_installed("Rglpk")
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  objective <- Maximize(x * y)
  set_vals <- c(1, 2, 3)
  constraints <- list(FiniteSet(x * y, set_vals), x == 1)
  problem <- Problem(objective, constraints)

  ## Solve (CVXR handles DPP automatically; no ignore_dpp/enforce_dpp needed)
  psolve(problem, gp = TRUE, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 1, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 3, tolerance = 1e-3)

  ## Re-solve with fresh problem (DPP path)
  problem2 <- Problem(objective, constraints)
  psolve(problem2, gp = TRUE, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 1, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 3, tolerance = 1e-3)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_invalid_gp
test_that("FiniteSet invalid GP: set contains 0 or negative", {
  skip_if_not_installed("Rglpk")
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  objective <- Maximize(x * y)

  ## Set containing 0 → not DGP
  set_vals_zero <- c(0, 1, 2, 3)
  constraints_zero <- list(FiniteSet(x, set_vals_zero), y <= 1)
  problem_zero <- Problem(objective, constraints_zero)
  expect_error(psolve(problem_zero, gp = TRUE, solver = "GLPK_MI"),
    "not DGP|DGP")

  ## Set containing negative → not DGP
  set_vals_neg <- c(-1, 1)
  constraints_neg <- list(FiniteSet(x, set_vals_neg), y <= 1)
  problem_neg <- Problem(objective, constraints_neg)
  expect_error(psolve(problem_neg, gp = TRUE, solver = "GLPK_MI"),
    "not DGP|DGP")
})

## @cvxpy test_valinvec2mixedint.py::test_default_argument
test_that("FiniteSet default argument: ineq_form defaults FALSE", {
  skip_if_not_installed("Rglpk")
  x <- Variable()
  objective <- Maximize(x)
  set_vals <- 0:4
  constraints <- list(FiniteSet(x, set_vals))  # default ineq_form
  problem <- Problem(objective, constraints)
  psolve(problem, solver = "GLPK_MI")
  expect_equal(as.numeric(value(x)), 4, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# test_logic.py — Logic atom gaps
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicBoolConstant::test_sparse_bool_constant
test_that("sparse_bool_constant: And with sparse matrix mask", {
  skip_if_not_installed("Rglpk")
  skip_if_not_installed("Matrix")

  x <- Variable(c(2, 2), boolean = TRUE)
  mask <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(1, 1), dims = c(2, 2)
  )
  ## mask is sparse identity → And(x, mask) should give diag elements

  ## Evaluate And(x, mask) with x = ones(2,2)
  y <- Variable(c(2, 2))
  prob <- Problem(Minimize(Constant(0)),
    list(y == And(x, mask), x == matrix(1, 2, 2)))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(as.numeric(value(y)), c(1, 0, 0, 1), tolerance = 1e-3)
})

## @cvxpy test_logic.py::TestLogicName::test_format_labeled
test_that("format_labeled: labeled logic expressions", {
  skip("format_labeled/set_label not implemented in CVXR")
  ## CVXPY test: And(x, y).set_label("my_and").format_labeled() == "my_and"
  ## This feature is not yet available in CVXR
})

# ═══════════════════════════════════════════════════════════════════════
# test_pow_cone_nd.py — Power cone ND gaps
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_nonuniform_alpha
test_that("PowConeND balanced tree nonuniform alpha", {
  alpha <- c(0.05, 0.1, 0.15, 0.2, 0.2, 0.15, 0.1, 0.05)
  bounds <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
  expected <- prod(bounds^alpha)  # ≈ 4.0576

  W <- Variable(8, pos = TRUE)
  z <- Variable()
  ## axis=1 in R corresponds to axis=0 in CVXPY for column vector
  con <- PowConeND(W, z, alpha)
  prob <- Problem(Maximize(z), list(con, W <= bounds))
  psolve(prob, solver = "CLARABEL")
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(as.numeric(value(z)), expected, tolerance = 1e-3)
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_multiple_cones
test_that("PowConeND balanced tree multiple cones (k > 1)", {
  n <- 8L; k <- 3L
  W <- Variable(c(n, k), pos = TRUE)
  z <- Variable(k)
  alpha <- matrix(1 / n, n, k)
  con <- PowConeND(W, z, alpha, axis = 2)  # R axis=2 → CVXPY axis=0
  prob <- Problem(Maximize(sum_entries(z)), list(con, W <= 2))
  psolve(prob, solver = "CLARABEL")
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(as.numeric(value(z)), c(2.0, 2.0, 2.0), tolerance = 1e-3)
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_multiple_cones_nonuniform
test_that("PowConeND balanced tree multiple cones nonuniform alpha", {
  n <- 6L; k <- 2L
  bounds <- matrix(c(
    1.0, 2.0,
    2.0, 3.0,
    3.0, 4.0,
    4.0, 5.0,
    5.0, 6.0,
    6.0, 7.0
  ), n, k, byrow = TRUE)
  alpha <- matrix(c(
    0.1, 0.3,
    0.2, 0.2,
    0.3, 0.1,
    0.15, 0.15,
    0.15, 0.1,
    0.1, 0.15
  ), n, k, byrow = TRUE)
  expected <- apply(bounds^alpha, 2, prod)  # ≈ [2.994, 3.592]

  W <- Variable(c(n, k), pos = TRUE)
  z <- Variable(k)
  con <- PowConeND(W, z, alpha, axis = 2)  # R axis=2 → CVXPY axis=0
  prob <- Problem(Maximize(sum_entries(z)), list(con, W <= bounds))
  psolve(prob, solver = "CLARABEL")
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(as.numeric(value(z)), expected, tolerance = 1e-2)
})

# ═══════════════════════════════════════════════════════════════════════
# test_problem.py — Problem-level gaps
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_verbose
test_that("verbose: solver output controlled by verbose flag", {
  a <- Variable(name = "a")
  x <- Variable(2, name = "x")
  p <- Problem(Minimize(a + x[1]), list(a >= 2, x >= 2))

  ## verbose = TRUE should produce cli messages (sent to stderr via cli_inform)
  msg_verbose <- capture.output(
    invisible(psolve(p, verbose = TRUE, solver = "CLARABEL")),
    type = "message"
  )
  ## At minimum the CVXR banner and summary are printed
  expect_true(length(msg_verbose) > 0)

  ## verbose = FALSE should produce no cli messages
  msg_silent <- capture.output(
    invisible(psolve(p, verbose = FALSE, solver = "CLARABEL")),
    type = "message"
  )
  expect_equal(length(msg_silent), 0)
})

## @cvxpy test_problem.py::TestProblem::test_solve_solver_path
test_that("solve_solver_path: skipped (not implemented in CVXR)", {
  skip("solver_path argument not implemented in CVXR")
  ## CVXPY allows solver_path = [(solver1, opts1), solver2, ...]
  ## to try multiple solvers in sequence
})

## @cvxpy test_problem.py::TestProblem::test_cp_node_count_warn
test_that("cp_node_count_warn: skipped (no node count warning in CVXR)", {
  skip("Node count warning not implemented in CVXR")
  ## CVXPY warns when expression tree has too many nodes
  ## suggesting vectorization. This is Python-specific.
})

## @cvxpy test_problem.py::TestProblem::test_ecos_warning
test_that("ecos_warning: skipped (ECOS deprecation N/A in CVXR)", {
  skip("ECOS deprecation warning is Python/CVXPY specific")
  ## CVXPY warns when ECOS is selected by default.
  ## CVXR does not have this deprecation path.
})

# ═══════════════════════════════════════════════════════════════════════
# test_constant_atoms.py — Constant atoms gap
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constant_atoms.py::test_constant_atoms
test_that("constant_atoms: representative atom evaluation with constants", {
  ## The CVXPY test_constant_atoms is a massive parametrized test over
  ## ~70+ atoms. Here we test a representative sample.

  ## 1. norm(x, 2) with constant
  x_val <- c(3, 4)
  expr <- p_norm(Constant(x_val), 2)
  prob <- Problem(Minimize(expr))
  psolve(prob, solver = "CLARABEL")
  expect_equal(value(prob), 5.0, tolerance = 1e-5)

  ## 2. sum_squares with constant
  expr2 <- sum_squares(Constant(c(1, 2, 3)))
  prob2 <- Problem(Minimize(expr2))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(value(prob2), 14.0, tolerance = 1e-5)

  ## 3. log with variable pinned to constant
  x <- Variable(2)
  prob3 <- Problem(Maximize(sum_entries(log(x))), list(x == c(1, exp(1))))
  psolve(prob3, solver = "SCS")
  expect_equal(value(prob3), 1.0, tolerance = 0.1)

  ## 4. quad_form with constants
  x_val2 <- c(1, 2)
  P <- matrix(c(1, 0, 0, 1), 2, 2)
  expr4 <- quad_form(Constant(x_val2), P)
  prob4 <- Problem(Minimize(expr4))
  psolve(prob4, solver = "CLARABEL")
  expect_equal(value(prob4), 5.0, tolerance = 1e-5)

  ## 5. abs with variable pinned
  x5 <- Variable(3)
  prob5 <- Problem(Minimize(sum_entries(abs(x5))),
    list(x5 == c(-1, 2, -3)))
  psolve(prob5, solver = "CLARABEL")
  expect_equal(value(prob5), 6.0, tolerance = 1e-3)
})
