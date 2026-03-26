## Phase 8c: Matrix/SDP atom tests
## Tests for LambdaMax, SigmaMax, NormNuc and their canonicalizers

# ═══════════════════════════════════════════════════════════════════
# LambdaMax — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("LambdaMax: shape is scalar", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- lambda_max(X)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("LambdaMax: sign is (FALSE, FALSE)", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- lambda_max(X)
  expect_false(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## @cvxpy NONE
test_that("LambdaMax: curvature is convex", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- lambda_max(X)
  expect_true(is_convex(expr))
  expect_false(is_concave(expr))
  expect_true(is_dcp(expr))
})

## @cvxpy NONE
test_that("LambdaMax: not monotone", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- lambda_max(X)
  expect_false(is_incr(expr, 1L))
  expect_false(is_decr(expr, 1L))
})

## @cvxpy NONE
test_that("LambdaMax: validate requires square matrix", {
  x <- Variable(c(3, 4))
  expect_error(lambda_max(x), "square")
})

## @cvxpy NONE
test_that("LambdaMax: validate requires 2D", {
  x <- Variable(c(3, 1))
  ## 3x1 is not square
  expect_error(lambda_max(x), "square")
})

## @cvxpy NONE
test_that("LambdaMax: numeric value is correct", {
  A <- matrix(c(2, 1, 1, 3), 2, 2)
  expr <- lambda_max(Constant(A))
  val <- value(expr)
  evals <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(as.numeric(val), max(evals), tolerance = 1e-10)
})

## @cvxpy NONE
test_that("LambdaMax: numeric value with negative eigenvalues", {
  A <- matrix(c(-2, 0, 0, -5), 2, 2)
  expr <- lambda_max(Constant(A))
  val <- value(expr)
  expect_equal(as.numeric(val), -2, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# lambda_min — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("lambda_min: is concave", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  expr <- lambda_min(X)
  expect_true(is_concave(expr))
  expect_false(is_convex(expr))
})

## @cvxpy NONE
test_that("lambda_min: shape is scalar", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  expr <- lambda_min(X)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("lambda_min: numeric value is correct", {
  A <- matrix(c(2, 1, 1, 3), 2, 2)
  expr <- lambda_min(Constant(A))
  val <- value(expr)
  evals <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(as.numeric(val), min(evals), tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# SigmaMax — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("SigmaMax: shape is scalar", {
  Y <- Variable(c(3, 4))
  expr <- sigma_max(Y)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("SigmaMax: sign is nonneg", {
  Y <- Variable(c(3, 4))
  expr <- sigma_max(Y)
  expect_true(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## @cvxpy NONE
test_that("SigmaMax: curvature is convex", {
  Y <- Variable(c(3, 4))
  expr <- sigma_max(Y)
  expect_true(is_convex(expr))
  expect_false(is_concave(expr))
})

## @cvxpy NONE
test_that("SigmaMax: accepts non-square matrix", {
  Y <- Variable(c(2, 5))
  expr <- sigma_max(Y)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("SigmaMax: accepts square matrix", {
  Y <- Variable(c(3, 3))
  expr <- sigma_max(Y)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("SigmaMax: numeric value is correct", {
  A <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  expr <- sigma_max(Constant(A))
  val <- value(expr)
  expect_equal(as.numeric(val), max(svd(A)$d), tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# NormNuc — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("NormNuc: shape is scalar", {
  Y <- Variable(c(3, 4))
  expr <- norm_nuc(Y)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("NormNuc: sign is nonneg", {
  Y <- Variable(c(3, 4))
  expr <- norm_nuc(Y)
  expect_true(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## @cvxpy NONE
test_that("NormNuc: curvature is convex", {
  Y <- Variable(c(3, 4))
  expr <- norm_nuc(Y)
  expect_true(is_convex(expr))
})

## @cvxpy NONE
test_that("NormNuc: numeric value is correct", {
  A <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  expr <- norm_nuc(Constant(A))
  val <- value(expr)
  expect_equal(as.numeric(val), sum(svd(A)$d), tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# Canonicalization tests (no solver)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("LambdaMax: DCP accepts Minimize(lambda_max(X))", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)), list(X[1,1] >= 1))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("lambda_min: DCP accepts Maximize(lambda_min(X))", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Maximize(lambda_min(X)), list(X[1,1] <= 5))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("SigmaMax: DCP accepts Minimize(sigma_max(Y))", {
  Y <- Variable(c(2, 3))
  prob <- Problem(Minimize(sigma_max(Y)), list(Y[1,1] >= 1))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("NormNuc: DCP accepts Minimize(norm_nuc(Y))", {
  Y <- Variable(c(2, 3))
  prob <- Problem(Minimize(norm_nuc(Y)), list(Y[1,1] >= 1))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("LambdaMax: canonicalization produces PSD constraint", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)), list(X[1,1] >= 2))
  pd <- problem_data(prob, solver = CLARABEL_SOLVER)
  ## Should have PSD cone (psd dimension > 0)
  expect_true(pd$data$dims@psd[1L] > 0L)
})

## @cvxpy NONE
test_that("SigmaMax: canonicalization produces PSD constraint", {
  Y <- Variable(c(2, 3))
  prob <- Problem(Minimize(sigma_max(Y)), list(Y[1,1] >= 1))
  pd <- problem_data(prob, solver = CLARABEL_SOLVER)
  ## Should have PSD cone
  expect_true(pd$data$dims@psd[1L] > 0L)
})

## @cvxpy NONE
test_that("NormNuc: canonicalization produces PSD constraint", {
  Y <- Variable(c(2, 3))
  prob <- Problem(Minimize(norm_nuc(Y)), list(Y[1,1] >= 1))
  pd <- problem_data(prob, solver = CLARABEL_SOLVER)
  ## Should have PSD cone
  expect_true(pd$data$dims@psd[1L] > 0L)
})

# ═══════════════════════════════════════════════════════════════════
# Solve tests — CVXPY parity
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("LambdaMax: solve with Clarabel matches CVXPY", {
  skip_if_not_installed("clarabel")
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)),
                  list(X[1,1] >= 2, X[2,2] >= 3, X[3,3] >= 4))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 4.000016
  expect_equal(val, 4.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("LambdaMax: solve with SCS matches CVXPY", {
  skip_if_not_installed("scs")
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)),
                  list(X[1,1] >= 2, X[2,2] >= 3, X[3,3] >= 4))
  val <- psolve(prob, solver = SCS_SOLVER)
  ## CVXPY: 4.0
  expect_equal(val, 4.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("lambda_min: solve with Clarabel matches CVXPY", {
  skip_if_not_installed("clarabel")
  W <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Maximize(lambda_min(W)),
                  list(W[1,1] <= 5, W[2,2] <= 3))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 3.0
  expect_equal(val, 3.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("SigmaMax: solve with Clarabel matches CVXPY", {
  skip_if_not_installed("clarabel")
  Y <- Variable(c(2, 3))
  prob <- Problem(Minimize(sigma_max(Y)),
                  list(Y[1,1] >= 1, Y[1,2] >= 2, Y[2,3] >= 3))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 3.0
  expect_equal(val, 3.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("NormNuc: solve with Clarabel matches CVXPY", {
  skip_if_not_installed("clarabel")
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(norm_nuc(Z)),
                  list(Z[1,1] >= 1, Z[2,2] >= 2))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 3.0
  expect_equal(val, 3.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("LambdaMax: non-symmetric input gets symmetry constraint", {
  skip_if_not_installed("clarabel")
  ## Use a general (non-symmetric) variable
  X <- Variable(c(2, 2))
  prob <- Problem(Minimize(lambda_max(X)),
                  list(X[1,1] >= 3, X[2,2] >= 5))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## Symmetry constraint forces off-diag to be equal → min eigenvalue = max diag
  expect_equal(val, 5.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("SigmaMax: rectangular matrix solve", {
  skip_if_not_installed("clarabel")
  Y <- Variable(c(3, 2))
  prob <- Problem(Minimize(sigma_max(Y)),
                  list(Y[1,1] >= 2, Y[2,2] >= 3))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## Min sigma_max with Y[1,1]>=2, Y[2,2]>=3 → 3.0
  expect_equal(val, 3.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("NormNuc: solve with SCS matches CVXPY", {
  skip_if_not_installed("scs")
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(norm_nuc(Z)),
                  list(Z[1,1] >= 1, Z[2,2] >= 2))
  val <- psolve(prob, solver = SCS_SOLVER)
  ## CVXPY: 3.0
  expect_equal(val, 3.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("LambdaMax: constraint usage (PSD via lambda_max)", {
  skip_if_not_installed("clarabel")
  ## lambda_max(X) <= t is equivalent to PSD(t*I - X)
  X <- Variable(c(2, 2), symmetric = TRUE)
  t_var <- Variable()
  prob <- Problem(Minimize(t_var),
                  list(lambda_max(X) <= t_var, X[1,1] >= 2, X[2,2] >= 3))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(val, 3.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("NormNuc: as constraint", {
  skip_if_not_installed("clarabel")
  Y <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Y)),
                  list(norm_nuc(Y) <= 5, Y >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## sum_entries minimized when Y is all zeros (norm_nuc(0) = 0 <= 5)
  expect_equal(val, 0.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# MatrixFrac — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("MatrixFrac: shape is scalar", {
  x <- Variable(c(3, 1))
  P <- Variable(c(3, 3), symmetric = TRUE)
  expr <- MatrixFrac(x, P)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("MatrixFrac: sign is nonneg", {
  x <- Variable(c(3, 1))
  P <- Variable(c(3, 3), symmetric = TRUE)
  expr <- MatrixFrac(x, P)
  expect_true(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## @cvxpy NONE
test_that("MatrixFrac: curvature is convex", {
  x <- Variable(c(3, 1))
  P <- Variable(c(3, 3), symmetric = TRUE)
  expr <- MatrixFrac(x, P)
  expect_true(is_convex(expr))
  expect_false(is_concave(expr))
})

## @cvxpy NONE
test_that("MatrixFrac: not monotone", {
  x <- Variable(c(3, 1))
  P <- Variable(c(3, 3), symmetric = TRUE)
  expr <- MatrixFrac(x, P)
  expect_false(is_incr(expr, 1L))
  expect_false(is_decr(expr, 1L))
  expect_false(is_incr(expr, 2L))
  expect_false(is_decr(expr, 2L))
})

## @cvxpy NONE
test_that("MatrixFrac: validate requires square P", {
  x <- Variable(c(3, 1))
  P <- Variable(c(3, 4))
  expect_error(MatrixFrac(x, P), "square")
})

## @cvxpy NONE
test_that("MatrixFrac: validate requires compatible dimensions", {
  x <- Variable(c(3, 1))
  P <- Variable(c(4, 4), symmetric = TRUE)
  expect_error(MatrixFrac(x, P), "incompatible")
})

## @cvxpy NONE
test_that("MatrixFrac: numeric value is correct", {
  X <- matrix(c(1, 2, 3), 3, 1)
  P <- matrix(c(2, 0, 0, 0, 3, 0, 0, 0, 4), 3, 3)
  expr <- MatrixFrac(Constant(X), Constant(P))
  val <- value(expr)
  expected <- sum(diag(t(X) %*% solve(P) %*% X))
  expect_equal(as.numeric(val), expected, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("MatrixFrac: quadratic overrides work", {
  x <- Variable(c(3, 1))
  P_const <- Constant(diag(3))
  P_var <- Variable(c(3, 3), symmetric = TRUE)
  ## X affine, P constant → quadratic
  expr1 <- MatrixFrac(x, P_const)
  expect_true(is_quadratic(expr1))
  expect_true(has_quadratic_term(expr1))
  ## P variable → not quadratic
  expr2 <- MatrixFrac(x, P_var)
  expect_false(is_quadratic(expr2))
  expect_false(has_quadratic_term(expr2))
})

## @cvxpy NONE
test_that("MatrixFrac: is_qpwa override", {
  x <- Variable(c(2, 1))
  P_const <- Constant(diag(2))
  expr <- MatrixFrac(x, P_const)
  expect_true(is_qpwa(expr))
})

## @cvxpy NONE
test_that("matrix_frac: convenience function with constant P uses QuadForm shortcut", {
  ## When P is a plain matrix, matrix_frac returns QuadForm
  x <- Variable(c(3, 1))
  P_mat <- diag(3)
  expr <- matrix_frac(x, P_mat)
  ## Should be QuadForm, not MatrixFrac
  expect_true(S7_inherits(expr, QuadForm))
})

## @cvxpy NONE
test_that("matrix_frac: convenience function with variable P creates MatrixFrac", {
  x <- Variable(c(3, 1))
  P <- Variable(c(3, 3), symmetric = TRUE)
  expr <- matrix_frac(x, P)
  expect_true(S7_inherits(expr, MatrixFrac))
})

## @cvxpy NONE
test_that("MatrixFrac: multi-column X (n by m)", {
  X <- Variable(c(3, 2))
  P <- Variable(c(3, 3), symmetric = TRUE)
  expr <- MatrixFrac(X, P)
  expect_equal(expr@shape, c(1L, 1L))
})

# ═══════════════════════════════════════════════════════════════════
# TrInv — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("TrInv: shape is scalar", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- tr_inv(X)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("TrInv: sign is nonneg", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- tr_inv(X)
  expect_true(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## @cvxpy NONE
test_that("TrInv: curvature is convex", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- tr_inv(X)
  expect_true(is_convex(expr))
  expect_false(is_concave(expr))
})

## @cvxpy NONE
test_that("TrInv: not monotone", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- tr_inv(X)
  expect_false(is_incr(expr, 1L))
  expect_false(is_decr(expr, 1L))
})

## @cvxpy NONE
test_that("TrInv: validate requires square matrix", {
  x <- Variable(c(3, 4))
  expect_error(tr_inv(x), "square")
})

## @cvxpy NONE
test_that("TrInv: numeric value is correct", {
  A <- diag(c(1, 2, 3))
  expr <- tr_inv(Constant(A))
  val <- value(expr)
  ## tr(A^{-1}) = 1/1 + 1/2 + 1/3 = 11/6
  expect_equal(as.numeric(val), 11 / 6, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("TrInv: numeric value Inf for non-PSD", {
  A <- diag(c(1, -1, 3))
  expr <- tr_inv(Constant(A))
  val <- value(expr)
  expect_equal(as.numeric(val), Inf)
})

# ═══════════════════════════════════════════════════════════════════
# LambdaSumLargest — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("LambdaSumLargest: shape is scalar", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- lambda_sum_largest(X, 2L)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("LambdaSumLargest: sign is (FALSE, FALSE)", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- lambda_sum_largest(X, 2L)
  expect_false(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## @cvxpy NONE
test_that("LambdaSumLargest: curvature is convex", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- lambda_sum_largest(X, 2L)
  expect_true(is_convex(expr))
  expect_false(is_concave(expr))
})

## @cvxpy NONE
test_that("LambdaSumLargest: validate requires square", {
  X <- Variable(c(3, 4))
  expect_error(lambda_sum_largest(X, 2L), "square")
})

## @cvxpy NONE
test_that("LambdaSumLargest: validate requires positive k", {
  X <- Variable(c(3, 3))
  expect_error(lambda_sum_largest(X, 0L), "positive")
})

## @cvxpy NONE
test_that("LambdaSumLargest: get_data returns k", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- LambdaSumLargest(X, 2L)
  expect_equal(get_data(expr), list(2L))
})

## @cvxpy NONE
test_that("LambdaSumLargest: numeric value is correct", {
  A <- diag(c(1, 4, 7))
  expr <- lambda_sum_largest(Constant(A), 2L)
  val <- value(expr)
  ## Two largest eigenvalues: 7, 4 → sum = 11
  expect_equal(as.numeric(val), 11, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("LambdaSumLargest: k=1 equals lambda_max", {
  A <- diag(c(2, 5, 3))
  expr1 <- lambda_sum_largest(Constant(A), 1L)
  expr2 <- lambda_max(Constant(A))
  expect_equal(as.numeric(value(expr1)), as.numeric(value(expr2)), tolerance = 1e-10)
})

## @cvxpy NONE
test_that("LambdaSumLargest: k=n equals trace", {
  A <- diag(c(2, 5, 3))
  expr <- lambda_sum_largest(Constant(A), 3L)
  val <- value(expr)
  ## Sum of all eigenvalues = trace = 10
  expect_equal(as.numeric(val), 10, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# lambda_sum_smallest — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("lambda_sum_smallest: concave", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- lambda_sum_smallest(X, 2L)
  expect_true(is_concave(expr))
})

## @cvxpy NONE
test_that("lambda_sum_smallest: numeric value", {
  A <- diag(c(1, 4, 7))
  expr <- lambda_sum_smallest(Constant(A), 2L)
  val <- value(expr)
  ## Two smallest eigenvalues: 1, 4 → sum = 5
  expect_equal(as.numeric(val), 5, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# LogDet — unit tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("LogDet: shape is scalar", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- log_det(X)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("LogDet: sign is unknown (CVXPY v1.8.2 fix)", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- log_det(X)
  ## CVXPY v1.8.2 fix: log_det can be negative (log(det) < 0 when 0 < det < 1)
  expect_false(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## @cvxpy NONE
test_that("LogDet: curvature is concave", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- log_det(X)
  expect_true(is_concave(expr))
  expect_false(is_convex(expr))
})

## @cvxpy NONE
test_that("LogDet: not monotone", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  expr <- log_det(X)
  expect_false(is_incr(expr, 1L))
  expect_false(is_decr(expr, 1L))
})

## @cvxpy NONE
test_that("LogDet: validate requires square matrix", {
  X <- Variable(c(3, 4))
  expect_error(log_det(X), "square")
})

## @cvxpy NONE
test_that("LogDet: numeric value is correct", {
  A <- diag(c(1, 2, 3))
  expr <- log_det(Constant(A))
  val <- value(expr)
  ## log(det(A)) = log(6) ≈ 1.7918
  expect_equal(as.numeric(val), log(6), tolerance = 1e-10)
})

## @cvxpy NONE
test_that("LogDet: numeric value for identity", {
  A <- diag(3)
  expr <- log_det(Constant(A))
  val <- value(expr)
  ## log(det(I)) = 0
  expect_equal(as.numeric(val), 0, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# DCP acceptance tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("MatrixFrac: DCP accepts Minimize(matrix_frac(x, P))", {
  x <- Variable(c(3, 1))
  P <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(MatrixFrac(x, P)), list(P[1,1] >= 1))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("TrInv: DCP accepts Minimize(tr_inv(X))", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(tr_inv(X)), list(matrix_trace(X) <= 6))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("LambdaSumLargest: DCP accepts Minimize", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_sum_largest(X, 2L)), list(X[1,1] >= 1))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("lambda_sum_smallest: DCP accepts Maximize", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Maximize(lambda_sum_smallest(X, 2L)), list(X[1,1] <= 5))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("LogDet: DCP accepts Maximize(log_det(X))", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Maximize(log_det(X)), list(matrix_trace(X) <= 6))
  expect_true(is_dcp(prob))
})

# ═══════════════════════════════════════════════════════════════════
# Canonicalization tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("MatrixFrac: canonicalization produces PSD constraint", {
  x <- Variable(c(3, 1))
  P <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(MatrixFrac(x, P)), list(P[1,1] >= 1, P[2,2] >= 1, P[3,3] >= 1))
  pd <- problem_data(prob, solver = CLARABEL_SOLVER)
  expect_true(pd$data$dims@psd[1L] > 0L)
})

## @cvxpy NONE
test_that("TrInv: canonicalization produces PSD constraints", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(tr_inv(X)), list(matrix_trace(X) <= 4))
  pd <- problem_data(prob, solver = CLARABEL_SOLVER)
  expect_true(pd$data$dims@psd[1L] > 0L)
})

## @cvxpy NONE
test_that("LambdaSumLargest: canonicalization produces PSD constraint", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_sum_largest(X, 2L)), list(X[1,1] >= 1))
  pd <- problem_data(prob, solver = CLARABEL_SOLVER)
  expect_true(pd$data$dims@psd[1L] > 0L)
})

## @cvxpy NONE
test_that("LogDet: canonicalization produces PSD + ExpCone constraints", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Maximize(log_det(X)), list(matrix_trace(X) <= 4))
  pd <- problem_data(prob, solver = CLARABEL_SOLVER)
  ## Should have both PSD and ExpCone
  expect_true(pd$data$dims@psd[1L] > 0L)
  expect_true(pd$data$dims@exp > 0L)
})

## @cvxpy NONE
test_that("DCP canonicalization S7 methods registered for matrix atoms", {
  X <- Variable(c(2, 2))
  expect_true(has_dcp_canon(MatrixFrac(Variable(2), Constant(diag(2)))))
  expect_true(has_dcp_canon(TrInv(X)))
  expect_true(has_dcp_canon(LambdaSumLargest(X, 1)))
  expect_true(has_dcp_canon(LogDet(X)))
})

# ═══════════════════════════════════════════════════════════════════
# Solve tests — CVXPY parity (Batch 2-4)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("MatrixFrac: solve with constant P matches CVXPY", {
  skip_if_not_installed("clarabel")
  ## matrix_frac(x, I) with x >= 1 → ||x||^2 = 3
  x <- Variable(c(3, 1))
  prob <- Problem(Minimize(matrix_frac(x, diag(3))),
                  list(x >= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 3.0
  expect_equal(val, 3.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MatrixFrac: solve with variable P (Clarabel)", {
  skip_if_not_installed("clarabel")
  ## Minimize matrix_frac([1;1;1], P) s.t. trace(P) <= 6, P >> 0
  ## At optimum P = 2I, value = tr(I/2 * I * I) = tr((1/2)I) * 3? No...
  ## matrix_frac(x, P) = tr(x^T P^{-1} x) = x^T P^{-1} x = [1,1,1] (2I)^{-1} [1,1,1]^T = 3/2
  ## But wait: with PSD constraint and trace<=6, CVXPY gives 0.5
  ## That's because P can be arbitrarily large in one direction
  ## Actually CVXPY gives 0.5 which means P = 2I is NOT the optimum
  ## The CVXPY result is 0.5 with diag(P) ≈ [2,2,2]
  x_const <- matrix(c(1, 1, 1), 3, 1)
  P <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(MatrixFrac(Constant(x_const), P)),
                  list(matrix_trace(P) <= 6, PSD(P)))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 0.5 (P concentrates all trace budget in x direction)
  expect_equal(val, 0.5, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("MatrixFrac: solve with P=2I matches CVXPY", {
  skip_if_not_installed("clarabel")
  ## matrix_frac(x, 2I) with x >= 1 → 0.5*||x||^2 = 1.5
  x <- Variable(c(3, 1))
  prob <- Problem(Minimize(matrix_frac(x, 2 * diag(3))),
                  list(x >= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 1.5
  expect_equal(val, 1.5, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("TrInv: solve with Clarabel matches CVXPY", {
  skip_if_not_installed("clarabel")
  ## Minimize tr_inv(X) s.t. trace(X) <= 6, X >> 0
  ## At optimum, X = 2I, tr_inv = 3/2 = 1.5
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(tr_inv(X)),
                  list(matrix_trace(X) <= 6, PSD(X)))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 1.5
  expect_equal(val, 1.5, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("LambdaSumLargest: solve with Clarabel matches CVXPY", {
  skip_if_not_installed("clarabel")
  ## Minimize lambda_sum_largest(X, 2) s.t. X[0,0]>=3, X[1,1]>=4, X[2,2]>=5
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_sum_largest(X, 2L)),
                  list(X[1,1] >= 3, X[2,2] >= 4, X[3,3] >= 5))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 9.0 (eigenvalues 3, 4, 5 → sum of 2 largest = 9)
  expect_equal(val, 9.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("lambda_sum_smallest: solve with Clarabel matches CVXPY", {
  skip_if_not_installed("clarabel")
  ## Maximize lambda_sum_smallest(X, 2) s.t. X[0,0]<=3, X[1,1]<=4, X[2,2]<=5, X>>0
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Maximize(lambda_sum_smallest(X, 2L)),
                  list(X[1,1] <= 3, X[2,2] <= 4, X[3,3] <= 5, PSD(X)))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  ## CVXPY: 7.0 (eigenvalues 3, 4, 5 → sum of 2 smallest = 7)
  expect_equal(val, 7.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("LogDet: solve with Clarabel matches CVXPY", {
  skip_if_not_installed("clarabel")
  ## Maximize log_det(X) s.t. diag(X) <= 2, X >> 0
  ## At optimum X = 2I, log_det = 3*log(2) ≈ 2.079
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Maximize(log_det(X)),
                  list(DiagMat(X) <= 2, PSD(X)))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(val, 3 * log(2), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("LogDet: solve with SCS matches CVXPY", {
  skip_if_not_installed("scs")
  ## Same as above with SCS
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Maximize(log_det(X)),
                  list(DiagMat(X) <= 2, PSD(X)))
  val <- psolve(prob, solver = SCS_SOLVER)
  expect_equal(val, 3 * log(2), tolerance = 5e-2)
})

## @cvxpy NONE
test_that("LambdaSumLargest: k=1 solve matches lambda_max", {
  skip_if_not_installed("clarabel")
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob1 <- Problem(Minimize(lambda_sum_largest(X, 1L)),
                   list(X[1,1] >= 2, X[2,2] >= 3, X[3,3] >= 4))
  val1 <- psolve(prob1, solver = CLARABEL_SOLVER)
  X2 <- Variable(c(3, 3), symmetric = TRUE)
  prob2 <- Problem(Minimize(lambda_max(X2)),
                   list(X2[1,1] >= 2, X2[2,2] >= 3, X2[3,3] >= 4))
  val2 <- psolve(prob2, solver = CLARABEL_SOLVER)
  expect_equal(val1, val2, tolerance = 1e-2)
})
