## Tests for Phase 9a: CRAN Cleanup + Pipeline Hardening

# ═══════════════════════════════════════════════════════════════
# is_qp / is_lp PSD attribute check
# ═══════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("is_qp returns FALSE for PSD variable", {
  x <- Variable(c(2, 2), PSD = TRUE)
  prob <- Problem(Minimize(sum_entries(x)))
  expect_false(is_qp(prob))
})

## @cvxpy NONE
test_that("is_lp returns FALSE for PSD variable", {
  x <- Variable(c(2, 2), PSD = TRUE)
  prob <- Problem(Minimize(sum_entries(x)))
  expect_false(is_lp(prob))
})

## @cvxpy NONE
test_that("is_qp returns FALSE for NSD variable", {
  x <- Variable(c(2, 2), NSD = TRUE)
  prob <- Problem(Minimize(sum_entries(x)))
  expect_false(is_qp(prob))
})

## @cvxpy NONE
test_that("is_qp still TRUE for regular QP (no regression)", {
  x <- Variable(3)
  prob <- Problem(Minimize(quad_form(x, diag(3))), list(x >= -1))
  expect_true(is_qp(prob))
})

## @cvxpy NONE
test_that("is_lp still TRUE for regular LP (no regression)", {
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  expect_true(is_lp(prob))
})

## @cvxpy NONE
test_that("is_qp still TRUE for LP with equality (no regression)", {
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0, sum_entries(x) == 10))
  expect_true(is_qp(prob))
  expect_true(is_lp(prob))
})

# ═══════════════════════════════════════════════════════════════
# NaN/Inf data validation
# ═══════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that(".check_finite catches NaN in dense vector", {
  expect_error(
    CVXR:::.check_finite(c(1, NaN, 3), "test"),
    "NaN"
  )
})

## @cvxpy NONE
test_that(".check_finite catches Inf in dense vector", {
  expect_error(
    CVXR:::.check_finite(c(1, Inf, 3), "test"),
    "Inf"
  )
})

## @cvxpy NONE
test_that(".check_finite catches NaN in sparse matrix", {
  m <- Matrix::sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(1, NaN), dims = c(3, 3))
  expect_error(
    CVXR:::.check_finite(m, "test"),
    "NaN"
  )
})

## @cvxpy NONE
test_that(".check_finite passes clean data", {
  expect_no_error(CVXR:::.check_finite(c(1, 2, 3), "test"))
  expect_no_error(CVXR:::.check_finite(NULL, "test"))
  expect_no_error(CVXR:::.check_finite(numeric(0), "test"))
})

## @cvxpy NONE
test_that(".validate_problem_data catches NaN in A matrix", {
  data <- list(A = Matrix::sparseMatrix(i = 1, j = 1, x = NaN, dims = c(2, 2)),
               b = c(0, 0), c = c(1, 0))
  expect_error(
    CVXR:::.validate_problem_data(data),
    "NaN"
  )
})

# ═══════════════════════════════════════════════════════════════
# warm_start parameter
# ═══════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("psolve accepts warm_start = FALSE (default behavior)", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER, warm_start = FALSE)
  expect_equal(val, 2, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("psolve accepts warm_start = TRUE without error", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER, warm_start = TRUE)
  expect_equal(val, 2, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════
# scalene
# ═══════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("scalene DCP: convex for positive alpha, beta", {
  x <- Variable(3)
  expr <- scalene(x, 2, 3)
  expect_true(is_convex(expr))
})

## @cvxpy NONE
test_that("scalene numeric value matches CVXPY", {
  ## CVXPY: scalene([-1, 2], 2, 3) = [3, 4]
  c_val <- Constant(c(-1, 2))
  expr <- scalene(c_val, 2, 3)
  result <- value(expr)
  expect_equal(as.numeric(result), c(3, 4), tolerance = 1e-6)
})

## @cvxpy NONE
test_that("scalene numeric all positive", {
  c_val <- Constant(c(1, 2))
  expr <- scalene(c_val, 2, 3)
  result <- value(expr)
  ## pos(1)=1, neg(1)=0 → 2*1+3*0=2; pos(2)=2, neg(2)=0 → 2*2+3*0=4
  expect_equal(as.numeric(result), c(2, 4), tolerance = 1e-6)
})

## @cvxpy NONE
test_that("scalene numeric all negative", {
  c_val <- Constant(c(-3, -1))
  expr <- scalene(c_val, 2, 3)
  result <- value(expr)
  ## pos(-3)=0, neg(-3)=3 → 0+9=9; pos(-1)=0, neg(-1)=1 → 0+3=3
  expect_equal(as.numeric(result), c(9, 3), tolerance = 1e-6)
})

# ═══════════════════════════════════════════════════════════════
# harmonic_mean
# ═══════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("harmonic_mean DCP: concave for nonneg x", {
  x <- Variable(2, nonneg = TRUE)
  expr <- harmonic_mean(x)
  expect_true(is_concave(expr))
})

## @cvxpy NONE
test_that("harmonic_mean numeric value matches CVXPY", {
  ## CVXPY: harmonic_mean([1, 4]) = 1.6
  c_val <- Constant(c(1, 4))
  expr <- harmonic_mean(c_val)
  expect_equal(scalar_value(value(expr)), 1.6, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("harmonic_mean solve parity with CVXPY", {
  ## Maximize harmonic_mean(x) s.t. x <= 4, sum(x) == 6
  ## CVXPY: optimal value ≈ 3.0, x ≈ [3, 3]
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Maximize(harmonic_mean(x)), list(x <= 4, sum_entries(x) == 6))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(val, 3.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════
# mixed_norm
# ═══════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("mixed_norm DCP: convex", {
  X <- Variable(c(2, 2))
  expr <- mixed_norm(X, 2, 1)
  expect_true(is_convex(expr))
})

## @cvxpy NONE
test_that("mixed_norm numeric value matches CVXPY", {
  ## CVXPY: mixed_norm([[3,4],[0,0]], 2, 1) = 5.0
  X <- Constant(matrix(c(3, 0, 4, 0), 2, 2))
  expr <- mixed_norm(X, 2, 1)
  expect_equal(scalar_value(value(expr)), 5.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("mixed_norm solve parity with CVXPY", {
  ## Minimize mixed_norm(X, 2, 1) s.t. X >= 1
  ## CVXPY: optimal value ≈ 2*sqrt(2) ≈ 2.8284
  X <- Variable(c(2, 2))
  prob <- Problem(Minimize(mixed_norm(X, 2, 1)), list(X >= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(val, 2 * sqrt(2), tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════
# cvar
# ═══════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cvar DCP: convex", {
  x <- Variable(4)
  expr <- cvar(x, 0.5)
  expect_true(is_convex(expr))
})

## @cvxpy NONE
test_that("cvar numeric value matches CVXPY", {
  ## CVXPY: cvar([1,2,3,4], 0.5) = 3.5
  c_val <- Constant(c(1, 2, 3, 4))
  expr <- cvar(c_val, 0.5)
  expect_equal(scalar_value(value(expr)), 3.5, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("cvar validation: beta < 0 errors", {
  x <- Variable(4)
  expect_error(cvar(x, -0.1), "beta")
})

## @cvxpy NONE
test_that("cvar validation: beta >= 1 errors", {
  x <- Variable(4)
  expect_error(cvar(x, 1.0), "beta")
  expect_error(cvar(x, 1.5), "beta")
})

## @cvxpy NONE
test_that("cvar validation: matrix input errors", {
  X <- Variable(c(2, 2))
  expect_error(cvar(X, 0.5), "vector")
})

## @cvxpy NONE
test_that("cvar works with row vector", {
  ## Row vector (1, 4) should also be accepted (min(shape) == 1)
  x <- Variable(c(1, 4))
  expr <- cvar(x, 0.5)
  expect_true(is_convex(expr))
})

## @cvxpy NONE
test_that("cvar beta=0 is valid (average of all elements)", {
  c_val <- Constant(c(1, 2, 3, 4))
  expr <- cvar(c_val, 0)
  ## k = (1-0)*4 = 4, sum_largest(x, 4)/4 = (1+2+3+4)/4 = 2.5
  expect_equal(scalar_value(value(expr)), 2.5, tolerance = 1e-4)
})
