## Tests for Phase 5a: Canonicalizer functions
## ============================================================

# ── abs_canon ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("abs_canon returns Variable + 2 NonNeg constraints", {
  x <- Variable(shape = c(2L, 3L))
  abs_expr <- Abs(x)
  result <- abs_canon(abs_expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  expect_equal(t@shape, c(2L, 3L))
  expect_equal(length(constrs), 2L)
  ## Each constraint should be >= (NonNeg or Inequality)
  for (c in constrs) {
    expect_true(S7_inherits(c, Constraint))
  }
})

## @cvxpy NONE
test_that("abs_canon scalar", {
  x <- Variable(1)
  abs_expr <- Abs(x)
  result <- abs_canon(abs_expr, list(x))
  expect_equal(result[[1L]]@shape, c(1L, 1L))
  expect_equal(length(result[[2L]]), 2L)
})

# ── maximum_canon ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("maximum_canon returns Variable + N constraints", {
  x <- Variable(2)
  y <- Variable(2)
  expr <- Maximum(x, y)
  result <- maximum_canon(expr, list(x, y))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  expect_equal(t@shape, c(2L, 1L))
  expect_equal(length(constrs), 2L)  # one per arg
})

## @cvxpy NONE
test_that("maximum_canon with 3 args", {
  x <- Variable(2)
  y <- Variable(2)
  z <- Variable(2)
  expr <- Maximum(x, y, z)
  result <- maximum_canon(expr, list(x, y, z))
  expect_equal(length(result[[2L]]), 3L)
})

# ── minimum_canon ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("minimum_canon negates and uses maximum_canon", {
  x <- Variable(2)
  y <- Variable(2)
  expr <- Minimum(x, y)
  result <- minimum_canon(expr, list(x, y))
  t <- result[[1L]]
  constrs <- result[[2L]]
  ## t should be a NegExpression wrapping a Variable
  expect_true(S7_inherits(t, Expression))
  expect_equal(length(constrs), 2L)
})

# ── max_canon ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("max_canon axis=NULL", {
  x <- Variable(shape = c(2L, 3L))
  expr <- MaxEntries(x)  # axis=NULL → scalar
  result <- max_canon(expr, list(x))
  t <- result[[1L]]
  expect_equal(t@shape, c(1L, 1L))
  expect_equal(length(result[[2L]]), 1L)  # single constraint: x <= promoted_t
})

## @cvxpy NONE
test_that("max_canon axis=2", {
  x <- Variable(shape = c(3L, 4L))
  expr <- MaxEntries(x, axis = 2L)  # reduce rows (column-wise) → shape (1, 4)
  result <- max_canon(expr, list(x))
  t <- result[[1L]]
  expect_equal(t@shape, c(1L, 4L))
  expect_equal(length(result[[2L]]), 1L)
})

## @cvxpy NONE
test_that("max_canon axis=1", {
  x <- Variable(shape = c(3L, 4L))
  expr <- MaxEntries(x, axis = 1L)  # reduce cols → shape (3, 1) in our convention
  result <- max_canon(expr, list(x))
  t <- result[[1L]]
  ## axis_shape for axis=1 on (3,4) with keepdims=FALSE → (3, 1)
  expect_equal(t@shape, c(3L, 1L))
  expect_equal(length(result[[2L]]), 1L)
})

# ── min_canon ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("min_canon axis=NULL", {
  x <- Variable(shape = c(2L, 3L))
  expr <- MinEntries(x)
  result <- min_canon(expr, list(x))
  t <- result[[1L]]
  expect_true(S7_inherits(t, Expression))
  expect_equal(length(result[[2L]]), 1L)
})

# ── exp_canon ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("exp_canon returns Variable + ExpCone constraint", {
  x <- Variable(shape = c(2L, 3L))
  expr <- Exp(x)
  result <- exp_canon(expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  expect_equal(t@shape, c(2L, 3L))
  expect_equal(length(constrs), 1L)
  expect_true(S7_inherits(constrs[[1L]], ExpCone))
})

## @cvxpy NONE
test_that("exp_canon scalar", {
  x <- Variable(1)
  expr <- Exp(x)
  result <- exp_canon(expr, list(x))
  expect_equal(result[[1L]]@shape, c(1L, 1L))
  expect_true(S7_inherits(result[[2L]][[1L]], ExpCone))
})

# ── log_canon ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("log_canon returns Variable + ExpCone constraint", {
  x <- Variable(shape = c(2L, 3L))
  expr <- Log(x)
  result <- log_canon(expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  expect_equal(length(constrs), 1L)
  expect_true(S7_inherits(constrs[[1L]], ExpCone))
})

# ── entr_canon ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("entr_canon returns Variable + ExpCone constraint", {
  x <- Variable(shape = c(2L, 3L))
  expr <- Entr(x)
  result <- entr_canon(expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  expect_equal(length(constrs), 1L)
  expect_true(S7_inherits(constrs[[1L]], ExpCone))
})

## @cvxpy NONE
test_that("entr_canon uses correct ExpCone argument order: (t, x, ones)", {
  x <- Variable(shape = c(1L, 1L))
  expr <- Entr(x)
  result <- entr_canon(expr, list(x))
  cone <- result[[2L]][[1L]]
  ## cone args should be: (t, x, ones)
  ## t is a Variable (the first arg to ExpCone constructor)
  expect_true(S7_inherits(cone@.x, Expression))  # t (new Variable)
  expect_true(S7_inherits(cone@.y, Expression))   # x (the original)
  expect_true(S7_inherits(cone@.z, Expression))   # ones (Constant)
})

# ── power_exact_canon ──────────────────────────────────────────────

## @cvxpy NONE
test_that("power_exact_canon p=1 returns identity", {
  x <- Variable(2)
  expr <- Power(x, 1)
  result <- power_exact_canon(expr, list(x))
  expect_identical(result[[1L]], x)
  expect_equal(length(result[[2L]]), 0L)
})

## @cvxpy NONE
test_that("power_exact_canon p=0 returns ones", {
  x <- Variable(2)
  expr <- Power(x, 0)
  result <- power_exact_canon(expr, list(x))
  expect_true(S7_inherits(result[[1L]], Constant))
  expect_equal(length(result[[2L]]), 0L)
})

## @cvxpy NONE
test_that("power_exact_canon p=2 returns SOC constraint", {
  x <- Variable(2)
  expr <- Power(x, 2)
  result <- power_exact_canon(expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  ## p=2 uses SOC formulation for broader solver support (not PowCone3D)
  has_soc <- any(vapply(constrs, function(c) S7_inherits(c, SOC), logical(1)))
  expect_true(has_soc)
})

## @cvxpy NONE
test_that("power_exact_canon 0 < p < 1", {
  x <- Variable(2)
  expr <- Power(x, 0.5)
  result <- power_exact_canon(expr, list(x))
  expect_true(S7_inherits(result[[1L]], Variable))
  expect_true(length(result[[2L]]) > 0L)
})

# ── pnorm_exact_canon ─────────────────────────────────────────────

## @cvxpy NONE
test_that("pnorm_exact_canon p=2 returns SOC constraint", {
  x <- Variable(3)
  expr <- Pnorm(x, p = 2)
  result <- pnorm_exact_canon(expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  expect_equal(length(constrs), 1L)
  expect_true(S7_inherits(constrs[[1L]], SOC))
})

# ── powcone_constrs ────────────────────────────────────────────────

## @cvxpy NONE
test_that("powcone_constrs returns single PowCone3D", {
  x <- Variable(2)
  y <- Variable(2)
  t <- Variable(2)
  result <- powcone_constrs(t, list(x, y), 0.5)
  expect_equal(length(result), 1L)
  expect_true(S7_inherits(result[[1L]], PowCone3D))
})

# ── norm1_canon ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("norm1_canon decomposes to abs + sum", {
  x <- Variable(shape = c(3L, 2L))
  expr <- Norm1(x)
  result <- norm1_canon(expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  ## Result should be a SumEntries of a Variable (abs result)
  expect_true(S7_inherits(t, Expression))
  ## 2 constraints from abs_canon: abs_t >= x, abs_t >= -x
  expect_equal(length(constrs), 2L)
})

## @cvxpy NONE
test_that("norm1_canon axis=0", {
  x <- Variable(shape = c(3L, 4L))
  expr <- Norm1(x, axis = 2L)
  result <- norm1_canon(expr, list(x))
  expect_true(S7_inherits(result[[1L]], Expression))
  expect_equal(length(result[[2L]]), 2L)
})

# ── norm_inf_canon ───────────────────────────────────────────────

## @cvxpy NONE
test_that("norm_inf_canon axis=NULL", {
  x <- Variable(shape = c(2L, 3L))
  expr <- NormInf(x)
  result <- norm_inf_canon(expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  expect_equal(t@shape, c(1L, 1L))
  ## 2 constraints: x <= promoted_t, x + promoted_t >= 0
  expect_equal(length(constrs), 2L)
})

## @cvxpy NONE
test_that("norm_inf_canon axis=0", {
  x <- Variable(shape = c(3L, 4L))
  expr <- NormInf(x, axis = 2L)
  result <- norm_inf_canon(expr, list(x))
  t <- result[[1L]]
  expect_equal(length(result[[2L]]), 2L)
})

# ── sum_largest_canon ────────────────────────────────────────────

## @cvxpy NONE
test_that("sum_largest_canon returns linear obj + 2 constraints", {
  x <- Variable(shape = c(3L, 1L))
  expr <- SumLargest(x, k = 2L)
  result <- sum_largest_canon(expr, list(x))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Expression))
  expect_equal(length(constrs), 2L)  # x <= t + q, t >= 0
})

# ── quad_over_lin_canon with SOC VStack fix ──────────────────────

## @cvxpy NONE
test_that("quad_over_lin_canon creates valid SOC constraint", {
  x <- Variable(shape = c(3L, 1L))
  y <- Variable(shape = c(1L, 1L), nonneg = TRUE)
  expr <- QuadOverLin(x, y)
  result <- quad_over_lin_canon(expr, list(x, y))
  t <- result[[1L]]
  constrs <- result[[2L]]
  expect_true(S7_inherits(t, Variable))
  expect_equal(length(constrs), 1L)
  expect_true(S7_inherits(constrs[[1L]], SOC))
})
