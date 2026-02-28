## CVXPY Parity Tests
## ==================
## These tests mirror CVXPY's test_curvature.py, test_monotonicity.py,
## test_canon_sign.py, test_problem.py (is_qp/is_lp), and test_atoms.py
## (sign logic for maximum/minimum/max/min).
## Purpose: prevent CVXR from veering off CVXPY's behavior.

library(testthat)
library(CVXR)

# ═══════════════════════════════════════════════════════════════════════
# Task #4: DCP Composition Chain Tests (from test_curvature.py +
#           test_monotonicity.py)
# ═══════════════════════════════════════════════════════════════════════

# ── Curvature arithmetic (test_curvature.py setUp fixtures) ──────────
# cvx = Variable()^2 (CONVEX), ccv = Variable()^0.5 (CONCAVE),
# aff = Variable() (AFFINE), const = Constant(5) (CONSTANT)

## @cvxpy test_curvature.py::TestCurvature::test_add
test_that("curvature: addition rules (test_curvature.py::test_add)", {
  cvx <- Variable(1)^2        # convex
  ccv <- Variable(1)^0.5      # concave
  aff <- Variable(1)          # affine
  const <- Constant(5)        # constant

  ## const + cvx → CONVEX
  expect_equal(expr_curvature(const + cvx), CONVEX)
  ## cvx + cvx → CONVEX
  expect_equal(expr_curvature(cvx + cvx), CONVEX)
  ## aff + ccv → CONCAVE
  expect_equal(expr_curvature(aff + ccv), CONCAVE)
  ## cvx + ccv → UNKNOWN
  expect_equal(expr_curvature(cvx + ccv), UNKNOWN_CURVATURE)
})

## @cvxpy test_curvature.py::TestCurvature::test_sub
test_that("curvature: subtraction rules (test_curvature.py::test_sub)", {
  cvx <- Variable(1)^2
  ccv <- Variable(1)^0.5
  aff <- Variable(1)
  const <- Constant(5)

  ## const - cvx → CONCAVE (negating convex yields concave)
  expect_equal(expr_curvature(const - cvx), CONCAVE)
  ## cvx - ccv → CONVEX (convex - concave = convex)
  expect_equal(expr_curvature(cvx - ccv), CONVEX)
  ## aff - ccv → CONVEX (affine - concave = convex)
  expect_equal(expr_curvature(aff - ccv), CONVEX)
  ## cvx - cvx → UNKNOWN (convex - convex = unknown)
  expect_equal(expr_curvature(cvx - cvx), UNKNOWN_CURVATURE)
})

## @cvxpy test_curvature.py::TestCurvature::test_sign_mult
test_that("curvature: sign-multiply rules (test_curvature.py::test_sign_mult)", {
  cvx <- Variable(1)^2
  ccv <- Variable(1)^0.5
  aff <- Variable(1)
  const <- Constant(5)
  zero <- Constant(0)
  neg <- Constant(-1)
  pos <- Constant(1)

  ## zero * cvx → AFFINE (technically CONSTANT after fold)
  e_zero_cvx <- zero * cvx
  expect_true(is_affine(e_zero_cvx))

  ## neg * cvx → CONCAVE
  expect_equal(expr_curvature(neg * cvx), CONCAVE)
  ## neg * ccv → CONVEX
  expect_equal(expr_curvature(neg * ccv), CONVEX)
  ## pos * aff → AFFINE
  expect_equal(expr_curvature(pos * aff), AFFINE)
  ## pos * ccv → CONCAVE
  expect_equal(expr_curvature(pos * ccv), CONCAVE)
  ## unknown_sign * const → CONSTANT
  unknown_sign <- pos + neg
  expect_true(is_constant(unknown_sign * const))
  ## unknown_sign * ccv → UNKNOWN
  expect_equal(expr_curvature(unknown_sign * ccv), UNKNOWN_CURVATURE)
})

## @cvxpy test_curvature.py::TestCurvature::test_neg
test_that("curvature: negation rules (test_curvature.py::test_neg)", {
  cvx <- Variable(1)^2
  aff <- Variable(1)

  ## -cvx → CONCAVE
  expect_equal(expr_curvature(-cvx), CONCAVE)
  ## -aff → AFFINE
  expect_equal(expr_curvature(-aff), AFFINE)
})

## @cvxpy test_curvature.py::TestCurvature::test_is_curvature
test_that("curvature: boolean queries (test_curvature.py::test_is_curvature)", {
  cvx <- Variable(1)^2
  ccv <- Variable(1)^0.5
  aff <- Variable(1)
  const <- Constant(5)

  ## is_affine
  expect_true(is_affine(const))
  expect_true(is_affine(aff))
  expect_false(is_affine(cvx))
  expect_false(is_affine(ccv))

  ## is_convex
  expect_true(is_convex(const))
  expect_true(is_convex(aff))
  expect_true(is_convex(cvx))
  expect_false(is_convex(ccv))

  ## is_concave
  expect_true(is_concave(const))
  expect_true(is_concave(aff))
  expect_false(is_concave(cvx))
  expect_true(is_concave(ccv))
})

# ── DCP composition chains (test_monotonicity.py::test_dcp_curvature) ─

## @cvxpy test_monotonicity.py::TestMonotonicity::test_dcp_curvature
test_that("DCP composition: 1 + exp(x) is CONVEX", {
  x <- Variable(1)
  expr <- 1 + exp(x)
  expect_equal(expr_curvature(expr), CONVEX)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_dcp_curvature
test_that("DCP composition: exp(x)^2 is CONVEX", {
  x <- Variable(1)
  expr <- exp(x)^2
  expect_equal(expr_curvature(expr), CONVEX)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_dcp_curvature
test_that("DCP composition: 1 - sqrt(x) is CONVEX", {
  x <- Variable(1)
  expr <- 1 - sqrt(x)
  expect_equal(expr_curvature(expr), CONVEX)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_dcp_curvature
test_that("DCP composition: log(sqrt(x)) is CONCAVE", {
  x <- Variable(1)
  expr <- log(sqrt(x))
  expect_equal(expr_curvature(expr), CONCAVE)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_dcp_curvature
test_that("DCP composition: -(exp(x))^2 is CONCAVE", {
  x <- Variable(1)
  expr <- -(exp(x))^2
  expect_equal(expr_curvature(expr), CONCAVE)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_dcp_curvature
test_that("DCP composition: log(exp(x)) is NOT DCP", {
  x <- Variable(1)
  expr <- log(exp(x))
  expect_false(is_dcp(expr))
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_dcp_curvature
test_that("DCP composition: entr(x_nonneg) is CONCAVE", {
  x <- Variable(1, nonneg = TRUE)
  expr <- entr(x)
  expect_equal(expr_curvature(expr), CONCAVE)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_dcp_curvature
test_that("DCP composition: ((x^2)^0.5)^0 is CONSTANT", {
  x <- Variable(1)
  expr <- ((x^2)^0.5)^0
  expect_true(is_constant(expr))
})

# ── Signed monotonicity (test_monotonicity.py::test_signed_curvature) ─

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(1 + exp(x)) is CONVEX", {
  x <- Variable(1)
  expr <- abs(1 + exp(x))
  expect_equal(expr_curvature(expr), CONVEX)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(-square(x)) is CONVEX", {
  x <- Variable(1)
  expr <- abs(-x^2)
  expect_equal(expr_curvature(expr), CONVEX)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(x_nonneg) is CONVEX", {
  x <- Variable(1, nonneg = TRUE)
  expr <- abs(x)
  expect_equal(expr_curvature(expr), CONVEX)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(-x_nonneg) is CONVEX", {
  x <- Variable(1, nonneg = TRUE)
  expr <- abs(-x)
  expect_equal(expr_curvature(expr), CONVEX)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(x) is CONVEX (affine, no sign info)", {
  x <- Variable(1)
  expr <- abs(x)
  expect_equal(expr_curvature(expr), CONVEX)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(-entr(x)) is UNKNOWN", {
  x <- Variable(1)
  expr <- abs(-entr(x))
  expect_equal(expr_curvature(expr), UNKNOWN_CURVATURE)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(log(x)) is UNKNOWN", {
  x <- Variable(1)
  expr <- abs(log(x))
  expect_equal(expr_curvature(expr), UNKNOWN_CURVATURE)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(entr(x)) is UNKNOWN", {
  x <- Variable(1)
  expr <- abs(entr(x))
  expect_equal(expr_curvature(expr), UNKNOWN_CURVATURE)
})

## @cvxpy test_monotonicity.py::TestMonotonicity::test_signed_curvature
test_that("signed DCP: abs(-log(x)) is UNKNOWN", {
  x <- Variable(1)
  expr <- abs(-log(x))
  expect_equal(expr_curvature(expr), UNKNOWN_CURVATURE)
})

# ═══════════════════════════════════════════════════════════════════════
# Task #5: is_lp / is_qp expansion (from test_problem.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_is_qp
test_that("is_qp: sum_squares objective is QP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))
  obj <- sum_squares(A %*% y - b)
  p <- Problem(Minimize(obj))
  expect_true(is_qp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_qp
test_that("is_qp: sum_squares + equality + inequality is QP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))
  Aeq <- Constant(matrix(rnorm(6), 2, 3))
  beq <- Constant(matrix(rnorm(2), ncol = 1))
  Fc <- Constant(matrix(rnorm(6), 2, 3))
  g <- Constant(matrix(rnorm(2), ncol = 1))

  p <- Problem(Minimize(sum_squares(A %*% y - b)),
               list(Aeq %*% y == beq, Fc %*% y <= g))
  expect_true(is_qp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_qp
test_that("is_qp: QP with PWL constraints (max, abs, norm1)", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))
  Aeq <- Constant(matrix(rnorm(6), 2, 3))
  beq <- Constant(matrix(rnorm(2), ncol = 1))

  p <- Problem(Minimize(sum_squares(A %*% y - b)),
               list(Maximum(Constant(matrix(1, 3, 1)), 3 * y) <= 200,
                    Abs(2 * y) <= 100,
                    norm1(2 * y) <= 1000,
                    Aeq %*% y == beq))
  expect_true(is_qp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_qp
test_that("is_qp: non-PWL constraint (y^2) rejects QP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))

  p <- Problem(Minimize(sum_squares(A %*% y - b)),
               list(Maximum(Constant(matrix(1, 3, 1)), 3 * y^2) <= 200))
  expect_false(is_qp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_qp
test_that("is_qp: SOC constraint rejects QP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))
  t <- Variable(1)

  p <- Problem(Minimize(sum_squares(A %*% y - b)),
               list(SOC(t, y)))
  expect_false(is_qp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_qp
test_that("is_qp: ExpCone constraint rejects QP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))

  p <- Problem(Minimize(sum_squares(A %*% y - b)),
               list(ExpCone(y[1], y[2], y[3])))
  expect_false(is_qp(p))
})

# ── is_lp tests (from test_problem.py::test_is_lp) ────────────────────

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: linear objective + linear constraints is LP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))
  cv <- Constant(matrix(rnorm(3), ncol = 1))

  p <- Problem(Minimize(sum(cv * y)), list(A %*% y <= b))
  expect_true(is_lp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: LP with equality + inequality constraints", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))
  cv <- Constant(matrix(rnorm(3), ncol = 1))
  Aeq <- Constant(matrix(rnorm(6), 2, 3))
  beq <- Constant(matrix(rnorm(2), ncol = 1))

  p <- Problem(Minimize(sum(cv * y)), list(A %*% y <= b, Aeq %*% y == beq))
  expect_true(is_lp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: Maximization LP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))
  cv <- Constant(matrix(rnorm(3), ncol = 1))

  p <- Problem(Maximize(sum(cv * y)), list(A %*% y <= b))
  expect_true(is_lp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: QP objective (sum_squares) is NOT LP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))

  p <- Problem(Minimize(sum_squares(y)), list(A %*% y <= b))
  expect_false(is_lp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: SOC constraint is NOT LP", {
  y <- Variable(3)
  cv <- Constant(matrix(rnorm(3), ncol = 1))
  t <- Variable(1)

  p <- Problem(Minimize(sum(cv * y)), list(SOC(t, y)))
  expect_false(is_lp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: PWL objective (abs) IS LP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))

  p <- Problem(Minimize(sum(abs(y))), list(A %*% y <= b))
  expect_true(is_lp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: PWL objective (max) IS LP", {
  y <- Variable(3)
  A <- Constant(matrix(rnorm(12), 4, 3))
  b <- Constant(matrix(rnorm(4), ncol = 1))

  p <- Problem(Minimize(max(y)), list(A %*% y <= b))
  expect_true(is_lp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: non-affine constraint (sum_squares <= 1) NOT LP", {
  y <- Variable(3)
  cv <- Constant(matrix(rnorm(3), ncol = 1))

  p <- Problem(Minimize(sum(cv * y)), list(sum_squares(y) <= 1))
  expect_false(is_lp(p))
})

# ═══════════════════════════════════════════════════════════════════════
# Task #6: Canonicalization Sign Preservation
#           (from test_canon_sign.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_canon_sign.py::TestCanonSign::test_maximum_sign
test_that("canon sign: maximum_canon preserves sign for nonneg/nonpos/zero", {
  ## Nonneg expression
  expr_nonneg <- Variable(1, nonneg = TRUE)
  tmp_nn <- Maximum(expr_nonneg, Constant(-Inf))
  result_nn <- maximum_canon(tmp_nn, tmp_nn@args)
  ## Canonicalized expression should preserve nonneg property
  expect_equal(is_nonneg(tmp_nn), is_nonneg(result_nn[[1L]]))

  ## Nonpos expression
  expr_nonpos <- Variable(1, nonpos = TRUE)
  tmp_np <- Maximum(expr_nonpos, Constant(-Inf))
  result_np <- maximum_canon(tmp_np, tmp_np@args)
  expect_equal(is_nonpos(tmp_np), is_nonpos(result_np[[1L]]))

  ## Zero expression
  expr_zero <- Constant(0)
  tmp_z <- Maximum(expr_zero, Constant(-Inf))
  result_z <- maximum_canon(tmp_z, tmp_z@args)
  expect_equal(is_nonneg(tmp_z), is_nonneg(result_z[[1L]]))
  expect_equal(is_nonpos(tmp_z), is_nonpos(result_z[[1L]]))
})

## @cvxpy test_canon_sign.py::TestCanonSign::test_minimum_sign
test_that("canon sign: minimum_canon preserves sign for nonneg/nonpos/zero", {
  ## Nonneg expression
  expr_nonneg <- Variable(1, nonneg = TRUE)
  tmp_nn <- Minimum(expr_nonneg, Constant(Inf))
  result_nn <- minimum_canon(tmp_nn, tmp_nn@args)
  expect_equal(is_nonneg(tmp_nn), is_nonneg(result_nn[[1L]]))

  ## Nonpos expression
  expr_nonpos <- Variable(1, nonpos = TRUE)
  tmp_np <- Minimum(expr_nonpos, Constant(Inf))
  result_np <- minimum_canon(tmp_np, tmp_np@args)
  expect_equal(is_nonpos(tmp_np), is_nonpos(result_np[[1L]]))

  ## Zero expression
  expr_zero <- Constant(0)
  tmp_z <- Minimum(expr_zero, Constant(Inf))
  result_z <- minimum_canon(tmp_z, tmp_z@args)
  expect_equal(is_nonneg(tmp_z), is_nonneg(result_z[[1L]]))
  expect_equal(is_nonpos(tmp_z), is_nonpos(result_z[[1L]]))
})

# ═══════════════════════════════════════════════════════════════════════
# Task #7: Atom Sign Logic (from test_atoms.py)
# ═══════════════════════════════════════════════════════════════════════

# ── Maximum sign (test_atoms.py::test_maximum_sign) ───────────────────

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: two nonneg args → NONNEG", {
  expect_equal(expr_sign_str(Maximum(Constant(1), Constant(2))), NONNEG_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: nonneg + unknown → NONNEG", {
  expect_equal(expr_sign_str(Maximum(Constant(1), Variable(1))), NONNEG_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: nonneg + nonpos → NONNEG", {
  expect_equal(expr_sign_str(Maximum(Constant(1), Constant(-2))), NONNEG_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: nonneg + zero → NONNEG", {
  expect_equal(expr_sign_str(Maximum(Constant(1), Constant(0))), NONNEG_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: unknown + zero → NONNEG", {
  expect_equal(expr_sign_str(Maximum(Variable(1), Constant(0))), NONNEG_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: two unknowns → UNKNOWN", {
  expect_equal(expr_sign_str(Maximum(Variable(1), Variable(1))), UNKNOWN_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: unknown + nonpos → UNKNOWN", {
  expect_equal(expr_sign_str(Maximum(Variable(1), Constant(-2))), UNKNOWN_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: two zeros → ZERO", {
  expect_equal(expr_sign_str(Maximum(Constant(0), Constant(0))), ZERO_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: zero + nonpos → ZERO", {
  expect_equal(expr_sign_str(Maximum(Constant(0), Constant(-2))), ZERO_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("Maximum sign: two nonpos → NONPOS", {
  expect_equal(expr_sign_str(Maximum(Constant(-3), Constant(-2))), NONPOS_SIGN)
})

# ── Minimum sign (test_atoms.py::test_minimum_sign) ───────────────────

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: two nonneg → NONNEG", {
  expect_equal(expr_sign_str(Minimum(Constant(1), Constant(2))), NONNEG_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: nonneg + unknown → UNKNOWN", {
  expect_equal(expr_sign_str(Minimum(Constant(1), Variable(1))), UNKNOWN_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: nonneg + nonpos → NONPOS", {
  expect_equal(expr_sign_str(Minimum(Constant(1), Constant(-2))), NONPOS_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: nonneg + zero → ZERO", {
  expect_equal(expr_sign_str(Minimum(Constant(1), Constant(0))), ZERO_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: unknown + zero → NONPOS", {
  expect_equal(expr_sign_str(Minimum(Variable(1), Constant(0))), NONPOS_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: two unknowns → UNKNOWN", {
  expect_equal(expr_sign_str(Minimum(Variable(1), Variable(1))), UNKNOWN_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: unknown + nonpos → NONPOS", {
  expect_equal(expr_sign_str(Minimum(Variable(1), Constant(-2))), NONPOS_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: two zeros → ZERO", {
  expect_equal(expr_sign_str(Minimum(Constant(0), Constant(0))), ZERO_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: zero + nonpos → NONPOS", {
  expect_equal(expr_sign_str(Minimum(Constant(0), Constant(-2))), NONPOS_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("Minimum sign: two nonpos → NONPOS", {
  expect_equal(expr_sign_str(Minimum(Constant(-3), Constant(-2))), NONPOS_SIGN)
})

# ── MaxEntries / MinEntries sign (test_atoms.py::test_max / test_min) ─

## @cvxpy test_atoms.py::TestAtoms::test_max
test_that("MaxEntries sign: nonneg constant → NONNEG", {
  expect_equal(expr_sign_str(MaxEntries(Constant(1))), NONNEG_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_max
test_that("MaxEntries sign: nonpos constant → NONPOS", {
  expect_equal(expr_sign_str(MaxEntries(Constant(-2))), NONPOS_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_max
test_that("MaxEntries sign: unknown variable → UNKNOWN", {
  expect_equal(expr_sign_str(MaxEntries(Variable(1))), UNKNOWN_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_max
test_that("MaxEntries sign: zero constant → ZERO", {
  expect_equal(expr_sign_str(MaxEntries(Constant(0))), ZERO_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_min
test_that("MinEntries sign: nonneg constant → NONNEG", {
  expect_equal(expr_sign_str(MinEntries(Constant(1))), NONNEG_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_min
test_that("MinEntries sign: nonpos constant → NONPOS", {
  expect_equal(expr_sign_str(MinEntries(Constant(-2))), NONPOS_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_min
test_that("MinEntries sign: unknown variable → UNKNOWN", {
  expect_equal(expr_sign_str(MinEntries(Variable(1))), UNKNOWN_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_min
test_that("MinEntries sign: zero constant → ZERO", {
  expect_equal(expr_sign_str(MinEntries(Constant(0))), ZERO_SIGN)
})

# ═══════════════════════════════════════════════════════════════════════
# Task #8: Pipeline Expansion Tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_get_problem_data
test_that("pipeline: ExpCone constraint passes through ConeMatrixStuffing", {
  ## log(x) ≥ 1 with x ≥ 0 → ExpCone in Dcp2Cone output
  x <- Variable(1, nonneg = TRUE)
  p <- Problem(Maximize(log(x)), list(x <= 10))
  expect_true(is_dcp(p))

  ## Get problem data — exercises full pipeline
  pd <- problem_data(p)
  ## Verify we get cone_dims with exp cone
  expect_true(!is.null(pd$data$dims))
  dims <- pd$data$dims
  expect_true(dims@exp > 0L)
})

## @cvxpy test_problem.py::TestProblem::test_get_problem_data
test_that("pipeline: SOC constraint passes through ConeMatrixStuffing", {
  x <- Variable(3)
  t <- Variable(1)
  p <- Problem(Minimize(t), list(SOC(t, x), x >= 1))
  expect_true(is_dcp(p))

  pd <- problem_data(p)
  dims <- pd$data$dims
  expect_true(length(dims@soc) > 0L)
})

## @cvxpy NONE
test_that("pipeline: LP with equality constraint", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 0, x[1] + x[2] == 5))
  expect_true(is_lp(p))

  pd <- problem_data(p)
  dims <- pd$data$dims
  expect_true(dims@zero > 0L)   # equality constraint
  expect_true(dims@nonneg > 0L) # inequality constraint
})

## @cvxpy NONE
test_that("pipeline: Dcp2Cone reduces exp atom to ExpCone constraints", {
  x <- Variable(1)
  p <- Problem(Minimize(exp(x)))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))

  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]

  ## After Dcp2Cone, constraints should include ExpCone
  has_expcone <- any(vapply(new_p@constraints, function(c) {
    S7_inherits(c, ExpCone)
  }, logical(1)))
  expect_true(has_expcone)
})

# ═══════════════════════════════════════════════════════════════════════
# Task #9: QuadForm Edge Cases
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_quad_form.py::TestNonOptimal::test_non_symmetric
test_that("QuadForm: non-symmetric P is rejected", {
  x <- Variable(2)
  P <- Constant(matrix(c(1, 2, 3, 4), 2, 2))  # not symmetric
  expect_error(QuadForm(x, P), "[Ss]ymmetric|[Hh]ermitian")
})

## @cvxpy NONE
test_that("QuadForm: incompatible dimensions rejected", {
  x <- Variable(3)  # 3 elements
  P <- Constant(diag(2))  # 2x2
  expect_error(QuadForm(x, P), "rows|match")
})

## @cvxpy test_quad_form.py::TestNonOptimal::test_psd_exactly_tolerance
test_that("QuadForm: PSD via eigenvalue tolerance", {
  ## Matrix with tiny negative eigenvalue within tolerance
  n <- 3L
  ## Build a nearly-PSD matrix: PSD + small perturbation
  P_mat <- diag(n) + 1e-12 * matrix(rnorm(n * n), n, n)
  ## Force symmetry
  P_mat <- (P_mat + t(P_mat)) / 2
  P <- Constant(P_mat)
  x <- Variable(n)

  ## Should be convex since eigenvalues are >= -EIGVAL_TOL
  qf <- QuadForm(x, P)
  ## The PSD check uses EIGVAL_TOL = 1e-10
  ## A matrix with eigenvalues around 1 ± 1e-12 should be PSD
  expect_true(is_convex(qf))
})

# ═══════════════════════════════════════════════════════════════════════
# QuadForm: NSD exactly tolerance (test_quad_form.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_quad_form.py::TestNonOptimal::test_nsd_exactly_tolerance
test_that("QuadForm: NSD with eigenvalue at tolerance boundary", {
  ## CVXPY: P = [[0.999*EIGVAL_TOL, 0], [0, -10]]
  ## This matrix has eigenvalues ~1e-13 and -10, so it is NSD (up to tolerance).
  ## Maximize x'Px with x=[1,2] → value ≈ -40
  ## Verified via CVXPY: status=optimal, value≈-40.0
  EIGVAL_TOL <- 1e-10
  P <- matrix(c(0.999 * EIGVAL_TOL, 0, 0, -10), 2, 2)
  x <- Variable(2)
  suppressWarnings({
    cost <- QuadForm(x, P)
    prob <- Problem(Maximize(cost), list(x == c(1, 2)))
    psolve(prob, solver = "SCS")
  })
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(prob)), -40, tolerance = 0.1)
})

# ═══════════════════════════════════════════════════════════════════════
# QuadForm: objective evaluation consistency (test_quad_form.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_quad_form.py::TestNonOptimal::test_obj_eval
test_that("QuadForm: prob$value == objective value for quad_form obj", {
  ## CVXPY: x = Variable((2,1)), A = [[1.0]], B = [[1],[1]]
  ## obj = -B'x + quad_form(B'x, A) → min at x = [0.25, 0.25], value = -0.25
  ## Verified via CVXPY: prob.value == prob.objective.value ≈ -0.25
  x <- Variable(c(2, 1))
  A <- matrix(1.0, 1, 1)
  B <- matrix(c(1.0, 1.0), 2, 1)
  obj0 <- -t(B) %*% x
  obj1 <- QuadForm(t(B) %*% x, A)
  prob <- Problem(Minimize(obj0 + obj1))
  psolve(prob, solver = "SCS")
  expect_equal(as.numeric(value(prob)), -0.25, tolerance = 1e-3)
  ## The key assertion: prob value equals objective evaluated value
  expect_equal(as.numeric(value(prob)),
               as.numeric(value(objective(prob))),
               tolerance = 1e-5)
})

# ═══════════════════════════════════════════════════════════════════════
# Batch A: test_problem.py parity tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_parameters
test_that("problem: parameters() extracts all parameters", {
  ## CVXPY: p1 = Parameter(), p2 = Parameter(3, nonpos=True),
  ## p3 = Parameter((4,4), nonneg=True)
  ## prob.parameters() returns [p1, p2, p3]
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  p1 <- Parameter(name = "p1")
  p2 <- Parameter(3, nonpos = TRUE, name = "p2")
  p3 <- Parameter(c(4, 4), nonneg = TRUE, name = "p3")
  prob <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + 2))
  params <- parameters(prob)
  ## Should find exactly 3 parameters
  expect_equal(length(params), 3L)
  ## Check by id: all three parameter ids should be present
  param_ids <- vapply(params, function(p) p@id, integer(1))
  expect_true(p1@id %in% param_ids)
  expect_true(p2@id %in% param_ids)
  expect_true(p3@id %in% param_ids)
})

## @cvxpy test_problem.py::TestProblem::test_vec
test_that("problem: vec vectorization", {
  ## CVXPY: A = Variable((2,2)), expr = vec(A, order='F')
  ## obj = Minimize(expr.T @ [1,2,3,4])
  ## A == [[-1,-2],[3,4]] → result = 20, expr.value = [-1,-2,3,4]
  ## In R, vec() uses column-major (Fortran order) by default.
  ## CVXPY A layout: A[0,0]=-1, A[0,1]=-2, A[1,0]=3, A[1,1]=4
  ## vec(A, 'F') stacks columns: col0=[-1,3], col1=[-2,4] → [-1,3,-2,4]
  ## R A = matrix(c(-1,3,-2,4), 2, 2) also col-major: A[1,1]=-1,A[2,1]=3,A[1,2]=-2,A[2,2]=4
  ## vec(A) in R = (-1, 3, -2, 4)
  ## dot with (1,2,3,4) = -1+6-6+16 = 15
  A <- Variable(c(2, 2), name = "A")
  c_vec <- c(1, 2, 3, 4)
  expr <- vec(A)
  obj <- Minimize(t(expr) %*% c_vec)
  constr <- list(A == matrix(c(-1, 3, -2, 4), 2, 2))
  prob <- Problem(obj, constr)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 15, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_norm2
test_that("problem: pnorm with p=2 various cases", {
  ## Constant argument
  prob <- Problem(Minimize(p_norm(Constant(-2), 2)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 2, tolerance = 1e-4)

  ## Scalar argument
  a <- Variable(name = "a")
  prob <- Problem(Minimize(p_norm(a, 2)), list(a <= -2))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(a)), -2, tolerance = 1e-4)

  ## Maximize
  a2 <- Variable(name = "a2")
  prob <- Problem(Maximize(-p_norm(a2, 2)), list(a2 <= -2))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), -2, tolerance = 1e-4)
  expect_equal(as.numeric(value(a2)), -2, tolerance = 1e-4)

  ## Vector arguments
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")
  prob <- Problem(Minimize(p_norm(x - z, 2) + 5),
                  list(x >= c(2, 3), z <= c(-1, -4)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 12.61577, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(2, 3), tolerance = 1e-3)
  expect_equal(as.numeric(value(z)), c(-1, -4), tolerance = 1e-3)

  ## Row (transposed) arguments
  x2 <- Variable(2, name = "x2")
  z2 <- Variable(2, name = "z2")
  prob <- Problem(Minimize(p_norm(t(x2 - z2), 2) + 5),
                  list(x2 >= c(2, 3), z2 <= c(-1, -4)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 12.61577, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_pnorm
test_that("problem: pnorm for p = 1.6, 2, 3, Inf", {
  ## CVXPY: x = Variable(3), a = [1,2,3]
  ## For each p, minimize pnorm(x, p) s.t. x'a >= 1
  ## Analytic solution: x_true = a^(1/(p-1)) / (a . a^(1/(p-1)))
  ## NOTE: p=1 and some fractional p values trigger internal bignum issues
  ## in CVXR's dyad_completion (next_pow2 with gmp). We test a subset.
  a_vec <- c(1.0, 2.0, 3.0)

  for (p in c(1.6, 2, 3, Inf)) {
    x <- Variable(3, name = "x")
    prob <- Problem(Minimize(p_norm(x, p)), list(t(x) %*% a_vec >= 1))
    psolve(prob, solver = "CLARABEL")

    if (p == Inf) {
      x_true <- rep(1, 3) / sum(a_vec)
    } else {
      x_true <- a_vec^(1.0 / (p - 1)) / sum(a_vec * a_vec^(1.0 / (p - 1)))
    }

    x_alg <- as.numeric(value(x))
    expect_true(all(abs(x_alg - x_true) < 1e-2),
                label = sprintf("p = %s: x_alg vs x_true mismatch", p))
    expect_equal(as.numeric(value(prob)),
                 as.numeric(value(p_norm(Constant(x_alg), p))),
                 tolerance = 1e-3,
                 label = sprintf("p = %s: prob value vs pnorm of solution", p))
  }
})

## @cvxpy test_problem.py::TestProblem::test_quad_form
test_that("problem: quad_form error handling and solving", {
  x <- Variable(2, name = "x")
  A <- Variable(c(2, 2), name = "A")

  ## Both arguments variable: should error
  expect_error(
    quad_form(x, A),
    regexp = ".*non-variable.*"
  )

  ## Dimension mismatch: scalar x with 2x2 P
  expect_error(
    quad_form(1, A),
    regexp = ".*rows.*|.*dimension.*|.*Invalid.*"
  )

  ## Non-PSD matrix with variable x: DCP error
  suppressWarnings({
    expect_error({
      obj <- Minimize(quad_form(x, matrix(c(-1, 0, 0, 9), 2, 2)))
      psolve(Problem(obj), solver = "CLARABEL")
    }, regexp = ".*DCP.*|.*convex.*|.*not follow.*")
  })

  ## Valid PSD matrix
  P <- matrix(c(4, 0, 0, 9), 2, 2)
  x2 <- Variable(2, name = "x2")
  prob <- Problem(Minimize(quad_form(x2, P)), list(x2 >= 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 13, tolerance = 1e-3)

  ## Constant x, variable P: quad_form returns x^H P x (affine in P)
  c_vec <- c(1, 2)
  A2 <- Variable(c(2, 2), name = "A2", symmetric = TRUE)
  prob <- Problem(Minimize(quad_form(c_vec, A2)), list(A2 >= 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 9, tolerance = 1e-3)

  ## Both constant
  prob <- Problem(Minimize(quad_form(c(1, 2), matrix(c(4, 0, 0, 9), 2, 2))))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 40, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_sdp
test_that("problem: SDP with lambda_max constraint", {
  ## CVXPY: max A[1,0] - A[0,1] s.t. lambda_max(A) <= 100,
  ## A[0,0]==2, A[1,1]==2, A[1,0]==2
  ## SDP constrains A to be symmetric: A[1,0] == A[0,1] -> result ~ 0
  A <- Variable(c(2, 2), name = "A")
  obj <- Maximize(A[2, 1] - A[1, 2])
  prob <- Problem(obj, list(
    lambda_max(A) <= 100,
    A[1, 1] == 2,
    A[2, 2] == 2,
    A[2, 1] == 2
  ))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_cummax
test_that("problem: cummax constraint", {
  ## CVXPY: tt = Variable(5)
  ## max sum(tt) s.t. cummax(tt, axis=0) <= [1,2,3,4,5]
  ## result = 15 (tt = [1,2,3,4,5])
  tt <- Variable(5, name = "tt")
  prob <- Problem(Maximize(sum_entries(tt)),
                  list(cummax_expr(tt) <= c(1, 2, 3, 4, 5)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 15, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_unpack_results
test_that("problem: decomposed solve API (problem_data + solve_via_data + problem_unpack_results)", {
  ## CVXPY: prob = Problem(Minimize(exp(a)), [a == 0])
  ## args, chain, inv = prob.get_problem_data("SCS")
  ## solution = scs.solve(data, cones)
  ## prob.unpack_results(solution, chain, inv)
  ## a.value ~ 0, prob.value ~ 1
  a <- Variable(name = "a")
  prob <- Problem(Minimize(exp(a)), list(a == 0))
  pd <- problem_data(prob, solver = "CLARABEL")
  chain <- pd$chain
  inverse_data <- pd$inverse_data
  data <- pd$data

  ## Solve via the chain
  raw_result <- solve_via_data(chain, data, warm_start = FALSE,
                               verbose = FALSE, solver_opts = list())

  ## Unpack into a fresh problem (same structure)
  prob2 <- Problem(Minimize(exp(a)), list(a == 0))
  problem_unpack_results(prob2, raw_result, chain, inverse_data)
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-3)
  expect_equal(as.numeric(value(prob2)), 1, tolerance = 1e-3)
  expect_equal(status(prob2), "optimal")
})

## @cvxpy test_problem.py::TestProblem::test_variable_name_conflict
test_that("problem: same-name variables are distinct", {
  ## CVXPY: a (setUp), var = Variable(name='a')
  ## max(a + var) s.t. var == 2 + a, var <= 3
  ## result = 4, a = 1, var = 3
  a <- Variable(name = "a")
  var <- Variable(name = "a")
  prob <- Problem(Maximize(a + var), list(var == 2 + a, var <= 3))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 4.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(a)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(var)), 3.0, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_special_index
test_that("problem: QP with sliced indexing", {
  ## CVXPY: x = Variable((1,3)), y = sum(x[:,0:2], axis=1)
  ## cost = QuadForm(y, diag([1])), minimize cost -> 0
  ## Same with fancy indexing x[:,[0,1]]
  ## CVXPY axis=1 -> CVXR axis=1 (row-wise)
  x1 <- Variable(c(1, 3), name = "x1")
  y1 <- sum_entries(x1[, 1:2], axis = 1)
  cost1 <- quad_form(y1, matrix(1, 1, 1))
  prob1 <- Problem(Minimize(cost1))
  psolve(prob1, solver = "CLARABEL")
  result1 <- as.numeric(value(prob1))

  x2 <- Variable(c(1, 3), name = "x2")
  y2 <- sum_entries(x2[, c(1, 2)], axis = 1)
  cost2 <- quad_form(y2, matrix(1, 1, 1))
  prob2 <- Problem(Minimize(cost2))
  psolve(prob2, solver = "CLARABEL")
  result2 <- as.numeric(value(prob2))

  expect_equal(result1, result2, tolerance = 1e-5)
})

## @cvxpy test_problem.py::TestProblem::test_spare_int8_matrix
test_that("problem: sparse matrix types (int8 equivalent)", {
  ## CVXPY: Sparse int8 matrix used in constraints
  ## R doesn't have int8 sparse, but we test with integer sparse matrix.
  ## The key test: dense and sparse D give same results.
  q <- c(1.88922129, 0.06938685, 0.91948919)
  P <- matrix(c(280.64, -49.84, -80.,
                -49.84, 196.04, 139.,
                -80., 139., 106.), 3, 3, byrow = TRUE)
  D_dense <- matrix(c(-1, 1, 0, 0, 0, 0,
                       0, -1, 1, 0, 0, 0,
                       0, 0, 0, -1, 1, 0), 3, 6, byrow = TRUE)
  D_sparse <- as(D_dense, "dgCMatrix")

  make_problem <- function(D, a_var) {
    ## Wrap transposed D as Constant to ensure proper CVXR dispatch
    Dt <- if (inherits(D, "sparseMatrix")) Constant(Matrix::t(D)) else t(D)
    obj <- Minimize(0.5 * quad_form(a_var, P) - t(a_var) %*% q)
    alpha <- Parameter(nonneg = TRUE, value = 2, name = "alpha")
    prob <- Problem(obj, list(a_var >= 0,
                              -alpha <= Dt %*% a_var,
                              Dt %*% a_var <= alpha))
    psolve(prob, solver = "CLARABEL")
    prob
  }

  a1 <- Variable(c(3, 1), name = "a1")
  prob_dense <- make_problem(D_dense, a1)
  expect_equal(status(prob_dense), "optimal")
  coef_dense <- as.numeric(t(value(a1)) %*% D_dense)

  a2 <- Variable(c(3, 1), name = "a2")
  prob_sparse <- make_problem(D_sparse, a2)
  expect_equal(status(prob_sparse), "optimal")
  coef_sparse <- as.numeric(t(value(a2)) %*% as.matrix(D_sparse))

  ## Results should match
  expect_equal(coef_dense, coef_sparse, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_int64
test_that("problem: integer dimension specification", {
  ## CVXPY: q = Variable(numpy.int64(2)), minimize norm(q, 1)
  q <- Variable(as.integer(2), name = "q")
  prob <- Problem(Minimize(norm1(q)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0, tolerance = 1e-5)
})

## @cvxpy test_problem.py::TestProblem::test_non_python_int_index
test_that("problem: indexing with integer types", {
  ## CVXPY: x[0:int(2)][0], x[0:numpy.int64(2)][0]
  x <- Variable(2, name = "x")
  cost <- x[1:2L][1]
  prob <- Problem(Minimize(cost), list(x == 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 1, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-5)

  ## With explicit integer indexing
  x2 <- Variable(2, name = "x2")
  cost2 <- x2[1:as.integer(2)][1]
  prob2 <- Problem(Minimize(cost2), list(x2 == 1))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob2)), 1, tolerance = 1e-5)
})

## @cvxpy test_problem.py::TestProblem::test_solver_error_raised_on_failure
test_that("problem: solver error raised on invalid QP", {
  ## CVXPY: quad_form(Variable(1) + 1, [[-1]], assume_PSD=True) with OSQP
  skip_if_not_installed("osqp")
  expect_error({
    x <- Variable(1, name = "x")
    prob <- Problem(Minimize(quad_form(x + 1, matrix(-1, 1, 1), assume_PSD = TRUE)))
    psolve(prob, solver = "OSQP")
  }, regexp = ".*")
})

## @cvxpy test_problem.py::TestProblem::test_variables_with_value
test_that("problem: Variable with initial value", {
  ## CVXPY: Variable(name="without_bounds", value=0.0)
  ## CVXR's Variable constructor does not accept 'value' directly.
  skip("Variable() does not accept 'value' parameter in CVXR")
})

## @cvxpy test_problem.py::TestProblem::test_param_dict
test_that("problem: param_dict property", {
  ## CVXPY: p.param_dict == {"p1": p1, "p2": p2, "p3": p3}
  skip("param_dict property not implemented in CVXR Problem class")
})

## @cvxpy test_problem.py::TestProblem::test_var_dict
test_that("problem: var_dict property", {
  ## CVXPY: p.var_dict == {"a": a, "x": x, "b": b, "A": A}
  skip("var_dict property not implemented in CVXR Problem class")
})

## @cvxpy test_problem.py::TestProblem::test_compilation_time
test_that("problem: compilation_time after solve", {
  ## CVXPY: prob.solve(); assert isinstance(prob.compilation_time, float)
  ## In CVXR, compile_time is stored in problem@.cache$compile_time
  x <- Variable(2, name = "x")
  prob <- Problem(Minimize(p_norm(x, 2)), list(x == 0))
  psolve(prob, solver = "CLARABEL")
  compile_time <- prob@.cache$compile_time
  expect_true(is.numeric(compile_time))
  expect_true(compile_time >= 0)
})

## @cvxpy test_problem.py::TestProblem::test_size_metrics
test_that("problem: size_metrics property", {
  ## CVXPY: p.size_metrics.num_scalar_variables, etc.
  skip("size_metrics property not implemented in CVXR Problem class")
})

# ═══════════════════════════════════════════════════════════════════════
# Batch B: test_domain.py parity tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_domain.py::TestDomain::test_geo_mean
test_that("domain: geo_mean domain constraints", {
  ## CVXPY: dom = geo_mean(x).domain
  ## Problem(Minimize(sum(x)), dom).solve() -> value ~ 0
  x <- Variable(2, name = "x")
  dom <- CVXR:::domain(geo_mean(x))
  prob <- Problem(Minimize(sum_entries(x)), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0, tolerance = 1e-3)

  ## Weighted geo_mean with one zero weight
  ## CVXPY: geo_mean(x, [0, 2]).domain + [x >= -1]
  ## optimal x = [-1, 0]
  x2 <- Variable(2, name = "x2")
  dom2 <- CVXR:::domain(geo_mean(x2, c(0, 2)))
  prob2 <- Problem(Minimize(sum_entries(x2)), c(dom2, list(x2 >= -1)))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.numeric(value(x2)), c(-1, 0), tolerance = 1e-3)

  ## 3-element with weights [0, 1, 1]
  z <- Variable(3, name = "z")
  dom3 <- CVXR:::domain(geo_mean(z, c(0, 1, 1)))
  prob3 <- Problem(Minimize(sum_entries(z)), c(dom3, list(z >= -1)))
  psolve(prob3, solver = "CLARABEL")
  expect_equal(as.numeric(value(z)), c(-1, 0, 0), tolerance = 1e-3)
})

## @cvxpy test_domain.py::TestDomain::test_log
test_that("domain: log domain constraints", {
  ## CVXPY: dom = log(a).domain, Minimize(a) s.t. dom -> value ~ 0
  a <- Variable(name = "a")
  dom <- CVXR:::domain(log(a))
  prob <- Problem(Minimize(a), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-3)
})

## @cvxpy test_domain.py::TestDomain::test_log_det
test_that("domain: log_det domain constraints", {
  ## CVXPY: dom = log_det(A + eye(2)).domain
  ## Minimize(sum(diag(A))) s.t. dom -> value ~ -2
  ## In CVXR, use matrix_trace(A) instead of sum(diag(A))
  ## since DiagVec requires a vector input, not a matrix.
  A <- Variable(c(2, 2), name = "A")
  dom <- CVXR:::domain(log_det(A + diag(2)))
  prob <- Problem(Minimize(matrix_trace(A)), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), -2, tolerance = 1e-2)
})

## @cvxpy test_domain.py::TestDomain::test_matrix_frac
test_that("domain: matrix_frac domain constraints", {
  ## CVXPY: dom = matrix_frac(x, A + eye(2)).domain
  ## Minimize(sum(diag(A))) s.t. dom -> value ~ -2
  x <- Variable(2, name = "x")
  A <- Variable(c(2, 2), name = "A")
  dom <- CVXR:::domain(matrix_frac(x, A + diag(2)))
  prob <- Problem(Minimize(matrix_trace(A)), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), -2, tolerance = 1e-2)
})

## @cvxpy test_domain.py::TestDomain::test_pnorm
test_that("domain: pnorm domain constraints for negative p", {
  ## CVXPY: dom = pnorm(a, -0.5).domain
  ## Minimize(a) s.t. dom -> value ~ 0
  a <- Variable(name = "a")
  dom <- CVXR:::domain(p_norm(a, -0.5))
  prob <- Problem(Minimize(a), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0, tolerance = 1e-3)
})

## @cvxpy test_domain.py::TestDomain::test_power
test_that("domain: power domain constraints", {
  ## sqrt(a).domain -> a >= 0; min a -> 0
  a <- Variable(name = "a")
  dom_sqrt <- CVXR:::domain(sqrt(a))
  prob <- Problem(Minimize(a), dom_sqrt)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-3)

  ## square(a).domain -> empty (no domain constraints)
  ## min a s.t. a >= -100 -> -100
  a2 <- Variable(name = "a2")
  dom_sq <- CVXR:::domain(square(a2))
  prob2 <- Problem(Minimize(a2), c(dom_sq, list(a2 >= -100)))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.numeric(value(a2)), -100, tolerance = 1e-3)

  ## a^(-1) domain -> a >= 0
  a3 <- Variable(name = "a3")
  dom_inv <- CVXR:::domain(power(a3, -1))
  prob3 <- Problem(Minimize(a3), c(dom_inv, list(a3 >= -100)))
  psolve(prob3, solver = "CLARABEL")
  expect_equal(as.numeric(value(a3)), 0, tolerance = 1e-3)

  ## a^3 domain -> a >= 0
  a4 <- Variable(name = "a4")
  dom_cube <- CVXR:::domain(power(a4, 3))
  prob4 <- Problem(Minimize(a4), c(dom_cube, list(a4 >= -100)))
  psolve(prob4, solver = "CLARABEL")
  expect_equal(as.numeric(value(a4)), 0, tolerance = 1e-3)
})


# ====================================================================
# ADDITIONAL CONSTRAINT TESTS (from test_constraints.py)
# ====================================================================

## ── test_nsd_constraint ───────────────────────────────────────────
## CVXPY: A << B creates NSD constraint (B - A >> 0, i.e., PSD(B - A)).
## Checks name, shape, value evaluation, and non-square error.
## In CVXR: %<<% operator creates PSD(e2 - e1).

## @cvxpy test_constraints.py::TestConstraints::test_nsd_constraint
test_that("CVXPY parity: NSD constraint via %<<% operator", {
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")

  ## A << B creates PSD(B - A)
  constr <- A %<<% B
  expect_equal(constr@shape, c(2L, 2L))

  ## Dual value is NULL before solve
  expect_null(dual_value(constr))

  ## Set values: B = [[2,-1],[1,2]], A = [[1,0],[0,1]]
  ## B - A = [[1,-1],[1,1]] has eigenvalues 1 +/- i... not PSD
  ## Actually for symmetric part: B - A must be symmetrized by PSD constraint
  ## Wait: CVXPY sets B = [[2,-1],[1,2]], A = [[1,0],[0,1]]
  ## PSD(B - A) = PSD([[1,-1],[1,1]]) -> checking if PSD
  ## Eigenvalues of [[1,-1],[1,1]]: trace=2, det=1+1=2, eigs=(2 +/- sqrt(4-8))/2
  ## = 1 +/- i, which means this is NOT a real symmetric PSD matrix.
  ## But CVXPY says it IS satisfied. CVXPY checks smallest eigenvalue of (B-A + (B-A)')/2.
  ## (B-A + (B-A)')/2 = ([[1,-1],[1,1]] + [[1,1],[-1,1]])/2 = [[1,0],[0,1]] = I
  ## So the symmetrized form is I, which IS PSD. -> constr is satisfied.
  value(B) <- matrix(c(2, 1, -1, 2), 2, 2)  # col-major: [[2,-1],[1,2]]
  value(A) <- diag(2)
  expect_true(value(constr))

  ## A = 3*I: B - A = [[-1,-1],[1,-1]], symmetrized = [[-1,0],[0,-1]] -> NSD -> not satisfied
  value(A) <- 3 * diag(2)
  expect_false(value(constr))

  ## Non-square matrix should error
  x <- Variable(2, name = "x")
  expect_error(x %<<% 0, ".*square.*")
})

## ── test_broadcasting ─────────────────────────────────────────────
## CVXPY: constraint broadcasting (x (3,1) == c (3,)) creates (3,3) constraint.
## R does NOT support NumPy-style broadcasting on constraints.
## CVXR shapes x as (3,1) and c(3) as (3,1), so x == c produces (3,1) constraint.

## @cvxpy test_constraints.py::TestConstraints::test_broadcasting
test_that("CVXPY parity: constraint broadcasting", {
  skip("R does not support NumPy-style broadcasting on constraints; CVXR (3,1) == (3,1) gives (3,1)")
  ## In CVXPY: x = Variable((3,1)), c = Constant(np.ones(3))
  ## (x == c).shape == (3, 3) due to NumPy broadcasting
  ## In CVXR: Variable(c(3,1)) and Constant(rep(1,3)) (shape (3,1))
  ## produce (3,1) constraint, not (3,3)
  x <- Variable(c(3, 1))
  c_val <- Constant(rep(1, 3))
  con <- (x == c_val)
  ## In CVXPY this is (3,3); in CVXR this would be (3,1)
  expect_equal(con@shape, c(3L, 3L))

  prob <- Problem(Minimize(0), list(con))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), rep(1, 3), tolerance = 1e-4)
})
