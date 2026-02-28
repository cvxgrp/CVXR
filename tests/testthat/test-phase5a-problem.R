## Tests for Phase 5a: Problem class, FlipObjective, Dcp2Cone
## ============================================================

# ── Problem class ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("Problem constructs with Minimize objective", {
  x <- Variable(2)
  obj <- Minimize(sum(x))
  p <- Problem(obj, list(x >= 0))
  expect_true(S7_inherits(p, Problem))
  expect_true(S7_inherits(p@objective, Minimize))
  expect_equal(length(p@constraints), 1L)
})

## @cvxpy NONE
test_that("Problem constructs with Maximize objective", {
  x <- Variable(2)
  obj <- Maximize(sum(x))
  p <- Problem(obj, list(x <= 1))
  expect_true(S7_inherits(p@objective, Maximize))
})

## @cvxpy NONE
test_that("Problem rejects non-Objective", {
  x <- Variable(2)
  expect_error(Problem(x), "Minimize.*Maximize")
})

## @cvxpy NONE
test_that("Problem rejects non-Constraint in constraints", {
  x <- Variable(2)
  obj <- Minimize(sum(x))
  expect_error(Problem(obj, list(x)), "Constraint")
})

## @cvxpy NONE
test_that("Problem constructs with empty constraints", {
  x <- Variable(1)
  p <- Problem(Minimize(x))
  expect_equal(length(p@constraints), 0L)
})

# ── is_dcp ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("is_dcp: convex minimize is DCP", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 0))
  expect_true(is_dcp(p))
})

## @cvxpy NONE
test_that("is_dcp: concave maximize is DCP", {
  x <- Variable(2)
  p <- Problem(Maximize(sum(x)), list(x <= 1))
  expect_true(is_dcp(p))
})

## @cvxpy NONE
test_that("is_dcp: non-DCP — minimize concave", {
  x <- Variable(2, nonneg = TRUE)
  ## -sum(x^2) is concave, so Minimize(-sum(x^2)) is NOT DCP
  ## Actually sum(x^2) is convex (via Power), not concave.
  ## Let's use a genuinely non-DCP case: maximize convex
  ## Maximize(sum_squares) needs QuadForm. Use simpler: Maximize(exp(x))
  ## Actually Maximize(exp(x)) on a scalar is non-DCP since exp is convex, not concave
  x1 <- Variable(1)
  p <- Problem(Maximize(exp(x1)))
  expect_false(is_dcp(p))
})

# ── variables / parameters / constants ─────────────────────────────

## @cvxpy NONE
test_that("variables collects from objective and constraints", {
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  p <- Problem(Minimize(sum(x)), list(y >= 0))
  vars_ <- variables(p)
  ids <- vapply(vars_, function(v) v@id, integer(1))
  expect_true(x@id %in% ids)
  expect_true(y@id %in% ids)
  expect_equal(length(vars_), 2L)
})

## @cvxpy NONE
test_that("parameters returns empty if none", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 0))
  expect_equal(length(parameters(p)), 0L)
})

## @cvxpy NONE
test_that("constants returns from objective + constraints", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 1))
  consts <- constants(p)
  ## There should be at least the constant 1
  expect_true(length(consts) > 0L)
})

# ── value / status ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("Problem value is NULL before solve", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)))
  expect_null(value(p))
})

## @cvxpy NONE
test_that("problem_status is NULL before solve", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)))
  expect_null(status(p))
})

# ── is_qp / is_lp ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("is_lp for linear problem", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 0, x <= 1))
  expect_true(is_lp(p))
  expect_true(is_qp(p))  ## LP is a special case of QP
})

## @cvxpy NONE
test_that("non-DCP problem is not QP", {
  x <- Variable(1)
  p <- Problem(Maximize(exp(x)))
  expect_false(is_qp(p))
})

# ── print ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Problem prints without error", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 0))
  expect_output(print(p), "Problem")
})

# ═══════════════════════════════════════════════════════════════════
# FlipObjective
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("FlipObjective flips Maximize to Minimize", {
  x <- Variable(1)
  p <- Problem(Maximize(x), list(x <= 5))
  fo <- FlipObjective()
  expect_true(reduction_accepts(fo, p))
  result <- reduction_apply(fo, p)
  new_p <- result[[1L]]
  expect_true(S7_inherits(new_p@objective, Minimize))
})

## @cvxpy NONE
test_that("FlipObjective flips Minimize to Maximize", {
  x <- Variable(1)
  p <- Problem(Minimize(x), list(x >= 0))
  fo <- FlipObjective()
  result <- reduction_apply(fo, p)
  new_p <- result[[1L]]
  expect_true(S7_inherits(new_p@objective, Maximize))
})

## @cvxpy NONE
test_that("FlipObjective preserves constraints", {
  x <- Variable(2)
  p <- Problem(Maximize(sum(x)), list(x >= 0, x <= 1))
  fo <- FlipObjective()
  result <- reduction_apply(fo, p)
  expect_equal(length(result[[1L]]@constraints), 2L)
})

## @cvxpy NONE
test_that("FlipObjective invert negates opt_val", {
  fo <- FlipObjective()
  sol <- Solution(status = "optimal", opt_val = 3.5)
  result <- reduction_invert(fo, sol, list())
  expect_equal(result@opt_val, -3.5)
})

# ═══════════════════════════════════════════════════════════════════
# Reduction / Canonicalization infrastructure
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Reduction base class errors on abstract methods", {
  r <- Reduction()
  p <- Problem(Minimize(Variable(1)))
  expect_error(reduction_accepts(r, p), "must implement")
  expect_error(reduction_apply(r, p), "must implement")
})

## @cvxpy NONE
test_that("InverseData stores variable mappings", {
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  p <- Problem(Minimize(sum(x)), list(y >= 0))
  inv <- InverseData(p)
  ## Should have 2 variables mapped
  expect_equal(length(inv@id_map), 2L)
  expect_equal(inv@x_length, 5L)  ## 2 + 3
})

## @cvxpy NONE
test_that("DCP canonicalization S7 generics are registered", {
  ## dcp_canonicalize is an S7 generic
  expect_true(inherits(dcp_canonicalize, "S7_generic"))

  ## has_dcp_canon returns TRUE for representative atoms that have canon methods
  x <- Variable(2)
  expect_true(has_dcp_canon(Abs(Variable())))
  expect_true(has_dcp_canon(Exp(Variable())))
  expect_true(has_dcp_canon(Log(Variable())))
  expect_true(has_dcp_canon(Entr(Variable())))
  expect_true(has_dcp_canon(Maximum(Variable(), Variable())))
  expect_true(has_dcp_canon(MaxEntries(Variable())))
  expect_true(has_dcp_canon(Power(Variable(), 2)))
  expect_true(has_dcp_canon(Pnorm(x)))
  expect_true(has_dcp_canon(Norm1(x)))
  expect_true(has_dcp_canon(NormInf(x)))
  expect_true(has_dcp_canon(SumLargest(x, 1)))
  expect_true(has_dcp_canon(QuadOverLin(Variable(), Variable())))

  ## has_dcp_canon returns FALSE for non-canon atoms (affine atoms, Variable)
  expect_false(has_dcp_canon(Variable()))
  expect_false(has_dcp_canon(AddExpression(list(Variable(), Variable()))))
})
