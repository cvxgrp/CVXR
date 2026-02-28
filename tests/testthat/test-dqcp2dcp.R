## Tests for Dqcp2Dcp reduction + bisection solver
## CVXPY SOURCE: tests/test_dqcp2dcp.py (selected parity tests)

## ── Bisection: ceil ──────────────────────────────────────────────

## @cvxpy NONE
test_that("DQCP bisection: minimize ceil(x) with x >= 0.5", {
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(ceiling(x)), list(x >= 0.5))

  ## Verify DQCP properties
  expect_true(is_dqcp(prob))
  expect_false(is_dcp(prob))

  ## Solve via bisection
  result <- psolve(prob, qcp = TRUE)
  expect_equal(as.numeric(result), 1.0, tolerance = 1e-3)
})

## ── Bisection: floor ─────────────────────────────────────────────

## @cvxpy NONE
test_that("DQCP bisection: minimize floor(x) with x >= 1.5", {
  x <- Variable()
  prob <- Problem(Minimize(floor(x)), list(x >= 1.5, x <= 5))

  expect_true(is_dqcp(prob))
  expect_false(is_dcp(prob))

  result <- psolve(prob, qcp = TRUE)
  expect_equal(as.numeric(result), 1.0, tolerance = 1e-3)
})

## ── Bisection: ratio (sublevel set) ─────────────────────────────

## @cvxpy NONE
test_that("DQCP bisection: minimize x/y ratio", {
  x <- Variable(nonneg = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x / y), list(x >= 1, y <= 2, x + y <= 5))

  expect_true(is_dqcp(prob))

  result <- psolve(prob, qcp = TRUE)
  expect_equal(as.numeric(result), 0.5, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 2.0, tolerance = 1e-2)
})

## ── Bisection: Maximize (product) ───────────────────────────────

## @cvxpy NONE
test_that("DQCP bisection: maximize x*y with x+y<=4", {
  x <- Variable(nonneg = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Maximize(x * y), list(x + y <= 4, x >= 0))

  expect_true(is_dqcp(prob))

  result <- psolve(prob, qcp = TRUE)
  expect_equal(as.numeric(result), 4.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 0.1)
  expect_equal(as.numeric(value(y)), 2.0, tolerance = 0.1)
})

## ── Bisection: infeasible ───────────────────────────────────────

## @cvxpy NONE
test_that("DQCP bisection: infeasible problem", {
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(ceiling(x)), list(x >= 5, x <= 2))

  expect_true(is_dqcp(prob))

  result <- psolve(prob, qcp = TRUE)
  expect_equal(as.numeric(result), Inf)
  expect_equal(status(prob), "infeasible")
})

## ── Error: gp and qcp both TRUE ─────────────────────────────────

## @cvxpy NONE
test_that("psolve errors when both gp and qcp are TRUE", {
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(x), list(x >= 1))
  expect_error(psolve(prob, gp = TRUE, qcp = TRUE), "gp.*qcp")
})

## ── Error: non-DQCP problem ─────────────────────────────────────

## @cvxpy NONE
test_that("psolve errors for non-DQCP with qcp=TRUE", {
  ## A problem that is neither DCP nor DQCP
  ## (for instance, maximizing a convex function without quasiconcave structure)
  x <- Variable()
  ## x^2 is convex, not quasiconcave, so Maximize(x^2) is not DQCP
  prob <- Problem(Maximize(square(x)), list(x >= -1, x <= 1))
  expect_false(is_dqcp(prob))
  expect_error(psolve(prob, qcp = TRUE), "not DQCP")
})

## ── DCP problem with qcp=TRUE goes through standard path ────────

## @cvxpy NONE
test_that("psolve with qcp=TRUE on DCP problem uses standard solver", {
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))

  ## DCP problem: qcp=TRUE should still work (standard path)
  result <- psolve(prob, qcp = TRUE)
  expect_equal(as.numeric(result), 5.0, tolerance = 1e-6)
})

## ── Inverse functions unit tests ────────────────────────────────

## @cvxpy NONE
test_that("DQCP invertibility checks", {
  x <- Variable(nonneg = TRUE)
  ## Ceil, Floor, Exp, Log are always invertible
  expect_true(.dqcp_invertible(Ceil(x)))
  expect_true(.dqcp_invertible(Floor(x)))
  expect_true(.dqcp_invertible(exp(x)))
  expect_true(.dqcp_invertible(log(x)))
  expect_true(.dqcp_invertible(log1p_atom(x)))
  expect_true(.dqcp_invertible(logistic(x)))
  expect_true(.dqcp_invertible(power(x, 2)))
  expect_true(.dqcp_invertible(abs(x)))

  ## Multiply with one constant arg is invertible
  c1 <- Constant(2)
  expect_true(.dqcp_invertible(c1 * x))
  ## Multiply with two variables is not
  y <- Variable()
  expect_false(.dqcp_invertible(x * y))
})

## ── Tighten functions ───────────────────────────────────────────

## @cvxpy NONE
test_that("DQCP tighten functions", {
  x <- Variable(nonneg = TRUE)

  ## Ceil → ceiling/floor tightening
  fns <- .tighten_fns(Ceil(x))
  expect_equal(fns$lower(1.3), 2)
  expect_equal(fns$upper(1.7), 1)

  ## Floor → ceiling/floor tightening
  fns <- .tighten_fns(Floor(x))
  expect_equal(fns$lower(1.3), 2)
  expect_equal(fns$upper(1.7), 1)

  ## Nonneg expression
  fns <- .tighten_fns(x)
  expect_equal(fns$lower(-1), 0)
  expect_equal(fns$lower(0.5), 0.5)
  expect_equal(fns$upper(1), 1)
})

## ── Sublevel/superlevel sets ────────────────────────────────────

## @cvxpy NONE
test_that("DQCP sublevel/superlevel set registries", {
  ## Multiply sublevel: nonneg * nonpos
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonpos = TRUE)
  expr <- x * y
  sub <- .dqcp_sublevel(expr, t = Constant(0))
  expect_true(length(sub) >= 1)

  ## Multiply superlevel: nonneg * nonneg
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  expr <- x * y
  sup <- .dqcp_superlevel(expr, t = Constant(0))
  expect_true(length(sup) >= 1)

  ## Ratio sublevel
  x <- Variable()
  y <- Variable(nonneg = TRUE)
  expr <- x / y
  sub <- .dqcp_sublevel(expr, t = Constant(1))
  expect_true(length(sub) >= 1)

  ## Ratio superlevel
  sup <- .dqcp_superlevel(expr, t = Constant(0))
  expect_true(length(sup) >= 1)
})

## ── Verbose output ──────────────────────────────────────────────

## @cvxpy NONE
test_that("DQCP bisection verbose output", {
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(ceiling(x)), list(x >= 0.5))
  expect_message(
    result <- psolve(prob, qcp = TRUE, verbose = TRUE),
    "bisection"
  )
  expect_equal(as.numeric(result), 1.0, tolerance = 1e-3)
})

## ── Bisection with explicit bounds ──────────────────────────────

## @cvxpy NONE
test_that("DQCP bisection with user-provided bounds", {
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(ceiling(x)), list(x >= 0.5))

  ## Provide explicit low/high to skip interval search
  result <- psolve(prob, qcp = TRUE, low = 0, high = 10)
  expect_equal(as.numeric(result), 1.0, tolerance = 1e-3)
})

## ── Dqcp2Dcp accepts ───────────────────────────────────────────

## @cvxpy NONE
test_that("Dqcp2Dcp accepts Minimize DQCP problems", {
  reducer <- Dqcp2Dcp()

  ## Accepts Minimize DQCP
  x <- Variable(nonneg = TRUE)
  prob_min <- Problem(Minimize(ceiling(x)), list(x >= 0.5))
  expect_true(reduction_accepts(reducer, prob_min))

  ## Rejects Maximize (FlipObjective handles that)
  prob_max <- Problem(Maximize(floor(x)), list(x <= 5))
  expect_false(reduction_accepts(reducer, prob_max))
})

## ── Exp/Log inversions ─────────────────────────────────────────

## @cvxpy NONE
test_that("DQCP bisection: minimize exp(x), x >= 1", {
  x <- Variable()
  prob <- Problem(Minimize(exp(x)), list(x >= 1))

  ## exp is quasiconvex (convex), so this is DCP
  ## But testing with qcp=TRUE should still work
  expect_true(is_dqcp(prob))
  result <- psolve(prob, qcp = TRUE)
  expect_equal(as.numeric(result), exp(1), tolerance = 1e-2)
})
