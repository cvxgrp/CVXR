## CVXPY parity tests for DQCP, perspective, and FiniteSet
## Phase 6: Comprehensive end-to-end parity verification
## CVXPY SOURCE: tests/test_dqcp.py, tests/test_perspective.py, tests/test_valinvec2mixedint.py
##
## All expected values verified via `uv run python` against CVXPY 1.8.1

## ── DQCP parity: ceil/floor ──────────────────────────────────────

## @cvxpy test_dqcp.py::TestDqcp::test_basic_with_interval
test_that("DQCP parity: minimize ceil(x), x in [12, 17]", {
  x <- Variable()
  prob <- Problem(Minimize(ceiling(x)), list(x >= 12, x <= 17))
  result <- psolve(prob, qcp = TRUE, low = 12, high = 17)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 12.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 12.0, tolerance = 0.1)
})

## @cvxpy test_dqcp.py::TestDqcp::test_basic_maximization_with_interval
test_that("DQCP parity: maximize ceil(x), x in [12, 17]", {
  x <- Variable()
  prob <- Problem(Maximize(ceiling(x)), list(x >= 12, x <= 17))
  result <- psolve(prob, qcp = TRUE)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 17.0, tolerance = 1e-3)
})

## @cvxpy test_dqcp.py::TestDqcp::test_basic_floor
test_that("DQCP parity: minimize floor(x), x in [11.8, 17]", {
  x <- Variable()
  prob <- Problem(Minimize(floor(x)), list(x >= 11.8, x <= 17))
  result <- psolve(prob, qcp = TRUE)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 11.0, tolerance = 1e-3)
  expect_true(as.numeric(value(x)) > 11.7)
})

## @cvxpy test_dqcp.py::TestDqcp::test_multiply_const
test_that("DQCP parity: 0.5 * ceil(x), x >= 10", {
  x <- Variable()
  prob <- Problem(Minimize(0.5 * ceiling(x)), list(x >= 10))
  result <- psolve(prob, qcp = TRUE)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 5.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 10.0, tolerance = 0.1)
})

## @cvxpy test_dqcp.py::TestDqcp::test_basic_maximum
test_that("DQCP parity: max(ceil(x), ceil(y))", {
  x <- Variable()
  y <- Variable()
  prob <- Problem(
    Minimize(max_elemwise(ceiling(x), ceiling(y))),
    list(x >= 12, x <= 17, y >= 17.4)
  )
  result <- psolve(prob, qcp = TRUE)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 18.0, tolerance = 1e-3)
  expect_true(as.numeric(value(y)) > 17.3)
})

## ── DQCP parity: products ────────────────────────────────────────

## @cvxpy test_dqcp.py::TestDqcp::test_basic_multiply_nonneg
test_that("DQCP parity: maximize x*y, nonneg", {
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  prob <- Problem(Maximize(x * y), list(x <= 12, y <= 6))
  result <- psolve(prob, qcp = TRUE)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 72.0, tolerance = 0.1)
  expect_equal(as.numeric(value(x)), 12.0, tolerance = 0.1)
  expect_equal(as.numeric(value(y)), 6.0, tolerance = 0.1)
})

## @cvxpy test_dqcp.py::TestDqcp::test_basic_multiply_nonpos
test_that("DQCP parity: maximize x*y, nonpos", {
  x <- Variable(nonpos = TRUE)
  y <- Variable(nonpos = TRUE)
  prob <- Problem(Maximize(x * y), list(x >= -12, y >= -6))
  result <- psolve(prob, qcp = TRUE)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 72.0, tolerance = 0.1)
  expect_equal(as.numeric(value(x)), -12.0, tolerance = 0.1)
  expect_equal(as.numeric(value(y)), -6.0, tolerance = 0.1)
})

## @cvxpy test_dqcp.py::TestDqcp::test_basic_multiply_qcvx
test_that("DQCP parity: minimize x*y, nonneg*nonpos (quasiconvex)", {
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonpos = TRUE)
  prob <- Problem(Minimize(x * y), list(x <= 7, y >= -6))
  result <- psolve(prob, qcp = TRUE)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), -42.0, tolerance = 0.1)
  expect_equal(as.numeric(value(x)), 7.0, tolerance = 0.1)
  expect_equal(as.numeric(value(y)), -6.0, tolerance = 0.1)
})

## @cvxpy test_dqcp.py::TestDqcp::test_concave_multiply
test_that("DQCP parity: maximize sqrt(x)*sqrt(y), nonneg", {
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  prob <- Problem(Maximize(sqrt(x) * sqrt(y)), list(x <= 4, y <= 9))
  result <- psolve(prob, qcp = TRUE)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 6.0, tolerance = 0.1)
  expect_equal(as.numeric(value(x)), 4.0, tolerance = 0.1)
  expect_equal(as.numeric(value(y)), 9.0, tolerance = 0.1)
})

## ── DQCP parity: ratio ──────────────────────────────────────────

## @cvxpy test_dqcp.py::TestDqcp::test_basic_ratio
test_that("DQCP parity: minimize x/y ratio", {
  x <- Variable()
  y <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(x / y), list(x == 12, y <= 6, y >= 1))
  result <- psolve(prob, qcp = TRUE, low = 0, high = 20)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 12.0, tolerance = 0.1)
  expect_equal(as.numeric(value(y)), 6.0, tolerance = 0.1)
})

## ── DQCP error messages ─────────────────────────────────────────

## @cvxpy NONE
test_that("DQCP suggestion in error message for non-DCP DQCP problem", {
  x <- Variable()
  prob <- Problem(Minimize(ceiling(x)), list(x >= 0.5))

  ## Not DCP, is DQCP -- should suggest qcp = TRUE
  expect_false(is_dcp(prob))
  expect_true(is_dqcp(prob))
  expect_error(psolve(prob), "qcp.*TRUE")
})

## @cvxpy NONE
test_that("DGP suggestion preferred over DQCP for DGP problems", {
  ## A DGP problem -- should suggest gp=TRUE not qcp=TRUE
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x * y), list(x >= 1, y >= 2))

  ## DGP problems are also DQCP, but DGP suggestion should come first
  expect_true(is_dgp(prob))
  expect_error(psolve(prob), "gp.*TRUE")
})

## @cvxpy NONE
test_that("non-DCP non-DQCP gives generic error", {
  ## Maximize(convex) is not quasiconcave → not DQCP
  x <- Variable()
  prob <- Problem(Maximize(square(x)), list(x >= -1, x <= 1))
  expect_false(is_dcp(prob))
  expect_false(is_dqcp(prob))
  expect_error(psolve(prob), "not DCP")
})

## ── Perspective parity ──────────────────────────────────────────

## @cvxpy test_perspective.py::test_exp
test_that("Perspective parity: minimize s*exp(x/s)", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- exp(x)
  prob <- Problem(Minimize(perspective(f, s)), list(s >= 1, x >= 1))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## CVXPY: 2.718282 (exp(1))
  expect_equal(opt_val, exp(1), tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(s)), 1.0, tolerance = 1e-2)
})

## @cvxpy test_perspective.py::test_quad_atom
test_that("Perspective parity: minimize s*(x/s)^2, s <= 1, x >= 2", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- square(x)
  prob <- Problem(Minimize(perspective(f, s)), list(s <= 1, x >= 2))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## CVXPY: 4.0 (x=2, s=1, perspective = 4/1 = 4)
  expect_equal(opt_val, 4.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(s)), 1.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("Perspective parity: scalar affine f", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(perspective(x - 1, s)), list(x >= 3.14, s <= 1))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## CVXPY: 2.14 = 3.14 - 1
  expect_equal(opt_val, 2.14, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Perspective parity: maximize s*log(x/s) (concave)", {
  skip_if_not_installed("clarabel")
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  prob <- Problem(
    Maximize(perspective(log(x), s)),
    list(x >= 1, x <= 2, s >= 1, s <= 2)
  )
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## CVXPY: 0.693147 (= log(2)), at x=2, s=1
  expect_equal(opt_val, log(2), tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(s)), 1.0, tolerance = 1e-2)
})

## @cvxpy test_perspective.py::test_parameter
test_that("Perspective parity: with Parameter", {
  skip_if_not_installed("clarabel")
  p <- Parameter(nonneg = TRUE)
  value(p) <- 99
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- p * square(x)
  prob <- Problem(Minimize(perspective(f, s)), list(s <= 1, x >= 2))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## CVXPY: 396 = 99 * 4
  expect_equal(opt_val, 396.0, tolerance = 1e-1)
})

## @cvxpy test_perspective.py::test_evaluate_persp
test_that("Perspective parity: numeric eval (x^2 + 3x - 5)", {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- square(x) + 3 * x - 5
  p <- perspective(f, s)

  ## (x=1, s=2): s * f(x/s) = 2 * ((1/2)^2 + 3*(1/2) - 5) = 2*(0.25+1.5-5) = -6.5
  value(x) <- 1.0
  value(s) <- 2.0
  expect_equal(as.numeric(value(p)), -6.5, tolerance = 1e-6)

  ## (x=5, s=0.25): s * f(x/s) = 0.25 * (20^2 + 3*20 - 5) = 0.25*(400+60-5) = 113.75
  value(x) <- 5.0
  value(s) <- 0.25
  expect_equal(as.numeric(value(p)), 113.75, tolerance = 1e-6)
})

## ── FiniteSet parity ─────────────────────────────────────────────

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_1
test_that("FiniteSet parity: minimize sum, set={1,3,5}, GLPK_MI", {
  skip_if_not_installed("Rglpk")
  x <- Variable(4)
  prob <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(1, 3, 5))))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  ## CVXPY: 4.0 (four 1's)
  expect_equal(opt_val, 4.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1, 1, 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("FiniteSet parity: ineq_form matches eq_form", {
  skip_if_not_installed("Rglpk")
  x <- Variable(4)
  prob_eq <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(1, 3, 5))))
  val_eq <- psolve(prob_eq, solver = "GLPK_MI")

  y <- Variable(4)
  prob_ineq <- Problem(Minimize(sum_entries(y)),
    list(FiniteSet(y, c(1, 3, 5), ineq_form = TRUE)))
  val_ineq <- psolve(prob_ineq, solver = "GLPK_MI")

  expect_equal(val_eq, val_ineq, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_11
test_that("FiniteSet parity: 2D variable, set={0,1} (binary)", {
  skip_if_not_installed("Rglpk")
  x <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(0, 1))))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  ## CVXPY: 0.0 (four 0's)
  expect_equal(opt_val, 0.0, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_independent_entries
test_that("FiniteSet parity: independent entry choices", {
  skip_if_not_installed("Rglpk")
  x <- Variable(3)
  prob <- Problem(
    Minimize(x[1] + x[2] - x[3]),
    list(FiniteSet(x, c(0, 5, 10)))
  )
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  ## x1=0, x2=0, x3=10 → objective = -10
  expect_equal(opt_val, -10.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("FiniteSet parity: Gurobi solver", {
  skip_if_not(requireNamespace("gurobi", quietly = TRUE))
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(1, 3, 5))))
  opt_val <- psolve(prob, solver = "GUROBI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 3.0, tolerance = 1e-4)
})
