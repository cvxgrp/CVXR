# Tests for the interval-arithmetic bounds module (CVXPY 1.9.0, #3080).
#
# These exercise utilities/bounds.R directly. Every expected value below was
# computed from CVXPY 1.9.0's cvxpy/utilities/bounds.py on the same inputs, so
# the R port is a verified mirror.
#
# CVXPY SOURCE: utilities/bounds.py

library(CVXR)

bp_lb <- function(b) as.numeric(b[[1L]])
bp_ub <- function(b) as.numeric(b[[2L]])

## @cvxpy NONE (verified against bounds.py)
test_that("constructors: unbounded / uniform / scalar", {
  u <- unbounded(c(2L, 1L))
  expect_equal(bp_lb(u), c(-Inf, -Inf))
  expect_equal(bp_ub(u), c(Inf, Inf))
  uf <- uniform_bounds(c(2L, 2L), -1, 3)
  expect_equal(bp_lb(uf), rep(-1, 4))
  expect_equal(bp_ub(uf), rep(3, 4))
  s <- scalar_bounds(2, 5)
  expect_equal(bp_lb(s), 2)
  expect_equal(bp_ub(s), 5)
})

## @cvxpy NONE
test_that("elementwise binary: add / neg / mul / div / scale", {
  expect_equal(bp_lb(add_bounds(c(-1, 0), c(2, 3), c(-3, 1), c(4, 5))), c(-4, 1))
  expect_equal(bp_ub(add_bounds(c(-1, 0), c(2, 3), c(-3, 1), c(4, 5))), c(6, 8))
  expect_equal(bp_lb(neg_bounds(c(-1, 2), c(3, 5))), c(-3, -5))

  m <- mul_bounds(c(-1, 0), c(2, 3), c(-3, 1), c(4, 5))
  expect_equal(bp_lb(m), c(-6, 0))
  expect_equal(bp_ub(m), c(8, 15))

  expect_equal(bp_lb(div_bounds(1, 2, -1, 4)), -Inf)   # divisor spans 0
  expect_equal(bp_lb(div_bounds(1, 2, 2, 4)), 0.25)
  expect_equal(bp_ub(div_bounds(1, 2, 2, 4)), 1)

  expect_equal(bp_lb(scale_bounds(c(1, 2), c(3, 4), -2)), c(-6, -8))
})

## @cvxpy NONE
test_that("abs: spanning / positive / negative intervals", {
  a <- abs_bounds(c(-2, 1, -5), c(3, 4, -2))
  expect_equal(bp_lb(a), c(0, 1, 2))
  expect_equal(bp_ub(a), c(3, 4, 5))
})

## @cvxpy NONE
test_that("power: even / odd / fractional / negative", {
  p2 <- power_bounds(c(-2, 1, -5), c(3, 4, -2), 2)
  expect_equal(bp_lb(p2), c(0, 1, 4))
  expect_equal(bp_ub(p2), c(9, 16, 25))
  p3 <- power_bounds(c(-2, 1), c(3, 4), 3)
  expect_equal(bp_lb(p3), c(-8, 1)); expect_equal(bp_ub(p3), c(27, 64))
  ph <- power_bounds(c(-1, 4), c(1, 9), 0.5)
  expect_equal(bp_lb(ph), c(-Inf, 2)); expect_equal(bp_ub(ph), c(Inf, 3))
  pn <- power_bounds(c(1, -2), c(2, 3), -1)
  expect_equal(bp_lb(pn), c(0.5, -Inf)); expect_equal(bp_ub(pn), c(1, Inf))
})

## @cvxpy NONE
test_that("monotone elementwise: exp / log / sqrt", {
  expect_equal(bp_ub(exp_bounds(c(0, 1), c(1, 2))), c(exp(1), exp(2)))
  lg <- log_bounds(c(-1, 1), c(0.5, exp(1)))
  expect_equal(bp_lb(lg), c(-Inf, 0)); expect_equal(bp_ub(lg), c(log(0.5), 1))
  sq <- sqrt_bounds(c(-1, 4), c(1, 9))
  expect_equal(bp_lb(sq), c(-Inf, 2)); expect_equal(bp_ub(sq), c(1, 3))
})

## @cvxpy NONE
test_that("axis reductions: sum (all/0/1) and max", {
  lb <- matrix(c(-1, 0, 2, 3), 2, 2); ub <- matrix(c(1, 2, 4, 5), 2, 2)
  expect_equal(bp_lb(sum_bounds(lb, ub, axis = NULL)), 4)
  expect_equal(bp_ub(sum_bounds(lb, ub, axis = NULL)), 12)
  expect_equal(bp_lb(sum_bounds(lb, ub, axis = 0)), c(-1, 5))   # per column
  expect_equal(bp_lb(sum_bounds(lb, ub, axis = 1)), c(1, 3))    # per row
  expect_equal(bp_lb(max_reduction_bounds(lb, ub, axis = 0)), c(0, 3))
})

## @cvxpy NONE
test_that("elementwise maximum / minimum over operands", {
  bl <- list(list(c(-1, 2), c(3, 4)), list(c(0, -5), c(2, 1)))
  expect_equal(bp_lb(maximum_bounds(bl)), c(0, 2))
  expect_equal(bp_ub(maximum_bounds(bl)), c(3, 4))
  expect_equal(bp_lb(minimum_bounds(bl)), c(-1, -5))
})

## @cvxpy NONE
test_that("matmul: constant (point) operand uses the exact split formula", {
  A <- matrix(c(1, 3, -2, 4), 2, 2)   # [[1,-2],[3,4]]
  mm <- matmul_bounds(A, A, matrix(c(-1, 0), 2, 1), matrix(c(1, 2), 2, 1))
  expect_equal(bp_lb(mm), c(-5, -3))
  expect_equal(bp_ub(mm), c(1, 11))
  ## both intervals -> unbounded
  mi <- matmul_bounds(matrix(c(0, 0, 0, 0), 2, 2), matrix(1, 2, 2),
                      matrix(c(-1, 0), 2, 1), matrix(c(1, 2), 2, 1))
  expect_true(all(is.infinite(bp_lb(mi))))
})

## @cvxpy NONE
test_that("norms: 1-norm and inf-norm; sign refinement", {
  expect_equal(bp_lb(norm1_bounds(c(-2, 1, -5), c(3, 4, -2))), 3)
  expect_equal(bp_ub(norm1_bounds(c(-2, 1, -5), c(3, 4, -2))), 12)
  expect_equal(bp_lb(norm_inf_bounds(c(-2, 1, -5), c(3, 4, -2))), 2)
  expect_equal(bp_ub(norm_inf_bounds(c(-2, 1, -5), c(3, 4, -2))), 5)
  r <- refine_bounds_from_sign(c(-2, 1), c(3, 4), TRUE, FALSE)
  expect_equal(bp_lb(r), c(0, 1))
})

## @cvxpy utilities/bounds.py::get_expr_bounds_if_supported
test_that("get_expr_bounds_if_supported filters useful auxiliary bounds", {
  x <- Variable(2, bounds = list(-2, 3))
  expr <- abs(x)

  no_bounds_solver <- SolverInfo(solver_name = "SCS", solver_supports_bounds = FALSE)
  expect_null(get_expr_bounds_if_supported(expr, no_bounds_solver))
  expect_null(get_expr_bounds_if_supported(expr, NULL))

  bounds_solver <- SolverInfo(solver_name = "HIGHS", solver_supports_bounds = TRUE)
  b <- get_expr_bounds_if_supported(expr, bounds_solver)
  expect_equal(as.numeric(b[[1L]]), c(0, 0))
  expect_equal(as.numeric(b[[2L]]), c(3, 3))

  y <- Variable(2, nonneg = TRUE)
  expect_null(get_expr_bounds_if_supported(y, bounds_solver))
})

## @cvxpy reductions/eliminate_pwl/canonicalizers/abs_canon.py
test_that("DCP auxiliary variables receive bounds when solver supports them", {
  x <- Variable(2, bounds = list(-2, 3))
  expr <- abs(x)
  ctx <- SolverInfo(solver_name = "HIGHS", solver_supports_bounds = TRUE)
  result <- abs_canon(expr, expr@args, solver_context = ctx)
  t <- result[[1L]]
  b <- get_bounds(t)
  expect_equal(as.numeric(b[[1L]]), c(0, 0))
  expect_equal(as.numeric(b[[2L]]), c(3, 3))
})

## @cvxpy expressions/variable.py::Variable.parameters
test_that("Variables report Parameters embedded in expression bounds", {
  lb <- Parameter(2, value = c(-1, -2))
  ub <- Parameter(2, value = c(3, 4))
  x <- Variable(2, bounds = list(lb, ub))
  expect_equal(vapply(parameters(x), function(p) p@id, integer(1)),
               c(lb@id, ub@id))
})

## @cvxpy expressions/variable.py::Variable.is_dpp
test_that("Variables validate expression bounds in DPP scope", {
  p <- Parameter(2, value = c(1, 2), nonneg = TRUE)
  x <- Variable(2, bounds = list(-p, p))
  expect_true(is_dpp(x))

  bad_bound <- square(p)
  y <- Variable(2, bounds = list(-bad_bound, bad_bound))
  expect_false(is_dpp(y))
})

## @cvxpy NONE
test_that("structural: reshape (F-order) / transpose / index", {
  lb <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3); ub <- lb + 10
  rs <- reshape_bounds(lb, ub, c(3L, 2L), order = "F")
  expect_equal(as.numeric(rs[[1L]]), 1:6)         # column-major preserved
  rc <- reshape_bounds(lb, ub, c(3L, 2L), order = "C")
  expect_equal(as.numeric(rc[[1L]]), c(1, 5, 4, 3, 2, 6))   # row-major (vs numpy 'C')
  tb <- transpose_bounds(lb, ub)
  expect_equal(dim(tb[[1L]]), c(3L, 2L))
  ib <- index_bounds(lb, ub, function(a) a[, 1L, drop = FALSE])
  expect_equal(as.numeric(ib[[1L]]), c(1, 2))
})
