## Tests for Category 2 atoms: cummax, dotsort, inv_prod, loggamma, log_normcdf

# ═══════════════════════════════════════════════════════════════════
# inv_prod
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("inv_prod: shape and curvature", {
  x <- Variable(3, nonneg = TRUE)
  e <- inv_prod(x)
  expect_equal(e@shape, c(1L, 1L))
  expect_true(is_convex(e))
})

## @cvxpy NONE
test_that("inv_prod: solves correctly with bounds", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Minimize(inv_prod(x)), list(x >= 1, x <= 2))
  val <- psolve(prob)
  ## 1/(2*2*2) = 0.125
  expect_equal(val, 0.125, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# loggamma
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("loggamma: returns expression", {
  x <- Variable(nonneg = TRUE)
  e <- loggamma(x)
  expect_true(S7::S7_inherits(e, Expression))
  expect_true(is_convex(e))
})

## @cvxpy NONE
test_that("loggamma: minimize over [1,3]", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(loggamma(x)), list(x >= 1, x <= 3))
  val <- psolve(prob)
  ## CVXPY reference: -0.126533 at x ≈ 1.342
  expect_equal(val, -0.1265, tolerance = 0.01)
})

# ═══════════════════════════════════════════════════════════════════
# log_normcdf
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("log_normcdf: shape preserved", {
  x <- Variable(3)
  e <- log_normcdf(x)
  expect_equal(e@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("log_normcdf: concave", {
  x <- Variable()
  e <- log_normcdf(x)
  expect_true(is_concave(e))
})

## @cvxpy NONE
test_that("log_normcdf: maximize at x=2", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Maximize(log_normcdf(x)), list(x >= -2, x <= 2))
  val <- psolve(prob)
  ## CVXPY reference: -0.023013 (log(pnorm(2)))
  expect_equal(val, -0.023013, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# cummax
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cummax_expr: shape preserved", {
  x <- Variable(4)
  e <- cummax_expr(x)
  expect_equal(e@shape, c(4L, 1L))
})

## @cvxpy NONE
test_that("cummax_expr: convex and increasing", {
  x <- Variable(3)
  e <- cummax_expr(x)
  expect_true(is_convex(e))
})

## @cvxpy NONE
test_that("cummax_expr: numeric value", {
  x <- Constant(matrix(c(1, 3, 2, 4), ncol = 1L))
  e <- cummax_expr(x)
  v <- numeric_value(e, list(matrix(c(1, 3, 2, 4), ncol = 1L)))
  expect_equal(as.numeric(v), c(1, 3, 3, 4))
})

## @cvxpy NONE
test_that("cummax_expr: axis=1 numeric", {
  m <- matrix(c(1, 4, 3, 2), nrow = 2, ncol = 2)
  e <- Cummax(Constant(m), axis = 1L)
  v <- numeric_value(e, list(m))
  ## cummax across rows: row 1: [1,3] → [1,3], row 2: [4,2] → [4,4]
  expect_equal(v[1, ], c(1, 3))
  expect_equal(v[2, ], c(4, 4))
})

## @cvxpy NONE
test_that("cummax_expr: has dcp canon", {
  expect_true(has_dcp_canon(Cummax(Variable(3))))
})

## @cvxpy NONE
test_that("cummax_expr: solves correctly", {
  skip_if_not_installed("clarabel")
  x <- Variable(4)
  prob <- Problem(Minimize(sum(cummax_expr(x))),
                  list(x == c(1, 3, 2, 4)))
  val <- psolve(prob)
  ## cummax([1,3,2,4]) = [1,3,3,4], sum = 11
  expect_equal(val, 11, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# dotsort
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("dotsort: shape is scalar", {
  x <- Variable(4)
  w <- c(1, 1, 1, 0)
  e <- dotsort(x, w)
  expect_equal(e@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("dotsort: convex and PWL", {
  x <- Variable(4)
  w <- c(1, 1, 1, 0)
  e <- dotsort(x, w)
  expect_true(is_convex(e))
  expect_true(is_pwl(e))
})

## @cvxpy NONE
test_that("dotsort: W must be constant", {
  x <- Variable(3)
  y <- Variable(3)
  expect_error(dotsort(x, y), "must be constant")
})

## @cvxpy NONE
test_that("dotsort: size(W) <= size(X)", {
  x <- Variable(2)
  w <- c(1, 2, 3)
  expect_error(dotsort(x, w), "less than or equal")
})

## @cvxpy NONE
test_that("dotsort: numeric value", {
  x <- Constant(c(-1, 0, 1, 2))
  w <- Constant(c(1, 1, 1, 0))
  e <- dotsort(x, w)
  v <- numeric_value(e, list(matrix(c(-1, 0, 1, 2), ncol = 1L),
                              matrix(c(1, 1, 1, 0), ncol = 1L)))
  ## sort(x)=[-1,0,1,2], sort(w)=[0,1,1,1], dot = 0+0+1+2 = 3
  expect_equal(as.numeric(v), 3)
})

## @cvxpy NONE
test_that("dotsort: has dcp canon", {
  expect_true(has_dcp_canon(Dotsort(Variable(3), Constant(c(3, 1, 2)))))
})

## @cvxpy NONE
test_that("dotsort: solves pinned values", {
  skip_if_not_installed("clarabel")
  x <- Variable(4)
  w <- c(1, 1, 1, 0)
  prob <- Problem(Minimize(dotsort(x, w)),
                  list(x == c(-1, 0, 1, 2)))
  val <- psolve(prob)
  ## sort(x)=[-1,0,1,2], sort(w)=[0,1,1,1], dot = 3
  expect_equal(val, 3, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("dotsort: weighted minimize", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  w <- c(3, 2, 1)
  prob <- Problem(Minimize(dotsort(x, w)),
                  list(x >= c(1, 2, 3), x <= c(1, 2, 3)))
  val <- psolve(prob)
  ## sort(x)=[1,2,3], sort(w)=[1,2,3], dot = 1+4+9 = 14
  expect_equal(val, 14, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("dotsort: sum_largest equivalence", {
  skip_if_not_installed("clarabel")
  x <- Variable(5)
  vals <- c(1, 5, 2, 4, 3)
  ## dotsort(x, [1,1,1]) should equal sum_largest(x, 3)
  prob1 <- Problem(Minimize(dotsort(x, c(1, 1, 1))),
                   list(x == vals))
  prob2 <- Problem(Minimize(sum_largest(x, 3)),
                   list(x == vals))
  v1 <- psolve(prob1)
  v2 <- psolve(prob2)
  expect_equal(v1, v2, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# S7 canon generics exist
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("dcp_canonicalize and has_dcp_canon are S7 generics", {
  expect_true(inherits(dcp_canonicalize, "S7_generic"))
  expect_true(inherits(has_dcp_canon, "S7_generic"))
})
