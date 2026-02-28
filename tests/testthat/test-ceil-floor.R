## Tests for ceil and floor atoms
## Phase 2: DQCP-only atoms (quasiconvex AND quasiconcave, not DCP)

## @cvxpy NONE
test_that("Ceil numeric evaluation", {
  x <- Variable()
  value(x) <- 2.3
  expect_equal(as.numeric(value(ceil_expr(x))), 3)

  value(x) <- -1.7
  expect_equal(as.numeric(value(ceil_expr(x))), -1)

  value(x) <- 5.0
  expect_equal(as.numeric(value(ceil_expr(x))), 5)
})

## @cvxpy NONE
test_that("Floor numeric evaluation", {
  x <- Variable()
  value(x) <- 2.7
  expect_equal(as.numeric(value(floor_expr(x))), 2)

  value(x) <- -1.3
  expect_equal(as.numeric(value(floor_expr(x))), -2)

  value(x) <- 5.0
  expect_equal(as.numeric(value(floor_expr(x))), 5)
})

## @cvxpy NONE
test_that("Ceil tolerance rounding", {
  ## CVXPY: ceil(np.around(2.9999, decimals=4)) = 3, not 4
  x <- Variable()
  value(x) <- 2.9999
  expect_equal(as.numeric(value(ceil_expr(x))), 3)
})

## @cvxpy NONE
test_that("ceil via Math handler (ceiling())", {
  x <- Variable()
  expr <- ceiling(x)
  expect_true(S7_inherits(expr, Ceil))

  value(x) <- 2.3
  expect_equal(as.numeric(value(expr)), 3)
})

## @cvxpy NONE
test_that("floor via Math handler (floor())", {
  x <- Variable()
  expr <- floor(x)
  expect_true(S7_inherits(expr, Floor))

  value(x) <- 2.7
  expect_equal(as.numeric(value(expr)), 2)
})

## @cvxpy NONE
test_that("Ceil sign preservation", {
  ## nonneg arg -> nonneg result
  x <- Variable(nonneg = TRUE)
  expr <- ceil_expr(x)
  expect_true(is_nonneg(expr))
  expect_false(is_nonpos(expr))

  ## nonpos arg -> nonpos result
  y <- Variable(nonpos = TRUE)
  expr2 <- ceil_expr(y)
  expect_false(is_nonneg(expr2))
  expect_true(is_nonpos(expr2))

  ## unknown sign
  z <- Variable()
  expr3 <- ceil_expr(z)
  expect_false(is_nonneg(expr3))
  expect_false(is_nonpos(expr3))
})

## @cvxpy NONE
test_that("Floor sign preservation", {
  x <- Variable(nonneg = TRUE)
  expr <- floor_expr(x)
  expect_true(is_nonneg(expr))

  y <- Variable(nonpos = TRUE)
  expr2 <- floor_expr(y)
  expect_true(is_nonpos(expr2))
})

## @cvxpy NONE
test_that("Ceil is NOT convex and NOT concave", {
  x <- Variable()
  expr <- ceil_expr(x)
  expect_false(is_atom_convex(expr))
  expect_false(is_atom_concave(expr))
  expect_false(is_convex(expr))
  expect_false(is_concave(expr))
})

## @cvxpy NONE
test_that("Floor is NOT convex and NOT concave", {
  x <- Variable()
  expr <- floor_expr(x)
  expect_false(is_atom_convex(expr))
  expect_false(is_atom_concave(expr))
  expect_false(is_convex(expr))
  expect_false(is_concave(expr))
})

## @cvxpy NONE
test_that("Ceil is quasiconvex AND quasiconcave", {
  x <- Variable()
  expr <- ceil_expr(x)
  expect_true(is_atom_quasiconvex(expr))
  expect_true(is_atom_quasiconcave(expr))
  expect_true(is_quasiconvex(expr))
  expect_true(is_quasiconcave(expr))
  expect_true(is_quasilinear(expr))
  expect_true(is_dqcp(expr))
})

## @cvxpy NONE
test_that("Floor is quasiconvex AND quasiconcave", {
  x <- Variable()
  expr <- floor_expr(x)
  expect_true(is_atom_quasiconvex(expr))
  expect_true(is_atom_quasiconcave(expr))
  expect_true(is_quasiconvex(expr))
  expect_true(is_quasiconcave(expr))
  expect_true(is_quasilinear(expr))
  expect_true(is_dqcp(expr))
})

## @cvxpy NONE
test_that("Ceil is increasing", {
  x <- Variable()
  expr <- Ceil(x)
  expect_true(is_incr(expr, 1L))
  expect_false(is_decr(expr, 1L))
})

## @cvxpy NONE
test_that("Floor is increasing", {
  x <- Variable()
  expr <- Floor(x)
  expect_true(is_incr(expr, 1L))
  expect_false(is_decr(expr, 1L))
})

## @cvxpy NONE
test_that("Ceil not DCP-compliant in DCP problem", {
  x <- Variable()
  prob <- Problem(Minimize(ceil_expr(x)), list(x >= 0.5))
  expect_false(is_dcp(prob))
  expect_true(is_dqcp(prob))
})

## @cvxpy NONE
test_that("Floor not DCP-compliant in DCP problem", {
  x <- Variable()
  prob <- Problem(Minimize(floor_expr(x)), list(x >= 0.5))
  expect_false(is_dcp(prob))
  expect_true(is_dqcp(prob))
})

## @cvxpy NONE
test_that("Vector ceil/floor", {
  x <- Variable(3)
  value(x) <- c(1.2, 2.7, -0.5)
  expect_equal(as.numeric(value(ceil_expr(x))), c(2, 3, 0))
  expect_equal(as.numeric(value(floor_expr(x))), c(1, 2, -1))
})

## @cvxpy NONE
test_that("Ceil/Floor names", {
  x <- Variable()
  expect_match(expr_name(ceil_expr(x)), "Ceil")
  expect_match(expr_name(floor_expr(x)), "Floor")
})

## @cvxpy NONE
test_that("Ceil of constant", {
  c <- Constant(2.3)
  expr <- ceil_expr(c)
  expect_equal(as.numeric(value(expr)), 3)
  ## Constants are always convex/concave/quasiconvex etc
  expect_true(is_quasiconvex(expr))
  expect_true(is_quasiconcave(expr))
})

# ══════════════════════════════════════════════════════════════════
# length_expr — last nonzero index
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("length_expr: curvature properties", {
  x <- Variable(5)
  expr <- length_expr(x)
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_false(is_quasiconcave(expr))
  expect_false(is_convex(expr))
  expect_false(is_concave(expr))
})

## @cvxpy NONE
test_that("length_expr: sign is always nonneg", {
  x <- Variable(5)
  expr <- length_expr(x)
  expect_true(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## @cvxpy NONE
test_that("length_expr: only accepts vectors", {
  X <- Variable(c(3, 3))
  expect_error(length_expr(X), "vectors")
})

## @cvxpy NONE
test_that("length_expr: numeric evaluation", {
  x <- Constant(c(1, 2, 0, 0, 0))
  expect_equal(as.numeric(value(length_expr(x))), 2)

  y <- Constant(c(0, 0, 0, 0, 3))
  expect_equal(as.numeric(value(length_expr(y))), 5)

  z <- Constant(c(0, 0, 0, 0, 0))
  expect_equal(as.numeric(value(length_expr(z))), 0)
})

## @cvxpy NONE
test_that("length_expr: monotonicity", {
  x <- Variable(5)
  ## increasing when arg is nonneg
  expect_true(is_incr(length_expr(abs(x)), 2L))
  expect_true(is_dqcp(length_expr(abs(x))))
  ## abs(x) - 1 is NOT nonneg, so length is not monotone → not DQCP
  expect_false(is_dqcp(length_expr(abs(x) - 1)))
  ## decreasing when arg is nonpos
  expect_true(is_decr(length_expr(-abs(x)), 2L))
})

## @cvxpy NONE
test_that("length_expr: solve basic DQCP problem", {
  x <- Variable(5)
  prob <- Problem(Minimize(length_expr(x)), list(x[1] == 2.0, x[2] == 1.0))
  r <- psolve(prob, qcp = TRUE)
  expect_equal(r, 2.0, tolerance = 1e-3)
  xv <- as.numeric(value(x))
  expect_equal(xv[1], 2.0, tolerance = 1e-3)
  expect_equal(xv[2], 1.0, tolerance = 1e-3)
  expect_equal(xv[3], 0.0, tolerance = 1e-3)
  expect_equal(xv[4], 0.0, tolerance = 1e-3)
  expect_equal(xv[5], 0.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("length_expr: min-length least squares", {
  set.seed(1)
  n <- 10
  A <- matrix(rnorm(n * n), n, n)
  x_star <- rnorm(n)
  b <- A %*% x_star
  epsilon <- 1e-2
  x <- Variable(n)
  mse <- sum_squares(A %*% x - b) / n
  prob <- Problem(Minimize(length_expr(x)), list(mse <= epsilon))
  expect_true(is_dqcp(prob))
  r <- psolve(prob, qcp = TRUE)
  ## R's rnorm(seed=1) gives different data than numpy; answer is 9 for R
  expect_equal(r, 9.0, tolerance = 1e-3)
})
