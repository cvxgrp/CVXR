## @cvxpy cvxpy/reductions/dqcp2dcp/inverse.py:56,92-98
## @cvxpy cvxpy/reductions/dqcp2dcp/sets.py:152-167,176,192
##
## DQCP dispatch keys on the ELEMENTWISE `multiply` atom, never on a matrix
## product. Upstream uses narrow tests at every site:
##
##   inverse.py:56   type(expr) == atoms.multiply          (exact type)
##   inverse.py:93   isinstance(expr, atoms.multiply)      (narrow class)
##   sets.py:176     SUBLEVEL_SETS[type(expr)]             (exact type)
##   sets.py:192     SUPERLEVEL_SETS[type(expr)]           (exact type)
##
## CVXR's Multiply EXTENDS MulExpression (binary_operators.R:190), so testing
## the parent swept genuine matrix products in. `.dqcp_inverse(A %*% x)` then
## returned `function(t) DivExpression(t, A)` -- elementwise division by the
## matrix, where the true inverse of a matrix product needs a solve. Likewise
## `.mul_sub` builds `y <= t * inv_pos(x)`, an elementwise relation that is not
## the sublevel set of `A %*% x`.
##
## Multiply is the only subclass of MulExpression, and DivExpression and
## AddExpression have none, so testing Multiply is exactly upstream's semantics.
##
## Ledger + paired CVXPY reproducer:
##   notes/audit/completeness_ledger.md finding 6
##   notes/audit/repro/03_dqcp_matmul.{py,R}

library(testthat)
library(CVXR)

test_that("Multiply is a subclass of MulExpression (the premise of the bug)", {
  x <- Variable(2)
  expect_true(inherits(x * c(1, 2), "CVXR::Multiply"))
  expect_true(inherits(x * c(1, 2), "CVXR::MulExpression"))
  ## ... and a matrix product is a MulExpression but NOT a Multiply.
  A <- Constant(matrix(c(1, 3, 2, 4), 2, 2))
  expect_true(inherits(A %*% x, "CVXR::MulExpression"))
  expect_false(inherits(A %*% x, "CVXR::Multiply"))
})

test_that("a matrix product is not DQCP-invertible", {
  A <- Constant(matrix(c(1, 3, 2, 4), 2, 2))
  x <- Variable(2)
  expect_false(.dqcp_invertible(A %*% x))
  expect_error(.dqcp_inverse(A %*% x), "Cannot compute inverse")
})

test_that("elementwise and scalar products stay DQCP-invertible", {
  ## The narrowing must not cost anything that worked. These are the forms
  ## CVXPY builds as `multiply`, verified against cvxpy 1.9.2.
  x <- Variable(2)
  expect_true(.dqcp_invertible(x * c(1, 2)))     # elementwise by a vector
  expect_true(.dqcp_invertible(2 * x))           # scalar on the left
  expect_true(.dqcp_invertible(x * 2))           # scalar on the right
  expect_true(.dqcp_invertible(x / 2))           # DivExpression, untouched
  expect_true(.dqcp_invertible(x + 1))           # AddExpression, untouched

  ## ... and each still yields a usable inverse.
  for (e in list(x * c(1, 2), 2 * x, x * 2)) {
    f <- .dqcp_inverse(e)
    expect_true(is.function(f))
    expect_true(inherits(f(Variable(2)), "CVXR::Expression"))
  }
})

test_that("a product of two variables is still not invertible", {
  ## Both arguments non-constant: `_non_const_idx()` has length 2.
  s <- Variable()
  expect_false(.dqcp_invertible(s * s))
})

test_that("sublevel/superlevel sets reject a matrix product", {
  A <- Constant(matrix(c(1, 3, 2, 4), 2, 2))
  x <- Variable(2)
  tt <- Variable()
  expect_error(.dqcp_sublevel(A %*% x, tt), "not yet supported in DQCP")
  expect_error(.dqcp_superlevel(A %*% x, tt), "not yet supported in DQCP")

  ## The elementwise forms still resolve, given the signs those handlers
  ## require. mul_sub needs OPPOSITE known signs and mul_sup needs MATCHING
  ## ones (sets.py mul_sub/mul_sup); an unsigned operand raises "Incorrect
  ## signs" in both libraries, so the operands here are declared nonneg.
  xp <- Variable(2, nonneg = TRUE)
  sub <- .dqcp_sublevel(xp * c(-1, -2), tt)     # nonneg * nonpos
  sup <- .dqcp_superlevel(xp * c(1, 2), tt)     # nonneg * nonneg
  expect_true(is.list(sub) && length(sub) == 1L)
  expect_true(is.list(sup) && length(sup) == 1L)
  expect_s3_class(sub[[1L]], "CVXR::Constraint")
  expect_s3_class(sup[[1L]], "CVXR::Constraint")

  ## Unsigned operands are rejected by the handler, not by the dispatch --
  ## a DIFFERENT error than the matrix-product one above, which is what shows
  ## the elementwise product still reaches .mul_sub at all.
  expect_error(.dqcp_sublevel(x * c(1, 2), tt), "Incorrect signs")
})

test_that("a ratio DQCP problem still solves after the narrowing", {
  ## End-to-end guard on the path that legitimately uses the Multiply branch:
  ## maximize a ratio, which bisection inverts through DivExpression/Multiply.
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Maximize(x / (x + 1)), list(x <= 4))
  val <- psolve(prob, qcp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 4 / 5, tolerance = 1e-4)
})
