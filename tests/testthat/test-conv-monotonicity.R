## @cvxpy cvxpy/atoms/affine/conv.py:64-70,104-112
##
## conv()'s monotonicity must follow the sign of the CONSTANT KERNEL, whichever
## argument holds it.
##
## Upstream reads args[0] unconditionally (conv.py:104-112). That is safe there
## only because validate_arguments (conv.py:64-70) rejects a non-constant first
## argument: "The first argument to conv must be constant." CVXR deliberately
## accepts the constant in EITHER position -- the constructor only requires that
## one of them be constant, and graph_implementation and is_dpp both locate it --
## so the monotonicity methods have to locate it too, the way Kron already does
## (kron.R:38-45, mirroring kron.py:77-87, for the atom where upstream itself
## permits either order).
##
## Before the fix, `conv(convex_expr, negative_kernel)` reported is_incr = TRUE,
## reading the sign of the VARIABLE rather than the kernel. Composition then
## called a concave expression convex, and is_dcp() accepted it for
## MINIMIZATION -- a soundness hole reachable only through CVXR's relaxation,
## since CVXPY cannot build that form at all.
##
## Ledger + paired CVXPY reproducer:
##   notes/audit/completeness_ledger.md finding 1
##   notes/audit/repro/01_conv_monotonicity.{py,R}

library(testthat)
library(CVXR)

test_that("conv monotonicity follows the kernel, not the variable", {
  y <- Variable(3, nonneg = TRUE)
  kpos <- c(1, 2)
  kneg <- c(-1, -2)

  ## Kernel first -- the only form CVXPY can express. Unchanged by the fix.
  expect_true(is_incr(conv(kpos, y), 1L))
  expect_false(is_decr(conv(kpos, y), 1L))
  expect_false(is_incr(conv(kneg, y), 1L))
  expect_true(is_decr(conv(kneg, y), 1L))

  ## Kernel SECOND -- CVXR-only. Must give the same answers: monotonicity is a
  ## property of the kernel, not of the argument's position.
  expect_true(is_incr(conv(y, kpos), 1L))
  expect_false(is_decr(conv(y, kpos), 1L))
  expect_false(is_incr(conv(y, kneg), 1L))
  expect_true(is_decr(conv(y, kneg), 1L))
})

test_that("conv monotonicity ignores the sign of the non-constant argument", {
  ## The pre-fix bug was reading the VARIABLE's sign. Pin that directly: a
  ## nonneg variable with a negative kernel must still be decreasing, and the
  ## answer must not change when the variable's declared sign changes.
  kneg <- c(-1, -2)
  for (v in list(Variable(3, nonneg = TRUE),
                 Variable(3, nonpos = TRUE),
                 Variable(3))) {
    expect_false(is_incr(conv(v, kneg), 1L))
    expect_true(is_decr(conv(v, kneg), 1L))
  }
})

test_that("conv(convex, negative kernel) is concave and is NOT minimizable", {
  ## The end-to-end consequence. Monotonicity only bites when the inner
  ## argument is non-affine, so feed conv a convex argument: with a negative
  ## kernel the composition is CONCAVE, and minimizing it is not DCP.
  y <- Variable(3)
  inner <- square(y)
  kneg <- c(-1, -2)

  for (e in list(conv(kneg, inner), conv(inner, kneg))) {
    expect_false(is_convex(e))
    expect_true(is_concave(e))
    expect_false(is_dcp(Problem(Minimize(sum_entries(e)))))
    expect_true(is_dcp(Problem(Maximize(sum_entries(e)))))
  }

  ## ... and with a POSITIVE kernel it is convex, minimizable, and not
  ## maximizable -- so the assertions above are discriminating, not vacuous.
  for (e in list(conv(c(1, 2), inner), conv(inner, c(1, 2)))) {
    expect_true(is_convex(e))
    expect_false(is_concave(e))
    expect_true(is_dcp(Problem(Minimize(sum_entries(e)))))
    expect_false(is_dcp(Problem(Maximize(sum_entries(e)))))
  }
})

test_that("conv with an unsigned kernel is monotone in neither direction", {
  y <- Variable(3)
  kmix <- c(-1, 2)
  expect_false(is_incr(conv(kmix, y), 1L))
  expect_false(is_decr(conv(kmix, y), 1L))
  expect_false(is_incr(conv(y, kmix), 1L))
  expect_false(is_decr(conv(y, kmix), 1L))
})

test_that("conv still solves correctly with the kernel in either position", {
  ## Guard against a fix that repairs the predicate but perturbs the numbers.
  ## graph_implementation already located the constant, so both orders must
  ## produce the same optimum.
  y <- Variable(3, nonneg = TRUE)
  k <- c(1, 2)
  target <- c(1, 3, 5, 2)

  v1 <- psolve(Problem(Minimize(sum_squares(conv(k, y) - target))),
               solver = "CLARABEL")
  y1 <- value(y)
  v2 <- psolve(Problem(Minimize(sum_squares(conv(y, k) - target))),
               solver = "CLARABEL")

  expect_equal(v1, v2, tolerance = 1e-6)
  expect_equal(y1, value(y), tolerance = 1e-6)
})
