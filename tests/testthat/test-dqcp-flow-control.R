## @cvxpy tests/test_dqcp.py
##
## Flow-control divergences in the DQCP lowering and the bisection driver:
## places where CVXR's transcription of an upstream branch carried an extra
## condition, or dropped one, and so refused a problem CVXPY solves.
##
## CVXPY SOURCE: dqcp2dcp.py:216-218 (case 2 is the plain `else` of case 1),
## bisection.py:74,88 (the zero-endpoint guard),
## logic.py:24-34 (_is_boolean_arg).

## ---------------------------------------------------------------------------
## Bisection: an endpoint sitting exactly at 0 must still be able to grow
##
## CVXPY writes `high = high * 2 if high != 0 else 1` (bisection.py:74) and
## `low = low * 2 if low != 0 else -1` (:88). CVXR used a bare `* 2`, which
## cannot move an endpoint off 0 -- the bracket never widens and the search
## burns all max_iters before failing with "Unable to find suitable interval".
## The defaults only land on 0 in the branch that is skipped, so this is
## reachable through an explicitly supplied low/high.
## ---------------------------------------------------------------------------

## The problem must be genuinely NON-DCP, or psolve() never enters the
## bisection path at all and the test measures nothing. `x / y` is
## quasiconvex and not convex; `inv_pos(x)` would NOT do -- it is DCP.
.ratio_prob <- function() {
  x <- Variable()
  y <- Variable(nonneg = TRUE)
  list(prob = Problem(Minimize(x / y), list(x >= -4, x <= -1, y >= 1, y <= 2)),
       x = x, y = y)
}

test_that("the ratio problem is DQCP and not DCP", {
  h <- .ratio_prob()
  expect_true(is_dqcp(h$prob))
  expect_false(is_dcp(h$prob))
})

test_that("an explicit zero bisection endpoint still brackets the solution", {
  ## Optimum is -4 (x = -4, y = 1), so the objective is NEGATIVE and `low = 0`
  ## is a legitimate starting endpoint the search must be able to move DOWN
  ## from. With the bare `low * 2` it could not: low stayed at 0 for all
  ## max_iters and the solve died. Measured pre-fix: "Unable to find suitable
  ## interval for bisection." CVXPY 1.9.2 returns -4 for all three spellings.
  h <- .ratio_prob()
  expect_equal(as.numeric(psolve(h$prob, qcp = TRUE, solver = "CLARABEL", low = 0)),
               -4, tolerance = 1e-4)

  h2 <- .ratio_prob()
  expect_equal(as.numeric(psolve(h2$prob, qcp = TRUE, solver = "CLARABEL", high = 0)),
               -4, tolerance = 1e-4)
})

test_that("the ordinary bisection path is unchanged", {
  h <- .ratio_prob()
  expect_equal(as.numeric(psolve(h$prob, qcp = TRUE, solver = "CLARABEL")),
               -4, tolerance = 1e-4)
  expect_equal(as.numeric(value(h$x)), -4, tolerance = 1e-3)
  expect_equal(as.numeric(value(h$y)), 1, tolerance = 1e-3)
})

## ---------------------------------------------------------------------------
## Dqcp2Dcp case 2 is the plain `else` of case 1
##
## Upstream reaches "constant <= quasiconcave" as the else-branch and only
## ASSERTS `rhs.is_quasiconcave()` (dqcp2dcp.py:216-218). CVXR made it a second
## `if` carrying an extra `&& !is_concave(rhs)`, so a quasiconcave-AND-concave
## RHS fell through to `cli_abort("Cannot canonicalize constraint in DQCP
## context.")`. It is reachable through the recursion: case 1 applies an inverse
## and re-enters with a CONSTANT lhs, for which `is_quasiconvex && !is_convex`
## is FALSE.
## ---------------------------------------------------------------------------

test_that("a concave quasiconcave RHS is lowered, not rejected", {
  ## sqrt is concave AND quasiconcave. After the ratio's inverse is applied the
  ## rewritten constraint arrives with a constant lhs and a concave rhs.
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  prob <- Problem(Maximize(sqrt(x)), list(sqrt(x) >= 1, x <= 9, y >= 0))
  expect_true(is_dqcp(prob))
  expect_no_error(psolve(prob, qcp = TRUE, solver = "CLARABEL"))
  expect_equal(as.numeric(psolve(prob, qcp = TRUE, solver = "CLARABEL")),
               3, tolerance = 1e-4)
})

test_that("a genuinely non-quasiconcave RHS is still rejected", {
  ## The abort must survive: only the extra conjunct was removed, not the check.
  expect_error(
    CVXR:::.dqcp_canonicalize_constraint(
      Inequality(Constant(1), square(Variable()))),
    "Cannot canonicalize")
})

## ---------------------------------------------------------------------------
## _is_boolean_arg accepts a partial boolean index
##
## Upstream tests `arg.attributes.get('boolean')` for TRUTHINESS (logic.py:29),
## so an index list qualifies. `isTRUE()` is FALSE for every index spelling, so
## CVXR RAISED on a program CVXPY builds -- and disagreed with its own
## `@.boolean_idx`, which the solver path and project() both read.
## ---------------------------------------------------------------------------

test_that("logic atoms accept every boolean spelling of a Variable", {
  whole <- Variable(2, boolean = TRUE)
  expect_no_error(And(whole, whole))
  flat <- Variable(2, boolean = c(1, 2))
  expect_no_error(And(flat, whole))
  partial <- Variable(2, boolean = c(1))
  expect_no_error(And(partial, whole))
  expect_no_error(Or(partial, whole))
  expect_no_error(Not(partial))
  mask <- Variable(2, boolean = c(TRUE, FALSE))
  expect_no_error(And(mask, whole))
})

test_that("logic atoms still reject a non-boolean argument", {
  plain <- Variable(2)
  whole <- Variable(2, boolean = TRUE)
  expect_error(And(plain, whole), "boolean")
})

## ---------------------------------------------------------------------------
## The Constant arm of _is_boolean_arg
##
## CVXR deliberately diverges from logic.py:31-33 here: upstream reads a flag
## recorded from the original numpy dtype, so `Constant(np.array([1, 1]))` is
## NOT boolean-valued while the dtype=bool version is. That split is a numpy
## artifact with no R counterpart, so CVXR asks what the value IS. These pin the
## two guards that rule needs, and the divergence itself.
## ---------------------------------------------------------------------------

test_that("a NaN or NA constant is not boolean, and says so", {
  b <- Variable(2, boolean = TRUE)
  ## `all(v == 0 | v == 1)` is NA for a NaN, and the caller does
  ## `if (!.is_boolean_arg(arg))`, so this used to die with R's internal
  ## "missing value where TRUE/FALSE needed" -- the same NA-into-`if` defect
  ## fixed in .validate_leaf_value.
  expect_false(CVXR:::.is_boolean_arg(Constant(c(NaN, 1))))
  expect_false(CVXR:::.is_boolean_arg(Constant(c(NA_real_, 1))))
  expect_error(And(b, Constant(c(NaN, 1))), "must be boolean")
  expect_error(And(b, Constant(c(NA_real_, 1))), "must be boolean")
})

test_that("a zero-size constant is not boolean by vacuous truth", {
  ## all(logical(0)) is TRUE, which reported an empty Constant as boolean.
  expect_false(CVXR:::.is_boolean_arg(Constant(numeric(0))))
})

test_that("0/1-valued constants are still accepted, in both R spellings", {
  b <- Variable(2, boolean = TRUE)
  expect_true(CVXR:::.is_boolean_arg(Constant(c(1, 0))))
  expect_true(CVXR:::.is_boolean_arg(Constant(c(TRUE, FALSE))))
  expect_no_error(And(b, Constant(c(1, 0))))
  ## and a genuinely non-0/1 constant is still rejected
  expect_false(CVXR:::.is_boolean_arg(Constant(c(0.5, 1))))
  expect_false(CVXR:::.is_boolean_arg(Constant(c(2, 0))))
  expect_error(And(b, Constant(c(0.5, 1))), "must be boolean")
})

test_that("a sparse 0/1 constant is boolean", {
  M <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L), x = c(1, 1),
                            dims = c(2L, 2L))
  expect_true(CVXR:::.is_boolean_arg(Constant(M)))
  M2 <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 0.5, dims = c(2L, 2L))
  expect_false(CVXR:::.is_boolean_arg(Constant(M2)))
})
