## @cvxpy tests/test_problem.py, tests/test_complex.py
##
## Three places where CVXR and CVXPY 1.9.2 disagreed about whether a problem is
## solvable at all: two where CVXPY succeeds and CVXR errored, one where CVXPY
## errors and CVXR let poisoned data reach the solver.
##
## CVXPY SOURCE: solving_chain.py:412-449 (_validate_problem_data),
## solving_chain.py:508-512 and :556-560 (the two entry points),
## complex2real.py:42-44 (accepts).

## ---------------------------------------------------------------------------
## Inf is allowed in the constraint RHS and in bounds
##
## CVXPY's `inf_allowed_keys = {s.B, s.G, s.LOWER_BOUNDS, s.UPPER_BOUNDS}`
## (solving_chain.py:427) exists because, in its own words at :419-420, "users
## sometimes use inf for unbounded constraints/variables". CVXR ran every key
## through one indiscriminate NaN-or-Inf check, so both problems below aborted
## with "Problem data 'b' contains Inf values" where CVXPY 1.9.2 returns 2.0 and
## -Inf respectively.
## ---------------------------------------------------------------------------

test_that("Inf in a constraint RHS is accepted, as upstream", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1, x <= Inf))
  expect_equal(as.numeric(psolve(prob, solver = "CLARABEL")), 2, tolerance = 1e-6)
  expect_equal(status(prob), "optimal")

  y <- Variable(2)
  prob2 <- Problem(Minimize(sum_entries(y)), list(y >= c(1, -Inf)))
  expect_equal(as.numeric(psolve(prob2, solver = "CLARABEL")), -Inf)
  expect_equal(status(prob2), "unbounded")
})

test_that("Inf in a variable bound is accepted", {
  z <- Variable(2, bounds = list(c(1, -Inf), c(Inf, Inf)))
  prob <- Problem(Minimize(sum_entries(z)), list(z >= -5))
  expect_equal(as.numeric(psolve(prob, solver = "CLARABEL")), -4, tolerance = 1e-6)
})

## Widening the Inf rule must NOT weaken the NaN rule -- NaN is still caught on
## every key, including the Inf-allowed ones (upstream :438-443 NaN-checks those).
test_that("NaN is still rejected everywhere, including Inf-allowed keys", {
  x <- Variable(2)
  expect_error(psolve(Problem(Minimize(sum_entries(x)), list(x >= c(1, NaN))),
                      solver = "CLARABEL"), "NaN")
})

## A NaN never reaches the chain from a Parameter, because assignment rejects
## it first -- but it must do so with the message the cascade exists to
## produce. In R every comparison against NaN is NA, so `if (!close_enough)`
## raised the internal "missing value where TRUE/FALSE needed" instead.
## `np.allclose` returns False for NaN (leaf.py:640-651), so CVXPY 1.9.2 says
## "Parameter value must be real." -- measured.
test_that("assigning NaN to a leaf reports the attribute, not an R internal", {
  p <- Parameter(2)
  expect_error({ value(p) <- c(NaN, 1) }, "must be real")
  q <- Parameter(2, nonneg = TRUE)
  expect_error({ value(q) <- c(NaN, 1) }, "must be nonnegative")
  v <- Variable(2)
  expect_error({ value(v) <- c(1, NaN) }, "must be real")
  ## Inf is NOT NaN and must still be assignable: the Inf - Inf path that the
  ## NaN mask exists for keeps working.
  w <- Parameter(2)
  expect_no_error({ value(w) <- c(Inf, 1) })
  expect_equal(as.numeric(value(w)), c(Inf, 1))
})

## Inf in the OBJECTIVE or the constraint MATRIX is still an error: those keys
## are not in inf_allowed_keys.
test_that("Inf outside the allowed keys is still rejected", {
  x <- Variable(2)
  q <- Parameter(2)
  value(q) <- c(Inf, 1)
  expect_error(psolve(Problem(Minimize(sum_entries(q * x)), list(x >= 1, x <= 2)),
                      solver = "CLARABEL"), "Inf")
})

## ---------------------------------------------------------------------------
## The decomposed API is the second entry point and must validate too
## (CVXPY solving_chain.py:556-560, whose comment names both entry points)
##
## Pre-fix, a NaN planted in `c` went straight through solve_via_data() to
## Clarabel, which returned status 10 and garbage that propagated to the caller.
## ---------------------------------------------------------------------------

test_that("solve_via_data() validates the data it is handed", {
  w <- Variable(2)
  pr <- Problem(Minimize(sum_entries(w)), list(w >= 1))
  pd <- problem_data(pr, solver = "CLARABEL")
  poisoned <- pd$data
  poisoned$c[1] <- NaN
  expect_error(solve_via_data(pd$chain, poisoned, FALSE, FALSE, list()), "NaN")
})

test_that("solve_via_data() still solves clean data", {
  w <- Variable(2)
  pr <- Problem(Minimize(sum_entries(w)), list(w >= 1))
  pd <- problem_data(pr, solver = "CLARABEL")
  raw <- solve_via_data(pd$chain, pd$data, FALSE, FALSE, list())
  soln <- problem_unpack_results(pr, raw, pd$chain, pd$inverse_data)
  expect_equal(status(pr), "optimal")
  expect_equal(as.numeric(value(w)), c(1, 1), tolerance = 1e-6)
})

## ---------------------------------------------------------------------------
## A complex Parameter must trigger Complex2Real
##
## `problem.constants()` collects only Constant leaves, so upstream's third
## term `problem.parameters()` (complex2real.py:42-44) is load-bearing. Without
## it the reduction was skipped, but `has_complex_params` still prepended
## EvalParams, turning the complex Parameter into a complex CONSTANT with no
## Complex2Real downstream -- so the solve died with
## "Inequality constraints cannot be complex." CVXPY 1.9.2: 4.123105626514879.
## ---------------------------------------------------------------------------

test_that("a problem whose only complex leaf is a Parameter is accepted", {
  x <- Variable(2)
  p <- Parameter(c(1, 2), complex = TRUE)
  value(p) <- matrix(c(1 + 2i, 3 - 1i), nrow = 1)
  prob <- Problem(Minimize(abs(p %*% x)), list(x >= 1))
  expect_true(CVXR:::complex2real_accepts(prob))
  expect_equal(as.numeric(psolve(prob, solver = "CLARABEL")),
               4.123105626514879, tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-5)
})

test_that("complex Variables and Constants still trigger Complex2Real", {
  z <- Variable(2, complex = TRUE)
  expect_true(CVXR:::complex2real_accepts(
    Problem(Minimize(sum_entries(abs(z))), list(Re(z) >= 1))))
  w <- Variable(2)
  expect_true(CVXR:::complex2real_accepts(
    Problem(Minimize(abs(Constant(matrix(c(1 + 1i, 2), nrow = 1)) %*% w)),
            list(w >= 1))))
  ## and a wholly real problem still does not
  v <- Variable(2)
  expect_false(CVXR:::complex2real_accepts(
    Problem(Minimize(sum_entries(v)), list(v >= 1))))
})
