## Bisection behavior when the solver FAILS at a query point -- neither
## reporting feasibility nor proving infeasibility.
##
## CVXPY SOURCE: bisection.py:106-155 (CVXPY 1.9.2).
##
## Upstream through 1.9.1 re-solved the identical failing subproblem until it
## raised "Max iters hit during bisection".  CVXR never had that infinite
## retry, but avoided it unsoundly: it treated the failure as infeasibility
## and moved `low`, the reported lower bound, on a point that had never been
## proven infeasible.  1.9.2 keeps `low` fixed and perturbs only the point the
## next query is drawn from.

## @cvxpy test_dqcp.py::TestDqcp::test_multiply_solver_failure_during_bisection
test_that("bisection survives a solver failure near the optimum", {
  skip(paste(
    "Scenario no longer reproducible after sqrt()/inv_pos() moved to the SOC",
    "path; see the block comment below."))
  ## WHY THIS IS SKIPPED, and why it was not simply re-pointed at another
  ## solver.  All figures below were measured, not estimated.
  ##
  ## The case needs a solver that FAILS at a query point -- neither feasible
  ## nor proven infeasible -- or it stops testing the thing it is named for.
  ## Under the PowCone3D canonicalization this scaling produced exactly one
  ## such failure, which the bisection then survived.  Once sqrt()/inv_pos()
  ## began canonicalizing to SOC (matching CVXPY: sqrt.py:19, inv_pos.py:7)
  ## the numerics shifted and CLARABEL now fails TWICE here, one of them in
  ## the interval finder, which aborts by design -- upstream's
  ## _find_bisection_interval raises on failure too (bisection.py:77-79), so
  ## this is not an unported branch.
  ##
  ## Retuning the scaling does not recover the scenario.  Measured over
  ## x + y <= {1e-2, 1e-3, 5e-4, 2e-4, 1e-4, 1e-5}: every bound gives either
  ## 0 failures and a clean solve, or 2 failures and an abort.  There is no
  ## window with one survivable failure.
  ##
  ## Re-pinning the solver would be a FALSE GREEN.  SCS and ECOS both pass
  ## this case -- in 23 probes with ZERO failures, so they never reach the
  ## code under test.  A test named "survives a solver failure" that never
  ## sees one is worse than a skipped test.
  ##
  ## CVXPY solves this problem (optimal_inaccurate, 1.897e-09) not through
  ## different logic but by drawing AlmostSolved -- a SOLUTION_PRESENT
  ## status -- on the probe where CLARABEL gives us InsufficientProgress.
  ##
  ## `accept_unknown` DOES NOT HELP HERE, contrary to the obvious guess (which
  ## an earlier version of this comment made).  bisection.py:32-42 `_solve`
  ## calls `problem.solve(solver=solver)` -- the solver ONLY -- so upstream
  ## forwards no solver options into subproblems, and neither does
  ## .bisect_solve().  The option cannot reach the stalling solve in EITHER
  ## library.  See test-clarabel-accept-unknown.R, which pins that limit.
  ##
  ## So there is nothing to wait for: restoring this test needs the underlying
  ## numerics to change (a CLARABEL that returns AlmostSolved here), not a
  ## CVXR feature.
  ##
  ## Badly scaled: optimal value is (5e-5)^2 = 2.5e-9 and the bisection
  ## tolerance is 1e-6, so the answer is 0 to that tolerance.
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  prob <- Problem(Maximize(x * y), list(x + y <= 1e-4))
  val <- psolve(prob, solver = "CLARABEL", qcp = TRUE)
  expect_true(status(prob) %in% SOLUTION_PRESENT)
  expect_equal(val, 2.5e-9, tolerance = 1e-6)
})

## @cvxpy NONE -- the well-behaved path must be untouched by the query_low
## change: with no solver failure, query_low tracks low exactly, so the
## bracket and the answer are identical to before.
test_that("ordinary DQCP bisection is unaffected", {
  x <- Variable()
  y <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(ceil_expr(x)), list(x >= 11.9, x <= 17))
  val <- psolve(prob, solver = "CLARABEL", qcp = TRUE)
  expect_equal(val, 12)

  ## A ratio objective, the canonical quasiconvex case.
  a <- Variable(nonneg = TRUE)
  b <- Variable(nonneg = TRUE)
  rp <- Problem(Minimize((a + 1) / (b + 1)), list(a >= 1, b <= 3, a <= 2))
  rv <- psolve(rp, solver = "CLARABEL", qcp = TRUE)
  expect_equal(rv, 0.5, tolerance = 1e-4)
})
