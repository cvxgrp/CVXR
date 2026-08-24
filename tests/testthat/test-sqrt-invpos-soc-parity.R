## @cvxpy tests/test_atoms.py, tests/test_conic_solvers.py
##
## `sqrt()` and `inv_pos()` must canonicalize through the SOC path, as CVXPY
## does -- sqrt.py:19 is `power(x, Fraction(1, 2))` and inv_pos.py:7 is
## `power(x, -1)`, both the `power()` WRAPPER, whose `approx = TRUE` default
## selects PowerApprox (SOC).  Building the `Power` class directly takes the
## exact path and emits PowCone3D.
##
## Why this file exists: CVXR built `Power` directly at both sites through
## 1.9.1, so every problem containing `sqrt()` demanded a power-cone solver.
## ECOS -- which has SOC and ExpCone but no power cone -- rejected such
## problems outright, and SCS returned `optimal_inaccurate` with a NaN
## objective and an out-of-domain iterate on the DQCP form.  Nothing failed
## loudly, because the default solver (CLARABEL) does support power cones, so
## the whole suite stayed green while CVXR handed solvers a different cone
## program than CVXPY for every square root in the package.
##
## The assertions below are therefore about WHICH CONE IS EMITTED, expressed
## as "a SOC-only solver can solve it". Testing the value alone would not have
## caught the bug -- CLARABEL returns the right value either way.

test_that("sqrt() canonicalizes to SOC, not PowCone3D", {
  skip_if_not("ECOS" %in% installed_solvers(), "ECOS not installed")

  ## The bug: this errored with
  ##   Solver "ECOS" does not support the required cone types: PowCone3D
  x <- Variable(nonneg = TRUE)
  p <- Problem(Maximize(sqrt(x)), list(x <= 4))
  expect_equal(psolve(p, solver = "ECOS"), 2, tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), 4, tolerance = 1e-6)

  ## p = 1/2 has denominator 2, so PowerApprox represents it with NO rational
  ## approximation error -- switching off the power cone costs no accuracy
  ## beyond the solver's own tolerance.  Guards against someone "fixing" this
  ## by raising max_denom, or reverting it on accuracy grounds.
  expect_equal(psolve(p, solver = "CLARABEL"), 2, tolerance = 1e-6)
})

test_that("inv_pos() canonicalizes to SOC, not PowCone3D", {
  skip_if_not("ECOS" %in% installed_solvers(), "ECOS not installed")

  ## min 1/x  s.t. x <= 2, x >= 0.5  ->  x* = 2, value 0.5
  x <- Variable(nonneg = TRUE)
  p <- Problem(Minimize(inv_pos(x)), list(x <= 2, x >= 0.5))
  expect_equal(psolve(p, solver = "ECOS"), 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), 2, tolerance = 1e-6)
})

test_that("the explicit exact path is still reachable", {
  ## `approx = FALSE` must keep emitting PowCone3D -- the SOC default is a
  ## default, not the removal of a capability.  CLARABEL supports both.
  x <- Variable(nonneg = TRUE)
  p_exact <- Problem(Maximize(power(x, 0.5, approx = FALSE)), list(x <= 4))
  expect_equal(psolve(p_exact, solver = "CLARABEL"), 2, tolerance = 1e-6)

  skip_if_not("ECOS" %in% installed_solvers(), "ECOS not installed")
  ## ...and it must still be REFUSED by a solver without power cones, rather
  ## than silently downgraded.
  expect_error(psolve(p_exact, solver = "ECOS"), "PowCone3D")
})

test_that("a sqrt-bearing DQCP problem reaches SOC-only solvers", {
  ## The concave-fractional example from the CVXPY docs, bounded away from the
  ## degenerate region.  Before the fix ECOS could not touch it and SCS
  ## returned NaN with x < 0 -- outside sqrt's domain -- while reporting
  ## `optimal_inaccurate`.  A wrong answer with a near-optimal status is worse
  ## than an error, which is why this case is pinned.
  truth <- 0.4288819424   # max sqrt(x)/exp(x) at x = 1/2

  for (slv in c("ECOS", "SCS")) {
    if (!(slv %in% installed_solvers())) next
    x <- Variable()
    p <- Problem(Maximize(sqrt(x) / exp(x)), list(x >= 0.05, x <= 5))
    val <- psolve(p, qcp = TRUE, solver = slv)
    expect_false(is.na(val), label = paste(slv, "objective is not NaN"))
    expect_equal(val, truth, tolerance = 1e-5)
    expect_true(as.numeric(value(x)) > 0,
                label = paste(slv, "iterate is in sqrt's domain"))
  }
})
