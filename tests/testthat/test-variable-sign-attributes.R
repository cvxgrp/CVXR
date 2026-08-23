## @cvxpy tests/test_constraints.py, tests/test_problem.py
##
## `Variable(pos = TRUE)` and `Variable(neg = TRUE)` must constrain the solve.
##
## Why this file exists: they did not. `CvxAttr2Constr` emitted a constraint for
## `nonneg`/`nonpos` only, so the `pos`/`neg` attribute reached `get_bounds()`
## (leaf.R already maps pos to a lower bound of 0) but nothing else: only a
## solver with BOUNDED_VARIABLES = TRUE reads those bounds, and of the conic
## solvers only HiGHS does. On CLARABEL/SCS/ECOS/OSQP/MOSEK the sign was
## silently dropped, so
##     minimize x  s.t.  x <= 5,  x = Variable(pos = TRUE)
## returned -Inf (UNBOUNDED) instead of 0. Measured identical in 1.8.2, 1.8.2.1,
## 1.8.2.9512, 1.9.1 (CRAN) and 1.9.1.9031 -- long-standing, not a regression.
##
## Nothing caught it because the attribute IS honored everywhere except the
## emitted cone program: `is_nonneg()` returns TRUE (leaf.py:287 does the same),
## so every sign, curvature and DCP assertion passed.
##
## CVXPY SOURCE: leaf.py:358-361 -- `nonneg or pos` => term >= 0, and
## `nonpos or neg` => term <= 0; cvx_attr2constr.py:32-52 lists all four sign
## attributes in CONVEX_ATTRIBUTES and BOUND_ATTRIBUTES.

## Tolerance note: 1e-5, not tighter, because SCS runs at eps_abs = eps_rel =
## 1e-5 and lands 4.6665e-06 from zero on these. That is not CVXR being loose --
## CVXPY's SCS returns the SAME 4.6665165683e-06 on the identical problems. A
## 1e-6 bar would fail upstream too. CLARABEL/ECOS reach ~1e-10 here.

.sign_solvers <- function() {
  intersect(c("CLARABEL", "SCS", "ECOS", "OSQP", "HIGHS", "MOSEK"),
            installed_solvers())
}

test_that("Variable(pos = TRUE) constrains the solve on every solver", {
  for (s in .sign_solvers()) {
    x <- Variable(pos = TRUE)
    prob <- Problem(Minimize(x), list(x <= 5))
    val <- psolve(prob, solver = s)
    expect_equal(as.numeric(val), 0, tolerance = 1e-5,
                 label = paste("objective under", s))
    expect_equal(status(prob), "optimal", label = paste("status under", s))
  }
})

test_that("Variable(neg = TRUE) constrains the solve on every solver", {
  for (s in .sign_solvers()) {
    y <- Variable(neg = TRUE)
    prob <- Problem(Maximize(y), list(y >= -5))
    val <- psolve(prob, solver = s)
    expect_equal(as.numeric(val), 0, tolerance = 1e-5,
                 label = paste("objective under", s))
    expect_equal(status(prob), "optimal", label = paste("status under", s))
  }
})

## nonneg/nonpos were always correct; they pin that the widened branches did not
## disturb the path that already worked.
test_that("Variable(nonneg/nonpos = TRUE) still constrain the solve", {
  for (s in .sign_solvers()) {
    z <- Variable(nonneg = TRUE)
    expect_equal(as.numeric(psolve(Problem(Minimize(z), list(z <= 5)), solver = s)),
                 0, tolerance = 1e-5, label = paste("nonneg under", s))
    w <- Variable(nonpos = TRUE)
    expect_equal(as.numeric(psolve(Problem(Maximize(w), list(w >= -5)), solver = s)),
                 0, tolerance = 1e-5, label = paste("nonpos under", s))
  }
})

test_that("pos/neg constrain vector and matrix variables", {
  x <- Variable(3, pos = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x <= 5))
  expect_equal(as.numeric(psolve(prob, solver = "CLARABEL")), 0, tolerance = 1e-6)
  expect_true(all(as.numeric(value(x)) >= -1e-6))

  M <- Variable(c(2, 2), pos = TRUE)
  prob2 <- Problem(Minimize(sum_entries(M)), list(M <= 5))
  expect_equal(as.numeric(psolve(prob2, solver = "CLARABEL")), 0, tolerance = 1e-6)
  expect_true(all(as.numeric(value(M)) >= -1e-6))
})

## The sign must actually bind, not merely be satisfiable: pushing the objective
## the other way must stop at the boundary rather than run away.
test_that("pos/neg bind when the objective pushes against them", {
  x <- Variable(pos = TRUE)
  ## Unconstrained below except by the attribute.
  expect_equal(as.numeric(psolve(Problem(Minimize(x), list(x <= 5)),
                                 solver = "CLARABEL")), 0, tolerance = 1e-6)
  ## An explicitly infeasible demand must report infeasible, not solve.
  x2 <- Variable(pos = TRUE)
  p2 <- Problem(Minimize(0), list(x2 <= -1))
  expect_equal(as.numeric(psolve(p2, solver = "CLARABEL")), Inf)
  expect_equal(status(p2), "infeasible")
})

## The DQCP nested-inverse lowering leans on this: with the sign dropped, the
## bisection's near-boundary subproblems were degenerate enough that Clarabel
## returned InsufficientProgress and the search stalled. See
## test-dqcp-nested-inverse.R.
test_that("pos = TRUE is honored inside a DQCP lowering", {
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(sqrt(inv_pos(x))), list(x <= 4, x >= 0.1))
  expect_true(is_dqcp(prob))
  val <- psolve(prob, qcp = TRUE, solver = "CLARABEL")
  expect_equal(as.numeric(val), 0.5, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 4, tolerance = 1e-4)
})
