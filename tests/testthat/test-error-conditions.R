## @cvxpy cvxpy/error.py:39-67
##
## CVXPY declares seven typed exception classes so callers can catch a specific
## failure mode:
##
##     class SolverError(Exception)     class DPPError(Exception)
##     class DCPError(Exception)        class DGPError(Exception)
##     class DNLPError(Exception)       class DQCPError(Exception)
##     class ParameterError(Exception)
##
## CVXR classed only two of them -- DCPError (binary_operators.R) and SolverError
## (problem.R) -- via `cli_abort(..., class = )`.  The other five aborted as plain
## rlang_errors, so a user could only catch them by matching the message TEXT,
## which is not an API.  Audit finding #3.
##
## Upstream raise sites this pins:
##     DPPError        solving_chain.py:201
##     DGPError        intermediate_chain.py:81
##     DQCPError       problem.py:1084
##     DNLPError       problem.py:1115
##     ParameterError  eval_params.py:15
##
## Each block below is written to FAIL against a build without the `class=`
## argument: `expect_error(class = )` fails both when no error is raised and when
## the error carries the wrong class, so it cannot pass vacuously.

test_that("a non-DPP problem under enforce_dpp raises a catchable DPPError", {
  p1 <- Parameter(); p2 <- Parameter()
  value(p1) <- 2; value(p2) <- 3
  x <- Variable()
  ## Parameters entering non-affinely (p1 * p2) is DCP but not DPP.
  prob <- Problem(Minimize(p1 * p2 * x), list(x >= 1))
  expect_false(is_dpp(prob))                       # guard on the guard
  expect_error(psolve(prob, enforce_dpp = TRUE), class = "DPPError")
})

test_that("gp = TRUE on a non-DGP problem raises a catchable DGPError", {
  y <- Variable()
  prob <- Problem(Minimize(y), list(y >= 1))       # DCP, not DGP
  expect_false(is_dgp(prob))
  expect_error(psolve(prob, gp = TRUE), class = "DGPError")
})

test_that("qcp = TRUE on a non-DQCP problem raises a catchable DQCPError", {
  z <- Variable()
  prob <- Problem(Minimize(-square(z)), list(z >= 1, z <= 2))
  expect_false(is_dcp(prob))                       # both must hold to reach
  expect_false(is_dqcp(prob))                      # the DQCP abort
  expect_error(psolve(prob, qcp = TRUE), class = "DQCPError")
})

test_that("nlp = TRUE on a non-DNLP problem raises a catchable DNLPError", {
  w <- Variable()
  ## Concave AND non-smooth under Minimize: not linearizable either way.
  prob <- Problem(Minimize(-abs(w)), list(w >= -2, w <= 2))
  expect_false(is_dnlp(prob))
  expect_error(psolve(prob, nlp = TRUE), class = "DNLPError")
})

## ParameterError is raised at TWO sites, on the two compilation paths --
## param_prob.R (DPP) and eval_params.R (non-DPP).  Both must carry the class, or
## whether a user can catch it would depend on which path the problem took.
##
## These are two separate test_that blocks ON PURPOSE.  A class-mismatched
## expect_error RETHROWS, which aborts the rest of its block -- so with both
## assertions in one block the second never ran, and the pre-fix run silently
## verified only one of the two sites.  Split, each site gets its own run and its
## own observed pre-fix failure.

test_that("an unspecified parameter raises ParameterError on the DPP path", {
  q <- Parameter(); t1 <- Variable()
  prob <- Problem(Minimize(q * t1), list(t1 >= 1))   # param-affine => DPP
  expect_true(is_dpp(prob))
  expect_error(psolve(prob), class = "ParameterError")
})

test_that("an unspecified parameter raises ParameterError on the non-DPP path", {
  q1 <- Parameter(); q2 <- Parameter(); t2 <- Variable()
  prob <- Problem(Minimize(q1 * q2 * t2), list(t2 >= 1))
  expect_false(is_dpp(prob))
  expect_error(psolve(prob), class = "ParameterError")
})

test_that("the typed conditions remain ordinary errors", {
  ## Adding a class must not break code that catches generically -- the class is
  ## prepended, so `error` and `condition` must still be in the class vector.
  y <- Variable()
  prob <- Problem(Minimize(y), list(y >= 1))
  cond <- tryCatch(psolve(prob, gp = TRUE), error = function(e) e)
  expect_s3_class(cond, "DGPError")
  expect_s3_class(cond, "error")
  expect_s3_class(cond, "condition")
  expect_true(nzchar(conditionMessage(cond)))
})

## NOT TESTED HERE, deliberately: the two conditions CVXR already classed.
## `a * b` for two Variables does NOT raise -- it builds a Multiply, and the
## DCPError abort lives in graph_implementation (binary_operators.R:133,291), so
## it fires during canonicalization, not construction.  Pinning it needs a solve
## and a check of which of the two "not DCP" sites wins, which is a separate
## question from this finding.  SolverError likewise needs an induced solver
## failure.  Left out rather than asserted loosely.
