## @cvxpy cvxpy/reductions/solvers/solving_chain.py:199-204
##
## CVXPY warns when a PARAMETERIZED problem is not DPP, and raises only under
## enforce_dpp:
##
##     if ignore_dpp or not is_dpp:
##         if not ignore_dpp and enforce_dpp: raise DPPError(DPP_ERROR_MSG)
##         if not ignore_dpp:                 warn(DPP_ERROR_MSG)
##
## CVXR raised but never warned, so a user whose re-solves were silently NOT
## being accelerated was never told -- the problem solves, and correctly, just
## without the DPP fast path.  Audit finding #5.
##
## The four states that matter, and why each is a separate block: the warning
## must fire when it should, STAY SILENT in the three cases where it should not
## (a DPP problem, a parameter-free problem, and ignore_dpp), and must not
## double up with the enforce_dpp abort.

## A parameterized problem that is NOT DPP: the parameters enter non-affinely.
non_dpp_problem <- function() {
  p1 <- Parameter(); p2 <- Parameter()
  value(p1) <- 2; value(p2) <- 3
  x <- Variable()
  Problem(Minimize(p1 * p2 * x), list(x >= 1))
}

## Parameter-affine, so DPP.
dpp_problem <- function() {
  p <- Parameter(); value(p) <- 2
  x <- Variable()
  Problem(Minimize(p * x), list(x >= 1))
}

test_that("a non-DPP parameterized problem warns", {
  prob <- non_dpp_problem()
  expect_false(is_dpp(prob))                       # guard on the guard
  expect_warning(psolve(prob, solver = "CLARABEL"), class = "DPPWarning")
})

test_that("the DPP warning names the consequence, not just the violation", {
  ## The point of the message is that RE-SOLVES will not be faster.  A warning
  ## that only said "not DPP" would not tell the user why they should care.
  w <- tryCatch(psolve(non_dpp_problem(), solver = "CLARABEL"),
                warning = function(w) w)
  expect_s3_class(w, "DPPWarning")
  expect_match(conditionMessage(w), "DPP")
  expect_match(conditionMessage(w), "faster")
})

test_that("the warning does not fire when the problem IS DPP", {
  prob <- dpp_problem()
  expect_true(is_dpp(prob))
  expect_no_warning(psolve(prob, solver = "CLARABEL"))
})

test_that("the warning does not fire when there are no parameters", {
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 1))
  expect_length(parameters(prob), 0L)
  expect_no_warning(psolve(prob, solver = "CLARABEL"))
})

test_that("ignore_dpp silences the warning", {
  ## The user has said they know; upstream suppresses both branches under
  ## ignore_dpp, and so does CVXR.
  prob <- non_dpp_problem()
  expect_false(is_dpp(prob))
  expect_no_warning(psolve(prob, solver = "CLARABEL", ignore_dpp = TRUE))
})

test_that("enforce_dpp still aborts, and does not also warn", {
  ## The abort returns before the warning, so exactly one signal is emitted.
  prob <- non_dpp_problem()
  expect_error(psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE),
               class = "DPPError")
  cond <- tryCatch(psolve(non_dpp_problem(), solver = "CLARABEL", enforce_dpp = TRUE),
                   warning = function(w) w, error = function(e) e)
  expect_s3_class(cond, "DPPError")
  expect_false(inherits(cond, "DPPWarning"))
})

test_that("the warning still solves the problem correctly", {
  ## A warning must not change the answer: min p1*p2*x s.t. x >= 1, p1*p2 = 6.
  prob <- non_dpp_problem()
  val <- suppressWarnings(psolve(prob, solver = "CLARABEL"))
  expect_equal(drop(val), 6, tolerance = 1e-6)
})
