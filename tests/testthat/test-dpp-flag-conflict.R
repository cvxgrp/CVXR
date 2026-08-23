## @cvxpy cvxpy/problems/problem.py:816-818
##
##     # Invalid DPP setting.
##     # Must be checked here to avoid cache issues.
##     if enforce_dpp and ignore_dpp:
##         raise DPPError("Cannot set enforce_dpp = True and ignore_dpp = True.")
##
## The flags are contradictory: enforce_dpp says "refuse unless this is DPP",
## ignore_dpp says "pretend it is not DPP".  CVXR had no such check, and the
## consequence differed by problem shape -- which is why each shape gets its own
## block below:
##
##   DPP problem        aborted "Problem does not follow DPP rules", a FALSE
##                      statement about a problem that IS DPP
##   non-DPP problem    same wrong message (right outcome, wrong reason)
##   parameter-free     the abort is gated on has_params, so the contradiction
##                      was SILENTLY ACCEPTED and the problem solved
##
## Upstream rejects all three identically; verified against CVXPY 1.9.2.
## Handoff 2026-08-17 section 5.2.

dpp_problem <- function() {           # parameter enters affinely => DPP
  p <- Parameter(); value(p) <- 2
  x <- Variable()
  Problem(Minimize(p * x), list(x >= 1))
}

non_dpp_problem <- function() {       # parameters enter non-affinely
  p1 <- Parameter(); p2 <- Parameter()
  value(p1) <- 2; value(p2) <- 3
  x <- Variable()
  Problem(Minimize(p1 * p2 * x), list(x >= 1))
}

param_free_problem <- function() {
  x <- Variable()
  Problem(Minimize(x), list(x >= 1))
}

## NOTE on what these first two blocks assert.  `class = "DPPError"` ALONE does
## not discriminate: before the guard existed both of these already raised a
## DPPError, just the wrong one ("Problem does not follow DPP rules").  So each
## must also pin the message, or it is a test that cannot fail -- the exact
## defect that let the sqrt bug survive two tests.

test_that("both flags TRUE is rejected on a DPP problem", {
  prob <- dpp_problem()
  expect_true(is_dpp(prob))           # guard on the guard
  expect_error(psolve(prob, solver = "CLARABEL",
                      enforce_dpp = TRUE, ignore_dpp = TRUE),
               class = "DPPError", regexp = "Cannot set")
})

test_that("both flags TRUE is rejected on a non-DPP problem", {
  prob <- non_dpp_problem()
  expect_false(is_dpp(prob))
  expect_error(psolve(prob, solver = "CLARABEL",
                      enforce_dpp = TRUE, ignore_dpp = TRUE),
               class = "DPPError", regexp = "Cannot set")
})

test_that("both flags TRUE is rejected on a parameter-free problem", {
  ## The case that previously SOLVED instead of erroring: the enforce_dpp abort
  ## is gated on has_params, so nothing caught the contradiction.
  prob <- param_free_problem()
  expect_length(parameters(prob), 0L)
  expect_error(psolve(prob, solver = "CLARABEL",
                      enforce_dpp = TRUE, ignore_dpp = TRUE),
               class = "DPPError")
})

test_that("the message names the conflict, not a DPP violation", {
  ## The old message asserted the problem did not follow DPP rules, which was
  ## false for a DPP problem.  Pin the distinction, not just the class.
  cond <- tryCatch(psolve(dpp_problem(), solver = "CLARABEL",
                          enforce_dpp = TRUE, ignore_dpp = TRUE),
                   error = function(e) e)
  expect_s3_class(cond, "DPPError")
  expect_match(conditionMessage(cond), "Cannot set")
  expect_false(grepl("does not follow DPP rules", conditionMessage(cond)))
})

test_that("either flag alone still works", {
  ## The guard must not break the legitimate single-flag uses.
  expect_equal(drop(psolve(dpp_problem(), solver = "CLARABEL",
                           enforce_dpp = TRUE)), 2, tolerance = 1e-6)
  expect_equal(drop(suppressWarnings(
    psolve(dpp_problem(), solver = "CLARABEL", ignore_dpp = TRUE))),
    2, tolerance = 1e-6)
  expect_equal(drop(psolve(param_free_problem(), solver = "CLARABEL",
                           enforce_dpp = TRUE)), 1, tolerance = 1e-6)
})

test_that("enforce_dpp alone still rejects a non-DPP problem, with ITS message", {
  ## The pre-existing behavior this change must not disturb.
  cond <- tryCatch(psolve(non_dpp_problem(), solver = "CLARABEL",
                          enforce_dpp = TRUE),
                   error = function(e) e)
  expect_s3_class(cond, "DPPError")
  expect_match(conditionMessage(cond), "does not follow DPP rules")
})

test_that("a rejected flag pair leaves the cache untouched", {
  ## Upstream's stated reason for checking BEFORE the cache key: "Must be
  ## checked here to avoid cache issues."  A rejected call must not disturb a
  ## chain already compiled, nor leave anything behind for the next call.
  prob <- dpp_problem()
  expect_equal(drop(psolve(prob, solver = "CLARABEL")), 2, tolerance = 1e-6)
  expect_error(psolve(prob, solver = "CLARABEL",
                      enforce_dpp = TRUE, ignore_dpp = TRUE),
               class = "DPPError")
  ## The good solve must still work, and still give the same answer.
  expect_equal(drop(psolve(prob, solver = "CLARABEL")), 2, tolerance = 1e-6)
})
