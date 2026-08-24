## @cvxpy tests/test_conic_solvers.py::TestCPLEX
##
## CPLEX status handling: a DECISION PROCEDURE, not a table.
##
## `cplex_conif.py:147-204` reads two secondary conditions --
## `model.solution.is_primal_feasible()` and `is_dual_feasible()` -- that decide
## the answer for eleven status codes. CVXR used a flat lookup, which cannot
## express "iteration limit, but a feasible point was found", so it answered
## that case wrongly: `psolve(prob, solver = "CPLEX", itlim = 1)` returned
##     Solver "CPLEX" failed.
## and discarded a usable solution. CVXPY returns `optimal_inaccurate` with it.
##
## The flat table was also mislabelled. Checked against
## CPLEX_Studio2211/cplex/include/ilcplex/cpxconst.h, the comments were shifted
## by one through the 10-13 block -- 10 was called "abort_user" but is
## ABORT_IT_LIM, 11 "abort_iteration_limit" but is ABORT_TIME_LIM, 12
## "abort_time_limit" but is ABORT_OBJ_LIM, 13 "abort_dettime_limit" but is
## ABORT_USER -- plus 6 called "feasible" (it is NUM_BEST; FEASIBLE is 23) and
## 115 called "MIP_infeasible_or_unbounded" (it is MIP_OPTIMAL_INFEAS; MIP
## INForUNBD is 119). The mappings had followed the labels.
##
## Most of these tests need no CPLEX licence: the ladder is a pure function of
## (solstat, pfeas, dfeas) and is asserted directly. That is deliberate -- a
## licence-gated test is how the mislabelling survived.

test_that("the status ladder resolves the limit codes on primal feasibility", {
  ## CVXPY group 2 (cplex_conif.py:160-171): pfeas ? OPTIMAL_INACCURATE : SOLVER_ERROR
  for (code in c(10L,   # CPX_STAT_ABORT_IT_LIM
                 11L,   # CPX_STAT_ABORT_TIME_LIM
                 12L,   # CPX_STAT_ABORT_OBJ_LIM
                 13L,   # CPX_STAT_ABORT_USER
                 21L,   # CPX_STAT_ABORT_PRIM_OBJ_LIM
                 22L,   # CPX_STAT_ABORT_DUAL_OBJ_LIM
                 24L,   # CPX_STAT_FIRSTORDER
                 25L,   # CPX_STAT_ABORT_DETTIME_LIM
                 126L)) {  # CPXMIP_ABORT_RELAXED
    expect_equal(CVXR:::.cplex_get_status(code, pfeas = TRUE, dfeas = FALSE),
                 "optimal_inaccurate", label = paste("code", code, "with a primal"))
    expect_equal(CVXR:::.cplex_get_status(code, pfeas = FALSE, dfeas = FALSE),
                 "solver_error", label = paste("code", code, "without a primal"))
  }
})

test_that("the feasible-but-not-proved codes resolve on dual feasibility", {
  ## CVXPY group 3 (cplex_conif.py:172-181): dfeas ? OPTIMAL : OPTIMAL_INACCURATE
  for (code in c(23L,    # CPX_STAT_FEASIBLE
                 104L,   # CPXMIP_SOL_LIM
                 105L,   # CPXMIP_NODE_LIM_FEAS
                 109L,   # CPXMIP_FAIL_FEAS
                 111L,   # CPXMIP_MEM_LIM_FEAS
                 116L,   # CPXMIP_FAIL_FEAS_NO_TREE
                 127L,   # CPXMIP_FEASIBLE -> feasible
                 128L)) {  # CPXMIP_POPULATESOL_LIM
    expect_equal(CVXR:::.cplex_get_status(code, pfeas = TRUE, dfeas = TRUE),
                 "optimal", label = paste("code", code, "dual feasible"))
    expect_equal(CVXR:::.cplex_get_status(code, pfeas = TRUE, dfeas = FALSE),
                 "optimal_inaccurate", label = paste("code", code, "not dual feasible"))
  }
})

test_that("the unconditional groups do not consult feasibility", {
  ## group 1 -- always SOLVER_ERROR (cplex_conif.py:154-159)
  for (code in c(6L, 106L, 110L, 112L, 117L)) {
    for (pf in c(TRUE, FALSE)) {
      expect_equal(CVXR:::.cplex_get_status(code, pfeas = pf, dfeas = pf),
                   "solver_error", label = paste("code", code, "pfeas", pf))
    }
  }
  ## group 4 -- always OPTIMAL (cplex_conif.py:183-188)
  for (code in c(1L, 5L, 102L, 129L, 130L)) {
    expect_equal(CVXR:::.cplex_get_status(code, pfeas = TRUE, dfeas = FALSE),
                 "optimal", label = paste("code", code))
  }
  ## group 5 -- always INFEASIBLE.
  ## Only code 3 can actually reach it: upstream's list also names
  ## optimal_relaxed_{sum,inf,quad} (cplex_conif.py:189-193), but
  ## _handle_solve_status intercepts every feasopt status and raises first
  ## (:110-123), so those three entries are DEAD -- as are the
  ## feasible_relaxed_* entries of the next branch (:196-199). The port keeps
  ## them for line-by-line correspondence; the test asserts what is reachable.
  expect_equal(CVXR:::.cplex_get_status(3L, pfeas = FALSE, dfeas = FALSE),
               "infeasible")
  for (code in c(120L, 121L, 122L, 123L, 124L, 125L)) {
    expect_error(CVXR:::.cplex_get_status(code, pfeas = TRUE, dfeas = TRUE),
                 "feasopt", label = paste("code", code))
  }
  expect_equal(CVXR:::.cplex_get_status(2L, TRUE, TRUE), "unbounded")
  expect_equal(CVXR:::.cplex_get_status(4L, FALSE, FALSE), "infeasible_or_unbounded")
})

test_that("MIP codes normalize onto their LP equivalents", {
  ## CVXPY SOURCE: cplex_conif.py:75-145 (_handle_solve_status).
  expect_equal(CVXR:::.cplex_normalize_solstat(101L), 1L)    # MIP_optimal -> optimal
  expect_equal(CVXR:::.cplex_normalize_solstat(103L), 3L)    # MIP_infeasible -> infeasible
  expect_equal(CVXR:::.cplex_normalize_solstat(107L), 11L)   # TIME_LIM_FEAS -> abort_time_limit
  expect_equal(CVXR:::.cplex_normalize_solstat(108L), 11L)   # TIME_LIM_INFEAS -> abort_time_limit
  expect_equal(CVXR:::.cplex_normalize_solstat(131L), 25L)   # DETTIME_LIM_FEAS -> abort_dettime
  expect_equal(CVXR:::.cplex_normalize_solstat(113L), 13L)   # MIP_ABORT_FEAS -> abort_user
  expect_equal(CVXR:::.cplex_normalize_solstat(115L), 5L)    # MIP_OPTIMAL_INFEAS -> optimal_infeasible
  expect_equal(CVXR:::.cplex_normalize_solstat(119L), 4L)    # MIP INForUNBD -> inf_or_unbd
  expect_equal(CVXR:::.cplex_normalize_solstat(118L), 2L)    # MIP_UNBOUNDED -> unbounded
  expect_equal(CVXR:::.cplex_normalize_solstat(127L), 23L)   # MIP_FEASIBLE -> feasible
  expect_equal(CVXR:::.cplex_normalize_solstat(41L), 6L)     # benders_num_best -> num_best
})

test_that("the codes the old flat table got wrong now match CVXPY", {
  ## 115 is the sharpest: CPXMIP_OPTIMAL_INFEAS means a solution IS available,
  ## with unscaled infeasibilities. The old table called it
  ## "MIP_infeasible_or_unbounded" and returned INFEASIBLE_OR_UNBOUNDED -- an
  ## assertion that the problem might have no solution at all.
  expect_equal(CVXR:::.cplex_get_status(115L, TRUE, FALSE), "optimal")
  ## 119 is the code that actually means MIP infeasible-or-unbounded, and was
  ## absent from the table entirely (so it fell through to solver_error).
  expect_equal(CVXR:::.cplex_get_status(119L, FALSE, FALSE), "infeasible_or_unbounded")
  ## 6 is NUM_BEST -- a solution exists but numerics prevented proving it.
  ## Upstream treats it as an error; the old table called it "feasible".
  expect_equal(CVXR:::.cplex_get_status(6L, TRUE, TRUE), "solver_error")
  ## 108 was INFEASIBLE_INACCURATE, an assertion about the problem; it is a
  ## time limit with no incumbent, so upstream says solver_error.
  expect_equal(CVXR:::.cplex_get_status(108L, pfeas = FALSE, dfeas = FALSE), "solver_error")
  expect_equal(CVXR:::.cplex_get_status(107L, pfeas = TRUE, dfeas = FALSE), "optimal_inaccurate")
})

test_that("a degenerate status never raises", {
  expect_equal(CVXR:::.cplex_get_status(NULL, TRUE, TRUE), "solver_error")
  expect_equal(CVXR:::.cplex_get_status(NA_integer_, TRUE, TRUE), "solver_error")
  expect_equal(CVXR:::.cplex_get_status(integer(0), TRUE, TRUE), "solver_error")
  expect_equal(CVXR:::.cplex_get_status(9999L, TRUE, TRUE), "solver_error")
})

test_that("feasopt and conflict-refiner statuses abort, as upstream asserts", {
  ## CVXR never calls feasopt or the conflict refiner, so one of these means the
  ## model came from somewhere else. cplex_conif.py:121-136 raises AssertionError.
  expect_error(CVXR:::.cplex_normalize_solstat(14L), "feasopt")
  expect_error(CVXR:::.cplex_normalize_solstat(120L), "feasopt")
  expect_error(CVXR:::.cplex_normalize_solstat(31L), "conflict refiner")
})

test_that(".cplex_solution_feasibility reads the Rcplex result", {
  ok <- list(xopt = c(1, 2), extra = list(lambda = c(0.5)))
  f <- CVXR:::.cplex_solution_feasibility(ok)
  expect_true(f$pfeas); expect_true(f$dfeas)
  ## a MIP has no meaningful dual (cplex_conif.py:149)
  f2 <- CVXR:::.cplex_solution_feasibility(ok, is_mip = TRUE)
  expect_true(f2$pfeas); expect_false(f2$dfeas)
  ## no solve happened
  f3 <- CVXR:::.cplex_solution_feasibility(list(xopt = NA, extra = list(lambda = NA)))
  expect_false(f3$pfeas); expect_false(f3$dfeas)
  f4 <- CVXR:::.cplex_solution_feasibility(list(xopt = NULL, extra = list()))
  expect_false(f4$pfeas); expect_false(f4$dfeas)
})

## The end-to-end case, which needs a licence.
test_that("an iteration-limited CPLEX solve returns its solution", {
  skip_if_not_installed("Rcplex")
  skip_if_not("CPLEX" %in% installed_solvers(), "CPLEX not available")
  set.seed(3)
  n <- 60; m <- 40
  A <- matrix(runif(m * n), m, n); bv <- rep(10, m); cv <- -runif(n)
  x <- Variable(n)
  prob <- Problem(Minimize(t(cv) %*% x), list(A %*% x <= bv, x >= 0, x <= 5))
  ## itlim = 1 makes CPLEX stop at CPX_STAT_ABORT_IT_LIM (10) WITH a feasible
  ## primal. Pre-fix this raised 'Solver "CPLEX" failed.'
  val <- suppressWarnings(psolve(prob, solver = "CPLEX", itlim = 1L))
  expect_true(is.numeric(val))
  expect_true(is.finite(as.numeric(val)))
  expect_equal(status(prob), "optimal_inaccurate")
  expect_false(anyNA(as.numeric(value(x))))

  ## and the unlimited solve is unchanged
  prob2 <- Problem(Minimize(t(cv) %*% x), list(A %*% x <= bv, x >= 0, x <= 5))
  expect_true(is.finite(as.numeric(psolve(prob2, solver = "CPLEX"))))
  expect_equal(status(prob2), "optimal")
})
