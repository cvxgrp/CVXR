## @cvxpy tests/test_conic_solvers.py, tests/test_qp_solvers.py
##
## Solver status maps: the table that decides whether a user is told "optimal",
## "infeasible", or "solver_error". A wrong entry here is a wrong ASSERTION
## ABOUT THE PROBLEM, and nothing downstream can catch it.
##
## Most of these solvers are not installed in most environments, so the tests
## that need one are skipped -- but the MAPS themselves are plain R lists and
## can be asserted directly, with no solver present. That is deliberate: the
## defect these pin (Gurobi's SUBOPTIMAL) had gone unnoticed precisely because
## every test touching it required a commercial license.
##
## CVXPY SOURCE: gurobi_conif.py:46-62, gurobi_qpif.py:36-52,
## mosek_conif.py:449-458, highs_conif.py:57-69.

## ---------------------------------------------------------------------------
## Gurobi
## ---------------------------------------------------------------------------

test_that("the Gurobi status map agrees with CVXPY's, key by key", {
  m <- CVXR:::GUROBI_STATUS_MAP
  ## CVXPY gurobi_conif.py:46-62 -- the integer codes with their enum names.
  ## R gurobi returns the NAME, so CVXR keys by name; the values must match.
  expect_equal(m[["OPTIMAL"]],         "optimal")            #  2
  expect_equal(m[["INFEASIBLE"]],      "infeasible")         #  3
  expect_equal(m[["INF_OR_UNBD"]],     "infeasible_or_unbounded")  # 4
  expect_equal(m[["UNBOUNDED"]],       "unbounded")          #  5
  expect_equal(m[["NUMERIC"]],         "solver_error")       # 12
  for (k in c("ITERATION_LIMIT", "NODE_LIMIT", "TIME_LIMIT", "SOLUTION_LIMIT",
              "INTERRUPTED", "INPROGRESS", "USER_OBJ_LIMIT", "WORK_LIMIT",
              "MEM_LIMIT")) {
    expect_equal(m[[k]], "user_limit", label = paste("Gurobi", k))
  }
})

test_that("Gurobi SUBOPTIMAL maps to user_limit, not optimal_inaccurate", {
  ## CVXPY gurobi_conif.py:57 and gurobi_qpif.py:47: `13: s.USER_LIMIT`.
  ##
  ## This was OPTIMAL_INACCURATE -- the ONE key of sixteen that disagreed --
  ## and the disagreement hid itself: reduction_invert demotes
  ## `USER_LIMIT with no primal solution` to INFEASIBLE_INACCURATE
  ## (gurobi_conif.py:301-302), and SUBOPTIMAL could never reach that rule.
  ## So a Gurobi run that stopped suboptimally with NO incumbent was reported
  ## "optimal_inaccurate": a claim that a nearly optimal solution exists.
  expect_equal(CVXR:::GUROBI_STATUS_MAP[["SUBOPTIMAL"]], "user_limit")
})

## ---------------------------------------------------------------------------
## A status lookup must never raise
##
## Upstream indexes its maps directly and would KeyError. CVXR guards with
## `if (is.null(status)) status <- SOLVER_ERROR`, but in R the failure happens
## BEFORE that guard: `map[[key]]` with `key` NULL or character(0) raises
## "subscript out of bounds". Ten sites repeated the idiom; .status_from_map()
## is now the one place it lives.
## ---------------------------------------------------------------------------

test_that(".status_from_map degrades to solver_error instead of raising", {
  m <- list("OPTIMAL" = "optimal", "7" = "optimal")
  expect_equal(CVXR:::.status_from_map(m, "OPTIMAL"), "optimal")
  expect_equal(CVXR:::.status_from_map(m, 7L), "optimal")        # coerced
  expect_equal(CVXR:::.status_from_map(m, "7"), "optimal")
  ## Every degenerate key returns the default rather than erroring.
  expect_equal(CVXR:::.status_from_map(m, NULL), "solver_error")
  expect_equal(CVXR:::.status_from_map(m, character(0)), "solver_error")
  expect_equal(CVXR:::.status_from_map(m, NA), "solver_error")
  expect_equal(CVXR:::.status_from_map(m, NA_character_), "solver_error")
  expect_equal(CVXR:::.status_from_map(m, "NOT_A_STATUS"), "solver_error")
  expect_equal(CVXR:::.status_from_map(m, c("a", "b")), "solver_error")
  ## and an explicit default is honored
  expect_equal(CVXR:::.status_from_map(m, NULL, default = "unbounded"), "unbounded")
})

## ---------------------------------------------------------------------------
## HiGHS: the map is keyed by ORDINAL, upstream keys by enum NAME
##
## This is the one map that would break SILENTLY if the R highs package
## renumbered HighsModelStatus -- the correspondence is positional and was
## asserted only by a trailing comment. Pin it against the package itself, so
## the next highs bump fails here rather than in a user's answer.
## ---------------------------------------------------------------------------

test_that("HiGHS model-status ordinals still mean what the map says", {
  skip_if_not_installed("highs")
  ## `status_message` is the enum name; solve a trivially optimal LP and check
  ## that the ordinal CVXR maps to OPTIMAL is the one HiGHS calls "Optimal".
  res <- highs::highs_solve(L = c(1, 1), lower = c(0, 0), upper = c(10, 10),
                            A = matrix(c(1, 1), 1, 2), lhs = 1, rhs = Inf)
  expect_equal(res$status_message, "Optimal")
  expect_equal(CVXR:::HIGHS_STATUS_MAP[[as.character(res$status)]], "optimal")

  ## Infeasible: lhs > rhs is unsatisfiable.
  inf_res <- highs::highs_solve(L = c(1, 1), lower = c(0, 0), upper = c(1, 1),
                                A = matrix(c(1, 1), 1, 2), lhs = 10, rhs = Inf)
  expect_equal(inf_res$status_message, "Infeasible")
  expect_equal(CVXR:::HIGHS_STATUS_MAP[[as.character(inf_res$status)]], "infeasible")

  ## The full ordinal -> meaning table CVXR relies on.
  expect_equal(CVXR:::HIGHS_STATUS_MAP[["7"]],  "optimal")
  expect_equal(CVXR:::HIGHS_STATUS_MAP[["8"]],  "infeasible")
  expect_equal(CVXR:::HIGHS_STATUS_MAP[["9"]],  "infeasible_or_unbounded")
  expect_equal(CVXR:::HIGHS_STATUS_MAP[["10"]], "unbounded")
  for (k in as.character(11:15)) {
    expect_equal(CVXR:::HIGHS_STATUS_MAP[[k]], "user_limit", label = paste("HiGHS", k))
  }
})

## ---------------------------------------------------------------------------
## MOSEK accept_unknown
##
## CVXPY mosek_conif.py:456-458 re-points solsta.unknown to OPTIMAL_INACCURATE
## when the option is set. CVXR hardcoded UNKNOWN -> SOLVER_ERROR and had no
## `accept_unknown` for MOSEK anywhere, so passing the option did NOTHING --
## silently, since an unrecognized psolve() argument is just forwarded.
## ---------------------------------------------------------------------------

test_that("the MOSEK map still reports UNKNOWN as an error by default", {
  expect_equal(CVXR:::MOSEK_STATUS_MAP[["UNKNOWN"]], "solver_error")
  expect_equal(CVXR:::MOSEK_ACCEPT_UNKNOWN, "accept_unknown")
})

test_that("MOSEK accept_unknown does not disturb an ordinary solve", {
  skip_if_not_installed("Rmosek")
  ## The option must be STRIPPED before the option list reaches Rmosek::mosek():
  ## an unrecognized name there is either rejected or silently dropped, and
  ## silently dropped is the worse of the two.
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  expect_equal(as.numeric(psolve(prob, solver = "MOSEK", accept_unknown = TRUE)),
               2, tolerance = 1e-6)
  expect_equal(status(prob), "optimal")
})
