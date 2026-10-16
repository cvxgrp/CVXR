## Bisection behavior when the solver FAILS at a query point -- neither
## reporting feasibility nor proving infeasibility.
##
## CVXPY SOURCE: bisection.py:103-157 (CVXPY 1.9.3).
##
## Upstream through 1.9.1 re-solved the identical failing subproblem until it
## raised "Max iters hit during bisection".  CVXR never had that infinite
## retry, but avoided it unsoundly: it treated the failure as infeasibility
## and moved `low`, the reported lower bound, on a point that had never been
## proven infeasible.  1.9.2 keeps `low` fixed and perturbs only the point the
## next query is drawn from.

## @cvxpy test_dqcp.py::TestDqcp::test_multiply_solver_failure_during_bisection
test_that("bisection survives a solver failure near the optimum", {
  ## Upstream's regression problem for #3407.  Traced sub-solve by sub-solve
  ## on CVXPY 1.9.3 and CVXR 1.9.2.9154 (claims audit decision 20): the
  ## feasibility check is optimal, t = -1 infeasible, 26 halvings down to
  ## t = -1.49e-8 all infeasible, CLARABEL FAILS at t = -7.45e-9 (the failure
  ## this test is named for), the perturbed query t = -3.73e-9 is
  ## optimal_inaccurate, and both libraries return 1.89734745333e-09 with the
  ## same x and y.  Through 9171 R made one extra sub-solve, an R-only seed
  ## solve at t = 0, which CLARABEL also fails (removed in 9172, decision 20
  ## side finding 3).  This test was skipped through 9153: the loop's R-only
  ## `is.null(soln)` abort fired on the t = -7.45e-9 failure.  Removed in 9154
  ## (C0723).
  ##
  ## Tolerance: upstream asserts `assertAlmostEqual(value, 2.5e-9, places=6)`,
  ## an ABSOLUTE test (|diff| < 5e-7); testthat's `tolerance` is relative, and
  ## 1.897e-9 vs 2.5e-9 is 24% off, so the comparison is written out as the
  ## absolute one.  Badly scaled on purpose: the optimum (5e-5)^2 = 2.5e-9 is
  ## below the bisection eps of 1e-6, so the value is 0 to that tolerance.
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  prob <- Problem(Maximize(x * y), list(x + y <= 1e-4))
  expect_warning(
    val <- psolve(prob, solver = "CLARABEL", qcp = TRUE),
    "Solver failed at iteration")                 # pytest.warns(RuntimeWarning)
  expect_true(status(prob) %in% SOLUTION_PRESENT)
  expect_lt(abs(val - 2.5e-9), 5e-7)
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

## A failure BEFORE any feasible query (claims audit C0723, decision 20, CVXR
## 1.9.2.9154).  Upstream's failure branch (bisection.py:142-156) never looks
## at whether a solution exists yet; through 9153 CVXR aborted with "Solver
## failed with status ..." when `soln` was NULL.  The tests below pass
## explicit `low` AND `high` so the interval search is skipped and the
## injected call numbers land on loop queries.
##
## No real solver fails deterministically at a chosen query point (see the
## skip above), so the failure is injected: `.bisect_solve` is replaced for the
## duration of one psolve() by a wrapper that reports SOLVER_ERROR on the
## listed call numbers exactly as the real function's error handler does, and
## solves normally otherwise.  Call 1 inside bisect() is the feasibility check,
## so call 2 is the first query of .bisect_loop.  The same intervention on
## CVXPY 1.9.3 (monkeypatched `bisection._solve`) gives: first query fails ->
## value 12, one RuntimeWarning at iteration 0, 5 solver calls; every query
## fails -> SolverError "Solver failed at every query point ..." after 55
## calls.
.failing_bisect_solve <- function(fail_calls) {
  real <- get(".bisect_solve", envir = asNamespace("CVXR"))
  counter <- new.env()
  counter$n <- 0L
  fake <- function(problem, solver) {
    counter$n <- counter$n + 1L
    if (counter$n %in% fail_calls) {
      problem@.cache$status <- CVXR:::SOLVER_ERROR
      problem@.cache$solution <- CVXR:::failure_solution(CVXR:::SOLVER_ERROR)
      return(invisible(NULL))
    }
    real(problem, solver)
  }
  list(fake = fake, counter = counter)
}

## @cvxpy NONE (control-flow check with an injected failure; upstream exercises
## bisection.py:142-156 only through test_multiply_solver_failure_during_bisection,
## skipped above because the real-solver scenario is not reproducible here)
test_that("a solver failure on the first bisection query perturbs instead of aborting", {
  x <- Variable()
  prob <- Problem(Minimize(ceil_expr(x)), list(x >= 11.9, x <= 17))
  inj <- .failing_bisect_solve(2L)
  expect_warning(
    val <- testthat::with_mocked_bindings(
      psolve(prob, solver = "CLARABEL", qcp = TRUE, low = 12, high = 17),
      .bisect_solve = inj$fake, .package = "CVXR"),
    "Solver failed at iteration 0, querying a point closer")
  expect_equal(val, 12)
  expect_true(status(prob) %in% SOLUTION_PRESENT)
  expect_equal(inj$counter$n, 5L)
})

## @cvxpy NONE (terminal case of the same branch: bisection.py:143-146)
test_that("a solver that fails at every query point stops with the upstream message", {
  a <- Variable(nonneg = TRUE)
  b <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize((a + 1) / (b + 1)), list(a >= 1, b <= 3, a <= 2))
  inj <- .failing_bisect_solve(2:1000)
  expect_error(
    suppressWarnings(testthat::with_mocked_bindings(
      psolve(prob, solver = "CLARABEL", qcp = TRUE, low = 0, high = 2),
      .bisect_solve = inj$fake, .package = "CVXR")),
    "Solver failed at every query point")
  ## The query point halves its distance to `high` on every failure and
  ## reaches it in finite steps, well inside max_iters = 100.
  expect_lt(inj$counter$n, 100L)
})

## Claims audit C0723 (CVXR 1.9.2.9167): one warning per failed query, in
## upstream's words (bisection.py:147-151), and SolverError as the class of both
## loop aborts (:144-146, :157). Measured on CVXPY 1.9.3 with its default warning
## filters: queries 2-4 failing give three RuntimeWarnings, "iteration 0", "1",
## "2", and the same three again when the problem is solved a second time in the
## same process, because the iteration number makes each message distinct.
## Through 9166 CVXR warned once per bisection, on the premise that Python shows
## one; these tests were seen failing on 9166 first.
.collect_bisect_warnings <- function(expr) {
  msgs <- character()
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, gsub("\\s+", " ", conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = msgs)
}

## @cvxpy NONE (injected failures; see the first-query test above)
test_that("every failed bisection query warns, in upstream's words", {
  x <- Variable()
  prob <- Problem(Minimize(ceil_expr(x)), list(x >= 11.9, x <= 17))
  inj <- .failing_bisect_solve(2:4)
  res <- .collect_bisect_warnings(testthat::with_mocked_bindings(
    psolve(prob, solver = "CLARABEL", qcp = TRUE, low = 12, high = 17),
    .bisect_solve = inj$fake, .package = "CVXR"))
  expect_equal(res$value, 12)
  expect_equal(inj$counter$n, 8L)
  expect_identical(res$warnings, sprintf(
    "Solver failed at iteration %d, querying a point closer to the feasible end of the interval.",
    0:2))
})

## @cvxpy NONE (bisection.py:143-146 and :157 raise error.SolverError)
test_that("the bisection loop aborts with class SolverError", {
  a <- Variable(nonneg = TRUE)
  b <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize((a + 1) / (b + 1)), list(a >= 1, b <= 3, a <= 2))
  inj <- .failing_bisect_solve(2:1000)
  expect_error(
    suppressWarnings(testthat::with_mocked_bindings(
      psolve(prob, solver = "CLARABEL", qcp = TRUE, low = 0, high = 2),
      .bisect_solve = inj$fake, .package = "CVXR")),
    "Solver failed at every query point", class = "SolverError")

  a <- Variable(nonneg = TRUE)
  b <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize((a + 1) / (b + 1)), list(a >= 1, b <= 3, a <= 2))
  expect_error(
    psolve(prob, solver = "CLARABEL", qcp = TRUE, low = 0, high = 2, max_iters = 3),
    "^Max iters hit during bisection\\.$", class = "SolverError")
})

## Claims audit C0723 side findings 1-3 (CVXR 1.9.2.9171): what a sub-solve
## lets through (bisection.py:32-43). Upstream ignores only warnings matching
## `.*DPP.*`, warns "Solver failed, attempting to continue bisection." on a
## SolverError and lets any other exception propagate. Through 9170 CVXR
## suppressed every warning and swallowed every error silently. These tests
## inject at psolve, BELOW .bisect_solve, so the real .bisect_solve runs; the
## same injection on CVXPY 1.9.3 (monkeypatched Problem.solve) is repro 31
## (`notes/audit/repro/31_bisection_subsolve.{R,py}`), whose numbers they assert.
.failing_psolve <- function(modes) {
  real <- CVXR::psolve
  counter <- new.env()
  counter$n <- 0L
  fake <- function(problem, ...) {
    counter$n <- counter$n + 1L
    mode <- unname(modes[as.character(counter$n)])
    if (identical(mode, "solvererror"))
      cli::cli_abort("injected SolverError", class = "SolverError")
    if (identical(mode, "othererror"))
      cli::cli_abort("injected non-solver error")
    if (identical(mode, "warn_inaccurate"))
      cli::cli_warn("Solution may be inaccurate. (injected)")
    if (identical(mode, "warn_dpp"))
      cli::cli_warn("injected: problem is not DPP", class = "DPPWarning")
    real(problem, ...)
  }
  list(fake = fake, counter = counter, real = real)
}

.ceil_prob <- function() {
  x <- Variable()
  Problem(Minimize(ceil_expr(x)), list(x >= 11.5))
}

## @cvxpy NONE (injected SolverError in a sub-solve: bisection.py:38-43)
test_that("a sub-solve SolverError warns that bisection continues, as upstream does", {
  prob <- .ceil_prob()
  inj <- .failing_psolve(c("2" = "solvererror", "3" = "solvererror", "4" = "solvererror"))
  res <- .collect_bisect_warnings(testthat::with_mocked_bindings(
    inj$real(prob, solver = "CLARABEL", qcp = TRUE, low = 12, high = 17),
    psolve = inj$fake, .package = "CVXR"))
  expect_equal(res$value, 12)
  expect_equal(inj$counter$n, 8L)
  ## CVXPY 1.9.3: six RuntimeWarnings, alternating; the repeated message is not
  ## collapsed, because _solve's catch_warnings resets the once-registry.
  expect_identical(res$warnings, as.vector(rbind(
    "Solver failed, attempting to continue bisection.",
    sprintf("Solver failed at iteration %d, querying a point closer to the feasible end of the interval.",
            0:2))))
})

## @cvxpy NONE (bisection.py:36-38 catches error.SolverError only)
test_that("a non-solver error in a sub-solve propagates out of the bisection", {
  prob <- .ceil_prob()
  inj <- .failing_psolve(c("2" = "othererror"))
  expect_error(
    testthat::with_mocked_bindings(
      inj$real(prob, solver = "CLARABEL", qcp = TRUE, low = 12, high = 17),
      psolve = inj$fake, .package = "CVXR"),
    "injected non-solver error")
  expect_equal(inj$counter$n, 2L)
})

## @cvxpy NONE (bisection.py:33-35 filters only messages matching .*DPP.*)
test_that("sub-solve warnings reach the user except the DPP one", {
  prob <- .ceil_prob()
  inj <- .failing_psolve(c("2" = "warn_inaccurate", "3" = "warn_dpp"))
  res <- .collect_bisect_warnings(testthat::with_mocked_bindings(
    inj$real(prob, solver = "CLARABEL", qcp = TRUE, low = 12, high = 17),
    psolve = inj$fake, .package = "CVXR"))
  expect_equal(res$value, 12)
  expect_identical(res$warnings, "Solution may be inaccurate. (injected)")
})

## @cvxpy NONE (bisection.py:78-79, :92-93 and :98-100 raise error.SolverError)
test_that("the interval search aborts with class SolverError, in upstream's words", {
  ## Call 2 is the solve at `high` (:70, raising at :78-79), call 3 the solve at
  ## `low` (:84, raising at :92-93).
  for (call in c("2", "3")) {
    prob <- .ceil_prob()
    inj <- .failing_psolve(stats::setNames("solvererror", call))
    err <- tryCatch(
      .collect_bisect_warnings(testthat::with_mocked_bindings(
        inj$real(prob, solver = "CLARABEL", qcp = TRUE),
        psolve = inj$fake, .package = "CVXR")),
      error = function(e) e)
    expect_s3_class(err, "SolverError")
    expect_identical(conditionMessage(err), "Solver failed with status solver_error")
    expect_equal(inj$counter$n, as.integer(call))
  }

  x <- Variable()
  prob <- Problem(Minimize(ceil_expr(x)))
  err <- tryCatch(psolve(prob, solver = "CLARABEL", qcp = TRUE), error = function(e) e)
  expect_s3_class(err, "SolverError")
  expect_identical(gsub("\\s+", " ", conditionMessage(err)),
                   "Unable to find suitable interval for bisection; your problem may be unbounded..")
})

## No solve at `high` between the interval search and the loop (claims audit
## decision 20 side finding 3, CVXR 1.9.2.9172).  Upstream's `_bisect` starts
## from `soln = None` (bisection.py:108) and `bisect` calls it straight after
## the interval search (:203-214), so a solution comes only from a feasible
## query.  When the optimum sits exactly at `high`, every midpoint is
## infeasible until the query rounds to `high` itself.  Sub-solve counts and
## outcomes measured on CVXPY 1.9.3 + CLARABEL (repro 32): 57, 58 and 26
## calls, and "Max iters hit" at max_iters = 30.  Through 9171 CVXR's R-only
## seed solve at `high` gave 25, 25 and 27 calls and an optimal answer at
## max_iters = 30.
## @cvxpy NONE (bisection.py:108 and :203-214; sub-solve counts are not asserted upstream)
test_that("bisection makes no seed solve at the upper bound", {
  count <- function(prob, ...) {
    inj <- .failing_bisect_solve(integer(0))
    val <- testthat::with_mocked_bindings(
      psolve(prob, solver = "CLARABEL", qcp = TRUE, ...),
      .bisect_solve = inj$fake, .package = "CVXR")
    list(val = val, n = inj$counter$n)
  }

  x <- Variable()
  r <- count(Problem(Minimize(sign(x)), list(-2 <= x, x <= -0.5)))
  expect_equal(r$val, -1)
  expect_equal(r$n, 57L)

  x <- Variable()
  r <- count(Problem(Minimize(sign(x)), list(0.5 <= x, x <= 2)))
  expect_equal(r$val, 1)
  expect_equal(r$n, 58L)

  x <- Variable()
  r <- count(Problem(Maximize(floor(x)), list(x <= 2.5)))
  expect_equal(r$val, 2)
  expect_equal(r$n, 26L)

  x <- Variable()
  prob <- Problem(Minimize(sign(x)), list(-2 <= x, x <= -0.5))
  err <- tryCatch(psolve(prob, solver = "CLARABEL", qcp = TRUE, max_iters = 30L),
                  error = function(e) e)
  expect_s3_class(err, "SolverError")
  expect_identical(conditionMessage(err), "Max iters hit during bisection.")
})
