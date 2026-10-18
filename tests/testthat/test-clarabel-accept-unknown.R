## @cvxpy cvxpy/reductions/solvers/conic_solvers/clarabel_conif.py:87,126-134,180-181
##
## CLARABEL's `accept_unknown` option: when set, and the solver returned BOTH a
## primal and a dual iterate, `InsufficientProgress` is read as
## OPTIMAL_INACCURATE rather than SOLVER_ERROR.
##
## Clarabel reports InsufficientProgress when it stalls, which on a
## near-degenerate problem can still leave a usable iterate. The default is to
## reject it; this option is how a caller asks for the weaker guarantee. It is
## opt-in precisely because "the solver gave up but here is a point" is not
## something to hand back silently.
##
## This was a parity gap dating to the ORIGINAL 1.8.2 port -- ACCEPT_UNKNOWN is
## in every CVXPY version CVXR has ever targeted, so it never showed up in an
## upstream delta, and all three parity passes were delta-driven.

test_that("accept_unknown does NOT reach DQCP bisection subproblems", {
  ## Pinning a limit rather than a capability, because the obvious assumption
  ## is wrong and cost real time: `accept_unknown` cannot change a DQCP solve.
  ##
  ## bisection.py:32-43 `_solve(problem, solver)` calls
  ## `problem.solve(solver=solver)` -- the SOLVER ONLY. Upstream forwards no
  ## solver options into the subproblems, and neither does CVXR's
  ## .bisect_solve(). So a stall inside bisection is SOLVER_ERROR however the
  ## outer psolve() was called. Do not "fix" .bisect_solve to forward options
  ## unless upstream does. (CVXPY goes further and REJECTS the option on the
  ## qcp path: `**kwargs` reach bisect() and `accept_unknown` is a TypeError,
  ## problem.py:1104-1105; R drops unknown kwargs there, problem.R:1073-1077.)
  ##
  ## Through CVXR 1.9.2.9153 this test expected an ERROR from this problem.
  ## That error was an R-only abort in the bisection loop (removed in 9154,
  ## claims audit C0723 / decision 20), not a consequence of the option: the
  ## problem itself solves, in CVXPY 1.9.3 and in CVXR, to optimal_inaccurate
  ## 1.897e-9 along an identical 30-query path (test-bisection-solver-failure.R
  ## runs it). What the option must NOT do is alter that path, so the pin is
  ## the comparison: every sub-solve's query point and status, and the final
  ## value and status, identical with and without the flag.
  trace_bisection <- function(...) {
    real <- get(".bisect_solve", envir = asNamespace("CVXR"))
    log <- new.env()
    log$rows <- list()
    tracer <- function(problem, solver) {
      ps <- parameters(problem)
      tv <- if (length(ps)) as.numeric(value(ps[[1L]])) else NA_real_
      real(problem, solver)
      log$rows[[length(log$rows) + 1L]] <- list(t = tv, status = status(problem))
      invisible(NULL)
    }
    x <- Variable(nonneg = TRUE)
    y <- Variable(nonneg = TRUE)
    prob <- Problem(Maximize(x * y), list(x + y <= 1e-4))
    val <- suppressWarnings(testthat::with_mocked_bindings(
      psolve(prob, solver = "CLARABEL", qcp = TRUE, ...),
      .bisect_solve = tracer, .package = "CVXR"))
    list(value = val, status = status(prob), trace = log$rows)
  }
  plain <- trace_bisection()
  flagged <- trace_bisection(accept_unknown = TRUE)
  expect_gt(length(plain$trace), 1L)
  expect_identical(flagged$trace, plain$trace)
  expect_identical(flagged$status, plain$status)
  expect_identical(flagged$value, plain$value)
  ## The flag did not promote anything: a stalled sub-solve stays solver_error
  ## in the flagged run exactly where it is in the plain one.
  expect_true(any(vapply(plain$trace, function(r) r$status, "") == "solver_error"))
})

test_that("the InsufficientProgress remap fires only when opted in", {
  ## Exercised through reduction_invert() with a synthetic solver result
  ## rather than a stalling problem: no small, portable problem makes CLARABEL
  ## return InsufficientProgress RELIABLY -- the one that does so inside DQCP
  ## bisection solves cleanly on its own, because the stall depends on the
  ## bisection's parameter value. A test that needs a solver to misbehave on
  ## cue is a flaky test; this pins the mapping itself, which is what was
  ## ported.
  ##
  ## status 11 = InsufficientProgress (clarabel::solver_status_descriptions()).
  mk_sol <- function(accept, status = 11L) {
    s <- list(x = c(1, 2), z = 0.5, s = 0.5, status = status,
              solve_time = 0.01, iterations = 7L, obj_val = 3)
    s[["accept_unknown"]] <- accept
    s
  }
  inv <- list(dims = CVXR:::ConeDims(list()), var_id = 1L,
              eq_constr = list(), other_constr = list())
  slv <- CVXR:::Clarabel_Solver()

  ## Default: a stall is an error.
  expect_equal(CVXR:::reduction_invert(slv, mk_sol(FALSE), inv)@status,
               "solver_error")
  ## Opted in, with x and z present: usable but flagged inaccurate.
  expect_equal(CVXR:::reduction_invert(slv, mk_sol(TRUE), inv)@status,
               "optimal_inaccurate")
  ## The remap touches ONE entry -- an infeasibility verdict is not softened
  ## just because the flag is on.
  expect_equal(CVXR:::reduction_invert(slv, mk_sol(TRUE, 3L), inv)@status,
               "infeasible")
})

test_that("accept_unknown is stripped before reaching clarabel_control()", {
  ## The option is CVXPY's, not Clarabel's. If it were passed through it would
  ## land in clarabel_control() as an unknown field -- which R accepts
  ## silently, so the option would be ignored rather than honored. A plain
  ## solve must be unaffected by its presence.
  x <- Variable(nonneg = TRUE)
  p <- Problem(Minimize(sum(x)), list(x >= 3))
  expect_equal(psolve(p, solver = "CLARABEL", accept_unknown = TRUE), 3,
               tolerance = 1e-6)
  expect_equal(psolve(p, solver = "CLARABEL"), 3, tolerance = 1e-6)
})

test_that("accept_unknown does not weaken ordinary statuses", {
  ## The remap touches exactly one entry. An infeasible problem must still be
  ## infeasible, not promoted because the flag happens to be on.
  x <- Variable(nonneg = TRUE)
  p <- Problem(Minimize(sum(x)), list(x >= 3, x <= 1))
  val <- suppressWarnings(
    psolve(p, solver = "CLARABEL", accept_unknown = TRUE))
  expect_equal(status(p), "infeasible")
  expect_true(is.infinite(val))
})
