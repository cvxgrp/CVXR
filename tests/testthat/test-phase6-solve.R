## Phase 6: Solver invocation end-to-end tests
## Tests psolve(), solve() S3 dispatch, solver stats, dual values,
## infeasibility, and CVXPY parity.

## @cvxpy NONE
test_that("Basic LP — Clarabel", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(val, 2, tolerance = 1e-5)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("Basic LP — SCS", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  val <- psolve(prob, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(val, 2, tolerance = 1e-4)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("LP with bounds", {
  x <- Variable(1)
  prob <- Problem(Minimize(2 * x), list(x >= 0, x <= 5))
  val <- psolve(prob, verbose = FALSE)
  expect_equal(val, 0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("LP with equality constraint", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x[1] + x[2] == 5, x >= 0))
  val <- psolve(prob, verbose = FALSE)
  expect_equal(val, 5, tolerance = 1e-5)
  expect_equal(sum(as.numeric(value(x))), 5, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("Maximize problem", {
  x <- Variable(2)
  prob <- Problem(Maximize(sum(x)), list(x <= 3))
  val <- psolve(prob, verbose = FALSE)
  expect_equal(val, 6, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(3, 3), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("SOCP — Clarabel", {
  x <- Variable(3)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))
  val <- psolve(prob, verbose = FALSE)
  expect_equal(val, sqrt(3), tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("SOCP — SCS", {
  x <- Variable(3)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))
  val <- psolve(prob, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(val, sqrt(3), tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Nonneg variable attribute", {
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Minimize(sum(x)), list(x <= 3))
  val <- psolve(prob, verbose = FALSE)
  expect_equal(val, 0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("Infeasible problem", {
  x <- Variable(1)
  prob <- Problem(Minimize(x), list(x >= 1, x <= -1))
  val <- psolve(prob, verbose = FALSE)
  expect_equal(status(prob), INFEASIBLE)
  expect_equal(value(prob), Inf)
})

## @cvxpy NONE
test_that("Dual values populated", {
  x <- Variable(2)
  constr <- list(x >= 1)
  prob <- Problem(Minimize(sum(x)), constr)
  psolve(prob, verbose = FALSE)
  dv <- dual_value(constr[[1]])
  expect_true(!is.null(dv))
  ## For min sum(x), x >= 1, dual should be approximately 1 for each variable
  expect_equal(as.numeric(dv), c(1, 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("SolverStats populated", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  psolve(prob, verbose = FALSE, solver = "CLARABEL")
  ss <- solver_stats(prob)
  expect_true(S7::S7_inherits(ss, SolverStats))
  expect_equal(ss@solver_name, CLARABEL_SOLVER)
  expect_true(!is.null(ss@solve_time))
  expect_true(ss@solve_time >= 0)
})

## @cvxpy NONE
test_that("solve() S3 dispatch works", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  result <- solve(prob, verbose = FALSE)
  expect_s3_class(result, "cvxr_result")
  expect_equal(result$value, 2, tolerance = 1e-5)
  expect_equal(result$status, OPTIMAL)
  expect_equal(status(prob), OPTIMAL)
})

## @cvxpy NONE
test_that("Non-DCP problem rejected", {
  x <- Variable(1)
  prob <- Problem(Maximize(exp(x)))
  expect_error(psolve(prob, verbose = FALSE), "not DCP")
})

## @cvxpy NONE
test_that("ExpCone problem — entropy", {
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Minimize(-sum(entr(x))), list(sum(x) == 1))
  val <- psolve(prob, verbose = FALSE)
  ## max entropy at uniform: -(-0.5*ln(0.5)*2) = ln(2) ≈ 0.6931
  expect_equal(val, -log(2), tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Log problem (ExpCone)", {
  x <- Variable(1)
  prob <- Problem(Maximize(log(x)), list(x <= 2))
  val <- psolve(prob, verbose = FALSE)
  expect_equal(val, log(2), tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 2, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("problem_solution accessor", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  psolve(prob, verbose = FALSE)
  sol <- solution(prob)
  expect_true(S7::S7_inherits(sol, Solution))
  expect_equal(sol@status, OPTIMAL)
})

## @cvxpy NONE
test_that("value(Problem) returns NULL before solve", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  expect_null(value(prob))
})

## @cvxpy NONE
test_that("Multiple constraints", {
  x <- Variable(1)
  y <- Variable(1)
  prob <- Problem(Minimize(x + y), list(x >= 2, y >= 3, x + y <= 10))
  val <- psolve(prob, verbose = FALSE)
  expect_equal(val, 5, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 2, tolerance = 1e-5)
  expect_equal(as.numeric(value(y)), 3, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("Solver options passed through", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  ## Pass max_iter to Clarabel — should still solve with 200 iters
  val <- psolve(prob, verbose = FALSE, solver = "CLARABEL", max_iter = 200L)
  expect_equal(val, 2, tolerance = 1e-5)
})
