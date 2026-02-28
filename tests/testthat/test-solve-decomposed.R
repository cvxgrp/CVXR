## Tests for the decomposed solve API:
## problem_data() -> solve_via_data() -> problem_unpack_results()

skip_if_not_installed("clarabel")

# ── Round-trip: decomposed solve matches normal solve ─────────────────────

## @cvxpy NONE
test_that("decomposed solve round-trip matches psolve()", {
  x <- Variable(5)
  prob <- Problem(Minimize(sum_squares(x - 1:5)), list(x >= 0))

  ## Normal solve
  psolve(prob, solver = "CLARABEL")
  val_normal <- value(prob)
  x_normal <- value(x)

  ## Decomposed solve (fresh problem to avoid cache leakage)
  x2 <- Variable(5)
  prob2 <- Problem(Minimize(sum_squares(x2 - 1:5)), list(x2 >= 0))
  pd <- problem_data(prob2, "CLARABEL")

  raw <- solve_via_data(pd$chain, pd$data, warm_start = FALSE,
                        verbose = FALSE, solver_opts = list())
  problem_unpack_results(prob2, raw, pd$chain, pd$inverse_data)

  expect_equal(status(prob2), "optimal")
  expect_equal(value(prob2), val_normal, tolerance = 1e-6)
  expect_equal(value(x2), x_normal, tolerance = 1e-6)
})

# ── solve_via_data on SolvingChain with problem arg ──────────────────────

## @cvxpy NONE
test_that("solve_via_data(chain, ..., problem=) manages solver_cache", {
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 2))
  pd <- problem_data(prob, "CLARABEL")

  ## Pass problem for solver_cache management
  raw <- solve_via_data(pd$chain, pd$data, problem = prob)
  problem_unpack_results(prob, raw, pd$chain, pd$inverse_data)

  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 6.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(2, 2, 2), tolerance = 1e-4)
})

# ── Multiple solves reusing compiled form ────────────────────────────────

## @cvxpy NONE
test_that("compiled form can be reused across multiple solve_via_data calls", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x >= 3))
  pd <- problem_data(prob, "CLARABEL")

  ## Solve twice with the same compiled data
  for (i in 1:3) {
    raw <- solve_via_data(pd$chain, pd$data)
    problem_unpack_results(prob, raw, pd$chain, pd$inverse_data)
    expect_equal(status(prob), "optimal")
    expect_equal(value(prob), 18.0, tolerance = 1e-4)
  }
})

# ── unpack_results backward-compat alias ─────────────────────────────────

## @cvxpy NONE
test_that("unpack_results() works as alias for problem_unpack_results()", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 5))
  pd <- problem_data(prob, "CLARABEL")

  raw <- solve_via_data(pd$chain, pd$data)
  expect_warning(
    unpack_results(prob, raw, pd$chain, pd$inverse_data),
    "deprecated"
  )
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 10.0, tolerance = 1e-4)
})

# ── LP problem ───────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("decomposed solve works for LP", {
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1, x <= 10))
  pd <- problem_data(prob, "CLARABEL")

  raw <- solve_via_data(pd$chain, pd$data)
  problem_unpack_results(prob, raw, pd$chain, pd$inverse_data)

  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 3.0, tolerance = 1e-4)
})

# ── Solver options forwarding ────────────────────────────────────────────

## @cvxpy NONE
test_that("solve_via_data forwards solver_opts", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))), list(x >= 0))
  pd <- problem_data(prob, "CLARABEL")

  ## Pass verbose as solver option (should not error)
  raw <- solve_via_data(pd$chain, pd$data, verbose = FALSE,
                        solver_opts = list(max_iter = 1000L))
  problem_unpack_results(prob, raw, pd$chain, pd$inverse_data)

  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-4)
})
