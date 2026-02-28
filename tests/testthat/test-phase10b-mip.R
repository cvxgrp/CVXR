## Phase 10b: MIP Tests (Boolean/Integer Variables via HiGHS)
## Tests for mixed-integer programming support: boolean LP, integer LP,
## MIQP rejection, solver compatibility validation, auto-selection.

skip_if_not_installed("highs")

# ── Boolean LP ───────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Boolean LP: knapsack problem", {
  ## max 3*x1 + 2*x2 + 5*x3 s.t. 2*x1 + x2 + 3*x3 <= 5, xi ∈ {0,1}
  x <- Variable(3, boolean = TRUE)
  prob <- Problem(Maximize(3 * x[1] + 2 * x[2] + 5 * x[3]),
                  list(2 * x[1] + x[2] + 3 * x[3] <= 5))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Optimal: x1=1, x2=0, x3=1 → value = 3+0+5 = 8
  expect_equal(result, 8.0, tolerance = 1e-5)
  xval <- round(as.numeric(value(x)))
  expect_equal(xval, c(1, 0, 1))
})

## @cvxpy NONE
test_that("Boolean LP: simple 0-1 selection", {
  x <- Variable(2, boolean = TRUE)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] >= 1))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 1.0, tolerance = 1e-5)
  ## One of x1,x2 should be 1, the other 0
  xval <- round(as.numeric(value(x)))
  expect_equal(sum(xval), 1)
})

# ── Integer LP ───────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Integer LP: simple integer variable", {
  ## min x s.t. x >= 1.5, x integer → x = 2
  x <- Variable(1, integer = TRUE)
  prob <- Problem(Minimize(x),
                  list(x >= 1.5))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 2.0, tolerance = 1e-5)
  expect_equal(round(as.numeric(value(x))), 2)
})

## @cvxpy NONE
test_that("Integer LP: assignment-type problem", {
  ## min x1 + x2 s.t. x1 + x2 >= 5, 0 <= xi, xi integer
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] >= 5, x >= 0))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 5.0, tolerance = 1e-5)
})

# ── Mixed continuous + boolean ───────────────────────────────────────────────

## @cvxpy NONE
test_that("Mixed continuous and boolean variables", {
  ## min y + z  s.t. y >= 2*z, z ∈ {0,1}, y >= 0.5
  y <- Variable(1)
  z <- Variable(1, boolean = TRUE)
  prob <- Problem(Minimize(y + z),
                  list(y >= 2 * z, y >= 0.5))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## z=0 → y=0.5, obj=0.5; z=1 → y=2, obj=3. Optimal: z=0, y=0.5
  expect_equal(result, 0.5, tolerance = 1e-5)
  expect_equal(round(as.numeric(value(z))), 0)
  expect_equal(as.numeric(value(y)), 0.5, tolerance = 1e-5)
})

# ── Mixed continuous + integer ───────────────────────────────────────────────

## @cvxpy NONE
test_that("Mixed continuous and integer variables", {
  ## min y + z  s.t. y + z >= 3.5, z integer, z >= 0, y >= 0
  y <- Variable(1)
  z <- Variable(1, integer = TRUE)
  prob <- Problem(Minimize(y + z),
                  list(y + z >= 3.5, z >= 0, y >= 0))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## z=4 (cheapest integer >= 3.5) or z=3, y=0.5 → obj=3.5. Optimal: z=3, y=0.5 or z=4, y=0
  ## Both give obj = 3.5
  expect_equal(result, 3.5, tolerance = 1e-5)
})

# ── MIP infeasible ───────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("MIP infeasible: boolean with contradictory constraints", {
  x <- Variable(1, boolean = TRUE)
  prob <- Problem(Minimize(x),
                  list(x >= 2))
  ## x ∈ {0,1} but x >= 2 → infeasible
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_true(status(prob) %in%
    c("infeasible", "infeasible_inaccurate", "infeasible_or_unbounded"))
})

# ── No dual variables for MIP ────────────────────────────────────────────────

## @cvxpy NONE
test_that("MIP problems do not produce dual variables", {
  x <- Variable(2, boolean = TRUE)
  con <- x[1] + x[2] >= 1
  prob <- Problem(Minimize(x[1] + x[2]), list(con))
  psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Dual values should be NULL or empty for MIP
  dv <- dual_value(con)
  expect_true(is.null(dv) || all(dv == 0))
})

# ── MIQP rejection ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("MIQP is rejected before solver call", {
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(sum_squares(x)),
                  list(x >= 1))
  expect_error(psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE),
               "MIQP")
})

## @cvxpy NONE
test_that("Boolean QP is rejected (MIQP)", {
  x <- Variable(2, boolean = TRUE)
  prob <- Problem(Minimize(sum_squares(x - c(0.5, 0.5))),
                  list(x[1] + x[2] >= 1))
  expect_error(psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE),
               "MIQP")
})

# ── Solver compatibility ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("MIP with SCS solver is rejected", {
  x <- Variable(2, boolean = TRUE)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] >= 1))
  expect_error(psolve(prob, solver = SCS_SOLVER, verbose = FALSE),
               "HIGHS")
})

## @cvxpy NONE
test_that("MIP with CLARABEL solver is rejected", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, boolean = TRUE)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] >= 1))
  expect_error(psolve(prob, solver = CLARABEL_SOLVER, verbose = FALSE),
               "HIGHS")
})

## @cvxpy NONE
test_that("MIP with OSQP solver is rejected", {
  skip_if_not_installed("osqp")
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x >= 1))
  expect_error(psolve(prob, solver = OSQP_SOLVER, verbose = FALSE),
               "HIGHS")
})

# ── Auto solver selection ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Auto-selection picks HiGHS for MIP", {
  x <- Variable(2, boolean = TRUE)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] >= 1))
  ## psolve without explicit solver → should auto-select HiGHS for MIP
  result <- psolve(prob, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 1.0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("Auto-selection picks Clarabel for non-MIP LP", {
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x >= 0, x[1] + x[2] >= 1))
  result <- psolve(prob, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 1.0, tolerance = 1e-5)
})

# ── is_mixed_integer query ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("is_mixed_integer detects boolean variables", {
  x <- Variable(2, boolean = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  expect_true(is_mixed_integer(prob))
})

## @cvxpy NONE
test_that("is_mixed_integer detects integer variables", {
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  expect_true(is_mixed_integer(prob))
})

## @cvxpy NONE
test_that("is_mixed_integer returns FALSE for continuous problems", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  expect_false(is_mixed_integer(prob))
})

## @cvxpy NONE
test_that("is_mixed_integer detects mixed continuous + boolean", {
  x <- Variable(2)
  y <- Variable(1, boolean = TRUE)
  prob <- Problem(Minimize(sum_entries(x) + y), list(x >= 0))
  expect_true(is_mixed_integer(prob))
})

# ── Larger MIP problems ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Boolean LP with multiple constraints", {
  ## Facility location-like problem
  n <- 5
  x <- Variable(n, boolean = TRUE)
  costs <- c(3, 7, 2, 5, 1)
  prob <- Problem(Minimize(sum(costs * x)),
                  list(sum_entries(x) >= 3))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  xval <- round(as.numeric(value(x)))
  ## Should pick the 3 cheapest: indices 5(1), 3(2), 1(3) → cost=1+2+3=6
  expect_equal(sum(xval), 3)
  expect_equal(result, 6.0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("Integer LP with bounds", {
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(-x[1] - 2 * x[2]),
                  list(x >= 0, x <= 10,
                       x[1] + x[2] <= 7,
                       2 * x[1] + x[2] <= 11))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## LP relaxation: x1=4, x2=3, obj=-10. Integer feasible too.
  xval <- round(as.numeric(value(x)))
  expect_true(all(xval >= 0))
  expect_true(all(xval <= 10))
  expect_true(xval[1] + xval[2] <= 7)
})
