## HiGHS Warm-Start Tests
## Tests for HiGHS warm-start support via the persistent-solver API
## (highs::hi_new_solver + hi_solver_set_solution + hi_solver_run).
##
## Covers both solver paths:
##   - HiGHS_Conic_Solver (LP / MILP)
##   - HiGHS_QP_Solver    (QP)
##
## Requires R highs >= 1.14 (persistent-solver API).  Earlier versions
## (1.12 on CRAN) lack hi_solver_set_solution().

skip_if_not_installed("highs", minimum_version = "1.14")

# -- HiGHS_QP_Solver path ----------------------------------------------------

## @cvxpy NONE
test_that("HiGHS QP warm-start: solve twice on same problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - 1)), list(x >= 0))

  r1 <- psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r1, 0.0, tolerance = 1e-4)

  r2 <- psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r2, 0.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS QP warm-start: parameter change between solves", {
  n <- 5L
  x <- Variable(n)
  b <- Parameter(n)

  prob <- Problem(Minimize(sum_squares(x - b)), list(x >= 0))

  ## First solve: b = (1,1,1,1,1)
  value(b) <- rep(1, n)
  psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), rep(1, n), tolerance = 1e-3)

  ## Second solve: b = (2,3,4,5,6)
  value(b) <- c(2, 3, 4, 5, 6)
  psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(2, 3, 4, 5, 6), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("HiGHS QP warm-start: cold and warm produce same solution", {
  x <- Variable(3)
  con <- list(x[1] + x[2] == 1, x[2] + x[3] == 2, x >= 0)

  ## Cold solve
  prob_cold <- Problem(Minimize(sum_squares(x)), con)
  r_cold <- psolve(prob_cold, solver = "HIGHS", warm_start = FALSE)
  x_cold <- as.numeric(value(x))

  ## Warm solve (first solve builds cache, second uses it)
  x2 <- Variable(3)
  con2 <- list(x2[1] + x2[2] == 1, x2[2] + x2[3] == 2, x2 >= 0)
  prob_warm <- Problem(Minimize(sum_squares(x2)), con2)
  psolve(prob_warm, solver = "HIGHS", warm_start = TRUE)
  r_warm <- psolve(prob_warm, solver = "HIGHS", warm_start = TRUE)
  x_warm <- as.numeric(value(x2))

  expect_equal(r_cold, r_warm, tolerance = 1e-4)
  expect_equal(x_cold, x_warm, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS QP warm-start preserves dual values", {
  x <- Variable(2)
  eq_con <- x[1] + x[2] == 1
  ineq_con <- x >= 0
  prob <- Problem(Minimize(sum_squares(x)), list(eq_con, ineq_con))

  ## First solve to warm the cache
  psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-3)

  ## Second solve uses cache; fresh Problem to avoid constraint cache leakage
  x2 <- Variable(2)
  eq_con2 <- x2[1] + x2[2] == 1
  ineq_con2 <- x2 >= 0
  prob2 <- Problem(Minimize(sum_squares(x2)), list(eq_con2, ineq_con2))
  psolve(prob2, solver = "HIGHS", warm_start = TRUE)
  psolve(prob2, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob2), "optimal")
  expect_equal(as.numeric(value(x2)), c(0.5, 0.5), tolerance = 1e-3)
  ## Equality dual is -1.0 for this problem
  expect_equal(as.numeric(dual_value(eq_con2)), -1.0, tolerance = 1e-2)
})

# -- HiGHS_Conic_Solver path (LP) --------------------------------------------

## @cvxpy NONE
test_that("HiGHS LP warm-start: solve twice on same problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  r1 <- psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r1, 2.0, tolerance = 1e-4)

  r2 <- psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r2, 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS LP warm-start: cold and warm produce same solution", {
  set.seed(7L)
  n <- 8L
  c_vec <- rnorm(n)
  x <- Variable(n)
  prob <- Problem(Minimize(t(c_vec) %*% x), list(sum_entries(x) == 1, x >= 0))

  r_cold <- psolve(prob, solver = "HIGHS", warm_start = FALSE)
  x_cold <- as.numeric(value(x))

  ## Fresh Problem so cold and warm don't share constraint dual cache
  x2 <- Variable(n)
  prob_w <- Problem(Minimize(t(c_vec) %*% x2), list(sum_entries(x2) == 1, x2 >= 0))
  psolve(prob_w, solver = "HIGHS", warm_start = TRUE)
  r_warm <- psolve(prob_w, solver = "HIGHS", warm_start = TRUE)

  expect_equal(r_cold, r_warm, tolerance = 1e-4)
  expect_equal(x_cold, as.numeric(value(x2)), tolerance = 1e-4)
})

# -- HiGHS_Conic_Solver path (MILP) ------------------------------------------

## @cvxpy NONE
test_that("HiGHS MILP warm-start: solve twice on same integer problem", {
  x <- Variable(3, integer = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 2, x <= 5))

  r1 <- psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r1, 6.0, tolerance = 1e-6)

  r2 <- psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r2, 6.0, tolerance = 1e-6)
})

# -- Cold after warm ---------------------------------------------------------

## @cvxpy NONE
test_that("HiGHS cold solve after warm solve creates fresh solver", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  r1 <- psolve(prob, solver = "HIGHS", warm_start = TRUE)
  expect_equal(r1, 2.0, tolerance = 1e-4)

  r2 <- psolve(prob, solver = "HIGHS", warm_start = FALSE)
  expect_equal(r2, 2.0, tolerance = 1e-4)
})

# -- Independent caches per Problem ------------------------------------------

## @cvxpy NONE
test_that("HiGHS warm-start: two Problems have independent caches", {
  x1 <- Variable(2)
  prob1 <- Problem(Minimize(sum_entries(x1)), list(x1 >= 1))

  x2 <- Variable(2)
  prob2 <- Problem(Minimize(sum_entries(x2)), list(x2 >= 5))

  r1 <- psolve(prob1, solver = "HIGHS", warm_start = TRUE)
  r2 <- psolve(prob2, solver = "HIGHS", warm_start = TRUE)

  expect_equal(r1, 2.0, tolerance = 1e-4)
  expect_equal(r2, 10.0, tolerance = 1e-4)

  r1b <- psolve(prob1, solver = "HIGHS", warm_start = TRUE)
  expect_equal(r1b, 2.0, tolerance = 1e-4)
})

# -- Dimension change fallback -----------------------------------------------

## @cvxpy NONE
test_that("HiGHS warm-start falls back on dimension change", {
  x <- Variable(2)
  prob1 <- Problem(Minimize(sum_squares(x)), list(x >= 0))
  psolve(prob1, solver = "HIGHS", warm_start = TRUE)

  ## Different size — the warm-path dimension check should skip set_solution
  ## and let the cold solve handle it.
  y <- Variable(4)
  prob2 <- Problem(Minimize(sum_squares(y)), list(y >= 0))
  val <- psolve(prob2, solver = "HIGHS", warm_start = TRUE)
  expect_true(is.finite(val))
  expect_equal(val, 0.0, tolerance = 1e-4)
})
