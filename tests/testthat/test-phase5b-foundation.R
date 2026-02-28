# Tests for Phase 5b: Solution, Chain, reduction utilities

library(CVXR)

# ═══════════════════════════════════════════════════════════════════
# Solution class
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Solution creation with all fields", {
  sol <- Solution(
    status = "optimal", opt_val = 42.0,
    primal_vars = list("1" = matrix(1:3)),
    dual_vars = list("2" = c(0.5, 0.5)),
    attr = list(solve_time = 0.01)
  )
  expect_equal(sol@status, "optimal")
  expect_equal(sol@opt_val, 42.0)
  expect_equal(length(sol@primal_vars), 1L)
  expect_equal(length(sol@dual_vars), 1L)
  expect_equal(sol@attr$solve_time, 0.01)
})

## @cvxpy NONE
test_that("Solution defaults", {
  sol <- Solution(status = "optimal")
  expect_equal(sol@opt_val, NA_real_)
  expect_equal(length(sol@primal_vars), 0L)
  expect_equal(length(sol@dual_vars), 0L)
})

## @cvxpy NONE
test_that("failure_solution for infeasible", {
  sol <- failure_solution("infeasible")
  expect_equal(sol@status, "infeasible")
  expect_equal(sol@opt_val, Inf)
  expect_equal(length(sol@primal_vars), 0L)
})

## @cvxpy NONE
test_that("failure_solution for unbounded", {
  sol <- failure_solution("unbounded")
  expect_equal(sol@status, "unbounded")
  expect_equal(sol@opt_val, -Inf)
})

## @cvxpy NONE
test_that("failure_solution for solver_error", {
  sol <- failure_solution("solver_error")
  expect_equal(sol@status, "solver_error")
  expect_true(is.na(sol@opt_val))
})

# ═══════════════════════════════════════════════════════════════════
# Reduction utilities
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("group_constraints groups by class name", {
  x <- Variable(3)
  cons <- list(
    x >= 0,    # Inequality -> NonNeg from lower_ineq_to_nonneg
    x == 1     # Equality
  )
  ## Direct instances
  z <- Zero(x - Constant(1))
  nn <- NonNeg(x)
  result <- group_constraints(list(z, nn))
  expect_equal(length(result[["Zero"]]), 1L)
  expect_equal(length(result[["NonNeg"]]), 1L)
  expect_equal(length(result[["SOC"]]), 0L)
})

## @cvxpy NONE
test_that("group_constraints ensures all cone types present", {
  result <- group_constraints(list())
  expected_types <- c("Zero", "NonNeg", "SOC", "PSD", "ExpCone",
                      "PowCone3D", "PowConeND")
  for (ct in expected_types) {
    expect_true(ct %in% names(result))
    expect_equal(length(result[[ct]]), 0L)
  }
})

## @cvxpy NONE
test_that("lower_equality converts Equality to Zero", {
  x <- Variable(2)
  eq <- (x == Constant(c(1, 2)))  # Equality
  z <- lower_equality(eq)
  expect_true(inherits(z, "CVXR::Zero"))
  expect_equal(z@id, eq@id)
})

## @cvxpy NONE
test_that("lower_ineq_to_nonneg converts Inequality to NonNeg", {
  x <- Variable(2)
  ineq <- (x <= Constant(c(3, 4)))  # Inequality: x <= c(3,4)
  nn <- lower_ineq_to_nonneg(ineq)
  expect_true(inherits(nn, "CVXR::NonNeg"))
  expect_equal(nn@id, ineq@id)
})

## @cvxpy NONE
test_that("nonpos2nonneg converts NonPos to NonNeg", {
  x <- Variable(2)
  np <- NonPos(x)
  nn <- nonpos2nonneg(np)
  expect_true(inherits(nn, "CVXR::NonNeg"))
  expect_equal(nn@id, np@id)
})

## @cvxpy NONE
test_that("are_args_affine works correctly", {
  x <- Variable(2)
  z <- Zero(x - Constant(c(1, 2)))
  expect_true(are_args_affine(list(z)))

  ## Non-affine arg: exp(x) is not affine
  nn <- NonNeg(exp(x))
  expect_false(are_args_affine(list(nn)))
})

## @cvxpy NONE
test_that("convex_attributes detects variable attributes", {
  x <- Variable(3, nonneg = TRUE)
  expect_true("nonneg" %in% convex_attributes(list(x)))

  y <- Variable(3)
  expect_equal(length(convex_attributes(list(y))), 0L)
})

# ═══════════════════════════════════════════════════════════════════
# Chain class
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Chain creation with empty reductions", {
  chain <- Chain()
  expect_equal(length(chain@reductions), 0L)
})

## @cvxpy NONE
test_that("Chain with FlipObjective", {
  x <- Variable(2)
  p <- Problem(Maximize(sum_entries(x)), list(x <= 1))
  chain <- Chain(reductions = list(FlipObjective()))
  result <- reduction_apply(chain, p)
  new_p <- result[[1L]]
  expect_true(inherits(new_p@objective, "CVXR::Minimize"))
})

## @cvxpy NONE
test_that("Chain apply and invert roundtrip", {
  x <- Variable(2)
  p <- Problem(Maximize(sum_entries(x)), list(x <= 1))
  chain <- Chain(reductions = list(FlipObjective()))

  result <- reduction_apply(chain, p)
  new_p <- result[[1L]]
  inv_data <- result[[2L]]

  ## Create a mock solution
  sol <- Solution(status = "optimal", opt_val = -5.0,
                  primal_vars = list(), dual_vars = list())
  recovered <- reduction_invert(chain, sol, inv_data)
  expect_equal(recovered@opt_val, 5.0)  # negated back
})

# ═══════════════════════════════════════════════════════════════════
# Settings constants
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Phase 5b settings constants exist", {
  expect_equal(SOLUTION_PRESENT, c("optimal", "optimal_inaccurate", "user_limit"))
  expect_true("solver_error" %in% ERROR_STATUS)
  expect_equal(SD_C, "c")
  expect_equal(SD_A, "A")
  expect_equal(SD_B, "b")
  expect_equal(SD_OFFSET, "offset")
})
