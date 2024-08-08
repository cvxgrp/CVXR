context("test_custom_solver")

CustomQPSolverCalled <- function() { stop("Custom QP solver called") }
CustomConicSolverCalled <- function() { stop("Custom conic solver called") }

CustomQPSolver <- setClass("CustomQPSolver", contains = "OSQP")
setMethod("name", "CustomQPSolver", function(x) { "CUSTOM_QP_SOLVER" })
setMethod("solve_via_data", "CustomQPSolver", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) { CustomQPSolverCalled() })

CustomConicSolver <- setClass("CustomConicSolver", contains = "SCS")
setMethod("name", "CustomConicSolver", function(x) { "CUSTOM_CONIC_SOLVER" })
setMethod("solve_via_data", "CustomConicSolver", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) { CustomConicSolverCalled() })

ConflictingCustomSolver <- setClass("ConflictingCustomSolver", contains = "OSQP")
setMethod("name", "ConflictingCustomSolver", function(x) { "OSQP" })

custom_qp_solver <- CustomQPSolver()
custom_conic_solver <- CustomConicSolver()

solve_example_qp <- function(solver) {
  x <- Variable()
  quadratic <- sum_squares(x)
  problem <- Problem(Minimize(quadratic))
  result <- solve(problem, solver = solver)
}

solve_example_mixed_integer_qp <- function(solver) {
  x <- Variable()
  z <- Variable(integer = TRUE)
  quadratic <- sum_squares(x + z)
  problem <- Problem(Minimize(quadratic))
  result <- solve(problem, solver = solver)
}

solve_example_socp <- function(solver) {
  x <- Variable(2)
  y <- Variable()
  quadratic <- sum_squares(x)
  problem <- Problem(Minimize(quadratic), list(SOC(y, x)))
  result <- solve(problem, solver = solver)
}

solve_example_mixed_integer_socp <- function(solver) {
  x <- Variable(2)
  y <- Variable()
  z <- Variable(integer = TRUE)
  quadratic <- sum_squares(x + z)
  problem <- Problem(Minimize(quadratic), list(SOC(y, x)))
  result <- solve(problem, solver = solver)
}

test_that("test custom continuous qp solver can solve continuous qp", {
  expect_error(solve_example_qp(solver = custom_qp_solver), "Custom QP solver called", fixed = TRUE)
})

test_that("test custom mip qp solver can solve mip qp", {
  custom_qp_solver@MIP_CAPABLE <- TRUE
  solve_example_mixed_integer_qp(solver = custom_qp_solver, "Custom QP solver called", fixed = TRUE)
})

test_that("test custom continuous qp solver cannot solve mip qp", {
  custom_qp_solver@MIP_CAPABLE <- FALSE
  solve_example_mixed_integer_qp(solver = custom_qp_solver)
})

test_that("test custom qp solver cannot solve socp", {
  expect_error(solve_example_socp(solver = custom_qp_solver))
})

test_that("test custom continuous conic solver can solve continuous socp", {
  expect_error(solve_example_socp(solver = custom_conic_solver), "Custom conic solver called", fixed = TRUE)
})

test_that("test custom mip conic solver can solve mip socp", {
  custom_conic_solver@MIP_CAPABLE <- TRUE
  supported_constraints <- custom_conic_solver@SUPPORTED_CONSTRAINTS
  custom_conic_solver@MI_SUPPORTED_CONSTRAINTS <- supported_constraints
  expect_error(solve_example_mixed_integer_socp(solver = custom_conic_solver), "Custom conic solver called", fixed = TRUE)
})

test_that("test custom conflicting solver fails", {
  expect_error(solve_example_qp(solver = ConflictingCustomSolver()))
})

