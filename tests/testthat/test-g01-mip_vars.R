context("test-g01-mip_vars")
TOL <- 1e-6

x_bool <- Variable(boolean = TRUE)
y_int <- Variable(integer = TRUE)
A_bool <- Variable(3, 2, boolean = TRUE)
B_int <- Variable(2, 3, integer = TRUE)

# Check for all installed QP solvers.
solvers <- CVXR:::MIP_SOLVERS

bool_prob <- function(solver) {
  # Bool in objective.
  obj <- Minimize(abs(x_bool - 0.2))
  p <- Problem(obj, list())
  result <- solve(p, solver = solver)
  expect_equal(result$value, 0.2, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = TOL)
  
  # Bool in constraints.
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(abs(x_bool) <= t))
  result <- solve(p, solver = solver)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = 1e-4)
  
  # Matrix Bool in objective.
  C <- cbind(c(0,1,0), c(1,1,1))
  obj <- Minimize(sum(abs(A_bool - C)))
  p <- Problem(obj, list())
  result <- solve(p, solver = solver)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(A_bool), C, tolerance = 1e-4)
  
  # Matrix Bool in constraint.
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(sum(abs(A_bool - C)) <= t))
  result <- solve(p, solver = solver)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(A_bool), C, tolerance = 1e-4)
}

int_prob <- function(solver) {
  # Int in objective.
  obj <- Minimize(abs(y_int - 0.2))
  p <- Problem(obj, list())
  result <- solve(p, solver = solver)
  expect_equal(result$value, 0.2, tolerance = TOL)
  expect_equal(result$getValue(y_int), 0, tolerance = TOL)
  
  # Infeasible integer problem.
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(y_int == 0.5, t >= 0))
  result <- solve(p, solver = solver)
  expect_true(result$status %in% INF_OR_UNB)
}

int_socp <- function(solver) {
  # Int in objective.
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(square(y_int - 0.2) <= t))
  result <- solve(p, solver = solver)
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(y_int), 0, tolerance = TOL)
}

bool_socp <- function(solver) {
  # Bool in objective.
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(square(x_bool - 0.2) <= t))
  result <- solve(p, solver = solver)
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = TOL)
}

test_that("test all solvers", {
  skip_on_cran()
  if(length(MIP_SOLVERS) == 0) {
    print("No mixed-integer solver is installed. Skipping test.")
    return()
  }
  
  for(solver in solvers) {
    bool_prob(solver)
    int_prob(solver)
    if(solver %in% c("CPLEX", "GUROBI", "MOSEK", "XPRESS")) {
      if(solver != "XPRESS")   # Issue #1815.
        bool_socp(solver)
      int_socp(solver)
    }
  }
})
