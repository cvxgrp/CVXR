context("test-g05-conic_solvers-scip")

test_that("test scip lp 0", {
  StandardTestLPs.test_lp_0(solver = "SCIP")
})

test_that("test scip lp 1", {
  StandardTestLPs.test_lp_1(solver = "SCIP")
})

test_that("test scip lp 2", {
  StandardTestLPs.test_lp_2(solver = "SCIP")
})

test_that("test scip lp 3", {
  StandardTestLPs.test_lp_3(solver = "SCIP")
})

test_that("test scip lp 4", {
  StandardTestLPs.test_lp_4(solver = "SCIP")
})

test_that("test scip socp 0", {
  StandardTestSOCPs.test_socp_0(solver = "SCIP")
})

test_that("test scip socp 1", {
  StandardTestSOCPs.test_socp_1(solver = "SCIP")
})

test_that("test scip socp 2", {
  StandardTestSOCPs.test_socp_2(solver = "SCIP")
})

test_that("test scip socp 3", {
  # Axis 1.
  StandardTestSOCPs.test_socp_3ax1(solver = "SCIP")
  # Axis 2.
  StandardTestSOCPs.test_socp_3ax2(solver = "SCIP")
})

test_that("test scip mi lp 0", {
  StandardTestLPs.test_mi_lp_0(solver = "SCIP")
})

test_that("test scip mi lp 1", {
  StandardTestLPs.test_mi_lp_1(solver = "SCIP")
})

test_that("test scip mi lp 2", {
  StandardTestLPs.test_mi_lp_2(solver = "SCIP")
})

test_that("test scip mi lp 3", {
  StandardTestLPs.test_mi_lp_3(solver = "SCIP")
})

test_that("test scip mi lp 5", {
  StandardTestLPs.test_mi_lp_5(solver = "SCIP")
})

test_that("test scip mi socp 1", {
  StandardTestSOCPs.test_mi_socp_1(solver = "SCIP", tolerance = 1e-3)
})

test_that("test scip mi socp 2", {
  StandardTestSOCPs.test_mi_socp_2(solver = "SCIP")
})

get_simple_problem <- function() {
  # Example problem that can be used within additional tests.
  x <- Variable()
  y <- Variable()
  constraints <- list(x >= 0, y >= 1, x + y <= 4)
  obj <- Maximize(x)
  prob <- Problem(obj, constraints)
  return(prob)
}

test_that("test scip test params - no params set", {
  prob <- get_simple_problem()
  result <- solve(prob, solver = "SCIP")
  # Important that passes without raising an error also check obj.
  expect_equal(result$value, 3)
})

test_that("test scip test params - valid params", {
  prob <- get_simple_problem()
  result <- solve(prob, solver = "SCIP", gp = FALSE)
  # Important that passes without raising an error also check obj.
  expect_equal(result$value, 3)
})

test_that("test scip test params - valid scip params", {
  prob <- get_simple_problem()
  result <- solve(prob, solver = "SCIP", scip_params = list('lp/fastmip' = 1, 'limits/gap' = 0.1))
  # Important that passes without raising an error also check obj.
  expect_equal(result$value, 3)
})

test_that("test scip test params - invalid params", {
  prob <- get_simple_problem()
  # Since an invalid NON-scip param is passed, an error is expected t be raised when calling solve.
  expect_error(solve(prob, solver = "SCIP", a = "what?"),
               "One or more solver params in ['a'] are not valid: 'Not a valid parameter name'")
})

test_that("test scip test params - invalid scip params", {
  prob <- get_simple_problem()
  # Since an invalid NON-scip param is passed, an error is expected t be raised when calling solve.
  expect_error(solve(prob, solver = "SCIP", scip_params = list(a = "what?")),
               "One or more solver params in ['a'] are not valid: 'Not a valid parameter name'")
})
