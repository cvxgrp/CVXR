context("test-g05-conic_solvers-pdlp")
TOL <- 1e-6

test_that("test pdlp lp 0", {
  StandardTestLPs.test_lp_0(solver = "PDLP")
})

test_that("test pdlp lp 1", {
  StandardTestLPs.test_lp_1(solver = "PDLP")
})

test_that("test pdlp lp 2", {
  StandardTestLPs.test_lp_2(solver = "PDLP")
})

test_that("test pdlp lp 3", {
  sth <- lp_3()
  result <- solve(sth@prob, solver = "PDLP")
  expect_equal(result$status, INFEASIBLE_OR_UNBOUNDED)
})

# We get the precise status when presolve is disabled.
test_that("test pdlp lp 3 no presolve", {
  params <- PrimalDualHybridGradientParams()
  params$presolve_options$use_glop <- FALSE
  StandardTestLPs.test_lp_3(solver = "PDLP", parameters_proto = params)
})

test_that("test pdlp lp 4", {
  sth <- lp_4()
  result <- solve(sth@prob, solver = "PDLP")
  expect_equal(result$status, INFEASIBLE_OR_UNBOUNDED)
})

test_that("test pdlp lp 4 no presolve", {
  params <- PrimalDualHybridGradientParams()
  params$presolve_options$use_glop <- FALSE
  StandardTestLPs.test_lp_4(solver = "PDLP", parameters_proto = params)
})

test_that("test pdlp lp 5", {
  StandardTestLPs.test_lp_5(solver = "PDLP")
})

test_that("test pdlp lp 6", {
  sth <- lp_6()
  result <- solve(sth@prob, solver = "PDLP")
  expect_equal(result$status, INFEASIBLE_OR_UNBOUNDED)
})

test_that("test pdlp bad parameters", {
  x <- Variable(1)
  prob <- Problem(Maximize(x), list(x <= 1))
  expect_error(solve(prob, solver = "PDLP", parameters_proto = "not a proto"))
})

test_that("test pdlp time limit", {
  sth <- lp_1()
  # Checks that the option doesn't error. A better test would be to solve a large
  # instance and check that the time limit is hit.
  result <- solve(sth, solver = "PDLP", time_limit_sec = 1.0)
})
