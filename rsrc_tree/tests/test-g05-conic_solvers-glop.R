context("test-g05-conic_solvers-glop")
TOL <- 1e-6

test_that("test glop lp 0", {
  StandardTestLPs.test_lp_0(solver = "GLOP")
})

test_that("test glop lp 1", {
  StandardTestLPs.test_lp_1(solver = "GLOP")
})

test_that("test glop lp 2", {
  StandardTestLPs.test_lp_2(solver = "GLOP")
})

test_that("test glop lp 3 no preprocessing", {
  params <- GlopParameters()
  params$use_preprocessing <- FALSE
  StandardTestLPs.test_lp_3(solver = "GLOP", parameters_proto = params)
})

# With preprocessing enabled, Glop internally detects
# INFEASIBLE_OR_UNBOUNDED. This status is translated to
# MPSOLVER_INFEASIBLE. See
# https://github.com/google/or-tools/blob/b37d9c786b69128f3505f15beca09e89bf078a89/ortools/linear_solver/glop_utils.cc#L25-L38.
test_that("test glop lp 3", {
  print("Skipping test. Known limitation of the GLOP interface")
  StandardTestLPs.test_lp_3(solver = "GLOP")
})

test_that("test glop lp 4", {
  StandardTestLPs.test_lp_4(solver = "GLOP")
})

test_that("test glop lp 5", {
  StandardTestLPs.test_lp_5(solver = "GLOP")
})

test_that("test glop lp 6 no preprocessing", {
  params <- GlopParameters()
  params$use_preprocessing <- FALSE
  StandardTestLPs.test_lp_6(solver = "GLOP", parameters_proto = params)
})

# Same issue as with test_glop_lp_3.
test_that("test glop lp 6", {
  print("Skipping test. Known limitation of the GLOP interface")
  StandardTestLPs.test_lp_6(solver = "GLOP")
})

test_that("test glop bad parameters", {
  x <- Variable(1)
  prob <- Problem(Maximize(x), list(x <= 1))
  expect_error(solve(prob, solver = "GLOP", parameters_proto = "not a proto"))
})

test_that("test glop time limit", {
  sth <- lp_1()
  # Checks that the option doesn't error. A better test would be to solve
  # a large instance and check that the time limit is hit.
  result <- solve(sth, solver = "GLOP", time_limit_sec = 1.0)
})
