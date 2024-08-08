context("test-g05-conic_solvers-copt")

test_that("test copt lp 0", {
  StandardTestLPs.test_lp_0(solver = "COPT")
})

test_that("test copt lp 1", {
  StandardTestLPs.test_lp_1(solver = "COPT")
})

test_that("test copt lp 2", {
  StandardTestLPs.test_lp_2(solver = "COPT")
})

test_that("test copt lp 3", {
  StandardTestLPs.test_lp_3(solver = "COPT")
})

test_that("test copt lp 4", {
  StandardTestLPs.test_lp_4(solver = "COPT")
})

test_that("test copt lp 5", {
  StandardTestLPs.test_lp_5(solver = "COPT")
})

test_that("test copt socp 0", {
  StandardTestSOCPs.test_socp_0(solver = "COPT")
})

test_that("test copt socp 1", {
  StandardTestSOCPs.test_socp_1(solver = "COPT", tolerance = 1e-3)
})

test_that("test copt socp 2", {
  StandardTestSOCPs.test_socp_2(solver = "COPT")
})

test_that("test copt socp 3", {
  # Axis 1.
  StandardTestSOCPs.test_socp_3ax1(solver = "COPT")
  # Axis 2.
  StandardTestSOCPs.test_socp_3ax2(solver = "COPT")
})

test_that("test copt mi lp 0", {
  StandardTestLPs.test_mi_lp_0(solver = "COPT")
})

test_that("test copt mi lp 1", {
  StandardTestLPs.test_mi_lp_1(solver = "COPT")
})

test_that("test copt mi lp 2", {
  StandardTestLPs.test_mi_lp_2(solver = "COPT")
})

test_that("test copt mi lp 3", {
  StandardTestLPs.test_mi_lp_3(solver = "COPT")
})

test_that("test copt mi lp 5", {
  StandardTestLPs.test_mi_lp_5(solver = "COPT")
})

test_that("test copt mi socp 1", {
  # COPT does not support MISOCP.
  expect_error(StandardTestSOCPs.test_mi_socp_1(solver = "COPT"), "do not support")
})

test_that("test copt sdp 1min", {
  StandardTestSDPs.test_sdp_1min(solver = "COPT")
})

test_that("test copt sdp 1max", {
  StandardTestSDPs.test_sdp_1max(solver = "COPT")
})

test_that("test copt sdp 2", {
  StandardTestSDPs.test_sdp_2(solver = "COPT")
})

test_that("test copt params", {
  n <- 10
  m <- 4
  set.seed(0)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  x <- matrix(rnorm(n))
  y <- A %*% x
  
  # Solve a simple basis pursuit problem for testing purposes.
  z <- Variable(n)
  objective <- Minimize(norm1(z))
  constraints <- list(A %*% z == y)
  problem <- Problem(objective, constraints)
  
  expect_error(solve(problem, solver = "COPT", invalid_kwarg = NULL))
  
  # Valid arg.
  result <- solve(problem, solver = "COPT", feastol = 1e-9)
})
