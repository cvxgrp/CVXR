context("test-g05-conic_solvers-sdpa")
TOl <- 1e-6

test_that("test sdpa lp 0", {
  StandardTestLPs.test_lp_0(solver = "SDPA")
})

test_that("test sdpa lp 1", {
  StandardTestLPs.test_lp_1(solver = "SDPA")
})

test_that("test sdpa lp 2", {
  StandardTestLPs.test_lp_2(solver = "SDPA")
})

test_that("test sdpa lp 3", {
  StandardTestLPs.test_lp_3(solver = "SDPA")
})

test_that("test sdpa lp 4", {
  StandardTestLPs.test_lp_4(solver = "SDPA")
})

test_that("test sdpa lp 5", {
  print("Skipping test. Known limitation of SDPA for degenerate LPs.")
  return()
  
  # This also tests the ability to pass solver options.
  StandardTestLPs.test_lp_5(solver = "SDPA", gammaStar = 0.86, epsilonDash = 8.0e-6, betaStar = 0.18, betaBar = 0.15)
})

test_that("test sdpa sdp 1", {
  # Minimization
  StandardTestSDPs.test_sdp_1min(solver = "SDPA")
  # Maximization
  StandardTestSDPs.test_sdp_1max(solver = "SDPA")
})

test_that("test sdpa sdp 2", {
  StandardTestSDPs.test_sdp_2(solver = "SDPA")
})
