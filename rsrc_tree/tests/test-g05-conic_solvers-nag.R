context("test-g05-conic_solvers-nag")

test_that("test nag lp 0", {
  StandardTestLPs.test_lp_0(solver = "NAG")
})

test_that("test nag lp 1", {
  StandardTestLPs.test_lp_1(solver = "NAG")
})

test_that("test nag lp 2", {
  StandardTestLPs.test_lp_2(solver = "NAG")
})

test_that("test nag lp 3", {
  StandardTestLPs.test_lp_3(solver = "NAG")
})

test_that("test nag lp 4", {
  StandardTestLPs.test_lp_4(solver = "NAG")
})

test_that("test nag lp 5", {
  StandardTestLPs.test_lp_5(solver = "NAG")
})

test_that("test nag socp 0", {
  StandardTestSOCPs.test_socp_0(solver = "NAG")
})

test_that("test nag socp 1", {
  StandardTestSOCPs.test_socp_1(solver = "NAG")
})

test_that("test nag socp 2", {
  StandardTestSOCPs.test_socp_2(solver = "NAG")
})

test_that("test nag socp 3", {
  # Axis 1.
  StandardTestSOCPs.test_socp_3ax1(solver = "NAG")
  # Axis 2.
  StandardTestSOCPs.test_socp_3ax2(solver = "NAG")
})
