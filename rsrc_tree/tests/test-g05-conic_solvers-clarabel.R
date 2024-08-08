context("test-g05-conic_solvers-clarabel")
TOL <- 1e-6

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test clarabel lp 0", {
  StandardTestLPs.test_lp_0(solver = "CLARABEL")
})

test_that("test clarabel lp 1", {
  StandardTestLPs.test_lp_1(solver = "CLARABEL")
})

test_that("test clarabel lp 2", {
  StandardTestLPs.test_lp_2(solver = "CLARABEL")
})

test_that("test clarabel lp 3", {
  StandardTestLPs.test_lp_3(solver = "CLARABEL")
})

test_that("test clarabel lp 4", {
  StandardTestLPs.test_lp_4(solver = "CLARABEL")
})

test_that("test clarabel lp 5", {
  StandardTestLPs.test_lp_5(solver = "CLARABEL")
})

test_that("test clarabel qp 0", {
  StandardTestQPs.test_qp_0(solver = "CLARABEL")
})

test_that("test clarabel qp 0 linear obj", {
  StandardTestQPs.test_qp_0(solver = "CLARABEL", use_quad_obj = FALSE)
})

test_that("test clarabel socp 0", {
  StandardTestSOCPs.test_socp_0(solver = "CLARABEL")
})

test_that("test clarabel socp 1", {
  StandardTestSOCPs.test_socp_1(solver = "CLARABEL")
})

test_that("test clarabel socp 2", {
  StandardTestSOCPs.test_socp_2(solver = "CLARABEL")
})

test_that("test clarabel socp 3", {
  # Axis 1
  StandardTestSOCPs.test_socp_3ax1(solver = "CLARABEL")
  # Axis 2
  StandardTestSOCPs.test_socp_3ax2(solver = "CLARABEL")
})

test_that("test clarabel expcone 1", {
  StandardTestECPs.test_expcone_1(solver = "CLARABEL")
})

test_that("test clarabel exp soc 1", {
  StandardTestMixedCPs.test_exp_soc_1(solver = "CLARABEL")
})

test_that("test clarabel pcp 0", {
  StandardTestSOCPs.test_socp_0(solver = "CLARABEL")
})

test_that("test clarabel pcp 1", {
  StandardTestSOCPs.test_socp_1(solver = "CLARABEL")
})

test_that("test clarabel pcp 2", {
  StandardTestSOCPs.test_socp_2(solver = "CLARABEL")
})
