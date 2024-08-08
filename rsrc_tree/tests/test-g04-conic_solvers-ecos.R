context("test-g04-conic_solvers-ecos")
TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test ECOS options", {
  # Test that all ECOS solver options work.
  # Test ecos
  # feastol, abstol, reltol, feastol_inacc,
  # abstol_inacc, and reltol_inacc for tolerance values
  # max_iters for the maximum number of iterations,
  EPS <- 1e-4
  
  prob <- Problem(Minimize(norm1(x) + 1.0), list(x == 0))
  for(i in 1:2) {
    result <- solve(prob, solver = "ECOS", feastol = EPS, abstol = EPS, reltol = EPS,
                    feastol_inacc = EPS, abstol_inacc = EPS, reltol_inacc = EPS,
                    max_iters = 20, verbose = TRUE, warm_start = TRUE)
  }
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
})

test_that("test ECOS lp 0", {
  StandardTestLPs.test_lp_0(solver = "ECOS")
})

test_that("test ECOS lp 1", {
  StandardTestLPs.test_lp_1(solver = "ECOS")
})

test_that("test ECOS lp 2", {
  StandardTestLPs.test_lp_2(solver = "ECOS")
})

test_that("test ECOS lp 3", {
  StandardTestLPs.test_lp_3(solver = "ECOS")
})

test_that("test ECOS lp 4", {
  StandardTestLPs.test_lp_4(solver = "ECOS")
})

test_that("test ECOS lp 5", {
  StandardTestLPs.test_lp_5(solver = "ECOS")
})

test_that("test ECOS socp 0", {
  StandardTestSOCPs.test_socp_0(solver = "ECOS")
})

test_that("test ECOS socp 1", {
  StandardTestSOCPs.test_socp_1(solver = "ECOS")
})

test_that("test ECOS socp 2", {
  StandardTestSOCPs.test_socp_2(solver = "ECOS")
})

test_that("test ECOS socp 3", {
  # Axis 1
  StandardTestSOCPs.test_socp_3ax1(solver = "ECOS")
  # Axis 2
  StandardTestSOCPs.test_socp3_ax2(solver = "ECOS")
})

test_that("test ECOS expcone 1", {
  StandardTestECPs.test_expcone_1(solver = "ECOS")
})

test_that("test ECOS exp soc 1", {
  StandardTestMixedCPs.test_exp_soc_1(solver = "ECOS")
})
