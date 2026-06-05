# Tests for NLP constraint-dual recovery (CVXR extension, CVXPY 1.9.0).
#
# CVXPY's UNO interface returns an empty dual map (uno_nlpif.py:94). CVXR
# recovers constraint duals from Uno's constraint_dual: the lowered-constraint
# duals are mapped back to the user's constraints (Zero from an equality keeps
# the sign; NonNeg from an inequality flips it). These tests pin the convention
# by cross-checking against the conic path (Clarabel), whose duals are the
# established reference.
#
# CVXPY SOURCE: reductions/solvers/nlp_solvers/uno_nlpif.py (invert)

library(CVXR)

## @cvxpy NONE
test_that("NLP equality dual matches the conic dual", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## min ||x||^2 s.t. sum(x) == 1.5  =>  dual = -1.5 (both paths).
  xc <- Variable(2); cc <- sum_entries(xc) == 1.5
  psolve(Problem(Minimize(sum_squares(xc)), list(cc)), solver = "CLARABEL")
  dconic <- as.numeric(dual_value(cc))

  xn <- Variable(2); cn <- sum_entries(xn) == 1.5
  psolve(Problem(Minimize(sum_squares(xn)), list(cn)), solver = "UNO", nlp = TRUE)
  dnlp <- as.numeric(dual_value(cn))

  expect_equal(dnlp, dconic, tolerance = 1e-2)
  expect_equal(dnlp, -1.5, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("NLP inequality dual (>=) matches the conic dual", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## min ||x||^2 s.t. sum(x) >= 1.5  (active)  =>  dual = +1.5.
  xc <- Variable(2); cc <- sum_entries(xc) >= 1.5
  psolve(Problem(Minimize(sum_squares(xc)), list(cc)), solver = "CLARABEL")
  dconic <- as.numeric(dual_value(cc))

  xn <- Variable(2); cn <- sum_entries(xn) >= 1.5
  psolve(Problem(Minimize(sum_squares(xn)), list(cn)), solver = "UNO", nlp = TRUE)
  dnlp <- as.numeric(dual_value(cn))

  expect_equal(dnlp, dconic, tolerance = 1e-2)
  expect_equal(dnlp, 1.5, tolerance = 1e-2)
  expect_true(dnlp > 0)            # NonNeg dual is nonnegative
})

## @cvxpy NONE
test_that("NLP inequality dual (<=) matches the conic dual", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## min ||x - (2,2)||^2 s.t. sum(x) <= 1.5  (active)  =>  positive dual.
  xc <- Variable(2); cc <- sum_entries(xc) <= 1.5
  psolve(Problem(Minimize(sum_squares(xc - c(2, 2))), list(cc)), solver = "CLARABEL")
  dconic <- as.numeric(dual_value(cc))

  xn <- Variable(2); cn <- sum_entries(xn) <= 1.5
  psolve(Problem(Minimize(sum_squares(xn - c(2, 2))), list(cn)), solver = "UNO", nlp = TRUE)
  dnlp <- as.numeric(dual_value(cn))

  expect_equal(dnlp, dconic, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("NLP recovers a vector-valued (multi-row) constraint dual", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## Elementwise equality x == c(1, 2): dual is length-2, matches conic.
  b <- c(1, 2)
  xc <- Variable(2); cc <- xc == b
  psolve(Problem(Minimize(sum_squares(xc)), list(cc)), solver = "CLARABEL")
  dconic <- as.numeric(dual_value(cc))

  xn <- Variable(2); cn <- xn == b
  psolve(Problem(Minimize(sum_squares(xn)), list(cn)), solver = "UNO", nlp = TRUE)
  dnlp <- as.numeric(dual_value(cn))

  expect_equal(length(dnlp), 2L)
  expect_equal(dnlp, dconic, tolerance = 1e-2)
})
