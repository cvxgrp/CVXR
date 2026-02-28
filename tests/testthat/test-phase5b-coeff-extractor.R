# Tests for Phase 5b: CoeffExtractor

library(CVXR)

# ═══════════════════════════════════════════════════════════════════
# CoeffExtractor basic extraction
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("CoeffExtractor construction from InverseData", {
  x <- Variable(3)
  p <- Problem(Minimize(sum_entries(x)))
  inv_data <- InverseData(p)

  ext <- CoeffExtractor(inv_data)
  expect_equal(ext@x_length, 3L)
  expect_equal(length(ext@id_to_col), 1L)
})

## @cvxpy NONE
test_that("coeff_affine extracts A and b for x + 3", {
  x <- Variable(1)
  p <- Problem(Minimize(x + 3))
  inv_data <- InverseData(p)
  ext <- CoeffExtractor(inv_data)

  result <- coeff_affine(ext, list(x + Constant(3)))
  A <- result$A
  b <- result$b

  ## A should be 1x1 with value 1 (coefficient of x)
  expect_equal(nrow(A), 1L)
  expect_equal(ncol(A), 1L)
  expect_equal(as.numeric(A[1, 1]), 1.0)

  ## b should be 3 (constant term)
  expect_equal(as.numeric(b), 3.0)
})

## @cvxpy NONE
test_that("coeff_affine extracts A and b for 2*x", {
  x <- Variable(1)
  p <- Problem(Minimize(2 * x))
  inv_data <- InverseData(p)
  ext <- CoeffExtractor(inv_data)

  result <- coeff_affine(ext, list(2 * x))
  A <- result$A
  b <- result$b

  expect_equal(as.numeric(A[1, 1]), 2.0)
  expect_equal(as.numeric(b), 0.0)
})

## @cvxpy NONE
test_that("coeff_affine handles vector variable", {
  x <- Variable(3)
  p <- Problem(Minimize(sum_entries(x)))
  inv_data <- InverseData(p)
  ext <- CoeffExtractor(inv_data)

  ## Extract coefficients for the sum expression
  result <- coeff_affine(ext, list(sum_entries(x)))
  A <- result$A
  b <- result$b

  ## sum(x) = 1*x1 + 1*x2 + 1*x3, so A is 1x3 with all ones
  expect_equal(nrow(A), 1L)
  expect_equal(ncol(A), 3L)
  expect_equal(as.numeric(Matrix::colSums(A)), rep(1.0, 3))
  expect_equal(as.numeric(b), 0.0)
})

## @cvxpy NONE
test_that("coeff_affine handles multiple expressions", {
  x <- Variable(2)
  p <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  inv_data <- InverseData(p)
  ext <- CoeffExtractor(inv_data)

  ## Extract for two expressions: x[1] and x[2]
  result <- coeff_affine(ext, list(x))
  A <- result$A
  b <- result$b

  ## x is a 2-vector, so A should be 2x2 identity
  expect_equal(nrow(A), 2L)
  expect_equal(ncol(A), 2L)
  expect_equal(as.numeric(b), c(0.0, 0.0))
})

## @cvxpy NONE
test_that("coeff_affine handles constant expression", {
  x <- Variable(1)
  p <- Problem(Minimize(x))
  inv_data <- InverseData(p)
  ext <- CoeffExtractor(inv_data)

  result <- coeff_affine(ext, list(Constant(5)))
  expect_equal(as.numeric(result$b), 5.0)
})
