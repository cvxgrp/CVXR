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

# ═══════════════════════════════════════════════════════════════════
# Block-structured quadratic extraction (CVXPY TestBlockQuadExtraction)
#
# `.extract_block_quad` is CVXR's `CoeffExtractor._extract_block_quad`
# (coeff_extractor.py:300-363) on the numeric path: triplets in, scaled
# triplets out, 1-based indices (CVXPY's are 0-based; +1 throughout).
# CVXPY's `c_part` has one column per parameter slice; the numeric path has
# the single constant slice, so `c_part` is a plain vector here.
# ═══════════════════════════════════════════════════════════════════

eye_coo <- function(n) data.frame(i = seq_len(n), j = seq_len(n), x = rep(1, n))

## @cvxpy test_coeff_extractor.py::TestBlockQuadExtraction::test_block_extraction_contiguous_identity
test_that("block extraction: contiguous blocks with identity P", {
  P_coo <- eye_coo(6L)
  block_indices <- list(1:3, 4:6)
  res <- .extract_block_quad(P_coo, c(1.0, 2.0), block_indices)
  expect_length(res$x, 6L)
  o <- order(res$i)
  expect_equal(res$i[o], 1:6)
  expect_equal(res$j[o], 1:6)
  expect_equal(res$x[o], c(1, 1, 1, 2, 2, 2))
})

## @cvxpy test_coeff_extractor.py::TestBlockQuadExtraction::test_block_extraction_noncontiguous
test_that("block extraction: non-contiguous blocks (row-wise reduction of a 3x2)", {
  P_coo <- eye_coo(6L)
  ## x of shape (3, 2) flattened column-major; one block per row.
  block_indices <- list(c(1L, 4L), c(2L, 5L), c(3L, 6L))
  res <- .extract_block_quad(P_coo, c(1.0, 2.0, 3.0), block_indices)
  expect_length(res$x, 6L)
  expect_equal(res$i, res$j)
  expected <- c(1, 2, 3, 1, 2, 3)
  expect_equal(res$x, expected[res$i])
})

## @cvxpy test_coeff_extractor.py::TestBlockQuadExtraction::test_block_extraction_sparse_blocks
test_that("block extraction: sparse (non-identity) P blocks", {
  P_coo <- data.frame(i = c(1, 1, 2, 2, 3, 4, 5, 6),
                      j = c(1, 2, 1, 2, 3, 4, 5, 6),
                      x = c(1.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0))
  res <- .extract_block_quad(P_coo, c(2.0, 3.0), list(1:3, 4:6))
  expect_length(res$x, 8L)
  for (k in seq_along(res$x)) {
    coef <- if (res$i[k] <= 3) 2.0 else 3.0
    orig <- P_coo$x[P_coo$i == res$i[k] & P_coo$j == res$j[k]]
    expect_equal(res$x[k], orig * coef)
  }
})

## @cvxpy test_coeff_extractor.py::TestBlockQuadExtraction::test_block_extraction_empty_blocks
test_that("block extraction: a block with no P entries contributes nothing", {
  P_coo <- data.frame(i = 1:3, j = 1:3, x = c(1, 1, 1))
  res <- .extract_block_quad(P_coo, c(1.0, 2.0), list(1:3, 4:6))
  expect_length(res$x, 3L)
  expect_equal(res$x, c(1, 1, 1))
})

## @cvxpy test_coeff_extractor.py::TestBlockQuadExtraction::test_block_extraction_zero_coefficients
test_that("block extraction: blocks with a zero coefficient are skipped", {
  res <- .extract_block_quad(eye_coo(4L), c(1.0, 0.0), list(1:2, 3:4))
  expect_length(res$x, 2L)
  expect_true(all(res$i <= 2))
})

## @cvxpy test_coeff_extractor.py::TestBlockQuadExtraction::test_block_extraction_all_empty
test_that("block extraction: all-zero coefficients give empty triplets", {
  res <- .extract_block_quad(eye_coo(4L), c(0.0, 0.0), list(1:2, 3:4))
  expect_length(res$x, 0L)
  expect_length(res$i, 0L)
  expect_length(res$j, 0L)
})

## @cvxpy test_coeff_extractor.py::TestBlockQuadExtraction::test_symbolic_quad_form_with_block_indices
test_that("SymbolicQuadForm stores block_indices and returns them from get_data", {
  x <- Variable(6)
  P <- diag(6)
  expr <- quad_form(x, P)
  block_indices <- list(1:3, 4:6)
  sqf <- SymbolicQuadForm(x, Constant(P), expr, block_indices = block_indices)
  expect_identical(sqf@block_indices, block_indices)
  data <- get_data(sqf)
  expect_length(data, 2L)
  expect_identical(data[[2L]], block_indices)
})

## @cvxpy test_coeff_extractor.py::TestBlockQuadExtraction::test_symbolic_quad_form_without_block_indices
test_that("SymbolicQuadForm without block_indices keeps NULL", {
  x <- Variable(2)
  P <- diag(2)
  sqf <- SymbolicQuadForm(x, Constant(P), quad_form(x, P))
  expect_null(sqf@block_indices)
  data <- get_data(sqf)
  expect_length(data, 2L)
  expect_null(data[[2L]])
})

## @cvxpy NONE
test_that("quad_over_lin with axis reaches the QP path with per-output blocks", {
  ## Before 1.9.2.9141 the quad canon returned NULL for any axis, so the atom
  ## went to the SOC canon and a QP-only solver refused the problem
  ## (claims-audit row C0670). The blocks are column-major flat indices.
  X <- Variable(c(3, 2))
  for (ax in c(2L, 1L)) {
    e <- quad_over_lin(X, 2, axis = ax)
    out <- quad_canonicalize(e, list(X, Constant(2)))
    expect_false(is.null(out))
    sqf <- out[[1L]]
    expect_true(inherits(sqf, "CVXR::SymbolicQuadForm"))
    if (ax == 2L) {
      expect_identical(sqf@block_indices, list(1:3, 4:6))
    } else {
      expect_identical(sqf@block_indices, list(c(1L, 4L), c(2L, 5L), c(3L, 6L)))
    }
    P <- value(.args(sqf)[[2L]])
    expect_equal(as.numeric(Matrix::diag(P)), rep(0.5, 6))
  }
  ## A non-constant denominator still goes to the cone canon (dcp2cone.py:235-237).
  y <- Variable(nonneg = TRUE)
  expect_null(quad_canonicalize(quad_over_lin(X, y, axis = 2L), list(X, y)))
})
