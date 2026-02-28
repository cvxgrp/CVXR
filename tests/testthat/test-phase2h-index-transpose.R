## Tests for Phase 2H: Index + Transpose Atoms
## Translated from CVXPY tests in test_problem.py

# ── Index construction ────────────────────────────────────────────────

## @cvxpy NONE
test_that("Index on column vector: x[1:2]", {
  x <- Variable(3)
  idx <- x[1:2]
  expect_s3_class(idx, "CVXR::Index")
  expect_equal(idx@shape, c(2L, 1L))
})

## @cvxpy NONE
test_that("Index on column vector: single element x[1]", {
  x <- Variable(3)
  idx <- x[1]
  expect_s3_class(idx, "CVXR::Index")
  expect_equal(idx@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("Index on matrix: x[1:2, ]", {
  x <- Variable(c(3, 4))
  idx <- x[1:2, ]
  expect_s3_class(idx, "CVXR::Index")
  expect_equal(idx@shape, c(2L, 4L))
})

## @cvxpy NONE
test_that("Index on matrix: x[, 2]", {
  x <- Variable(c(3, 4))
  idx <- x[, 2]
  expect_s3_class(idx, "CVXR::Index")
  expect_equal(idx@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Index on matrix: x[1, 2]", {
  x <- Variable(c(3, 4))
  idx <- x[1, 2]
  expect_s3_class(idx, "CVXR::Index")
  expect_equal(idx@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("Index on matrix: x[1:2, 3:4]", {
  x <- Variable(c(3, 4))
  idx <- x[1:2, 3:4]
  expect_s3_class(idx, "CVXR::Index")
  expect_equal(idx@shape, c(2L, 2L))
})

## @cvxpy NONE
test_that("Index on matrix: full slice x[,]", {
  x <- Variable(c(3, 4))
  result <- x[, ]
  ## This should return x itself (no-op)
  expect_s3_class(result, "CVXR::Expression")
})

# ── Index numeric value ───────────────────────────────────────────────

## @cvxpy NONE
test_that("Index numeric_value on Constant vector", {
  c1 <- Constant(matrix(c(1, 2, 3), 3, 1))
  idx <- c1[1:2]
  v <- value(idx)
  expect_equal(v, matrix(c(1, 2), 2, 1))
})

## @cvxpy NONE
test_that("Index numeric_value on Constant matrix", {
  m <- matrix(1:12, 3, 4)
  c1 <- Constant(m)
  idx <- c1[1:2, 3:4]
  v <- value(idx)
  expect_equal(v, m[1:2, 3:4, drop = FALSE])
})

## @cvxpy NONE
test_that("Index single element from Constant", {
  c1 <- Constant(matrix(c(10, 20, 30), 3, 1))
  idx <- c1[2]
  v <- value(idx)
  expect_equal(as.numeric(v), 20)
})

## @cvxpy NONE
test_that("Index on column slice of Constant matrix", {
  m <- matrix(1:6, 2, 3)
  c1 <- Constant(m)
  idx <- c1[, 2]
  v <- value(idx)
  expect_equal(v, matrix(c(3, 4), 2, 1))
})

# ── Index on expressions ──────────────────────────────────────────────

## @cvxpy NONE
test_that("Index on expression (Variable + Constant)", {
  x <- Variable(3)
  expr <- x + Constant(matrix(c(10, 20, 30), 3, 1))
  idx <- expr[1]
  expect_s3_class(idx, "CVXR::Index")
  expect_equal(idx@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("Index preserves DCP (affine stays affine)", {
  x <- Variable(3)
  idx <- x[1:2]
  expect_true(is_affine(idx))
  expect_true(is_convex(idx))
  expect_true(is_concave(idx))
})

# ── Index with reference semantics (value propagation) ────────────────

## @cvxpy NONE
test_that("Index value propagates from Variable", {
  x <- Variable(3)
  idx <- x[2]
  value(x) <- matrix(c(10, 20, 30), 3, 1)
  expect_equal(as.numeric(value(idx)), 20)
})

## @cvxpy NONE
test_that("Index on row slice propagates", {
  x <- Variable(c(3, 2))
  idx <- x[1:2, ]
  value(x) <- matrix(1:6, 3, 2)
  expect_equal(value(idx), matrix(c(1, 2, 4, 5), 2, 2))
})

# ── Index error cases ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("Index out of bounds", {
  x <- Variable(3)
  expect_error(x[5], "out of bounds")
})

## @cvxpy NONE
test_that("Index on matrix with single subscript errors", {
  x <- Variable(c(3, 4))
  expect_error(x[5], "not supported")
})

# ── Index naming ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Index has meaningful name", {
  x <- Variable(3, name = "x")
  idx <- x[1:2]
  nm <- expr_name(idx)
  expect_true(grepl("x", nm))
})

# ── Transpose construction ────────────────────────────────────────────

## @cvxpy NONE
test_that("Transpose reverses shape", {
  x <- Variable(c(3, 4))
  tx <- t(x)
  expect_s3_class(tx, "CVXR::Transpose")
  expect_equal(tx@shape, c(4L, 3L))
})

## @cvxpy NONE
test_that("Transpose of column vector gives row vector", {
  x <- Variable(3)
  tx <- t(x)
  expect_equal(tx@shape, c(1L, 3L))
})

## @cvxpy NONE
test_that("Transpose of row vector gives column vector", {
  x <- Variable(c(1, 3))
  tx <- t(x)
  expect_equal(tx@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Transpose of scalar is scalar", {
  x <- Variable(1)
  tx <- t(x)
  expect_equal(tx@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("Transpose of square matrix", {
  x <- Variable(c(3, 3))
  tx <- t(x)
  expect_equal(tx@shape, c(3L, 3L))
})

# ── Transpose numeric value ──────────────────────────────────────────

## @cvxpy NONE
test_that("Transpose numeric_value on Constant matrix", {
  m <- matrix(1:6, 2, 3)
  c1 <- Constant(m)
  tc <- t(c1)
  v <- value(tc)
  expect_equal(v, t(m))
})

## @cvxpy NONE
test_that("Transpose numeric_value on Constant vector", {
  c1 <- Constant(matrix(c(1, 2, 3), 3, 1))
  tc <- t(c1)
  v <- value(tc)
  expect_equal(v, matrix(c(1, 2, 3), 1, 3))
})

# ── Transpose DCP ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Transpose of affine is affine", {
  x <- Variable(c(3, 4))
  tx <- t(x)
  expect_true(is_affine(tx))
  expect_true(is_convex(tx))
  expect_true(is_concave(tx))
})

# ── Transpose symmetry ──────────────────────────────────────────────

## @cvxpy NONE
test_that("Transpose of symmetric is symmetric", {
  x <- Variable(c(3, 3), symmetric = TRUE)
  tx <- t(x)
  expect_true(is_symmetric(tx))
})

# ── Transpose name ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("Transpose has meaningful name", {
  x <- Variable(c(3, 4), name = "x")
  tx <- t(x)
  nm <- expr_name(tx)
  expect_equal(nm, "t(x)")
})

# ── Transpose reference semantics ────────────────────────────────────

## @cvxpy NONE
test_that("Transpose value propagates from Variable", {
  x <- Variable(c(2, 3))
  tx <- t(x)
  value(x) <- matrix(1:6, 2, 3)
  expect_equal(value(tx), t(matrix(1:6, 2, 3)))
})

# ── Composition: Transpose of Index ──────────────────────────────────

## @cvxpy NONE
test_that("t(x[1:2, ]) gives correct shape", {
  x <- Variable(c(3, 4))
  result <- t(x[1:2, ])
  expect_equal(result@shape, c(4L, 2L))
})

# ── Composition: Index of Transpose ──────────────────────────────────

## @cvxpy NONE
test_that("t(x)[1:2, ] gives correct shape", {
  x <- Variable(c(3, 4))
  tx <- t(x)
  result <- tx[1:2, ]
  expect_equal(result@shape, c(2L, 3L))
})

# ── Composition: Index in constraints ─────────────────────────────────

## @cvxpy NONE
test_that("x[1] == 5 creates Equality constraint", {
  x <- Variable(3)
  constr <- x[1] == 5
  expect_s3_class(constr, "CVXR::Equality")
  expect_equal(constr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("x[1:2] <= 10 creates Inequality constraint", {
  x <- Variable(3)
  constr <- x[1:2] <= 10
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(2L, 1L))
})

## @cvxpy NONE
test_that("x[1] >= 0 creates Inequality constraint", {
  x <- Variable(3)
  constr <- x[1] >= 0
  expect_s3_class(constr, "CVXR::Inequality")
})

# ── Composition: t(x) %*% y is scalar for vectors ────────────────────

## @cvxpy NONE
test_that("t(x) %*% y gives scalar for vectors", {
  x <- Variable(3)
  y <- Variable(3)
  product <- t(x) %*% y
  expect_s3_class(product, "CVXR::MulExpression")
  expect_equal(product@shape, c(1L, 1L))
})

# ── Objective with Index ──────────────────────────────────────────────

## @cvxpy NONE
test_that("Minimize(x[1]) works", {
  x <- Variable(3)
  obj <- Minimize(x[1])
  expect_s3_class(obj, "CVXR::Minimize")
  expect_true(is_dcp(obj))
})

# ── Objective with t(x) %*% y ────────────────────────────────────────

## @cvxpy NONE
test_that("Minimize(t(x) %*% c) for inner product with constant", {
  x <- Variable(3)
  c_vec <- Constant(matrix(c(1, 2, 3), 3, 1))
  obj <- Minimize(t(x) %*% c_vec)
  expect_s3_class(obj, "CVXR::Minimize")
  expect_equal(obj@args[[1L]]@shape, c(1L, 1L))
  expect_true(is_dcp(obj))
})

# ── End-to-end smoke test from plan ──────────────────────────────────

## @cvxpy NONE
test_that("Phase 2 end-to-end smoke test", {
  x <- Variable(3)
  y <- Variable(3)
  c_vec <- Constant(matrix(c(1, 1, 1), 3, 1))
  expr <- x + 2 * y                  # AddExpression, shape (3,1)
  neg_expr <- -x                     # NegExpression
  constr <- list(x >= 0, x <= 10, x[1] == 5)
  ## t(x) %*% c is affine (constant on one side), hence DCP for Minimize
  obj <- Minimize(t(x) %*% c_vec)

  expect_equal(expr@shape, c(3L, 1L))
  expect_s3_class(neg_expr, "CVXR::NegExpression")
  expect_s3_class(constr[[1]], "CVXR::Inequality")
  expect_s3_class(constr[[2]], "CVXR::Inequality")
  expect_s3_class(constr[[3]], "CVXR::Equality")
  expect_true(is_dcp(obj))
  expect_true(all(vapply(constr, is_dcp, logical(1))))
})

# ── Logical indexing ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("Logical indexing on vector", {
  x <- Variable(3)
  idx <- x[c(TRUE, FALSE, TRUE)]
  expect_equal(idx@shape, c(2L, 1L))
})
