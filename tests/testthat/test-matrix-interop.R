## Tests for Matrix package interoperability with CVXR
##
## Matrix objects (dgCMatrix, dgeMatrix, etc.) are S4 classes. S4 dispatch
## preempts S3/S7, so raw Matrix objects cannot be used directly with CVXR
## operators. Use as_cvxr_expr() to wrap Matrix objects as CVXR Constants.
##
## Base R matrix/numeric objects work natively (no wrapping needed).

skip_if_not_installed("Matrix")

# ═══════════════════════════════════════════════════════════════════════
# Setup
# ═══════════════════════════════════════════════════════════════════════

A_sparse <- Matrix::sparseMatrix(
  i = c(1, 2, 3, 1, 2, 3),
  j = c(1, 1, 1, 2, 2, 2),
  x = c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
  dims = c(3, 3)
)
A_dense <- Matrix::Matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
A_diag <- Matrix::Diagonal(3, x = c(1.0, 2.0, 3.0))


# ═══════════════════════════════════════════════════════════════════════
# 1. as_cvxr_expr() wraps Matrix objects as Constant
# ═══════════════════════════════════════════════════════════════════════

test_that("as_cvxr_expr wraps sparse Matrix", {
  expr <- as_cvxr_expr(A_sparse)
  expect_true(S7_inherits(expr, Constant))
  expect_equal(expr@shape, c(3L, 3L))
})

test_that("as_cvxr_expr wraps dense Matrix", {
  expr <- as_cvxr_expr(A_dense)
  expect_true(S7_inherits(expr, Constant))
})

test_that("as_cvxr_expr wraps diagonal Matrix", {
  expr <- as_cvxr_expr(A_diag)
  expect_true(S7_inherits(expr, Constant))
})

test_that("as_cvxr_expr wraps sparseVector", {
  sv <- Matrix::sparseVector(x = c(1.0, 2.0), i = c(1L, 3L), length = 3L)
  expr <- as_cvxr_expr(sv)
  expect_true(S7_inherits(expr, Constant))
  expect_equal(expr@shape, c(3L, 1L))
})

test_that("as_cvxr_expr passes through Expression unchanged", {
  x <- Variable(3)
  expect_identical(as_cvxr_expr(x), x)
})


# ═══════════════════════════════════════════════════════════════════════
# 2. Arithmetic with wrapped Matrix objects
# ═══════════════════════════════════════════════════════════════════════

test_that("wrapped sparse + Variable", {
  x <- Variable(c(3, 3))
  expr <- as_cvxr_expr(A_sparse) + x
  expect_true(S7_inherits(expr, Expression))
  expect_equal(expr@shape, c(3L, 3L))
})

test_that("Variable - wrapped sparse", {
  x <- Variable(c(3, 3))
  expr <- x - as_cvxr_expr(A_sparse)
  expect_true(S7_inherits(expr, Expression))
})

test_that("wrapped sparse * Variable", {
  x <- Variable(c(3, 3))
  expr <- as_cvxr_expr(A_sparse) * x
  expect_true(S7_inherits(expr, Expression))
})

test_that("wrapped dense + Variable", {
  x <- Variable(c(3, 3))
  expr <- as_cvxr_expr(A_dense) + x
  expect_true(S7_inherits(expr, Expression))
})

test_that("wrapped diagonal - Variable", {
  x <- Variable(c(3, 3))
  expr <- as_cvxr_expr(A_diag) - x
  expect_true(S7_inherits(expr, Expression))
})


# ═══════════════════════════════════════════════════════════════════════
# 3. Comparison constraints with wrapped Matrix objects
# ═══════════════════════════════════════════════════════════════════════

test_that("wrapped sparse >= Variable", {
  x <- Variable(c(3, 3))
  constr <- as_cvxr_expr(A_sparse) >= x
  expect_true(S7_inherits(constr, Inequality))
})

test_that("Variable <= wrapped sparse", {
  x <- Variable(c(3, 3))
  constr <- x <= as_cvxr_expr(A_sparse)
  expect_true(S7_inherits(constr, Inequality))
})

test_that("wrapped sparse == Variable", {
  x <- Variable(c(3, 3))
  constr <- as_cvxr_expr(A_sparse) == x
  expect_true(S7_inherits(constr, Equality))
})


# ═══════════════════════════════════════════════════════════════════════
# 4. Matrix multiply with wrapped objects
# ═══════════════════════════════════════════════════════════════════════

test_that("wrapped sparse %*% Variable", {
  x <- Variable(3)
  expr <- as_cvxr_expr(A_sparse) %*% x
  expect_true(S7_inherits(expr, Expression))
  expect_equal(expr@shape, c(3L, 1L))
})

test_that("Variable %*% wrapped sparse", {
  x <- Variable(c(1, 3))
  expr <- x %*% as_cvxr_expr(A_sparse)
  expect_true(S7_inherits(expr, Expression))
  expect_equal(expr@shape, c(1L, 3L))
})

test_that("wrapped dense %*% Variable", {
  x <- Variable(3)
  expr <- as_cvxr_expr(A_dense) %*% x
  expect_true(S7_inherits(expr, Expression))
  expect_equal(expr@shape, c(3L, 1L))
})

test_that("wrapped diagonal %*% Variable", {
  x <- Variable(3)
  expr <- as_cvxr_expr(A_diag) %*% x
  expect_true(S7_inherits(expr, Expression))
  expect_equal(expr@shape, c(3L, 1L))
})

test_that("t(wrapped sparse) %*% Variable (crossprod equivalent)", {
  x <- Variable(3)
  expr <- t(as_cvxr_expr(A_sparse)) %*% x
  expect_true(S7_inherits(expr, Expression))
  expect_equal(expr@shape, c(3L, 1L))
})


# ═══════════════════════════════════════════════════════════════════════
# 5. Base matrix works natively (no wrapping needed)
# ═══════════════════════════════════════════════════════════════════════

test_that("base matrix %*% Variable (native)", {
  x <- Variable(3)
  D <- matrix(1:9, 3, 3)
  expr <- D %*% x
  expect_true(S7_inherits(expr, Expression))
  expect_equal(expr@shape, c(3L, 1L))
})

test_that("base matrix + Variable (native)", {
  x <- Variable(c(3, 3))
  D <- matrix(1:9, 3, 3)
  expr <- D + x
  expect_true(S7_inherits(expr, Expression))
})

test_that("numeric scalar * Variable (native)", {
  x <- Variable(3)
  expr <- 2.0 * x
  expect_true(S7_inherits(expr, Expression))
})


# ═══════════════════════════════════════════════════════════════════════
# 6. Pure Matrix operations are NOT broken
# ═══════════════════════════════════════════════════════════════════════

test_that("pure Matrix: sparse %*% sparse unaffected", {
  result <- A_sparse %*% A_sparse
  expect_true(inherits(result, "Matrix"))
})

test_that("pure Matrix: sparse + sparse unaffected", {
  result <- A_sparse + A_sparse
  expect_true(inherits(result, "Matrix"))
})

test_that("pure Matrix: sparse >= sparse unaffected", {
  result <- A_sparse >= A_sparse
  expect_true(inherits(result, "Matrix"))
})


# ═══════════════════════════════════════════════════════════════════════
# 7. End-to-end solve with sparse data
# ═══════════════════════════════════════════════════════════════════════

test_that("end-to-end: sparse regression", {
  set.seed(42)
  n <- 10; p <- 3
  A_sp <- Matrix::sparseMatrix(
    i = rep(1:n, each = p),
    j = rep(1:p, times = n),
    x = rnorm(n * p),
    dims = c(n, p)
  )
  b <- rnorm(n)
  beta <- Variable(p)

  prob <- Problem(Minimize(sum_squares(as_cvxr_expr(A_sp) %*% beta - b)))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")

  ## Compare with dense equivalent
  beta2 <- Variable(p)
  prob2 <- Problem(Minimize(sum_squares(as.matrix(A_sp) %*% beta2 - b)))
  result2 <- psolve(prob2)
  expect_equal(value(beta), value(beta2), tolerance = 1e-5)
})

test_that("end-to-end: sparse constraint", {
  set.seed(123)
  n <- 5
  A_sp <- Matrix::sparseMatrix(
    i = 1:n, j = 1:n, x = rep(1.0, n), dims = c(n, n)
  )
  x <- Variable(n)

  prob <- Problem(Minimize(sum_squares(x)),
                  list(as_cvxr_expr(A_sp) %*% x >= 1))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_true(all(value(x) >= 1 - 1e-5))
})

test_that("end-to-end: mixed sparse + dense", {
  set.seed(7)
  n <- 5
  A_sp <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1.0, dims = c(n, n))
  D <- matrix(rnorm(n * n), n, n)
  x <- Variable(n)
  y <- Variable(n)

  prob <- Problem(
    Minimize(sum_squares(as_cvxr_expr(A_sp) %*% x) + sum_squares(D %*% y)),
    list(x + y == 1)
  )
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
})
