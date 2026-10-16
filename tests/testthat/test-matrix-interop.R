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


# ═══════════════════════════════════════════════════════════════════════
# 8. Constant() normalizes its value (intf_convert)
# ═══════════════════════════════════════════════════════════════════════
## CVXPY's Constant.__init__ (constant.py:51-54) sends every sparse input
## through SparseArrayInterface.const_to_matrix (sparse_matrix_interface.py:
## 40-51): any 2-D sparse format becomes a double csc_array, bool and int are
## cast to double, explicit zeros are kept. Dense inputs go through
## NDArrayInterface.const_to_matrix the same way (shape kept, cast to double).
## CVXR's intf_convert is the twin; the rest of CVXR assumes a general double
## dgCMatrix (the C++ bridge reads @i/@p/@x, intf_sign reads @x, the quad_form
## P tensor is built from the stored value). Matrix's class zoo is R-specific,
## so these tests have no upstream twin. Expected solve values were measured
## against CVXPY 1.9.3 with the scipy formats (dia, csr, coo, bool csc).

## General, symmetric, triangular and diagonal 3 x 3 matrices, all nonsingular,
## built so Matrix() does not re-class them.
.cvt_M <- matrix(c(2, 1, 0, 0, 2, 1, 1, 0, 2), 3, 3)   # general
.cvt_S <- .cvt_M + t(.cvt_M)                            # symmetric
.cvt_L <- .cvt_M * lower.tri(.cvt_M, diag = TRUE)       # lower triangular
.cvt_d <- c(1, 2, 3)

## Every sparseMatrix class a user might hand in, keyed by the class the
## constructor is expected to produce (checked, so a Matrix change is visible).
.cvt_sparse_inputs <- function() {
  dgC <- as(Matrix::Matrix(.cvt_M, sparse = TRUE), "CsparseMatrix")
  list(
    dgCMatrix = dgC,
    dgTMatrix = as(dgC, "TsparseMatrix"),
    dgRMatrix = as(dgC, "RsparseMatrix"),
    dsCMatrix = Matrix::Matrix(.cvt_S, sparse = TRUE),
    dtCMatrix = Matrix::Matrix(.cvt_L, sparse = TRUE),
    ddiMatrix = Matrix::Diagonal(x = .cvt_d),
    ngCMatrix = as(dgC, "nMatrix"),
    lgCMatrix = as(dgC, "lMatrix"),
    ltCMatrix = as(Matrix::Matrix(.cvt_L, sparse = TRUE), "lMatrix")
  )
}

## @cvxpy NONE
test_that("Constant: every sparseMatrix class becomes a general double dgCMatrix", {
  require_solver("CLARABEL")
  ## Deprecated Matrix coercions (as(<lMatrix>, "dgCMatrix")) become errors.
  old_opt <- options(Matrix.warnDeprecatedCoerce = 2L)
  on.exit(options(old_opt), add = TRUE)
  inputs <- .cvt_sparse_inputs()
  for (nm in names(inputs)) {
    A <- inputs[[nm]]
    expect_s4_class(A, nm)
    dense <- unname(as.matrix(A) + 0)           # pattern/logical -> 0/1 numeric
    cst <- Constant(A)
    v <- value(cst)
    expect_s4_class(v, "dgCMatrix")
    expect_type(v@x, "double")
    expect_equal(as.matrix(v), dense, info = nm)
    expect_true(cst@.sparse, info = nm)
    ## The constant must also solve correctly through the compiled path.
    y <- Variable(3)
    prob <- Problem(Minimize(sum(y)), list(as_cvxr_expr(A) %*% y == 1))
    val <- psolve(prob, solver = "CLARABEL")
    expect_equal(status(prob), "optimal", info = nm)
    expect_equal(val, sum(solve(dense, rep(1, 3))), tolerance = 1e-6, info = nm)
  }
})

## @cvxpy NONE
test_that("Constant: a unit diagonal materializes its ones (sign and DCP)", {
  I3 <- Matrix::Diagonal(3)                      # ddiMatrix with diag = "U", empty @x
  cst <- Constant(I3)
  v <- value(cst)
  expect_s4_class(v, "dgCMatrix")
  expect_equal(length(v@x), 3L)
  expect_equal(as.matrix(v), diag(3))
  expect_true(is_nonneg(cst))
  expect_false(is_nonpos(cst))                   # was TRUE: empty @x read as all zeros
  ## CVXPY 1.9.3: Problem(Maximize(sum(eye(3) @ square(y)))).is_dcp() is False.
  y <- Variable(3)
  expect_false(is_dcp(Problem(Maximize(sum(as_cvxr_expr(I3) %*% square(y))))))
})

## @cvxpy NONE
test_that("Constant: a sparseVector becomes a double sparse column", {
  sv_l <- Matrix::sparseVector(x = c(TRUE, TRUE), i = c(1L, 3L), length = 4L)
  cst <- Constant(sv_l)
  v <- value(cst)
  expect_equal(cst@shape, c(4L, 1L))
  expect_s4_class(v, "dgCMatrix")
  expect_type(v@x, "double")
  expect_equal(as.numeric(as.matrix(v)), c(1, 0, 1, 0))
  expect_true(cst@.sparse)
  sv_d <- Matrix::sparseVector(x = c(1.5, 2), i = c(1L, 3L), length = 4L)
  expect_equal(as.numeric(as.matrix(value(Constant(sv_d)))), c(1.5, 0, 2, 0))
})

## @cvxpy test_expressions.py::TestExpressions::test_constants
test_that("Constant: dense inputs keep their shape and are stored as double", {
  ## A logical matrix lost its dim (as.numeric drops it) and became a column.
  M <- matrix(c(TRUE, FALSE, TRUE, TRUE), 2, 2)
  cst <- Constant(M)
  expect_equal(cst@shape, c(2L, 2L))
  expect_type(value(cst), "double")
  expect_equal(value(cst), matrix(c(1, 0, 1, 1), 2, 2))
  x <- Variable(c(2, 2))
  expect_equal((x + M)@shape, c(2L, 2L))
  ## Integer storage is cast to double, as CVXPY casts int to float64.
  expect_type(value(Constant(matrix(1:4, 2, 2))), "double")
  expect_type(value(Constant(1:3)), "double")
  ## Scalars, vectors and 1-d arrays become a column; complex is preserved.
  expect_equal(Constant(2L)@shape, c(1L, 1L))
  expect_equal(Constant(c(TRUE, FALSE, TRUE))@shape, c(3L, 1L))
  expect_equal(Constant(array(1:3))@shape, c(3L, 1L))
  cz <- Constant(matrix(c(1i, 2), 1, 2))
  expect_equal(cz@shape, c(1L, 2L))
  expect_type(value(cz), "complex")
  ## Character never becomes a constant: a character matrix used to turn
  ## into NAs with a warning (storage.mode<-), not an error.
  expect_error(Constant("a"), class = "rlang_error")
  expect_error(Constant(matrix("a", 1, 1)), class = "rlang_error")
  expect_error(Constant(array(1:8, c(2, 2, 2))), class = "rlang_error")
})

## @cvxpy NONE
test_that("Constant: explicit zeros in a dgCMatrix are kept", {
  ## CVXPY keeps them (csc_array(value) does not eliminate zeros; nnz stays 2).
  A0 <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L), x = c(0, 1), dims = c(2L, 2L))
  expect_equal(length(A0@x), 2L)
  expect_equal(length(value(Constant(A0))@x), 2L)
})
