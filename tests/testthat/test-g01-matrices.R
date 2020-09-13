context("test-g01-matrices")

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

assert_expression <- function(expr, dim) {
  skip_on_cran()
  expect_true(is(expr, "Expression") || is(expr, "Constraint"))
  expect_equal(dim(expr), dim)
}

test_that("Test R vectors", {
  skip_on_cran()
  ## Vector
  v <- 1:2
  assert_expression(x + v, c(2, 1))
  assert_expression(v + x, c(2, 1))
  assert_expression(x - v, c(2, 1))
  assert_expression(v - x, c(2, 1))
  assert_expression(x <= v, c(2, 1))
  assert_expression(v <= x, c(2, 1))
  assert_expression(x == v, c(2, 1))
  assert_expression(v == x, c(2, 1))

  ## Matrix
  Amat <- matrix(1:8, nrow = 4, ncol = 2)
  assert_expression(Amat %*% x, c(4, 1))

  ## PSD inequalities
  Amat <- matrix(1, nrow = 2, ncol = 2)
  assert_expression(Amat %<<% A, c(2, 2))
  assert_expression(Amat %>>% A, c(2, 2))
})

test_that("Test R matrices", {
  skip_on_cran()
  ## Vector
  v <- matrix(1:2, nrow = 2, ncol = 1)
  assert_expression(x + v, c(2, 1))
  assert_expression(v + v + x, c(2, 1))
  assert_expression(x - v, c(2, 1))
  assert_expression(v - v - x, c(2, 1))
  assert_expression(x <= v, c(2, 1))
  assert_expression(v <= x, c(2, 1))
  assert_expression(x == v, c(2, 1))
  assert_expression(v == x, c(2, 1))

  ## Matrix
  Amat <- matrix(1:8, nrow = 4, ncol = 2)
  assert_expression(Amat %*% x, c(4, 1))
  assert_expression((t(Amat) %*% Amat) %*% x, c(2, 1))

  ## PSD inequalities
  Amat <- matrix(rep(1, 4), nrow = 2, ncol = 2)
  assert_expression(Amat %<<% A, c(2, 2))
  assert_expression(Amat %>>% A, c(2, 2))
})

test_that("Test R scalars", {
  skip_on_cran()
  v <- 2.0
  assert_expression(x + v, c(2, 1))
  assert_expression(v + x, c(2, 1))
  assert_expression(v * x, c(2, 1))
  assert_expression(x - v, c(2, 1))
  assert_expression(v - v - x, c(2, 1))
  assert_expression(x <= v, c(2, 1))
  assert_expression(v <= x, c(2, 1))
  assert_expression(x == v, c(2, 1))
  assert_expression(v == x, c(2, 1))

  ## PSD inequalities
  assert_expression(v %<<% A, c(2, 2))
  assert_expression(v %>>% A, c(2, 2))
})

test_that("Test sparseMatrix objects from the Matrix library", {
  skip_on_cran()
    ##require(Matrix)

    ## Constants
    A <- matrix(1:8, nrow = 4, ncol = 2)
    A <- Matrix::Matrix(A, sparse = TRUE)
    A <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = rep(1, 2))
    Aidx <- Matrix::Matrix(A[1,], sparse = TRUE)
    expect_equal(dim(Aidx), c(2, 1))
    expect_equal(Aidx[1,1], 1)
    expect_equal(Aidx[2,1], 0)

    ## Linear ops
    var <- Variable(4, 2)
    A <- matrix(1:8, nrow = 4, ncol = 2)
    A <- Matrix::Matrix(A, sparse = TRUE)
    B <- cbind(A, A)
    assert_expression(var + A, c(4, 2))
    assert_expression(A + var, c(4, 2))
    assert_expression(B %*% var, c(4, 2))
    assert_expression(var - A, c(4, 2))
    assert_expression(A - A - var, c(4, 2))
})
