a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, naem = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

assert_expression <- function(expr, size) {
  expect_true(is(expr, "Expression") || is(expr, "Constraint"))
  expect_equal(size(expr), size)
}

test_that("Test R vectors", {
  # Vector
  v <- 1:2
  w <- Variable(1, 2, name = "w")
  assert_expression(w + v, c(1, 2))
  assert_expression(v + w, c(1, 2))
  assert_expression(w - v, c(1, 2))
  assert_expression(v - w, c(1, 2))
  assert_expression(w <= v, c(1, 2))
  assert_expression(v <= w, c(1, 2))
  assert_expression(w == v, c(1, 2))
  assert_expression(v == w, c(1, 2))
  
  # Matrix
  Amat <- matrix(1:4, nrow = 4, ncol = 1)
  assert_expression(Amat * w, c(4, 1))
})

test_that("Test R matrices", {
  # Vector
  v <- matrix(1:2, nrow = 2, ncol = 1)
  assert_expression(x + v, c(2, 1))
  assert_expression(v + v + x, c(2, 1))
  assert_expression(x - v, c(2, 1))
  assert_expression(v - v - x, c(2, 1))
  assert_expression(x <= v, c(2, 1))
  assert_expression(v <= x, c(2, 1))
  assert_expression(x == v, c(2, 1))
  assert_expression(v == x, c(2, 1))
  
  # Matrix
  Amat <- matrix(1:8, nrow = 4, ncol = 2)
  assert_expression(Amat*x, c(4, 1))
  assert_expression((t(Amat) %*% Amat) * x, c(2, 1))
  
  # PSD inequalities
  Amat <- matrix(rep(1, 4), nrow = 2, ncol = 2)
  # assert_expression(Amat << A, c(2, 2))
  # assert_expression(Amat >> A, c(2, 2))
})

test_that("Test R scalars", {
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
  
  # PSD inequalities
  # assert_expression(v << A, c(2, 2))
  # assert_expression(v >> A, c(2, 2))
})

test_that("Test sparseMatrix objects from the Matrix library", {
  # Constants
  A <- matrix(1:8, nrow = 4, ncol = 2)
  A <- Matrix(A, sparse = TRUE)
  A <- Matrix(i = 1:2, j = 1:2, x = rep(1, 2))
  Aidx <- A[1,]
  Aidx <- A[1:2,]
  # expect_equal(dim(Aidx), c(1, 2))
  # expect_equal(Aidx[1,1], 1)
  # expect_equal(Aidx[1,2], 0)
  
  # Linear ops
  var <- Variable(4, 2)
  A <- matrix(1:8, nrow = 4, ncol = 2)
  A <- Matrix(A, sparse = TRUE)
  B <- rbind(A, A)
  assert_expression(var + A, c(4, 2))
  assert_expression(A + var, c(4, 2))
  assert_expression(B * var, c(4, 2))
  assert_expression(var - A, c(4, 2))
  assert_expression(A - A - var, c(4, 2))
})
