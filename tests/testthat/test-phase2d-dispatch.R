## Tests for Phase 2D: Arithmetic Operator Dispatch
## Verifies that R operators (+, -, *, /, %*%) produce the correct
## CVXR expression types with correct shapes and values.

# ═══════════════════════════════════════════════════════════════════
# Addition (+)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Addition: Variable + Variable → AddExpression", {
  x <- Variable(3); y <- Variable(3)
  result <- x + y
  expect_true(S7_inherits(result, AddExpression))
  expect_equal(result@shape, c(3L, 1L))
  expect_length(result@args, 2L)
})

## @cvxpy NONE
test_that("Addition: Variable + scalar promotes", {
  x <- Variable(3)
  result <- x + 1
  expect_true(S7_inherits(result, AddExpression))
  expect_equal(result@shape, c(3L, 1L))
  expect_true(is_constant(result@args[[2L]]))
})

## @cvxpy NONE
test_that("Addition: scalar + Variable (reverse dispatch)", {
  x <- Variable(3)
  result <- 1 + x
  expect_true(S7_inherits(result, AddExpression))
  expect_equal(result@shape, c(3L, 1L))
  expect_true(is_constant(result@args[[1L]]))
})

## @cvxpy NONE
test_that("Addition: matrix + Variable", {
  x <- Variable(3)
  A <- matrix(1:3, 3, 1)
  result <- A + x
  expect_true(S7_inherits(result, AddExpression))
  expect_equal(result@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Addition: Constant + Constant", {
  a <- Constant(5); b <- Constant(3)
  result <- a + b
  expect_true(S7_inherits(result, AddExpression))
  expect_equal(drop(value(result)), 8)
})

# ═══════════════════════════════════════════════════════════════════
# Unary + and -
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Unary minus: -Variable → NegExpression", {
  x <- Variable(3)
  result <- -x
  expect_true(S7_inherits(result, NegExpression))
  expect_equal(result@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Unary minus: -Constant evaluates correctly", {
  a <- Constant(matrix(1:3, 3, 1))
  result <- -a
  expect_equal(as.matrix(value(result)), matrix(c(-1L, -2L, -3L), 3, 1))
})

## @cvxpy NONE
test_that("Unary plus: +Variable → identity", {
  x <- Variable(3)
  result <- +x
  expect_true(S7_inherits(result, Expression))
  expect_identical(result@id, x@id)
})

## @cvxpy NONE
test_that("Double negation: -(-x) → NegExpression(NegExpression)", {
  x <- Variable(3)
  result <- -(-x)
  expect_true(S7_inherits(result, NegExpression))
  expect_true(S7_inherits(result@args[[1L]], NegExpression))
})

# ═══════════════════════════════════════════════════════════════════
# Subtraction (-)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Subtraction: Variable - Variable → Add(x, Neg(y))", {
  x <- Variable(3); y <- Variable(3)
  result <- x - y
  expect_true(S7_inherits(result, AddExpression))
  expect_length(result@args, 2L)
  expect_true(S7_inherits(result@args[[2L]], NegExpression))
})

## @cvxpy NONE
test_that("Subtraction: Variable - scalar", {
  x <- Variable(3)
  result <- x - 1
  expect_true(S7_inherits(result, AddExpression))
  expect_true(S7_inherits(result@args[[2L]], NegExpression))
})

## @cvxpy NONE
test_that("Subtraction: scalar - Variable", {
  x <- Variable(3)
  result <- 1 - x
  expect_true(S7_inherits(result, AddExpression))
  ## 1 - x → Add(Constant(1), Neg(x))
  expect_true(is_constant(result@args[[1L]]))
  expect_true(S7_inherits(result@args[[2L]], NegExpression))
})

# ═══════════════════════════════════════════════════════════════════
# Elementwise multiplication (*)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Multiplication: Variable * scalar → Multiply", {
  x <- Variable(3)
  result <- x * 2
  expect_true(S7_inherits(result, Multiply))
  expect_equal(result@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Multiplication: scalar * Variable → Multiply", {
  x <- Variable(3)
  result <- 2 * x
  expect_true(S7_inherits(result, Multiply))
  expect_equal(result@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Multiplication: Variable * Variable → Multiply (not DCP but constructs)", {
  x <- Variable(3); y <- Variable(3)
  result <- x * y
  expect_true(S7_inherits(result, Multiply))
  ## Not DCP (product of two non-constants), but construction should succeed
  expect_false(is_dcp(result))
})

# ═══════════════════════════════════════════════════════════════════
# Division (/)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Division: Variable / scalar → DivExpression", {
  x <- Variable(3)
  result <- x / 2
  expect_true(S7_inherits(result, DivExpression))
  expect_equal(result@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Division: scalar / Variable → DivExpression (not DCP)", {
  x <- Variable(3)
  result <- 1 / x
  expect_true(S7_inherits(result, DivExpression))
  ## Not DCP because denominator is not constant
  expect_false(is_dcp(result))
})

# ═══════════════════════════════════════════════════════════════════
# Matrix multiplication (%*%)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Matmul: matrix %*% Variable → MulExpression", {
  A <- matrix(1:6, 3, 2)
  x <- Variable(c(2, 1))
  result <- A %*% x
  expect_true(S7_inherits(result, MulExpression))
  expect_equal(result@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Matmul: Variable %*% matrix → MulExpression", {
  x <- Variable(c(3, 4))
  A <- matrix(1:8, 4, 2)
  result <- x %*% A
  expect_true(S7_inherits(result, MulExpression))
  expect_equal(result@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("Matmul: Variable %*% Variable → MulExpression", {
  x <- Variable(c(3, 4)); y <- Variable(c(4, 2))
  result <- x %*% y
  expect_true(S7_inherits(result, MulExpression))
  expect_equal(result@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("Matmul: scalar %*% Variable → MulExpression (scalar mult)", {
  x <- Variable(3)
  result <- 5 %*% x
  expect_true(S7_inherits(result, MulExpression))
  expect_equal(result@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Matmul: Constant %*% Constant evaluates", {
  A <- Constant(matrix(1:6, 2, 3))
  b <- Constant(matrix(1:3, 3, 1))
  result <- A %*% b
  expected <- matrix(1:6, 2, 3) %*% matrix(1:3, 3, 1)
  expect_equal(as.matrix(value(result)), expected)
})

# ═══════════════════════════════════════════════════════════════════
# Nested flattening (AddExpression)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Flattening: (x + y) + z → 3 args", {
  x <- Variable(3); y <- Variable(3); z <- Variable(3)
  result <- (x + y) + z
  expect_true(S7_inherits(result, AddExpression))
  expect_length(result@args, 3L)
})

## @cvxpy NONE
test_that("Flattening: x + (y + z) → 3 args", {
  x <- Variable(3); y <- Variable(3); z <- Variable(3)
  result <- x + (y + z)
  expect_true(S7_inherits(result, AddExpression))
  expect_length(result@args, 3L)
})

## @cvxpy NONE
test_that("Flattening: (x + y) + (z + w) → 4 args", {
  x <- Variable(3); y <- Variable(3)
  z <- Variable(3); w <- Variable(3)
  result <- (x + y) + (z + w)
  expect_true(S7_inherits(result, AddExpression))
  expect_length(result@args, 4L)
})

## @cvxpy NONE
test_that("Subtraction does NOT flatten: (x - y) - z structure", {
  x <- Variable(3); y <- Variable(3); z <- Variable(3)
  result <- (x - y) - z
  ## x - y → Add(x, Neg(y)); then - z → Add(Add(x,Neg(y)), Neg(z))
  ## BUT AddExpression flattens: Add(x, Neg(y), Neg(z)) — 3 args
  expect_true(S7_inherits(result, AddExpression))
  expect_length(result@args, 3L)
  expect_true(S7_inherits(result@args[[2L]], NegExpression))
  expect_true(S7_inherits(result@args[[3L]], NegExpression))
})

# ═══════════════════════════════════════════════════════════════════
# Numeric evaluation with operators
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Numeric evaluation: vector arithmetic", {
  a <- Constant(matrix(1:3, 3, 1))
  b <- Constant(matrix(4:6, 3, 1))

  ## Addition
  expect_equal(as.matrix(value(a + b)), matrix(c(5L, 7L, 9L), 3, 1))

  ## Subtraction
  expect_equal(as.matrix(value(a - b)), matrix(c(-3L, -3L, -3L), 3, 1))

  ## Elementwise multiplication
  expect_equal(as.matrix(value(a * b)), matrix(c(4L, 10L, 18L), 3, 1))

  ## Negation
  expect_equal(as.matrix(value(-a)), matrix(c(-1L, -2L, -3L), 3, 1))

  ## Scalar multiplication
  expect_equal(as.matrix(value(a * 2)), matrix(c(2L, 4L, 6L), 3, 1))

  ## Division
  expect_equal(as.matrix(value(b / 2)), matrix(c(2.0, 2.5, 3.0), 3, 1))
})

## @cvxpy NONE
test_that("Numeric evaluation: matrix multiplication", {
  A <- Constant(matrix(1:6, 2, 3))
  b <- Constant(matrix(1:3, 3, 1))
  expected <- matrix(1:6, 2, 3) %*% matrix(1:3, 3, 1)
  expect_equal(as.matrix(value(A %*% b)), expected)
})

## @cvxpy NONE
test_that("Numeric evaluation: scalar arithmetic", {
  a <- Constant(10); b <- Constant(3)
  expect_equal(drop(value(a + b)), 13)
  expect_equal(drop(value(a - b)), 7)
  expect_equal(drop(value(a * b)), 30)
  expect_equal(drop(value(a / b)), 10 / 3)
})

# ═══════════════════════════════════════════════════════════════════
# Complex expressions
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Complex: 2*x + 3*y builds correctly", {
  x <- Variable(3); y <- Variable(3)
  result <- 2 * x + 3 * y
  expect_true(S7_inherits(result, AddExpression))
  expect_equal(result@shape, c(3L, 1L))
  expect_length(result@args, 2L)
  expect_true(S7_inherits(result@args[[1L]], Multiply))
  expect_true(S7_inherits(result@args[[2L]], Multiply))
})

## @cvxpy NONE
test_that("Complex: (x + y) / 2 builds correctly", {
  x <- Variable(3); y <- Variable(3)
  result <- (x + y) / 2
  expect_true(S7_inherits(result, DivExpression))
  expect_true(S7_inherits(result@args[[1L]], AddExpression))
})

## @cvxpy NONE
test_that("Complex: x - 2*y builds correctly", {
  x <- Variable(3); y <- Variable(3)
  result <- x - 2 * y
  expect_true(S7_inherits(result, AddExpression))
  expect_length(result@args, 2L)
  ## Second arg is NegExpression(Multiply(2, y))
  expect_true(S7_inherits(result@args[[2L]], NegExpression))
  expect_true(S7_inherits(result@args[[2L]]@args[[1L]], Multiply))
})

## @cvxpy NONE
test_that("Complex: numeric evaluation of 2*a + 3*b - c", {
  a <- Constant(matrix(c(1, 2, 3), 3, 1))
  b <- Constant(matrix(c(4, 5, 6), 3, 1))
  cc <- Constant(matrix(c(10, 10, 10), 3, 1))
  result <- 2 * a + 3 * b - cc
  ## 2*[1,2,3] + 3*[4,5,6] - [10,10,10] = [2,4,6] + [12,15,18] - [10,10,10]
  ## = [4, 9, 14]
  expect_equal(as.matrix(value(result)), matrix(c(4, 9, 14), 3, 1))
})

# ═══════════════════════════════════════════════════════════════════
# DCP properties through operators
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DCP: affine expression is_convex, is_concave, is_affine", {
  x <- Variable(3)
  result <- 2 * x + 1
  expect_true(is_convex(result))
  expect_true(is_concave(result))
  expect_true(is_affine(result))
})

## @cvxpy NONE
test_that("DCP: -x is affine", {
  x <- Variable(3)
  result <- -x
  expect_true(is_affine(result))
})

## @cvxpy NONE
test_that("DCP: x / constant is affine", {
  x <- Variable(3)
  result <- x / 2
  expect_true(is_affine(result))
})

## @cvxpy NONE
test_that("DCP: x + y is affine", {
  x <- Variable(3); y <- Variable(3)
  result <- x + y
  expect_true(is_affine(result))
})

## @cvxpy NONE
test_that("DCP: x - y is affine", {
  x <- Variable(3); y <- Variable(3)
  result <- x - y
  expect_true(is_affine(result))
})

# ═══════════════════════════════════════════════════════════════════
# expr_name via operators
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("expr_name: addition", {
  x <- Variable(3, name = "x"); y <- Variable(3, name = "y")
  expect_equal(expr_name(x + y), "x + y")
})

## @cvxpy NONE
test_that("expr_name: negation", {
  x <- Variable(3, name = "x")
  expect_equal(expr_name(-x), "-x")
})

## @cvxpy NONE
test_that("expr_name: subtraction", {
  x <- Variable(3, name = "x"); y <- Variable(3, name = "y")
  ## x - y → Add(x, Neg(y)) → "x + -y"
  result <- x - y
  expect_equal(expr_name(result), "x + -y")
})

# ═══════════════════════════════════════════════════════════════════
# Interaction with Matrix package
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Matrix objects: use as.matrix() for operator interop", {
  skip_if_not_installed("Matrix")
  ## Matrix package S4 dispatch intercepts Ops before our S3 handler,
  ## for BOTH directions when a Matrix object is involved.
  ## Workaround: coerce Matrix to base matrix via as.matrix().
  x <- Variable(3)
  A <- Matrix::Matrix(1:3, 3, 1)

  ## as.matrix(Matrix) + Expression works
  result <- as.matrix(A) + x
  expect_true(S7_inherits(result, AddExpression))
  expect_equal(result@shape, c(3L, 1L))

  ## Expression + as.matrix(Matrix) works
  result2 <- x + as.matrix(A)
  expect_true(S7_inherits(result2, AddExpression))
  expect_equal(result2@shape, c(3L, 1L))

  ## as.matrix(Matrix) * Expression works
  result3 <- as.matrix(A) * x
  expect_true(S7_inherits(result3, Multiply))
})

# ═══════════════════════════════════════════════════════════════════
# Edge cases
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Parameter in arithmetic", {
  p <- Parameter(3)
  x <- Variable(3)
  result <- p * x
  expect_true(S7_inherits(result, Multiply))
  ## Parameter is constant in the DCP sense
  expect_true(is_dcp(result))
})

## @cvxpy NONE
test_that("Chained operations: x + y + z + w", {
  x <- Variable(3); y <- Variable(3)
  z <- Variable(3); w <- Variable(3)
  result <- x + y + z + w
  ## Should flatten completely
  expect_true(S7_inherits(result, AddExpression))
  expect_length(result@args, 4L)
})

## @cvxpy NONE
test_that("Mixed shapes: scalar + vector broadcasts", {
  x <- Variable(3)
  result <- x + 5
  expect_equal(result@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Mixed shapes: scalar * matrix", {
  x <- Variable(c(3, 4))
  result <- 2 * x
  expect_true(S7_inherits(result, Multiply))
  expect_equal(result@shape, c(3L, 4L))
})

## @cvxpy NONE
test_that("Unsupported operator gives clear error", {
  x <- Variable(3); y <- Variable(3)
  ## ^: now routed to power() in Phase 3
  p <- x ^ 2
  expect_s3_class(p, "CVXR::Power")
  expect_equal(p@shape, c(3L, 1L))
  ## %%, %/% are in Ops group but not supported
  expect_error(x %% y, "not yet implemented")
  expect_error(x %/% y, "not yet implemented|not meaningful for factors")
})
