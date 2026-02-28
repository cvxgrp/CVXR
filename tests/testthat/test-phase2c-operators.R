## Tests for Phase 2C: Concrete Affine Operators
##
## Tests for: Promote, NegExpression, AddExpression,
##             BinaryOperator, MulExpression, Multiply, DivExpression

# ═══════════════════════════════════════════════════════════════════
# Promote
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Promote class exists and inherits correctly", {
  expect_true(inherits(Promote, "S7_class"))
  p <- Promote(Constant(1), c(3L, 1L))
  expect_true(S7_inherits(p, Promote))
  expect_true(S7_inherits(p, AffAtom))
  expect_true(S7_inherits(p, Atom))
})

## @cvxpy NONE
test_that("Promote shape is target shape", {
  p <- Promote(Constant(5), c(3L, 4L))
  expect_equal(p@shape, c(3L, 4L))
})

## @cvxpy NONE
test_that("Promote value broadcasts scalar", {
  p <- Promote(Constant(5), c(3L, 2L))
  v <- value(p)
  expect_equal(dim(v), c(3L, 2L))
  expect_true(all(v == 5))
})

## @cvxpy NONE
test_that("Promote is affine", {
  x <- Variable(1)
  p <- Promote(x, c(3L, 1L))
  expect_true(is_affine(p))
  expect_true(is_dcp(p))
})

## @cvxpy NONE
test_that("Promote sign propagates", {
  x <- Variable(1, nonneg = TRUE)
  p <- Promote(x, c(3L, 1L))
  expect_true(is_nonneg(p))
  expect_false(is_nonpos(p))
})

## @cvxpy NONE
test_that("Promote is_symmetric for square shapes", {
  p <- Promote(Constant(1), c(3L, 3L))
  expect_true(is_symmetric(p))
  p2 <- Promote(Constant(1), c(3L, 2L))
  expect_false(is_symmetric(p2))
})

## @cvxpy NONE
test_that("Promote get_data returns shape", {
  p <- Promote(Constant(1), c(3L, 4L))
  d <- get_data(p)
  expect_equal(d, list(c(3L, 4L)))
})

## @cvxpy NONE
test_that("Promote graph_implementation creates promote LinOp", {
  p <- Promote(Constant(1), c(3L, 1L))
  cf <- canonical_form(p)
  expect_true(is.list(cf))
  expect_equal(length(cf), 2L)
})

## @cvxpy NONE
test_that("cvxr_promote promotes scalar to vector", {
  x <- Variable(1)
  p <- cvxr_promote(x, c(3L, 1L))
  expect_true(S7_inherits(p, Promote))
  expect_equal(p@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("cvxr_promote returns expr unchanged if same shape", {
  x <- Variable(c(3, 1))
  p <- cvxr_promote(x, c(3L, 1L))
  expect_identical(p, x)
})

## @cvxpy NONE
test_that("cvxr_promote errors for non-scalar", {
  x <- Variable(c(3, 1))
  expect_error(cvxr_promote(x, c(4L, 1L)), "Only scalars")
})

## @cvxpy NONE
test_that("broadcast_args promotes scalar lhs", {
  x <- Constant(5)
  y <- Variable(c(3, 1))
  result <- broadcast_args(x, y)
  expect_true(S7_inherits(result[[1L]], Promote))
  expect_identical(result[[2L]], y)
})

## @cvxpy NONE
test_that("broadcast_args promotes scalar rhs", {
  x <- Variable(c(3, 1))
  y <- Constant(5)
  result <- broadcast_args(x, y)
  expect_identical(result[[1L]], x)
  expect_true(S7_inherits(result[[2L]], Promote))
})

## @cvxpy NONE
test_that("broadcast_args leaves same-shape args unchanged", {
  x <- Variable(c(3, 1))
  y <- Variable(c(3, 1))
  result <- broadcast_args(x, y)
  expect_identical(result[[1L]], x)
  expect_identical(result[[2L]], y)
})


# ═══════════════════════════════════════════════════════════════════
# NegExpression
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("NegExpression class exists and inherits correctly", {
  n <- NegExpression(Variable(3))
  expect_true(S7_inherits(n, NegExpression))
  expect_true(S7_inherits(n, AffAtom))
})

## @cvxpy NONE
test_that("NegExpression shape matches arg shape", {
  x <- Variable(c(3, 4))
  n <- NegExpression(x)
  expect_equal(n@shape, c(3L, 4L))
})

## @cvxpy NONE
test_that("NegExpression flips sign", {
  x <- Variable(1, nonneg = TRUE)
  n <- NegExpression(x)
  expect_false(is_nonneg(n))
  expect_true(is_nonpos(n))
})

## @cvxpy NONE
test_that("NegExpression of nonpositive becomes nonneg", {
  x <- Variable(1, nonpos = TRUE)
  n <- NegExpression(x)
  expect_true(is_nonneg(n))
  expect_false(is_nonpos(n))
})

## @cvxpy NONE
test_that("NegExpression is affine", {
  x <- Variable(3)
  n <- NegExpression(x)
  expect_true(is_affine(n))
  expect_true(is_dcp(n))
})

## @cvxpy NONE
test_that("NegExpression is decreasing", {
  n <- NegExpression(Variable(1))
  expect_false(is_incr(n, 1L))
  expect_true(is_decr(n, 1L))
})

## @cvxpy NONE
test_that("NegExpression value", {
  n <- NegExpression(Constant(5))
  expect_equal(drop(value(n)), -5)
})

## @cvxpy NONE
test_that("NegExpression of vector constant", {
  v <- matrix(c(1, 2, 3), ncol = 1)
  n <- NegExpression(Constant(v))
  expect_equal(as.numeric(value(n)), c(-1, -2, -3))
})

## @cvxpy NONE
test_that("NegExpression is_symmetric delegates to arg", {
  x <- Variable(c(3, 3), symmetric = TRUE)
  n <- NegExpression(x)
  expect_true(is_symmetric(n))
})

## @cvxpy NONE
test_that("NegExpression graph_implementation produces neg LinOp", {
  n <- NegExpression(Constant(5))
  cf <- canonical_form(n)
  expect_true(is.list(cf))
  expect_equal(length(cf), 2L)
})

## @cvxpy NONE
test_that("NegExpression expr_name has minus prefix", {
  x <- Variable(1, name = "x")
  n <- NegExpression(x)
  expect_true(grepl("^-", expr_name(n)))
})


# ═══════════════════════════════════════════════════════════════════
# AddExpression
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("AddExpression class exists and inherits correctly", {
  x <- Variable(3)
  y <- Variable(3)
  a <- AddExpression(list(x, y))
  expect_true(S7_inherits(a, AddExpression))
  expect_true(S7_inherits(a, AffAtom))
})

## @cvxpy NONE
test_that("AddExpression shape from sum_shapes", {
  x <- Variable(c(3, 1))
  y <- Variable(c(3, 1))
  a <- AddExpression(list(x, y))
  expect_equal(a@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("AddExpression broadcasts (1,3) + (3,3)", {
  x <- Variable(c(1, 3))
  y <- Variable(c(3, 3))
  a <- AddExpression(list(x, y))
  expect_equal(a@shape, c(3L, 3L))
})

## @cvxpy NONE
test_that("AddExpression flattens nested AddExpressions", {
  x <- Variable(3)
  y <- Variable(3)
  z <- Variable(3)
  inner <- AddExpression(list(x, y))
  outer <- AddExpression(list(inner, z))
  ## Should have 3 args, not 2 (inner was flattened)
  expect_equal(length(outer@args), 3L)
})

## @cvxpy NONE
test_that("AddExpression deep flattening", {
  a <- Variable(3)
  b <- Variable(3)
  c <- Variable(3)
  d <- Variable(3)
  ab <- AddExpression(list(a, b))
  cd <- AddExpression(list(c, d))
  abcd <- AddExpression(list(ab, cd))
  expect_equal(length(abcd@args), 4L)
})

## @cvxpy NONE
test_that("AddExpression is affine", {
  x <- Variable(3)
  y <- Variable(3)
  a <- AddExpression(list(x, y))
  expect_true(is_affine(a))
  expect_true(is_dcp(a))
})

## @cvxpy NONE
test_that("AddExpression sign: both nonneg", {
  x <- Variable(1, nonneg = TRUE)
  y <- Variable(1, nonneg = TRUE)
  a <- AddExpression(list(x, y))
  expect_true(is_nonneg(a))
  expect_false(is_nonpos(a))
})

## @cvxpy NONE
test_that("AddExpression sign: mixed → unknown", {
  x <- Variable(1, nonneg = TRUE)
  y <- Variable(1, nonpos = TRUE)
  a <- AddExpression(list(x, y))
  expect_false(is_nonneg(a))
  expect_false(is_nonpos(a))
})

## @cvxpy NONE
test_that("AddExpression value", {
  a <- AddExpression(list(Constant(3), Constant(4)))
  expect_equal(drop(value(a)), 7)
})

## @cvxpy NONE
test_that("AddExpression of vector constants", {
  v1 <- matrix(c(1, 2, 3), ncol = 1)
  v2 <- matrix(c(10, 20, 30), ncol = 1)
  a <- AddExpression(list(Constant(v1), Constant(v2)))
  expect_equal(as.numeric(value(a)), c(11, 22, 33))
})

## @cvxpy NONE
test_that("AddExpression with variable → NULL value", {
  a <- AddExpression(list(Variable(3), Constant(matrix(1:3, ncol = 1))))
  expect_null(value(a))
})

## @cvxpy NONE
test_that("AddExpression is_symmetric", {
  x <- Variable(c(3, 3), symmetric = TRUE)
  y <- Variable(c(3, 3), symmetric = TRUE)
  a <- AddExpression(list(x, y))
  expect_true(is_symmetric(a))
})

## @cvxpy NONE
test_that("AddExpression is_symmetric FALSE for non-square", {
  x <- Variable(c(3, 2))
  y <- Variable(c(3, 2))
  a <- AddExpression(list(x, y))
  expect_false(is_symmetric(a))
})

## @cvxpy NONE
test_that("AddExpression expr_name with + separator", {
  x <- Variable(1, name = "x")
  y <- Variable(1, name = "y")
  a <- AddExpression(list(x, y))
  name <- expr_name(a)
  expect_true(grepl("\\+", name))
})

## @cvxpy NONE
test_that("AddExpression three args", {
  a <- AddExpression(list(Constant(1), Constant(2), Constant(3)))
  expect_equal(length(a@args), 3L)
  expect_equal(drop(value(a)), 6)
})

## @cvxpy NONE
test_that("AddExpression graph_implementation", {
  a <- AddExpression(list(Constant(1), Constant(2)))
  cf <- canonical_form(a)
  expect_true(is.list(cf))
  expect_equal(length(cf), 2L)
})

## @cvxpy NONE
test_that("AddExpression expr_copy preserves flattening", {
  x <- Variable(3)
  y <- Variable(3)
  a <- AddExpression(list(x, y))
  copy <- expr_copy(a)
  expect_equal(length(copy@args), 2L)
  expect_true(S7_inherits(copy, AddExpression))
})


# ═══════════════════════════════════════════════════════════════════
# MulExpression (matrix multiplication)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("MulExpression class exists and inherits correctly", {
  m <- MulExpression(Constant(matrix(1:6, 3, 2)), Variable(c(2, 1)))
  expect_true(S7_inherits(m, MulExpression))
  expect_true(S7_inherits(m, BinaryOperator))
  expect_true(S7_inherits(m, AffAtom))
})

## @cvxpy NONE
test_that("MulExpression shape from mul_shapes", {
  lhs <- Constant(matrix(1:6, 3, 2))  ## 3×2
  rhs <- Variable(c(2, 1))            ## 2×1
  m <- MulExpression(lhs, rhs)
  expect_equal(m@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("MulExpression shape (3,4) x (4,2) → (3,2)", {
  lhs <- Variable(c(3, 4))
  rhs <- Constant(matrix(1:8, 4, 2))
  m <- MulExpression(lhs, rhs)
  expect_equal(m@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("MulExpression scalar × vector", {
  lhs <- Constant(2)       ## 1×1
  rhs <- Variable(c(3, 1)) ## 3×1
  m <- MulExpression(lhs, rhs)
  expect_equal(m@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("MulExpression incompatible dimensions error", {
  expect_error(MulExpression(Variable(c(3, 4)), Variable(c(5, 2))),
               "Incompatible")
})

## @cvxpy NONE
test_that("MulExpression with one constant is affine", {
  lhs <- Constant(matrix(1:6, 3, 2))
  rhs <- Variable(c(2, 1))
  m <- MulExpression(lhs, rhs)
  expect_true(is_atom_convex(m))
  expect_true(is_atom_concave(m))
  expect_true(is_affine(m))
  expect_true(is_dcp(m))
})

## @cvxpy NONE
test_that("MulExpression with two non-constants is NOT convex", {
  m <- MulExpression(Variable(c(1, 3)), Variable(c(3, 1)))
  expect_false(is_atom_convex(m))
  expect_false(is_atom_concave(m))
})

## @cvxpy NONE
test_that("MulExpression monotonicity: left const nonneg", {
  lhs <- Constant(matrix(c(1, 2, 3, 4), 2, 2))  ## nonneg constant
  rhs <- Variable(c(2, 1))
  m <- MulExpression(lhs, rhs)
  ## is_incr(1): checks args[[2]] (rhs) nonneg? rhs is variable → FALSE
  expect_false(is_incr(m, 1L))
  ## is_incr(2): checks args[[1]] (lhs) nonneg? Yes, constant positive
  expect_true(is_incr(m, 2L))
})

## @cvxpy NONE
test_that("MulExpression sign: nonneg const * nonneg var", {
  lhs <- Constant(2)
  rhs <- Variable(1, nonneg = TRUE)
  m <- MulExpression(lhs, rhs)
  expect_true(is_nonneg(m))
})

## @cvxpy NONE
test_that("MulExpression value: const * const", {
  lhs <- Constant(matrix(1:4, 2, 2))
  rhs <- Constant(matrix(c(1, 0, 0, 1), 2, 2))
  m <- MulExpression(lhs, rhs)
  expect_equal(as.matrix(value(m)), matrix(1:4, 2, 2))
})

## @cvxpy NONE
test_that("MulExpression value: scalar * vector", {
  m <- MulExpression(Constant(3), Constant(matrix(c(1, 2), ncol = 1)))
  expect_equal(as.numeric(value(m)), c(3, 6))
})

## @cvxpy NONE
test_that("MulExpression graph_implementation: left constant", {
  m <- MulExpression(Constant(matrix(1:4, 2, 2)), Variable(c(2, 1)))
  cf <- canonical_form(m)
  expect_true(is.list(cf))
})

## @cvxpy NONE
test_that("MulExpression graph_implementation: right constant", {
  m <- MulExpression(Variable(c(1, 2)), Constant(matrix(1:4, 2, 2)))
  cf <- canonical_form(m)
  expect_true(is.list(cf))
})


# ═══════════════════════════════════════════════════════════════════
# Multiply (elementwise multiplication)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Multiply class exists and inherits correctly", {
  m <- Multiply(Constant(2), Variable(3))
  expect_true(S7_inherits(m, Multiply))
  expect_true(S7_inherits(m, MulExpression))
})

## @cvxpy NONE
test_that("Multiply shape from broadcast", {
  lhs <- Constant(matrix(1:3, ncol = 1))  ## 3×1
  rhs <- Variable(c(3, 1))                ## 3×1
  m <- Multiply(lhs, rhs)
  expect_equal(m@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Multiply scalar * vector broadcasts", {
  lhs <- Constant(2)       ## scalar
  rhs <- Variable(c(3, 1)) ## 3×1
  m <- Multiply(lhs, rhs)
  expect_equal(m@shape, c(3L, 1L))
  ## LHS should be promoted
  expect_true(S7_inherits(m@args[[1L]], Promote))
})

## @cvxpy NONE
test_that("Multiply with one constant is affine", {
  m <- Multiply(Constant(3), Variable(c(3, 1)))
  expect_true(is_affine(m))
  expect_true(is_dcp(m))
})

## @cvxpy NONE
test_that("Multiply value: elementwise", {
  v1 <- matrix(c(1, 2, 3), ncol = 1)
  v2 <- matrix(c(10, 20, 30), ncol = 1)
  m <- Multiply(Constant(v1), Constant(v2))
  expect_equal(as.numeric(value(m)), c(10, 40, 90))
})

## @cvxpy NONE
test_that("Multiply value: scalar * vector", {
  m <- Multiply(Constant(5), Constant(matrix(c(1, 2, 3), ncol = 1)))
  expect_equal(as.numeric(value(m)), c(5, 10, 15))
})

## @cvxpy NONE
test_that("Multiply sign: nonneg * nonneg", {
  lhs <- Variable(1, nonneg = TRUE)
  rhs <- Variable(1, nonneg = TRUE)
  m <- Multiply(lhs, rhs)
  expect_true(is_nonneg(m))
})

## @cvxpy NONE
test_that("Multiply sign: nonneg * nonpos → nonpos", {
  lhs <- Variable(1, nonneg = TRUE)
  rhs <- Variable(1, nonpos = TRUE)
  m <- Multiply(lhs, rhs)
  expect_true(is_nonpos(m))
})

## @cvxpy NONE
test_that("Multiply PSD propagation", {
  x <- Variable(c(3, 3), PSD = TRUE)
  y <- Variable(c(3, 3), PSD = TRUE)
  m <- Multiply(x, y)
  expect_true(is_psd(m))
})

## @cvxpy NONE
test_that("Multiply graph_implementation: left constant", {
  m <- Multiply(Constant(matrix(1:3, ncol = 1)), Variable(c(3, 1)))
  cf <- canonical_form(m)
  expect_true(is.list(cf))
})


# ═══════════════════════════════════════════════════════════════════
# DivExpression
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DivExpression class exists and inherits correctly", {
  d <- DivExpression(Variable(3), Constant(2))
  expect_true(S7_inherits(d, DivExpression))
  expect_true(S7_inherits(d, BinaryOperator))
})

## @cvxpy NONE
test_that("DivExpression shape is numerator shape", {
  d <- DivExpression(Variable(c(3, 2)), Constant(2))
  expect_equal(d@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("DivExpression scalar / scalar", {
  d <- DivExpression(Variable(1), Constant(2))
  expect_equal(d@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("DivExpression broadcasts scalar denominator", {
  d <- DivExpression(Variable(c(3, 1)), Constant(2))
  expect_equal(d@shape, c(3L, 1L))
  ## Denominator should be promoted
  expect_true(S7_inherits(d@args[[2L]], Promote))
})

## @cvxpy NONE
test_that("DivExpression with constant denominator is affine", {
  d <- DivExpression(Variable(3), Constant(2))
  expect_true(is_atom_convex(d))
  expect_true(is_atom_concave(d))
  expect_true(is_affine(d))
  expect_true(is_dcp(d))
})

## @cvxpy NONE
test_that("DivExpression with non-constant denominator is NOT convex", {
  d <- DivExpression(Variable(1), Variable(1))
  expect_false(is_atom_convex(d))
  expect_false(is_dcp(d))
})

## @cvxpy NONE
test_that("DivExpression monotonicity: nonneg denominator", {
  d <- DivExpression(Variable(1), Constant(2))  ## const is nonneg
  ## is_incr(1): denominator nonneg → TRUE
  expect_true(is_incr(d, 1L))
  ## is_decr(1): denominator nonpos → FALSE
  expect_false(is_decr(d, 1L))
})

## @cvxpy NONE
test_that("DivExpression monotonicity: nonpos denominator", {
  d <- DivExpression(Variable(1), Constant(-2))
  ## is_incr(1): denominator nonneg → FALSE
  expect_false(is_incr(d, 1L))
  ## is_decr(1): denominator nonpos → TRUE
  expect_true(is_decr(d, 1L))
})

## @cvxpy NONE
test_that("DivExpression sign: nonneg / nonneg → nonneg", {
  d <- DivExpression(Variable(1, nonneg = TRUE), Constant(2))
  expect_true(is_nonneg(d))
})

## @cvxpy NONE
test_that("DivExpression sign: nonneg / nonpos → nonpos", {
  d <- DivExpression(Variable(1, nonneg = TRUE), Constant(-2))
  expect_true(is_nonpos(d))
})

## @cvxpy NONE
test_that("DivExpression value", {
  d <- DivExpression(Constant(10), Constant(2))
  expect_equal(drop(value(d)), 5)
})

## @cvxpy NONE
test_that("DivExpression vector / scalar value", {
  v <- matrix(c(10, 20, 30), ncol = 1)
  d <- DivExpression(Constant(v), Constant(10))
  expect_equal(as.numeric(value(d)), c(1, 2, 3))
})

## @cvxpy NONE
test_that("DivExpression graph_implementation", {
  d <- DivExpression(Variable(3), Constant(2))
  cf <- canonical_form(d)
  expect_true(is.list(cf))
})


# ═══════════════════════════════════════════════════════════════════
# BinaryOperator sign (mul_sign rules)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("BinaryOperator sign: nonneg * nonneg", {
  m <- Multiply(Variable(1, nonneg = TRUE), Constant(3))
  expect_true(is_nonneg(m))
  expect_false(is_nonpos(m))
})

## @cvxpy NONE
test_that("BinaryOperator sign: nonneg * nonpos → nonpos", {
  m <- Multiply(Variable(1, nonneg = TRUE), Constant(-3))
  expect_false(is_nonneg(m))
  expect_true(is_nonpos(m))
})

## @cvxpy NONE
test_that("BinaryOperator sign: nonpos * nonpos → nonneg", {
  m <- Multiply(Variable(1, nonpos = TRUE), Constant(-3))
  expect_true(is_nonneg(m))
  expect_false(is_nonpos(m))
})

## @cvxpy NONE
test_that("BinaryOperator sign: zero * anything → zero", {
  m <- Multiply(Constant(0), Variable(3))
  expect_true(is_nonneg(m))
  expect_true(is_nonpos(m))
  expect_true(is_zero(m))
})


# ═══════════════════════════════════════════════════════════════════
# Complex propagation through BinaryOperator
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("BinaryOperator is_complex: complex * real → complex", {
  x <- Variable(1, complex = TRUE)
  y <- Constant(2)
  m <- Multiply(x, y)
  expect_true(is_complex(m))
})

## @cvxpy NONE
test_that("BinaryOperator is_imag: imag * real → imag", {
  x <- Variable(1, imag = TRUE)
  y <- Constant(2)
  m <- Multiply(x, y)
  expect_true(is_imag(m))
})

## @cvxpy NONE
test_that("BinaryOperator is_complex: imag * imag → NOT complex (real)", {
  x <- Variable(1, imag = TRUE)
  y <- Variable(1, imag = TRUE)
  m <- Multiply(x, y)
  expect_false(is_complex(m))
})


# ═══════════════════════════════════════════════════════════════════
# Expression tree construction
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Nested: Neg(Add(x, y)) is affine", {
  x <- Variable(3)
  y <- Variable(3)
  a <- AddExpression(list(x, y))
  n <- NegExpression(a)
  expect_true(is_affine(n))
})

## @cvxpy NONE
test_that("Nested: Add(Neg(x), y) is affine", {
  x <- Variable(3)
  y <- Variable(3)
  n <- NegExpression(x)
  a <- AddExpression(list(n, y))
  expect_true(is_affine(a))
})

## @cvxpy NONE
test_that("Nested: Mul(const, Add(x, y)) is affine", {
  x <- Variable(c(3, 1))
  y <- Variable(c(3, 1))
  a <- AddExpression(list(x, y))
  c_mat <- Constant(matrix(1:9, 3, 3))
  m <- MulExpression(c_mat, a)
  expect_true(is_affine(m))
})

## @cvxpy NONE
test_that("Variables collected from nested expression", {
  x <- Variable(c(3, 1))
  y <- Variable(c(3, 1))
  a <- AddExpression(list(x, y))
  n <- NegExpression(a)
  vars <- variables(n)
  expect_equal(length(vars), 2L)
})

## @cvxpy NONE
test_that("Constants collected from nested expression", {
  x <- Variable(c(3, 1))
  c1 <- Constant(matrix(1:3, ncol = 1))
  m <- Multiply(c1, x)
  consts <- constants(m)
  expect_equal(length(consts), 1L)
})


# ═══════════════════════════════════════════════════════════════════
# expr_name formatting
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("MulExpression expr_name contains %*%", {
  m <- MulExpression(Constant(matrix(1, 1, 1)), Variable(1))
  expect_true(grepl("%\\*%", expr_name(m)))
})

## @cvxpy NONE
test_that("Multiply expr_name contains *", {
  m <- Multiply(Constant(2), Variable(1))
  expect_true(grepl("\\*", expr_name(m)))
})

## @cvxpy NONE
test_that("DivExpression expr_name contains /", {
  d <- DivExpression(Variable(1), Constant(2))
  expect_true(grepl("/", expr_name(d)))
})

## @cvxpy NONE
test_that("AddExpression expr_name contains +", {
  a <- AddExpression(list(Variable(1), Variable(1)))
  expect_true(grepl("\\+", expr_name(a)))
})
