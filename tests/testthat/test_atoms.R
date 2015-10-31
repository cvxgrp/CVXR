a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test the NormInf class", {
  exp <- x + y
  atom <- NormInf(exp)
  
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_true(is_convex(atom))
  expect_true(is_concave(-atom))
  expect_equal(curvature(NormInf(atom)), Curvature.CONVEX)
  expect_equal(curvature(NormInf(-atom)), Curvature.CONVEX)
})

test_that("test the Norm1 class", {
  exp <- x + y
  atom <- Norm1(exp)
  
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(curvature(Norm1(atom)), Curvature.CONVEX)
  expect_equal(curvature(Norm1(-atom)), Curvature.CONVEX)
})

test_that("test the Norm2 class", {
  exp <- x + y
  atom <- Norm2(exp)
  
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(curvature(Norm2(atom)), Curvature.CONVEX)
  expect_equal(curvature(Norm2(-atom)), Curvature.CONVEX)
})

test_that("test the Power class", {
  for(size in list(c(1,1), c(3,1), c(2,3))) {
    x_pow <- Variable(size[1], size[2])
    y_pow <- Variable(size[1], size[2])
    exp <- x_pow + y_pow
    
    for(p in c(0, 1, 2, 3, 2.7, 0.67, -1, -2.3, 4/5)) {
      atom <- Power(exp, p)
      expect_equal(size(atom), size)
      
      if(p > 1 || p < 0)
        expect_equal(curvature(atom), Curvature.CONVEX)
      else if(p == 1)
        expect_equal(curvature(atom), Curvature.AFFINE)
      else if(p == 0)
        expect_equal(curvature(atom), Curvature.CONSTANT)
      else
        expect_equal(curvature(atom), Curvature.CONCAVE)
      
      if(p != 1)
        expect_equal(sign(atom), Sign.POSITIVE)
    }
  }
})

test_that("test the HarmonicMean class", {
  atom <- HarmonicMean(x)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
})

test_that("test the Pnorm class", {
  atom <- Pnorm(x, p = 1.5)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = 1)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = 2)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = Inf)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = 0.5)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = 0.7)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = -0.1)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = -1)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = -1.3)
  expect_equal(size(atom), c(1,1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
})

test_that("test the QuadOverLin class", {
  atom <- QuadOverLin(Square(x), a)
  expect_equal(curvature(atom), Curvature.CONVEX)
  
  atom <- QuadOverLin(-Square(x), a)
  expect_equal(curvature(atom), Curvature.CONVEX)
  
  atom <- QuadOverLin(Sqrt(x), a)
  expect_equal(curvature(atom), Curvature.UNKNOWN)
  expect_false(is_dcp(atom))
  
  expect_error(QuadOverLin(x, x))
})

test_that("test the sign for MaxEntries", {
  expect_equal(sign(MaxEntries(1)), Sign.POSITIVE)
  expect_equal(sign(MaxEntries(-2)), Sign.NEGATIVE)
  expect_equal(sign(MaxEntries(Variable())), Sign.UNKNOWN)
  expect_equal(sign(MaxEntries(0)), Sign.ZERO)
})

test_that("test the sign for MinEntries", {
  expect_equal(sign(MinEntries(1)), Sign.POSITIVE)
  expect_equal(sign(MinEntries(-2)), Sign.NEGATIVE)
  expect_equal(sign(MinEntries(Variable())), Sign.UNKNOWN)
  expect_equal(sign(MinEntries(0)), Sign.ZERO)
})

test_that("test sign logic for MaxElemwise", {
  expect_equal(sign(MaxElemwise(1, 2)), Sign.POSITIVE)
  expect_equal(sign(MaxElemwise(1, Variable())), Sign.POSITIVE)
  expect_equal(sign(MaxElemwise(1, -2)), Sign.POSITIVE)
  expect_equal(sign(MaxElemwise(1, 0)), Sign.POSITIVE)
  
  expect_equal(sign(MaxElemwise(Variable(), 0)), Sign.POSITIVE)
  expect_equal(sign(MaxElemwise(Variable(), Variable())), Sign.UNKNOWN)
  expect_equal(sign(MaxElemwise(Variable(), -2)), Sign.UNKNOWN)
  
  expect_equal(sign(MaxElemwise(0, 0)), Sign.ZERO)
  expect_equal(sign(MaxElemwise(0, -2)), Sign.ZERO)
  
  expect_equal(sign(MaxElemwise(-3, -2)), Sign.NEGATIVE)
  
  expect_equal(sign(MaxElemwise(-2, Variable(), 0, -1, Variable(), -1)), Sign.POSITIVE)
  
  expect_equal(sign(MaxElemwise(1, Variable(2))), Sign.POSITIVE)
  expect_equal(size(MaxElemwise(1, Variable(2))), c(2,1))
})

test_that("test sign logic for MinElemwise", {
  expect_equal(sign(MinElemwise(1, 2)), Sign.POSITIVE)
  expect_equal(sign(MinElemwise(1, Variable())), Sign.UNKNOWN)
  expect_equal(sign(MinElemwise(1, -2)), Sign.NEGATIVE)
  expect_equal(sign(MinElemwise(1, 0)), Sign.ZERO)
  
  expect_equal(sign(MinElemwise(Variable(), 0)), Sign.NEGATIVE)
  expect_equal(sign(MinElemwise(Variable(), Variable())), Sign.UNKNOWN)
  expect_equal(sign(MinElemwise(Variable(), -2)), Sign.NEGATIVE)
  
  expect_equal(sign(MinElemwise(0, 0)), Sign.ZERO)
  expect_equal(sign(MinElemwise(0, -2)), Sign.NEGATIVE)
  
  expect_equal(sign(MinElemwise(-3, -2)), Sign.NEGATIVE)
  
  expect_equal(sign(MinElemwise(-2, Variable(), 0, -1, Variable(), 1)), Sign.NEGATIVE)
  
  expect_equal(sign(MinElemwise(-1, Variable(2))), Sign.NEGATIVE)
  expect_equal(size(MinElemwise(-1, Variable(2))), c(2, 1))
})

test_that("test the SumEntries class", {
  expect_equal(sign(SumEntries(1)), Sign.POSITIVE)
  expect_equal(sign(SumEntries(c(1, -1))), Sign.UNKNOWN)
  expect_equal(curvature(SumEntries(c(1, -1))), Curvature.CONSTANT)
  expect_equal(sign(SumEntries(Variable(2))), Sign.UNKNOWN)
  expect_equal(size(SumEntries(Variable(2))), c(1,1))
  expect_equal(curvature(SumEntries(Variable(2))), Curvature.AFFINE)
  
  # Mixed curvature
  expect_equal(curvature(SumEntries( c(1, -1) * Square(Variable(2)) )), Curvature.UNKNOWN)
})

test_that("test the Reshape class", {
  expr <- Reshape(A, 4, 1)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(4,1))
  
  expr <- Reshape(expr, 2, 2)
  expect_equal(size(expr), c(2,2))
  
  expr <- Reshape(Square(x), 1, 2)
  expect_equal(sign(expr), Sign.POSITIVE)
  expect_equal(curvature(expr), Curvature.CONVEX)
  expect_equal(size(expr), c(1,2))
  
  expect_error(Reshape(C, 5, 4))
})

test_that("test the Diag class", {
  expr <- diag(x)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(2,2))
  
  expr <- diag(A)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(2,1))
  
  expect_error(diag(C))
})

test_that("test the Trace class", {
  expr <- Trace(A)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(1,1))
  
  expect_error(Trace(C))
})

test_that("test the SumLargest class", {
  expect_error(SumLargest(x, -1))
  # expect_error(LambdaSumLargest(x, 2.4))   # TODO: Currently unimplemented
  # expect_error(LambdaSumLargest(Variable(2,2), 2.4))
})

test_that("test the SumSmallest class", {
  expect_error(SumSmallest(x, -1))
  # expect_error(LambdaSumSmallest(Variable(2,2), 2.4))   TODO: Currently unimplemented
})

test_that("test the Conv class", {
  a <- matrix(1, nrow = 3, ncol = 1)
  b <- Parameter(2, sign = "positive")
  expr <- Conv(a, b)
  expect_true(is_positive(expr))
  expect_equal(size(expr), c(4,1))
  
  b <- Parameter(2, sign = "negative")
  expr <- Conv(a, b)
  expect_true(is_negative(expr))
  expect_error(Conv(x, -1))
})

test_that("test the Kron class", {
  a <- matrix(1, nrow = 3, ncol = 2)
  b <- Parameter(2, sign = "positive")
  expr <- Kron(a, b)
  expect_true(is_positive(expr))
  expect_equal(size(expr), c(6,2))
  
  b <- Parameter(2, sign = "negative")
  expr <- Kron(a, b)
  expect_true(is_negative(expr))
  expect_error(Kron(x, -1))
})
