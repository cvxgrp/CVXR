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
