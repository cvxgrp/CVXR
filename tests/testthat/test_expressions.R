a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test the Variable class", {
  x <- Variable(2, name = "x")
  y <- Variable()
  expect_equal(size(x), c(2,1))
  expect_equal(size(y), c(1,1))
  expect_equal(curvature(x), Curvature.AFFINE)
  expect_equal(canonical_form(x)[[1]]$size, c(2,1))
  expect_equal(canonical_form(x)[[2]], list())
})

test_that("test transposing variables", {
  var <- t(a)
  expect_equal(size(var), c(1,1))
  
  var <- t(x)
  expect_equal(size(var), c(1,2))
  
  var <- t(C)
  expect_equal(size(var), c(2,3))
  
  var <- t(t(x))
  expect_equal(size(var), c(2,1))
})

test_that("test the Constant class", {
  c <- Constant(2)
  expect_equal(size(c), c(1,1))
  expect_equal(curvature(c), Curvature.CONSTANT)
  expect_equal(sign(c), Sign.POSITIVE)
  expect_equal(sign(Constant(-2)), Sign.NEGATIVE)
  expect_equal(sign(Constant(0)), Sign.ZERO)
  expect_equal(canonical_form(c)[[1]]$size, c(1,1))
  expect_equal(canonical_form(c)[[2]], list())
})

test_that("test the NegExpression class", {
  exp <- -x
  expect_equal(curvature(exp), Curvature.AFFINE)
  expect_true(is_affine(exp))
  expect_equal(sign(exp), Sign.UNKNOWN)
  expect_true(!is_positive(exp))
  # expect_equal(size(canonical_form(exp)[[1]]), c(2,1))   # TODO: Need to get canonical form working
  # expect_equal(canonical_form(exp)[[2]], list())
  expect_equal(size(exp), size(x))
  
  exp <- -C
  expect_equal(curvature(exp), Curvature.AFFINE)
  expect_equal(size(exp), c(3,2))
})
