context("test-g01-dimensions")

sum_dims <- CVXR:::sum_dims
mul_dims <- CVXR:::mul_dims

test_that("test addition of matching dimensions", {
  expect_equal(sum_dims(list(c(3,4), c(3,4))), c(3,4))
  expect_equal(sum_dims(lapply(1:5, function(i) { c(3,4) })), c(3,4))
})

test_that("test broadcasting of dimensions during addition", {
  # Broadcasting with scalars is permitted.
  expect_equal(sum_dims(list(c(3,4), c(1,1))), c(3,4))
  expect_equal(sum_dims(list(c(1,1), c(3,4))), c(3,4))
  
  # All other types of broadcasting is not permitted.
  expect_error(sum_dims(list(c(4,2), c(4,1))))
  expect_error(sum_dims(list(c(4,1), c(4,2))))
})

test_that("test addition of incompatible dimensions raises an error", {
  expect_error(sum_dims(list(c(4,2), c(4,1))))
})

test_that("test multiplication by scalars raises an error", {
  expect_error(mul_dims(c(1,1), c(5,9)))
  expect_error(mul_dims(c(5,9), c(1,1)))
})

test_that("test multiplication of dimensions", {
  # Test multiplication where at least one of the shapes is >= 2D.
  expect_equal(mul_dims(c(5,9), c(9,2)), c(5,2))
  expect_equal(mul_dims(c(3,5,9), c(3,9,2)), c(3,5,2))
  
  expect_error(mul_dims(c(5,3), c(9,2)))
  expect_error(mul_dims(c(3,5,9), c(4,9,2)))
})