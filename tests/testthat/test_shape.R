test_that("test shape size", {
  expect_equal(size(Shape(1,3)), c(1,3))
  expect_equal(size(Shape(2,1)), c(2,1))
})

test_that("test shape addition", {
  expect_equal(size(Shape(3,4) + Shape(3,4)), c(3,4))
  expect_error(Shape(1,3) + Shape(4,3))
  expect_equal(size(Shape(3,4) + Shape(1,1)), c(3,4))
  expect_equal(size(Shape(1,1) + Shape(3,4)), c(3,4))
})

test_that("test shape multiplication", {
  expect_equal(size(Shape(5,9) * Shape(9,2)), c(5,2))
  expect_error(Shape(5,3) * Shape(9,2))
  expect_equal(size(Shape(3,4) * Shape(1,1)), c(3,4))
  expect_equal(size(Shape(1,1) * Shape(3,4)), c(3,4))
})