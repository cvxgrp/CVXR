test_that("test shape addition", {
  expect_equal(sum_shapes(list(c(3,4), c(3,4))), c(3,4))
  expect_error(sum_shapes(list(c(1,3), c(4,3))))
  
  # Promotion
  expect_equal(sum_shapes(list(c(3,4), c(1,1))), c(3,4))
  expect_equal(sum_shapes(list(c(1,1), c(3,4))), c(3,4))
})

test_that("test shape multiplication", {
  expect_equal(mul_shapes(c(5,9), c(9,2)), c(5,2))
  expect_error(mul_shapes(c(5,3), c(9,2)))
  
  # Promotion
  expect_equal(mul_shapes(c(3,4), c(1,1)), c(3,4))
  expect_equal(mul_shapes(c(1,1), c(3,4)), c(3,4))
})