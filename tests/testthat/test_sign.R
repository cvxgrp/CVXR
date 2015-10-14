test_that("test sign addition", {
  expect_equal(Sign.POSITIVE + Sign.NEGATIVE, Sign.UNKNOWN)
  expect_equal(Sign.NEGATIVE + Sign.ZERO, Sign.NEGATIVE)
  expect_equal(Sign.POSITIVE + Sign.POSITIVE, Sign.POSITIVE)
  expect_equal(Sign.UNKNOWN + Sign.ZERO, Sign.UNKNOWN)
})

test_that("test sign subtraction", {
  expect_equal(Sign.POSITIVE - Sign.NEGATIVE, Sign.POSITIVE)
  expect_equal(Sign.NEGATIVE - Sign.ZERO, Sign.NEGATIVE)
  expect_equal(Sign.POSITIVE - Sign.POSITIVE, Sign.UNKNOWN)
})

test_that("test sign multiplication", {
  expect_equal(Sign.ZERO * Sign.POSITIVE, Sign.ZERO)
  expect_equal(Sign.UNKNOWN * Sign.POSITIVE, Sign.UNKNOWN)
  expect_equal(Sign.POSITIVE * Sign.NEGATIVE, Sign.NEGATIVE)
  expect_equal(Sign.POSITIVE * Sign.POSITIVE, Sign.POSITIVE)
  expect_equal(Sign.POSITIVE * Sign("positive"), Sign.POSITIVE)
  expect_equal(Sign.NEGATIVE * Sign.NEGATIVE, Sign.POSITIVE)
  expect_equal(Sign.ZERO * Sign.UNKNOWN, Sign.ZERO)
})

test_that("test sign negation", {
  expect_equal(-Sign.ZERO, Sign.ZERO)
  expect_equal(-Sign.POSITIVE, Sign.NEGATIVE)
})

test_that("test if sign is positive, negative, or zero", {
  expect_true(is_positive(Sign.POSITIVE))
  expect_false(is_positive(Sign.NEGATIVE))
  expect_false(is_positive(Sign.UNKNOWN))
  expect_true(is_positive(Sign.ZERO))
  
  expect_false(is_negative(Sign.POSITIVE))
  expect_true(is_negative(Sign.NEGATIVE))
  expect_false(is_negative(Sign.UNKNOWN))
  expect_true(is_negative(Sign.ZERO))
  
  expect_true(is_zero(Sign.ZERO))
  expect_false(is_zero(Sign.NEGATIVE))
  
  expect_true(is_unknown(Sign.UNKNOWN))
})
