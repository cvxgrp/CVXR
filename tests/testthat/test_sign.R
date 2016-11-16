pos <- Constant(1)
neg <- Constant(-1)
zero <- Constant(0)
unknown <- Variable()

test_that("test sign addition", {
  expect_equal(sign(pos + neg), sign(unknown))
  expect_equal(sign(neg + zero), sign(neg))
  expect_equal(sign(pos + pos), sign(pos))
  expect_equal(sign(unknown + zero), sign(unknown))
})

test_that("test sign subtraction", {
  expect_equal(sign(pos - neg), sign(pos))
  expect_equal(sign(neg - zero), sign(neg))
  expect_equal(sign(pos - pos), sign(unknown))
})

test_that("test sign multiplication", {
  expect_equal(sign(zero %*% pos), sign(zero))
  expect_equal(sign(unknown %*% pos), sign(unknown))
  expect_equal(sign(pos %*% neg), sign(neg))
  expect_equal(sign(pos %*% pos), sign(pos))
  expect_equal(sign(neg %*% neg), sign(pos))
  expect_equal(sign(zero %*% unknown), sign(zero))
})

test_that("test sign negation", {
  expect_equal(sign(-zero), sign(zero))
  expect_equal(sign(-pos), sign(neg))
})

test_that("test if sign is positive, negative, or zero", {
  expect_true(is_positive(pos))
  expect_false(is_positive(neg))
  expect_false(is_positive(unknown))
  expect_true(is_positive(zero))
  
  expect_false(is_negative(pos))
  expect_true(is_negative(neg))
  expect_false(is_negative(unknown))
  expect_true(is_negative(zero))
  
  expect_true(is_zero(zero))
  expect_false(is_zero(neg))
  
  expect_false(is_positive(unknown) || is_negative(unknown))
})
