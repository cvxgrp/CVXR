context("test-g01-sign")

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
  expect_equal(sign(zero * pos), sign(zero))
  expect_equal(sign(unknown * pos), sign(unknown))
  expect_equal(sign(pos * neg), sign(neg))
  expect_equal(sign(pos * pos), sign(pos))
  expect_equal(sign(neg * neg), sign(pos))
  expect_equal(sign(zero * unknown), sign(zero))
})

test_that("test sign negation", {
  expect_equal(sign(-zero), sign(zero))
  expect_equal(sign(-pos), sign(neg))
})

test_that("test if sign is nonnegative or nonpositive", {
  expect_true(is_nonneg(pos))
  expect_false(is_nonneg(neg))
  expect_false(is_nonneg(unknown))
  expect_true(is_nonneg(zero))
  
  expect_false(is_nonpos(pos))
  expect_true(is_nonpos(neg))
  expect_false(is_nonpos(unknown))
  expect_true(is_nonpos(zero))
  
  expect_true(is_zero(zero))
  expect_false(is_zero(neg))
  
  expect_false(is_nonneg(unknown) || is_nonpos(unknown))
})
