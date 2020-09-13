context("test-g01-curvature")

cvx <- Variable()^2
ccv <- Variable()^0.5
aff <- Variable()
const <- Constant(5)
unknown_curv <- log(Variable()^3)

pos <- Constant(1)
neg <- Constant(-1)
zero <- Constant(0)
unknown_sign <- Parameter()

test_that("test curvature addition", {
  skip_on_cran()
  expect_equal(curvature(const + cvx), curvature(cvx))
  expect_equal(curvature(unknown_curv + ccv), curvature(unknown_curv))
  expect_equal(curvature(cvx + ccv), curvature(unknown_curv))
  expect_equal(curvature(cvx + cvx), curvature(cvx))
  expect_equal(curvature(aff + ccv), curvature(ccv))
})

test_that("test curvature subtraction", {
  skip_on_cran()
  expect_equal(curvature(const - cvx), curvature(ccv))
  expect_equal(curvature(unknown_curv - ccv), curvature(unknown_curv))
  expect_equal(curvature(cvx - ccv), curvature(cvx))
  expect_equal(curvature(cvx - cvx), curvature(unknown_curv))
  expect_equal(curvature(aff - ccv), curvature(cvx))
})

test_that("test multiplication of sign and curvature", {
  skip_on_cran()
  expect_equal(curvature(zero * cvx), curvature(aff))
  expect_equal(curvature(neg * cvx), curvature(ccv))
  expect_equal(curvature(neg * ccv), curvature(cvx))
  expect_equal(curvature(neg * unknown_curv), curvature(unknown_curv))
  expect_equal(curvature(pos * aff), curvature(aff))
  expect_equal(curvature(pos * ccv), curvature(ccv))
  expect_equal(curvature(unknown_sign * const), curvature(const))
  expect_equal(curvature(unknown_sign * ccv), curvature(unknown_curv))
})

test_that("test curvature negation", {
  skip_on_cran()
  expect_equal(curvature(-cvx), curvature(ccv))
  expect_equal(curvature(-aff), curvature(aff))
})

test_that("test if curvature is affine, convex, or concave", {
  skip_on_cran()
  expect_true(is_affine(const))
  expect_true(is_affine(aff))
  expect_false(is_affine(cvx))
  expect_false(is_affine(ccv))
  expect_false(is_affine(unknown_curv))

  expect_true(is_convex(const))
  expect_true(is_convex(aff))
  expect_true(is_convex(cvx))
  expect_false(is_convex(ccv))
  expect_false(is_convex(unknown_curv))

  expect_true(is_concave(const))
  expect_true(is_concave(aff))
  expect_false(is_concave(cvx))
  expect_true(is_concave(ccv))
  expect_false(is_concave(unknown_curv))
})
