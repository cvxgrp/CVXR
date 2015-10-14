test_that("test curvature addition", {
  expect_equal(Curvature.CONSTANT + Curvature.CONVEX, Curvature.CONVEX)
  expect_equal(Curvature.UNKNOWN + Curvature.CONCAVE, Curvature.UNKNOWN)
  expect_equal(Curvature.CONVEX + Curvature.CONCAVE, Curvature.UNKNOWN)
  expect_equal(Curvature.CONVEX + Curvature.CONVEX, Curvature.CONVEX)
  expect_equal(Curvature.AFFINE + Curvature.CONCAVE, Curvature.CONCAVE)
})

test_that("test curvature subtraction", {
  expect_equal(Curvature.CONSTANT - Curvature.CONVEX, Curvature.CONCAVE)
  expect_equal(Curvature.UNKNOWN - Curvature.CONCAVE, Curvature.UNKNOWN)
  expect_equal(Curvature.CONVEX - Curvature.CONCAVE, Curvature.CONVEX)
  expect_equal(Curvature.CONVEX - Curvature.CONVEX, Curvature.UNKNOWN)
  expect_equal(Curvature.AFFINE - Curvature.CONCAVE, Curvature.CONVEX)
})

test_that("test multiplication of sign and curvature", {
  expect_equal(Sign.ZERO * Curvature.CONVEX, Curvature.CONSTANT)
  expect_equal(Sign.NEGATIVE * Curvature.CONVEX, Curvature.CONCAVE)
  expect_equal(Sign.NEGATIVE * Curvature.CONCAVE, Curvature.CONVEX)
  expect_equal(Sign.POSITIVE * Curvature.AFFINE, Curvature.AFFINE)
  expect_equal(Sign.POSITIVE * Curvature.CONCAVE, Curvature.CONCAVE)
  expect_equal(Sign.UNKNOWN * Curvature.CONSTANT, Curvature.CONSTANT)
  expect_equal(Sign.UNKNOWN * Curvature.CONCAVE, Curvature.UNKNOWN)
})

test_that("test curvature negation", {
  expect_equal(-Curvature.CONVEX, Curvature.CONCAVE)
  expect_equal(-Curvature.AFFINE, Curvature.AFFINE)
})

test_that("test if curvature is affine, convex, or concave", {
  expect_true(is_affine(Curvature.CONSTANT))
  expect_true(is_affine(Curvature.AFFINE))
  expect_false(is_affine(Curvature.CONVEX))
  expect_false(is_affine(Curvature.CONCAVE))
  expect_false(is_affine(Curvature.UNKNOWN))
  
  expect_true(is_convex(Curvature.CONSTANT))
  expect_true(is_convex(Curvature.AFFINE))
  expect_true(is_convex(Curvature.CONVEX))
  expect_false(is_convex(Curvature.CONCAVE))
  expect_false(is_convex(Curvature.UNKNOWN))
  
  expect_true(is_concave(Curvature.CONSTANT))
  expect_true(is_concave(Curvature.AFFINE))
  expect_false(is_concave(Curvature.CONVEX))
  expect_true(is_concave(Curvature.CONCAVE))
  expect_false(is_concave(Curvature.UNKNOWN))
})
