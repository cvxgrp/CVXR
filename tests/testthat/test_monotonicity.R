test_that("test DCP curvature function", {
  expect_equal(dcp_curvature(INCREASING, Curvature.AFFINE, Sign.POSITIVE, Curvature.CONVEX), Curvature.CONVEX)
  expect_equal(dcp_curvature(NONMONOTONIC, Curvature.AFFINE, Sign.POSITIVE, Curvature.AFFINE), Curvature.AFFINE)
  expect_equal(dcp_curvature(DECREASING, Curvature.UNKNOWN, Sign.POSITIVE, Curvature.CONSTANT), Curvature.CONSTANT)
  expect_equal(dcp_curvature(INCREASING, Curvature.CONVEX, Sign.POSITIVE, Curvature.CONVEX), Curvature.CONVEX)
  expect_equal(dcp_curvature(DECREASING, Curvature.CONVEX, Sign.POSITIVE, Curvature.CONCAVE), Curvature.CONVEX)
  expect_equal(dcp_curvature(INCREASING, Curvature.CONCAVE, Sign.POSITIVE, Curvature.CONCAVE), Curvature.CONCAVE)
  expect_equal(dcp_curvature(DECREASING, Curvature.CONCAVE, Sign.POSITIVE, Curvature.CONVEX), Curvature.CONCAVE)
  expect_equal(dcp_curvature(INCREASING, Curvature.CONCAVE, Sign.POSITIVE, Curvature.CONVEX), Curvature.UNKNOWN)
  expect_equal(dcp_curvature(NONMONOTONIC, Curvature.CONCAVE, Sign.POSITIVE, Curvature.AFFINE), Curvature.CONCAVE)
  expect_equal(dcp_curvature(NONMONOTONIC, Curvature.CONSTANT, Sign.POSITIVE, Curvature.UNKNOWN), Curvature.UNKNOWN)
})

test_that("test signed curvature", {
  # Convex argument
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.POSITIVE, Curvature.CONVEX), Curvature.CONVEX)
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.NEGATIVE, Curvature.CONVEX), Curvature.UNKNOWN)
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.UNKNOWN, Curvature.CONVEX), Curvature.UNKNOWN)
  
  # Concave argument
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.POSITIVE, Curvature.CONCAVE), Curvature.UNKNOWN)
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.NEGATIVE, Curvature.CONCAVE), Curvature.CONVEX)
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.UNKNOWN, Curvature.CONCAVE), Curvature.UNKNOWN)
  
  # Affine argument
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.POSITIVE, Curvature.AFFINE), Curvature.CONVEX)
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.NEGATIVE, Curvature.AFFINE), Curvature.CONVEX)
  expect_equal(dcp_curvature(SIGNED, Curvature.CONVEX, Sign.UNKNOWN, Curvature.AFFINE), Curvature.CONVEX)
})
