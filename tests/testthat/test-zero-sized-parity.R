## CVXPY 1.9 zero-sized parity audit.
##
## CVXR rejects zero dimensions for leaves and reshape targets because its
## public expression model is positive-dimensional and 2D. The all-false
## logical-indexing cases below are the supported zero-length edge that falls
## out of R indexing.

## @cvxpy test_zero_sized.py::TestZeroSizedFromBooleanIndex::test_all_false_boolean_index_1d test_zero_sized.py::TestZeroSizedFromBooleanIndex::test_all_false_boolean_index_2d_rows test_zero_sized.py::TestZeroSizedFromBooleanIndex::test_all_false_boolean_index_2d_cols
test_that("all-false logical indexing can create zero-length selections", {
  x <- Variable(5)
  expect_equal(x[rep(FALSE, 5), ]@shape, c(0L, 1L))

  X <- Variable(c(3, 4))
  expect_equal(X[rep(FALSE, 3), ]@shape, c(0L, 4L))
  expect_equal(X[, rep(FALSE, 4)]@shape, c(3L, 0L))
})

## @cvxpy test_zero_sized.py::TestZeroSizedFromBooleanIndex::test_sum_of_all_false_boolean_index
test_that("sum of all-false logical index is a scalar zero contribution", {
  x <- Variable(5)
  expr <- sum_entries(x[rep(FALSE, 5), ])
  expect_equal(expr@shape, c(1L, 1L))

  prob <- Problem(Minimize(sum_entries(x) + expr), list(x >= 1))
  expect_equal(psolve(prob, solver = "CLARABEL"), 5, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), rep(1, 5), tolerance = 1e-5)
})
