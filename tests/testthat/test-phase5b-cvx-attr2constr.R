# Tests for Phase 5b: CvxAttr2Constr

library(CVXR)

# ===================================================================
# CvxAttr2Constr
# ===================================================================

## @cvxpy NONE
test_that("CvxAttr2Constr passes through non-attributed variables", {
  x <- Variable(3)
  p <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  ca <- CvxAttr2Constr()
  result <- reduction_apply(ca, p)
  ## When no convex attributes, returns problem unchanged
  new_p <- result[[1]]
  expect_equal(length(new_p@constraints), length(p@constraints))
})

## @cvxpy NONE
test_that("CvxAttr2Constr strips nonneg attribute and adds NonNeg constraint", {
  x <- Variable(2, nonneg = TRUE)
  p <- Problem(Minimize(sum_entries(x)), list(x <= 3))
  ca <- CvxAttr2Constr()
  result <- reduction_apply(ca, p)
  new_p <- result[[1]]

  ## Should add NonNeg constraint for x >= 0 plus the original x <= 3
  ## New problem should have: NonNeg + lowered(x <= 3)
  expect_true(length(new_p@constraints) > length(p@constraints))

  ## Check that at least one constraint is NonNeg
  cls_names <- vapply(new_p@constraints, function(c) sub("^.*::", "", class(c)[[1]]),
                      character(1))
  expect_true("NonNeg" %in% cls_names)
})

## @cvxpy NONE
test_that("CvxAttr2Constr nonneg produces correct solver data", {
  x <- Variable(2, nonneg = TRUE)
  p <- Problem(Minimize(sum_entries(x)), list(x <= 3))
  data <- problem_data(p, CLARABEL_SOLVER)$data

  ## CVXPY reference:
  ## c = [1, 1]
  ## A = [[-1, 0], [0, -1], [1, 0], [0, 1]]
  ## b = [0, 0, 3, 3]
  ## dims: zero=0, nonneg=4
  expect_equal(data$c, c(1, 1))
  expect_equal(as.numeric(data$b), c(0, 0, 3, 3))
  expect_equal(data$dims@nonneg, 4L)
  expect_equal(data$dims@zero, 0L)
})

## @cvxpy NONE
test_that("CvxAttr2Constr nonpos adds NonPos -> NonNeg", {
  x <- Variable(2, nonpos = TRUE)
  p <- Problem(Maximize(sum_entries(x)), list(x >= -3))
  ca <- CvxAttr2Constr()
  result <- reduction_apply(ca, p)
  new_p <- result[[1]]

  cls_names <- vapply(new_p@constraints, function(c) sub("^.*::", "", class(c)[[1]]),
                      character(1))
  ## Should have NonPos (from attribute) and possibly others
  expect_true("NonPos" %in% cls_names || "NonNeg" %in% cls_names)
})

## @cvxpy NONE
test_that("convex_attributes utility", {
  x <- Variable(3, nonneg = TRUE)
  y <- Variable(2)
  expect_true("nonneg" %in% convex_attributes(list(x)))
  expect_equal(length(convex_attributes(list(y))), 0L)
  expect_true("nonneg" %in% convex_attributes(list(x, y)))
})
