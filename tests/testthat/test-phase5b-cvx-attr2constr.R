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

# ===================================================================
# reduce_bounds parameter (CVXPY 1.9.0 -- V19-PARTIAL slice)
# CVXPY SOURCE: cvx_attr2constr.py:108-119,168-171
# ===================================================================

## @cvxpy NONE
test_that("CvxAttr2Constr defaults to reduce_bounds = TRUE", {
  ## R DEVIATION (documented in cvx_attr2constr.R header): default is TRUE
  ## for direct reduction calls; solver chains pass an explicit value.
  expect_true(CvxAttr2Constr()@reduce_bounds)
  expect_false(CvxAttr2Constr(reduce_bounds = FALSE)@reduce_bounds)
})

## @cvxpy NONE
test_that("reduce_bounds = TRUE lowers box bounds to constraints", {
  x <- Variable(2, bounds = list(-1, 3))
  p <- Problem(Minimize(sum_entries(x)), list())
  new_p <- reduction_apply(CvxAttr2Constr(reduce_bounds = TRUE), p)[[1]]
  ## bounds become >= and <= inequality constraints
  expect_equal(length(new_p@constraints), 2L)
  ## and the reduced variable no longer carries the bounds attribute
  v <- variables(new_p)[[1]]
  expect_null(v@attributes$bounds)
})

## @cvxpy NONE
test_that("reduce_bounds = FALSE keeps bounds on the variable, emits no constraints", {
  x <- Variable(2, bounds = list(-1, 3))
  p <- Problem(Minimize(sum_entries(x)), list())
  new_p <- reduction_apply(CvxAttr2Constr(reduce_bounds = FALSE), p)[[1]]
  ## No bound constraints added ...
  expect_equal(length(new_p@constraints), 0L)
  ## ... and the variable still carries the bounds attribute, readable via
  ## get_bounds() (the prior V19 leaf slice the NLP path consumes).
  v <- variables(new_p)[[1]]
  expect_false(is.null(v@attributes$bounds))
  gb <- get_bounds(v)
  ## get_bounds returns shaped (matrix) bounds as of Wave 2 (#3080).
  expect_equal(as.numeric(gb[[1]]), c(-1, -1))
  expect_equal(as.numeric(gb[[2]]), c(3, 3))
})

## @cvxpy NONE
test_that("reduce_bounds = FALSE preserves sign attribute instead of lowering it", {
  x <- Variable(2, nonneg = TRUE)
  p <- Problem(Minimize(sum_entries(x)), list())
  new_p <- reduction_apply(CvxAttr2Constr(reduce_bounds = FALSE), p)[[1]]
  expect_equal(length(new_p@constraints), 0L)
  v <- variables(new_p)[[1]]
  expect_true(isTRUE(v@attributes$nonneg))
})

# ===================================================================
# Dimension-reducing parameter attributes (CVXPY 1.9.0 -- 2D CVXR slice)
# CVXPY SOURCE: cvx_attr2constr.py:102-164,278-290,303-318
# ===================================================================

## @cvxpy reductions/cvx_attr2constr.py::build_dim_reduced_expression
test_that("PSD Parameters are reduced to upper-triangle Parameters", {
  P <- Parameter(c(2, 2), PSD = TRUE, value = matrix(c(2, 1, 1, 3), 2, 2))
  X <- Variable(c(2, 2))
  p <- Problem(Minimize(sum_entries(X + P)), list())

  ca <- CvxAttr2Constr()
  new_p <- reduction_apply(ca, p)[[1]]
  params <- parameters(new_p)

  expect_equal(length(params), 1L)
  expect_equal(params[[1L]]@shape, c(3L, 1L))
  expect_equal(as.numeric(value(params[[1L]])), c(2, 1, 3))
  expect_true(is_psd(new_p@objective@args[[1L]]@args[[1L]]@args[[2L]]))
})

## @cvxpy reductions/cvx_attr2constr.py::update_parameters
test_that("CvxAttr2Constr refreshes reduced Parameter values", {
  P <- Parameter(c(2, 2), PSD = TRUE, value = matrix(c(2, 1, 1, 3), 2, 2))
  X <- Variable(c(2, 2))
  p <- Problem(Minimize(sum_entries(X + P)), list())

  ca <- CvxAttr2Constr()
  reduction_apply(ca, p)
  reduced <- get(as.character(P@id), envir = ca@.cache$parameters)
  expect_equal(as.numeric(value(reduced)), c(2, 1, 3))

  value(P) <- matrix(c(5, 2, 2, 7), 2, 2)
  update_parameters(ca, p)
  expect_equal(as.numeric(value(reduced)), c(5, 2, 7))

  ## Dict-in/dict-out (#3147 part A): keyed by leaf id, returns a dict.
  pid <- as.character(P@id); rid <- as.character(reduced@id)
  pf_in <- list(); pf_in[[pid]] <- matrix(c(1, 4, 4, 9), 2, 2)
  fwd <- param_forward(ca, pf_in)
  expect_equal(as.numeric(fwd[[rid]]), c(1, 4, 9))

  pb_in <- list(); pb_in[[rid]] <- c(1, 4, 9)
  bwd <- param_backward(ca, pb_in)
  expect_equal(as.matrix(bwd[[pid]]), matrix(c(1, 4, 4, 9), 2, 2))
})

## @cvxpy reductions/cvx_attr2constr.py::lower_value
test_that("diagonal Parameters are reduced to diagonal vectors", {
  D <- Parameter(c(3, 3), diag = TRUE, value = diag(1:3))
  X <- Variable(c(3, 3))
  p <- Problem(Minimize(sum_entries(X + D)), list())

  ca <- CvxAttr2Constr()
  new_p <- reduction_apply(ca, p)[[1]]
  params <- parameters(new_p)

  expect_equal(length(params), 1L)
  expect_equal(params[[1L]]@shape, c(3L, 1L))
  expect_equal(as.numeric(value(params[[1L]])), 1:3)
  expect_true(S7_inherits(new_p@objective@args[[1L]]@args[[1L]]@args[[2L]], DiagVec))
})
