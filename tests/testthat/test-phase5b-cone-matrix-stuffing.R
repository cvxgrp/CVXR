# Tests for Phase 5b: ConeMatrixStuffing + ConeDims

library(CVXR)

# ===================================================================
# ConeDims
# ===================================================================

## @cvxpy NONE
test_that("ConeDims from empty constraint map", {
  cmap <- group_constraints(list())
  cd <- ConeDims(cmap)
  expect_equal(cd@zero, 0L)
  expect_equal(cd@nonneg, 0L)
  expect_equal(cd@exp, 0L)
  expect_equal(length(cd@soc), 0L)
  expect_equal(length(cd@psd), 0L)
})

## @cvxpy NONE
test_that("ConeDims with NonNeg constraints", {
  x <- Variable(3)
  nn <- NonNeg(x)
  cmap <- group_constraints(list(nn))
  cd <- ConeDims(cmap)
  expect_equal(cd@zero, 0L)
  expect_equal(cd@nonneg, 3L)
})

## @cvxpy NONE
test_that("ConeDims with Zero + NonNeg", {
  x <- Variable(2)
  z <- Zero(x - Constant(c(1, 2)))
  nn <- NonNeg(x)
  cmap <- group_constraints(list(z, nn))
  cd <- ConeDims(cmap)
  expect_equal(cd@zero, 2L)
  expect_equal(cd@nonneg, 2L)
})

## @cvxpy NONE
test_that("ConeDims with SOC", {
  t_var <- Variable(1)
  x <- Variable(3)
  soc <- SOC(t_var, x)
  cmap <- group_constraints(list(soc))
  cd <- ConeDims(cmap)
  expect_equal(cd@soc, 4L)  # 1 (t) + 3 (x) = 4
})

## @cvxpy NONE
test_that("dims_to_solver_dict", {
  x <- Variable(3)
  nn <- NonNeg(x)
  cmap <- group_constraints(list(nn))
  cd <- ConeDims(cmap)
  d <- dims_to_solver_dict(cd)
  expect_equal(d$f, 0L)
  expect_equal(d$l, 3L)
  expect_equal(d$ep, 0L)
})

# ===================================================================
# ConeMatrixStuffing
# ===================================================================

## @cvxpy NONE
test_that("ConeMatrixStuffing accepts affine problem", {
  x <- Variable(3)
  p <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  cms <- ConeMatrixStuffing()
  expect_true(reduction_accepts(cms, p))
})

## @cvxpy NONE
test_that("ConeMatrixStuffing rejects non-affine", {
  x <- Variable(1)
  ## exp(x) is not affine, so direct ConeMatrixStuffing should reject
  ## (unless Dcp2Cone has already transformed it)
  p <- Problem(Minimize(exp(x)))
  cms <- ConeMatrixStuffing()
  expect_false(reduction_accepts(cms, p))
})

## @cvxpy NONE
test_that("ConeMatrixStuffing: LP min sum(x) s.t. x >= 1", {
  x <- Variable(3)
  ## After Dcp2Cone, this produces NonNeg(x - 1)
  p_dcp <- reduction_apply(Dcp2Cone(), Problem(Minimize(sum_entries(x)), list(x >= 1)))[[1]]
  cms <- ConeMatrixStuffing()
  result <- reduction_apply(cms, p_dcp)
  data <- result[[1]]

  ## ConeMatrixStuffing returns PRE-negation data: A*x + b >= 0 for NonNeg
  ## For NonNeg(x - 1): A = I, b = -1 (A*x + b = x - 1 >= 0)
  ## ConicSolver later negates A for solver convention.
  expect_equal(data$c, c(1, 1, 1))
  expect_equal(as.numeric(data$b), c(-1, -1, -1))
  A <- as.matrix(data$A)
  expect_equal(A, diag(3))
  expect_equal(data$dims@nonneg, 3L)
  expect_equal(data$dims@zero, 0L)
})

## @cvxpy NONE
test_that("ConeMatrixStuffing: LP min 2*x s.t. x>=0, x<=5", {
  x <- Variable(1)
  p <- Problem(Minimize(2 * x), list(x >= 0, x <= 5))
  p_dcp <- reduction_apply(Dcp2Cone(), p)[[1]]
  result <- reduction_apply(ConeMatrixStuffing(), p_dcp)
  data <- result[[1]]

  ## Pre-negation: NonNeg(x) has A=I, b=0; NonNeg(5-x) has A=-I, b=5
  expect_equal(data$c, 2.0)
  expect_equal(as.numeric(data$b), c(0, 5))
  A <- as.matrix(data$A)
  expect_equal(A[1, 1], 1)    # x >= 0: A*x + 0 >= 0
  expect_equal(A[2, 1], -1)   # x <= 5: -x + 5 >= 0
  expect_equal(data$dims@nonneg, 2L)
})

## @cvxpy NONE
test_that("ConeMatrixStuffing: LP with equality", {
  xy <- Variable(2)
  p <- Problem(Minimize(sum_entries(xy)), list(xy[1] + xy[2] == 5, xy >= 0))
  p_dcp <- reduction_apply(Dcp2Cone(), p)[[1]]
  result <- reduction_apply(ConeMatrixStuffing(), p_dcp)
  data <- result[[1]]

  ## Pre-negation: Zero(xy[1]+xy[2]-5) has A=[1,1], b=-5
  ## NonNeg(xy) has A=I, b=0
  expect_equal(data$c, c(1, 1))
  expect_equal(as.numeric(data$b), c(-5, 0, 0))
  expect_equal(data$dims@zero, 1L)
  expect_equal(data$dims@nonneg, 2L)
})

## @cvxpy NONE
test_that("ConeMatrixStuffing: SOC constraint", {
  xv <- Variable(1)
  yz <- Variable(2)
  p <- Problem(Minimize(xv), list(SOC(xv, yz)))
  p_dcp <- reduction_apply(Dcp2Cone(), p)[[1]]
  result <- reduction_apply(ConeMatrixStuffing(), p_dcp)
  data <- result[[1]]

  ## CVXPY: c=[1,0,0], A=-I_3, b=[0,0,0], dims={soc:[3]}
  expect_equal(length(data$c), 3L)
  expect_equal(data$c[1], 1)
  expect_equal(sum(abs(data$c[2:3])), 0)
  expect_equal(as.numeric(data$b), c(0, 0, 0))
  expect_equal(data$dims@soc, 3L)
})

## @cvxpy NONE
test_that("ConeMatrixStuffing constant offset", {
  x <- Variable(1)
  p <- Problem(Minimize(x + 10))
  p_dcp <- reduction_apply(Dcp2Cone(), p)[[1]]
  result <- reduction_apply(ConeMatrixStuffing(), p_dcp)
  data <- result[[1]]

  expect_equal(data$c, 1.0)
  expect_equal(as.numeric(data$offset), 10.0)
})
