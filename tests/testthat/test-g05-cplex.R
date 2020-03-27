context("test-g05-cplex")

tol <- 1e-4

CPLEX_AVAILABLE <- "CPLEX" %in% installed_solvers()
if (CPLEX_AVAILABLE) library(Rcplex)

test_that("test a simple LP with CPLEX", {
  # Example lpex1.m in CPLEX
  skip_on_cran()
  skip_if_not(CPLEX_AVAILABLE, "Skipping CPLEX test as it is not available.!")
  cvec <- c(1,2,3)
  Amat <- matrix(c(-1,1,1,-1,3,-1),byrow=TRUE,nc=3)
  bvec <- c(20, -30)
  ub <- c(40, Inf, Inf)
  cplex <- Rcplex(cvec, Amat, bvec, ub=ub, objsense="max", sense=c('L','G'))

  xvar <- Variable(3, nonneg=T)
  obj <- Maximize(t(as.matrix(cvec)) %*% xvar)
  constraints <- list(-xvar[1] + xvar[2] + xvar[3] <= 20,
                      xvar[1] - 3 * xvar[2] + xvar[3] <= 30,
                      xvar[1] <= 40)
  prob <- Problem(obj, constraints)
  cvxr <- solve(prob, solver = "CPLEX")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, cplex$obj, tolerance = tol)
  expect_equal(cvxr$getValue(xvar), matrix(cplex$xopt), tolerance = tol)
})

test_that("test a simple QP with CPLEX", {
  # Example qpex1.m in CPLEX
  skip_on_cran()
  skip_if_not(CPLEX_AVAILABLE, "Skipping CPLEX test as it is not available.!")
  cvec <- c(1,2,3)
  Qmat <- matrix(c(-33, 6, 0, 6,-22, 11.5, 0, 11.5, -11),
                 byrow=TRUE, nc=3)
  Amat <- matrix(c(-1, 1, 1, 1, -3, 1),
                 byrow=TRUE, nc=3)
  bvec <- c(20, 30)
  ub <- c(40, Inf, Inf)
  cplex <- Rcplex(cvec, Amat, bvec, Qmat, ub=ub, objsense="max")

  xvar <- Variable(3, nonneg=T)
  obj <- Maximize(t(as.matrix(cvec)) %*% xvar + .5 * quad_form(xvar, Qmat))
  constraints <- list(-xvar[1] + xvar[2] + xvar[3] <= 20,
                      xvar[1] - 3 * xvar[2] + xvar[3] <= 30,
                      xvar[1] <= 40)
  prob <- Problem(obj, constraints)
  cvxr <- solve(prob, solver = "CPLEX")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, cplex$obj, tolerance = tol)
  expect_equal(cvxr$getValue(xvar), matrix(cplex$xopt), tolerance = tol)
})

test_that("test a mixed integer linear program with CPLEX", {
  # Example mipex1.m in CPLEX
  skip_on_cran()
  skip_if_not(CPLEX_AVAILABLE, "Skipping CPLEX test as it is not available.!")
  cvec <- c(1,2,3,1)
  Amat <- matrix(c(-1, 1, 1, 10, 1, -3,
                   1, 0, 0, 1, 0, -3.5), byrow = T, nc = 4)
  bvec <- c(20, 30, 0)
  lb <- c(0, 0, 0, 2)
  ub <- c(40, Inf, Inf, 3)
  vtype <- c(rep("C", 3), "I")
  cplex <- Rcplex(cvec, Amat, bvec, lb=lb, ub=ub, sense=c("L","L","E"),
                objsense="max", vtype=vtype)

  xnonint <- Variable(3, nonneg=T)
  xint <- Variable(integer = T)
  xvar <- rbind(xnonint, xint)
  obj <- Maximize(t(as.matrix(cvec)) %*% xvar)
  constraints <- list(-xvar[1] + xvar[2] + xvar[3] + 10 * xvar[4] <= 20,
                      xvar[1] - 3 * xvar[2] + xvar[3] <= 30,
                      xvar[2] - 3.5 * xvar[4] == 0,
                      xvar[1] <= 40,
                      xvar[4] >= 2,
                      xvar[4] <= 3)
  prob <- Problem(obj, constraints)
  cvxr <- solve(prob, solver = "CPLEX")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, cplex$obj, tolerance = tol)
  expect_equal(cvxr$getValue(xvar), matrix(cplex$xopt), tolerance = tol)
})

test_that("test a mixed integer quadratic program with CPLEX", {
  # Example miqpex1.m in CPLEX
  skip_on_cran()
  skip_if_not(CPLEX_AVAILABLE, "Skipping CPLEX test as it is not available.!")
  cvec <- c(1,2,3,1)
  Qmat <- matrix(c(-33, 6, 0, 0, 6,-22, 11.5, 0,
                   0, 11.5, -11,0, 0, 0, 0, 0),
                 byrow=TRUE, nc=4)
  Amat <- matrix(c(-1, 1, 1, 10,
                   1, -3, 1, 0,
                   0, 1, 0,-3.5),
                 byrow=TRUE, nc=4)
  bvec <- c(20, 30, 0)
  ub <- c(40, Inf, Inf, 3)
  vtype <- c(rep("C",3), "I")
  cplex <- Rcplex(cvec, Amat, bvec, Qmat=Qmat,
                  ub=ub, sense=c("L","L","E"),
                objsense="max", vtype=vtype)

  xnonint <- Variable(3, nonneg=T)
  xint <- Variable(integer = T)
  xvar <- rbind(xnonint, xint)
  obj <- Maximize(t(as.matrix(cvec)) %*% xvar + .5 * quad_form(xvar, Qmat))
  constraints <- list(-xvar[1] + xvar[2] + xvar[3] + 10 * xvar[4] <= 20,
                      xvar[1] - 3 * xvar[2] + xvar[3] <= 30,
                      xvar[2] - 3.5 * xvar[4] == 0,
                      xvar[1] <= 40,
                      xvar[4] >= 0,
                      xvar[4] <= 3)
  prob <- Problem(obj, constraints)
  cvxr <- solve(prob, solver = "CPLEX")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, cplex$obj, tolerance = tol)
  expect_equal(cvxr$getValue(xvar), matrix(cplex$xopt), tolerance = tol)
})
