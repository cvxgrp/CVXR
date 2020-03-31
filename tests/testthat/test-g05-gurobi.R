context("test-g05-gurobi")

GUROBI_AVAILABLE <- "GUROBI" %in% installed_solvers()
if (GUROBI_AVAILABLE) {
    library(gurobi)
    library(Matrix)
}

test_that("test a simple mixed integer program for GUROBI", {
  skip_on_cran()
  skip_if_not(GUROBI_AVAILABLE, "Skipping GUROBI test as it is not available.!")
  # Example from
  # https://www.gurobi.com/documentation/9.0/examples/mip_r.html
  model <- list()
  model$A          <- matrix(c(1,2,3,1,1,0), nrow=2, ncol=3, byrow=T)
  model$obj        <- c(1,1,2)
  model$modelsense <- 'max'
  model$rhs        <- c(4,1)
  model$sense      <- c('<', '>')
  model$vtype      <- 'B'

  params <- list(OutputFlag=0)

  gurobiResult <- gurobi(model, params)

  xvar <- Variable(3, boolean=T)
  obj <- Maximize(xvar[1] + xvar[2] + 2 * xvar[3])
  constraint <- list(xvar[1] + 2*xvar[2] + 3*xvar[3] <= 4,
                     xvar[1] + xvar[2] >= 1)
  prob <- Problem(obj, constraint)
  cvxr <- solve(prob, solver="GUROBI")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, gurobiResult$objval, tolerance = 1e-4)
  expect_equal(cvxr$getValue(xvar), matrix(gurobiResult$x), tolerance = 1e-4)

})

test_that("test a simple QP for GUROBI",{
  skip_on_cran()
  skip_if_not(GUROBI_AVAILABLE, "Skipping GUROBI test as it is not available.!")
  # Example from
  # https://www.gurobi.com/documentation/9.0/examples/qp_r.html

  model <- list()

  model$A     <- matrix(c(1,2,3,1,1,0), nrow=2, byrow=T)
  model$Q     <- matrix(c(1,0.5,0,0.5,1,0.5,0,0.5,1), nrow=3, byrow=T)
  model$obj   <- c(2,0,0)
  model$rhs   <- c(4,1)
  model$sense <- c('>', '>')

  result <- gurobi(model)
  gurobiResult <- gurobi(model)

  xvar <- Variable(3, nonneg=T)
  Pmat <- matrix(c(1,0.5,0,0.5,1,0.5,0,0.5,1), nrow=3, byrow=T)
  obj <- Minimize(quad_form(xvar, Pmat) + 2 * xvar[1])
  constraint <- list(xvar[1] + 2*xvar[2] + 3*xvar[3] >= 4,
                     xvar[1] + xvar[2] >= 1)
  prob <- Problem(obj, constraint)
  cvxr <- solve(prob, solver="GUROBI")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, gurobiResult$objval, tolerance = 1e-4)
  expect_equal(cvxr$getValue(xvar), matrix(gurobiResult$x), tolerance = 1e-4)

})

# SOC/PSD constraints not supported by GUROBI
# test_that("test a simple QCP for GUROBI",{
#   skip_on_cran()
#   skip_if_not(GUROBI_AVAILABLE, "Skipping GUROBI test as it is not available.!")
#   # Example from
#   # https://www.gurobi.com/documentation/9.0/examples/qcp_r.html
#
#   model <- list()
#
#   model$A          <- matrix(c(1,1,1), nrow=1, byrow=T)
#   model$modelsense <- 'max'
#   model$obj        <- c(1,0,0)
#   model$rhs        <- c(1)
#   model$sense      <- c('=')
#
#   # First quadratic constraint: x^2 + y^2 - z^2 <= 0
#   qc1 <- list()
#   qc1$Qc <- spMatrix(3, 3, c(1, 2, 3), c(1, 2, 3), c(1.0, 1.0, -1.0))
#   qc1$rhs <- 0.0
#
#   # Second quadratic constraint: x^2 - yz <= 0
#   qc2 <- list()
#   qc2$Qc <- spMatrix(3, 3, c(1, 2), c(1, 3), c(1.0, -1.0))
#   qc2$rhs <- 0.0
#
#   model$quadcon <- list(qc1, qc2)
#
#   gurobiResult <- gurobi(model)
#
#   xvar <- Variable(3, nonneg=T)
#   obj <- Maximize(xvar[1])
#   constraint <- list( xvar[1] + xvar[2] + xvar[3] == 1,
#                       CVXR:::SOC(xvar[3], xvar[1:2]),
#                       CVXR:::PSDConstraint(rbind(cbind(xvar[2], xvar[1]), cbind(xvar[1], xvar[3])))
#                       )
#   prob <- Problem(obj, constraint)
#   cvxr <- solve(prob, solver="GUROBI")
# })

# SOC/PSD constraints not supported by GUROBI
# test_that("test a piecewise linear problem for GUROBI", {
#   skip_on_cran()
#   skip_if_not(GUROBI_AVAILABLE, "Skipping GUROBI test as it is not available.!")
#   # Example from
#   # https://www.gurobi.com/documentation/9.0/examples/piecewise_r.html#subsubsection:piecewise.R
#   model <- list()
#
#   model$A     <- matrix(c(1,2,3,1,1,0), nrow=2, byrow=T)
#   model$obj   <- c(0,-1,0)
#   model$ub    <- c(1,1,1)
#   model$rhs   <- c(4,1)
#   model$sense <- c('<', '>')
#
#   # Uniformly spaced points in [0.0, 1.0]
#   u <- seq(from=0, to=1, by=0.01)
#
#   # First piecewise-linear function: f(x) = exp(-x)
#   pwl1     <- list()
#   pwl1$var <- 1
#   pwl1$x   <- u
#   pwl1$y   <- sapply(u, function(x) exp(-x))
#
#   # Second piecewise-linear function: g(z) = 2 z^2 - 4 z
#   pwl2     <- list()
#   pwl2$var <- 3
#   pwl2$x   <- u
#   pwl2$y   <- sapply(u, function(z) 2 * z * z - 4 * z)
#
#   model$pwlobj <- list(pwl1, pwl2)
#
#   gurobiResult <- gurobi(model)
#
#   xvar <- Variable(3)
#   obj <- Minimize(CVXR:::Exp(xvar[1]) - xvar[2] + 2 *square(xvar[3]) - 4 * xvar[3])
#   constraint <- list(xvar[1] + 2 * xvar[2] + 3 * xvar[3] <= 4,
#                      xvar[1] + xvar[2] >= 1,
#                      xvar[1] <= 1,
#                      xvar[2] <= 1,
#                      xvar[3] <= 1)
#
#   prob <- Problem(obj, constraint)
#   cvxr <- solve(prob, solver="GUROBI")
#
# })
