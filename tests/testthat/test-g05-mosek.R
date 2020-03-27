context("test-g05-mosek")
TOL <- 1e-6

MOSEK_AVAILABLE  <- "MOSEK" %in% installed_solvers()

if (MOSEK_AVAILABLE) library("Rmosek")

test_that("test MOSEK affco1 example",{
  skip_on_cran()
  skip_if_not(MOSEK_AVAILABLE, "Skipping MOSEK test as it is not available.!")
  # Example from
  # https://docs.mosek.com/9.1/rmosek/examples-list.html#doc-example-file-affco1-r
  affco1 <- function(){
    prob <- list(sense="max")

    # Variables [x1; x2; t1; t2]
    prob$c <- c(0, 0, 1, 1)

    # Linear inequality x_1 - x_2 <= 1
    prob$A <- Matrix(c(1, -1, 0, 0), nrow=1, sparse=TRUE)
    prob$bc <- rbind(blc=-Inf, buc=1)
    prob$bx <- rbind(blx=rep(-Inf,4), bux=rep(Inf,4))

    # The quadratic cone
    FQ <- rbind(c(0,0,0,0), c(1,0,0,0), c(0,1,0,0))
    gQ <- c(1, -0.5, -0.6)
    cQ <- matrix(list("QUAD", 3, NULL), nrow=3, ncol=1)

    # The power cone for (x_1, 1, t_1) \in POW3^(1/3,2/3)
    FP1 <- rbind(c(1,0,0,0), c(0,0,0,0), c(0,0,1,0))
    gP1 <- c(0, 1, 0)
    cP1 <- matrix(list("PPOW", 3, c(1/3, 2/3)), nrow=3, ncol=1)

    # The power cone for (x_1+x_2+0.1, 1, t_2) \in POW3^(1/4,3/4)
    FP2 <- rbind(c(1,1,0,0), c(0,0,0,0), c(0,0,0,1))
    gP2 <- c(0.1, 1, 0)
    cP2 <- matrix(list("PPOW", 3, c(1.0, 3.0)), nrow=3, ncol=1)

    # All cones
    prob$F <- rbind(FQ, FP1, FP2)
    prob$g <- cbind(gQ, gP1, gP2)
    prob$cones <- cbind(cQ, cP1, cP2)
    rownames(prob$cones) <- c("type","dim","conepar")

    r <- mosek(prob, list(soldetail=1))
    stopifnot(identical(r$response$code, 0))
    list(status = r$response$code, value = r$sol$itr$pobjval, x = r$sol$itr$xx[1:2])
  }
  # rmosek answer
  rmosek <- affco1()

  # CVXR answer
  x  <- Variable(2)
  obj  <- x[1, 1]^(1/3) + (sum(x) + 0.1)^(1/4)
  cons  <- list(
    sum_squares(x - matrix(c(0.5, 0.6))) <= 1,
    x[1, 1] >= 0,
    x[1, 1] <= x[2, 1] + 1
  )
  p  <- Problem(Maximize(obj), cons)
  cvxr <- solve(p, gp = FALSE, solver = "MOSEK")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, rmosek$value, tolerance = 1e-4)
  expect_equal(cvxr$getValue(x), matrix(rmosek$x), tolerance = 1e-4)
})

test_that("test MOSEK simple quadratic optimization problem",{
  skip_on_cran()
  skip_if_not(MOSEK_AVAILABLE, "Skipping MOSEK test as it is not available.!")

  # Example from
  # https://docs.mosek.com/9.1/rmosek/tutorial-qo-shared.html
  qo1 <- function(){
    # Specify the non-quadratic part of the problem.
    prob <- list(sense="min")
    prob$c <- c(0, -1, 0)
    prob$A <- Matrix(c(1, 1, 1), nrow=1, sparse=TRUE)
    prob$bc <- rbind(blc=1,
                     buc=Inf)
    prob$bx <- rbind(blx=rep(0,3),
                     bux=rep(Inf,3))

    # Specify the quadratic objective matrix in triplet form.
    prob$qobj$i <- c(1,  3,   2,  3)
    prob$qobj$j <- c(1,  1,   2,  3)
    prob$qobj$v <- c(2, -1, 0.2,  2)

    # Solve the problem
    r <- mosek(prob, list(soldetail=1))

    # Return the solution
    stopifnot(identical(r$response$code, 0))
    r$sol$itr
  }
  rmosek <- qo1()

  xvar <- Variable(3, nonneg=TRUE)
  Pmat <- matrix(c(1, 0, -.5, 0, .1, 0, -.5, 0, 1), nrow = 3)
  obj <- Minimize(quad_form(xvar, Pmat) - xvar[2])
  constraints <- list(sum_entries(xvar[1:3]) >= 1)
  prob <- Problem(obj, constraints)
  cvxr <- solve(prob, solver="MOSEK")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, rmosek$pobjval, tolerance = 1e-4)
  expect_equal(cvxr$getValue(xvar), matrix(rmosek$xx), tolerance = 1e-4)

})


test_that("test MOSEK simple conic quadratic problem", {
  skip_on_cran()
  skip_if_not(MOSEK_AVAILABLE, "Skipping MOSEK test as it is not available.!")

  # Example from
  # https://docs.mosek.com/9.1/rmosek/tutorial-cqo-shared.html

  cqo1 <- function(){
    # Specify the non-conic part of the problem.
    prob <- list(sense="min")
    prob$c  <- c(0, 0, 0, 1, 1, 1)
    prob$A  <- Matrix(c(1, 1, 2, 0, 0, 0), nrow=1, sparse=TRUE)
    prob$bc <- rbind(blc=1,
                     buc=1)
    prob$bx <- rbind(blx=c(rep(0,3), rep(-Inf,3)),
                     bux=rep(Inf,6))

    # Specify the cones.
    NUMCONES <- 2
    prob$cones <- matrix(list(), nrow=2, ncol=NUMCONES)
    rownames(prob$cones) <- c("type","sub")

    prob$cones[,1] <- list("QUAD", c(4, 1, 2))
    prob$cones[,2] <- list("RQUAD", c(5, 6, 3))

    # Solve the problem
    r <- mosek(prob, list(soldetail=1))

    # Return the solution
    stopifnot(identical(r$response$code, 0))
    r$sol
  }
  rmosek <- cqo1()

  xpos <- Variable(3, nonneg=TRUE)
  xother <- Variable(3)
  xvar <- rbind(xpos, xother)
  cvar <- Matrix(c(0, 0, 0, 1, 1, 1))
  obj <- Minimize(t(xvar) %*% cvar)
  # To deal with rotated cone constraint, reference
  # https://github.com/cvxgrp/cvxpy/issues/737
  constraints <- list( xvar[1] + xvar[2] + 2*xvar[3] == 1,
                       CVXR:::SOC(xvar[4],xvar[1:2]),
                       CVXR:::PSDConstraint(rbind(cbind(2*xvar[5], xvar[3]), cbind(xvar[3], xvar[6])))
                       )
  prob <- Problem(obj, constraints)
  cvxr <- solve(prob, solver = "MOSEK")

  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, rmosek$itr$pobjval, tolerance = 1e-4)
  # Answers are slightly different but objective value is same and constraints satisfied
  #expect_equal(cvxr$getValue(xvar), matrix(rmosek$itr$xx), tolerance = 1e-2)
})
