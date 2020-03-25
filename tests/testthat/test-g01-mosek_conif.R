context("test-g01-mosek_conif")
library("Rmosek")
TOL <- 1e-6

test_that("test MOSEK SOCP", {
  skip_on_cran()
  # Formulate the following SOCP with CVXR
  #    min 3 * x[1] + 2 * x[2] + x[3]
  #       s.t. p_norm(x,2) <= y[1]
  #            p_norm(x,2) <= y[2]
  #            x[1] + x[2] + 3*x[3] >= 1.0
  #            y <= 5
  # and solve with MOSEK and ECOS. Compare MOSEK and ECOS primal and dual solutions.

  if("MOSEK" %in% installed_solvers()) {
    x <- Variable(3)
    y <- Variable(2)
    constraints <- list(p_norm(x, 2) <= y[1],
                        p_norm(x, 2) <= y[2],
                        x[1] + x[2] + 3 * x[3] >= 1.0,
                        y <= 5)
    obj <- Minimize(3 * x[1] + 2 * x[2] + x[3])
    prob <- Problem(obj, constraints)
    result_mosek <- solve(prob, solver = "MOSEK")
    x_mosek <- result_mosek$getValue(x)
    duals_mosek <- lapply(constraints, function(c) { result_mosek$getDualValue(c) })

    result_ecos <- solve(prob, solver = "ECOS")
    x_ecos <- result_ecos$getValue(x)
    duals_ecos <- lapply(constraints, function(c) { result_ecos$getDualValue(c) })

    expect_equal(x_mosek, x_ecos, tolerance = TOL)
    expect_equal(length(duals_ecos), length(duals_mosek))
    for(i in seq_along(duals_mosek))
      expect_equal(duals_mosek[[i]], duals_ecos[[i]], tolerance = 1e-4)
  }
})

test_that("test MOSEK SDP", {
  skip_on_cran()
  # Solve "Example 8.3" from Convex Optimization by Boyd & Vandenberghe.
  # Verify (1) optimal objective values, (2) that the dual variable to the PSD constraint
  # belongs to the correct cone (i.e. the dual variable is itself PSD), and (3) that
  # complementary slackness holds with the PSD primal variable and its dual variable.

  if("MOSEK" %in% installed_solvers()) {
    # This is an example from Convex Optimization by B&V.
    # Example 8.3 (page 408 in the 19th printing).
    rho <- Variable(4, 4)
    constraints <- list(0.6 <= rho[1, 2], rho[1, 2] <= 0.9,
                        0.8 <= rho[1, 3], rho[1, 3] <= 0.9,
                        0.5 <= rho[2, 4], rho[2, 4] <= 0.7,
                       -0.8 <= rho[3, 4], rho[3, 4] <= -0.4,
                        rho[1, 1] == 1, rho[2, 2] == 1, rho[3, 3] == 1, rho[4, 4] == 1,
                        rho %>>% 0)

    # First, minimize rho[1, 4]
    obj <- Minimize(rho[1, 4])
    prob <- Problem(obj, constraints)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(round(result$value, 2), -0.39)
    Y <- result$getDualValue(constraints[[length(constraints)]])
    eigs <- eigen(Y, only.values = TRUE)
    expect_true(all(eigs$value >= 0))
    complementary_slackness <- sum(diag(result$getValue(rho) %*% Y))
    expect_equal(complementary_slackness, 0.0, tolereance = TOL)

    # Then, maximize rho[1, 4]
    obj <- Maximize(rho[1, 4])
    prob <- Problem(obj, constraints)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(round(result$value, 2), 0.23)
    Y <- result$getDualValue(constraints[[length(constraints)]])
    eigs <- eigen(Y, only.values = TRUE)
    expect_true(all(eigs$value >= 0))
    complementary_slackness <- sum(diag(result$getValue(rho) %*% Y))
    expect_equal(complementary_slackness, 0.0, tolerance = TOL)
  }
})

# test_that("test MOSEK exponential cone", {
#   # Formulate the following exponential cone problem with cvxpy
#   #   min   3 * x[1] + 2 * x[2] + x[3]
#   #     s.t.  0.1 <= x[1] + x[2] + x[3] <= 1
#   #           x >= 0
#   #           x[1] >= x[2] * exp(x[3] / x[2])
#   # and solve with MOSEK and ECOS. Ensure that MOSEK and ECOS have the same
#   # primal and dual solutions.
#   #
#   # Note that the exponential cone constraint can be rewritten in terms of the
#   # relative entropy cone. The correspondence is as follows:
#   #     x[1] >= x[2] * exp(x[3] / x[2])
#   #       iff
#   #     x[2] * log(x[2] / x[1]) + x[3] <= 0.
#
#   if("MOSEK" %in% installed_solvers()) {
#     # TODO: Only run if exponential cone is supported.
#     # Formulate and solve the problem with CVXR
#     x <- Variable(3, 1)
#     constraints <- list(sum(x) <= 1.0, sum(x) >= 0.1, x >= 0.01,
#                         kl_div(x[2], x[1]) + x[2] - x[1] + x[3] <= 0)
#     obj <- Minimize(3 * x[1] + 2 * x[2] + x[1])
#     prob <- Problem(obj, constraints)
#     result_mosek <- solve(prob, solver = "MOSEK")
#     val_mosek = result_mosek$value
#     x_mosek = result_mosek$getValue(x)
#     duals_mosek = lapply(constraints, function(c) { result_mosek$getDualValue(c) })
#
#     result_ecos <- solve(prob, solver = "ECOS")
#     val_ecos = result_ecos$value
#     x_ecos = result_ecos$getValue(x)
#     duals_ecos = lapply(constraints, function(c) { result_ecos$getDualValue(c) })
#
#     # Verify results
#     expect_equal(val_mosek, val_ecos, tolerance = TOL)
#     expect_equal(x_mosek, x_ecos, tolerance = 1e-4)
#     expect_equal(length(duals_ecos), length(duals_mosek))
#     for(i in seq_along(duals_mosek))
#       expect_equal(duals_mosek[[i]], duals_ecos[[i]], tolerance = 1e-4)
#   }
# })

test_that("test MOSEK MI SOCP", {
  skip_on_cran()
  # Formulate the following mixed-integer SOCP with CVXR
  #   min 3 * x[1] + 2 * x[2] + x[3] + y[1] + 2 * y[2]
  #     s.t. p_norm(x,2) <= y[1]
  #          p_norm(x,2) <= y[2]
  #          x[1] + x[2] + 3*x[3] >= 0.1
  #          y <= 5, y integer.
  # and solve with MOSEK.

  if("MOSEK" %in% installed_solvers()) {
    x <- Variable(3)
    y <- Variable(2, integer = TRUE)
    constraints <- list(p_norm(x, 2) <= y[1],
                        p_norm(x, 2) <= y[2],
                        x[1] + x[2] + 3 * x[3] >= 0.1,
                        y <= 5)
    obj <- Minimize(3 * x[1] + 2 * x[2] + x[3] + y[1] + 2 * y[2])
    prob <- Problem(obj, constraints)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$getValue(y), as.matrix(c(1,1)), tolerance = TOL)
  }
})

test_that("test MOSEK LP solution selection", {
  skip_on_cran()
  if("MOSEK" %in% installed_solvers()) {
    # Problem from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp

    x <- Variable(2)
    objective <- Minimize(-4 * x[1] - 5 * x[2])
    constraints <- list(2 * x[1] + x[2] <= 3,
                        x[1] + 2 * x[2] <= 3,
                        x[1] >= 0, x[2] >= 0)
    prob <- Problem(objective, constraints)

    # Default solution (interior point)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1,1)), tolerance = 1e-5)

    # NOTE: Rmosek does not provide the bfs option.
    # Basic feasible solution
    # result <- solve(prob, solver = "MOSEK", bfs = TRUE)
    # expect_equal(result$value, -9, tolerance = 1e-10)
    # expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = 1e-10)
  }
})

test_that("test MOSEK affco1 example",{
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
                       norm2(xvar[1:2]) <= xvar[4],
                       CVXR:::PSDConstraint(rbind(cbind(2*xvar[5], xvar[3]), cbind(xvar[3], xvar[6])))
                       )
  prob <- Problem(obj, constraints)
  cvxr <- solve(prob, solver = "MOSEK")
  
  expect_equal(cvxr$status, "optimal")
  expect_equal(cvxr$value, rmosek$itr$pobjval, tolerance = 1e-4)
  # Answers are slightly different but objective value is same and constraints satisfied
  #expect_equal(cvxr$getValue(xvar), matrix(rmosek$itr$xx), tolerance = 1e-2)
})