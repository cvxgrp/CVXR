context("test_conic_solvers")

INSTALLED_SOLVERS <- installed_solvers()

#######################
#                     #
#   Setup Functions   #
#                     #
#######################

test_ecos_setup <- function() {
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  TOL <- 1e-6
  
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  c <- Variable(name = "c")
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  z <- Variable(2, name = "z")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(TOL = TOL, a = a, b = b, c = c, x = x, y = y, z = z, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

test_scs_setup <- function() {
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  TOL <- 1e-2
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(TOL = TOL, x = x, y = y, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

test_clarabel_setup <- function() {
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  TOL <- 1e-6
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(TOL = TOL, x = x, y = y, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

test_mosek_setup <- function() {
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
}

test_cvxopt_setup <- function() {
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  TOL <- 1e-6
  
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  c <- Variable(name = "c")
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  z <- Variable(2, name = "z")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(TOL = TOL, a = a, b = b, c = c, x = x, y = y, z = z, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

test_cbc_setup <- function() {
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  TOL <- 1e-6
  
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  c <- Variable(name = "c")
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  z <- Variable(2, name = "z")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(TOL = TOL, a = a, b = b, c = c, x = x, y = y, z = z, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

test_cplex_setup <- function() {
  if("CPLEX" %in% INSTALLED_SOLVERS)
    require(Rcplex)
  
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  TOL <- 1e-4
  
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  c <- Variable(name = "c")
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  z <- Variable(2, name = "z")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(TOL = TOL, a = a, b = b, c = c, x = x, y = y, z = z, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

test_gurobi_setup <- function() {
  if("GUROBI" %in% INSTALLED_SOLVERS)
    require(gurobi)
  
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  TOL <- 1e-5
  
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  c <- Variable(name = "c")
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  z <- Variable(2, name = "z")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(TOL = TOL, a = a, b = b, c = c, x = x, y = y, z = z, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

test_xpress_setup <- function() {
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  c <- Variable(name = "c")
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  z <- Variable(2, name = "z")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(a = a, b = b, c = c, x = x, y = y, z = z, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

test_all_solvers_setup <- function() {
  rm(list = ls(), envir = .GlobalEnv)   # Clear global environment.
  
  TOL <- 1e-5
  
  a <- Variable(name = "a")
  b <- Variable(name = "b")
  c <- Variable(name = "c")
  
  x <- Variable(2, name = "x")
  y <- Variable(3, name = "y")
  z <- Variable(2, name = "z")
  
  A <- Variable(2, 2, name = "A")
  B <- Variable(2, 2, name = "B")
  C <- Variable(3, 2, name = "C")
  
  parms <- list(TOL = TOL, a = a, b = b, c = c, x = x, y = y, z = z, A = A, B = B, C = C)
  
  # Save to global environment.
  for(nam in names(parms))
    assign(nam, parms[[nam]], envir = .GlobalEnv)
}

#######################
#                     #
#     ECOS Tests      #
#                     #
#######################

test_ecos_setup()

test_that("test ECOS options", {
  # Test that all ECOS solver options work.
  # Test ecos
  # feastol, abstol, reltol, feastol_inacc,
  # abstol_inacc, and reltol_inacc for tolerance values
  # max_iters for the maximum number of iterations,
  EPS <- 1e-4
  
  prob <- Problem(Minimize(norm1(x) + 1.0), list(x == 0))
  for(i in 1:2) {
    result <- solve(prob, solver = "ECOS", feastol = EPS, abstol = EPS, reltol = EPS,
                    feastol_inacc = EPS, abstol_inacc = EPS, reltol_inacc = EPS,
                    max_iters = 20, verbose = TRUE, warm_start = TRUE)
  }
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
})

test_that("test ECOS lp 0", {
  StandardTestLPs.test_lp_0(solver = "ECOS")
})

test_that("test ECOS lp 1", {
  StandardTestLPs.test_lp_1(solver = "ECOS")
})

test_that("test ECOS lp 2", {
  StandardTestLPs.test_lp_2(solver = "ECOS")
})

test_that("test ECOS lp 3", {
  StandardTestLPs.test_lp_3(solver = "ECOS")
})

test_that("test ECOS lp 4", {
  StandardTestLPs.test_lp_4(solver = "ECOS")
})

test_that("test ECOS lp 5", {
  StandardTestLPs.test_lp_5(solver = "ECOS")
})

test_that("test ECOS socp 0", {
  StandardTestSOCPs.test_socp_0(solver = "ECOS")
})

test_that("test ECOS socp 1", {
  StandardTestSOCPs.test_socp_1(solver = "ECOS")
})

test_that("test ECOS socp 2", {
  StandardTestSOCPs.test_socp_2(solver = "ECOS")
})

test_that("test ECOS socp 3", {
  # Axis 1
  StandardTestSOCPs.test_socp_3ax1(solver = "ECOS")
  # Axis 2
  StandardTestSOCPs.test_socp3_ax2(solver = "ECOS")
})

test_that("test ECOS expcone 1", {
  StandardTestECPs.test_expcone_1(solver = "ECOS")
})

test_that("test ECOS exp soc 1", {
  StandardTestMixedCPs.test_exp_soc_1(solver = "ECOS")
})

#######################
#                     #
#      SCS Tests      #
#                     #
#######################

test_scs_setup()

test_that("test that SCS retry doesn't trigger a crash", {
  n_sec <- 20
  set.seed(315)
  mu <- runif(n_sec)
  random_mat <- matrix(runif(n_sec^2), nrow = n_sec, ncol = n_sec)
  C <- random_mat %*% t(random_mat)
  
  xv <- Variable(n_sec)
  prob <- Problem(Minimize(QuadForm(xv, C)), list(sum(xv) == 1, 0 <= xv, xv <= 1, t(xv) %*% mu >= max(mu) - 1e-6))
  result <- solve(prob, solver = "SCS")
  expect_true(result$status %in% c(OPTIMAL, OPTIMAL_INACCURATE))
})

test_that("test that all SCS solver options work", {
  # Test SCS
  # MAX_ITERS, EPS, ALPHA, UNDET_TOL, VERBOSE, and NORMALIZE.
  # If opts is missing, then the algorithm uses default settings.
  # USE_INDIRECT = True
  EPS <- 1e-4
  
  xv <- Variable(2, name = "x")
  prob <- Problem(Minimize(norm1(x) + 1.0), list(x == 0))
  for(i in 1:2) {
    result <- solve(prob, solver = "SCS", max_iters = 50, eps = EPS, alpha = 1.2, verbose = TRUE, normalize = TRUE, use_indirect = FALSE)
  }
  expect_equal(result$value, 1.0, tolerance = 1e-2)
  expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = 1e-2)
})

test_that("test log problem", {
  # Log in objective.
  obj <- Maximimize(sum(log(x)))
  constr <- list(x <= c(1, exp(1)))
  p <- Problem(obj, constr)
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 1, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, exp(1))), tolerance = TOL)
  
  # Log in constraint.
  obj <- Minimize(sum(x))
  constr <- list(log(x) >= 0, x <= c(1,1))
  p <- Problem(obj, constr)
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
  
  # Index into log.
  obj <- Maximize(log(x)[2])
  constr <- list(x <= c(1, exp(1)))
  p <- Problem(obj, constr)
  result <- solve(p, solver = "SCS")
})

test_that("test sigma_max", {
  const <- Constant(rbind(c(1,2,3), c(4,5,6)))
  constr <- list(C == const)
  prob <- Problem(Minimize(norm(C, "2")), constr)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, value(norm(const, "2")), tolerance = TOL)
  expect_equal(result$getValue(C), value(const), tolerance = TOL)
})

test_that("test sdp var", {
  const <- Constant(rbind(c(1,2,3), c(4,5,6), c(7,8,9)))
  X <- Variable(3, 3, PSD = TRUE)
  prob <- Problem(Minimize(0), list(X == const))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, INFEASIBLE)
})

test_that("test complex matrices", {
  # Complex-valued matrix.
  set.seed(0)
  K <- matrix(runif(2*2), nrow = 2, ncol = 2) + 1i*matrix(runif(2*2), nrow = 2, ncol = 2)   # Example matrix
  n1 <- svd(K, only.values = TRUE)$values   # Trace norm of K.
  
  # Dual Problem.
  X <- Variable(2, 2, complex = TRUE)
  Y <- Variable(2, 2, complex = TRUE)
  
  # X, Y >= 0 so trace is real.
  objective <- Minimize(Re(0.5*matrix_trace(X) + 0.5*matrix_trace(Y)))
  constraints <- list(bmat(list(list(X, Conj(-K)), list(-K, Y))) %>>% 0, X %>>% 0, Y %>>% 0)
  problem <- Problem(objective, constraints)
  
  result <- solve(problem, solver = "SCS")
  sol_scs <- result$value
  expect_equal(dim(result$getDualValue(constraints[[1]])), c(4,4))
  expect_equal(dim(result$getDualValue(constraints[[2]])), c(2,2))
  expect_equal(dim(result$getDualValue(constraints[[3]])), c(2,2))
  expect_equal(sol_scs, n1, tolerance = TOL)
})

test_that("test a problem with entr", {
  for(n in c(5, 10, 25)) {
    print(n)
    x <- Variable(n)
    obj <- Minimize(sum(exp(x)))
    p <- Problem(obj, list(sum(x) == 1))
    result <- solve(p, solver = "SCS")
    expect_equal(result$getValue(x), matrix(rep(1/n, n)), tolerance = TOL)
  }
})

test_that("test a problem with log", {
  for(n in c(5, 10, 25)) {
    print(n)
    x <- Variable(n)
    obj <- Maximize(sum(log(x)))
    p <- Problem(obj, list(sum(x) == 1))
    result <- solve(p, solver = "SCS")
    expect_equal(result$getValue(x), matrix(rep(1/n, n)), tolerance = TOL)
  }
})

test_that("test a problem with log twice", {
  n <- 5
  x <- Variable(n)
  obj <- Maximize(sum(log(x)))
  p <- Problem(obj, list(sum(x) == 1))
  result1 <- solve(p, solver = "SCS")
  first_value <- result1$getValue(x)
  expect_equal(first_value, matrix(rep(1/n, n )), tolerance = tOL)
  
  result2 <- solve(p, solver = "SCS")
  second_value <- result2$getValue(x)
  expect_equal(first_value, second_value, tolerance = TOL)
})

test_that("test warm start", {
  x <- Variable(10)
  obj <- Minimize(sum(exp(x)))
  prob <- Problem(obj, list(sum(x) == 1))
  result <- solve(prob, solver = "SCS")
  time <- result$solver_stats$solve_time
  result2 <- solve(prob, solver = "SCS", warm_start = TRUE)
  time2 <- result$solver_stats$solver_time
  expect_equal(result2$value, result$value, tolerance = 1e-2)
  print(time > time2)
})

test_that("test warm start in diffcp", {
  if(!require(diffcp, quietly = TRUE)) {
    print("Skipping test, diffcp not installed")
    return()
  }
  
  x <- Variable(10)
  obj <- Minimize(sum(exp(x)))
  prob <- Problem(obj, list(sum(x) == 1))
  result <- solve(prob, solver = "DIFFCP")
  result2 <- solve(prob, solver = "DIFFCP", warm_start = TRUE)
  expect_equal(result2$value, result$value, tolerance = 1e-2)
})

test_that("test PSD constraint", {
  s <- Variable(2, 2)
  obj <- Maximize(min_elemwise(s[1,2], 10))
  const <- list(s %>>% 0, diag(s) == c(1,1))
  prob <- Problem(obj, const)
  result <- solve(prob, solver = "SCS")
  r <- result$value
  s <- result$getValue(s)
  
  print(residual(const[[1]]))
  print(paste("value", r))
  print(paste("s", s))
  
  eigs <- eigen(s + t(s), only.values = TRUE)$values
  print(paste("eigs", eigs))
  expect_true(all(eigs >= 0))
})

test_that("test SCS canonicalization with a quadratic objective", {
  # Note: Only relevant for SCS >= 3.0.0.
  require(scs)
  # TODO: Get scs version and only run if version >= 3.0.0.
  x <- Variable(2)
  expr <- sum_squares(x)
  constr <- list(x >= 1)
  prob <- Problem(Minimize(expr), constr)
  data <- get_problem_data(prob, solver = "SCS")
  expect_equal(as.matrix(data[[1]]$P), 2*diag(2), tolerance = TOL)
  solution1 <- solve(prob, solver = "SCS")
  
  # When use_quad_obj = FALSE, the quadratic objective is canonicalized to a SOC constraint.
  prob - Problem(Minimize(expr), constr)
  solver_opts <- list(use_quad_obj = FALSE)
  data <- get_problem_data(prob, solver = "SCS", solver_opts = solver_opts)
  expect_false("P" %in% names(data[[1]]))
  solution2 <- do.call("solve", c(list(a = prob, solver = "SCS"), solver_opts))
  
  expect_true(is.allclose(solution1$value, solution2$value))
  
  # Check that there is no P for non-quadratic objectives.
  expr <- norm1(x)
  prob <- Problem(Minimize(expr), constr)
  data <- get_problem_data(prob, solver = "SCS")
  expect_false("P" %in% names(data[[1]]))
})

test_that("test SCS lp 3", {
  StandardTestLPs.test_lp_3(solver = "SCS")
})

test_that("test SCS lp 4", {
  StandardTestLPs.test_lp_4(solver = "SCS")
})

test_that("test SCS lp 5", {
  StandardTestLPs.test_lp_5(solver = "SCS", eps = 1e-6)
})

test_that("test SCS socp 1", {
  # Axis 1
  StandardTestSOCPs.test_socp_3ax1(solver = "SCS")
  # Axis 2
  StandardTestSOCPs.test_socp_3ax2(solver = "SCS")
})

test_that("test SCS sdp 1min", {
  StandardTestSDPs.test_sdp_1min(solver = "SCS")
})

test_that("test SCS sdp 2", {
  StandardTestSDPs.test_sdp_2(solver = "SCS", eps = 1e-5)
})

test_that("test SCS expcone 1", {
  StandardTestECPs.test_expcone_1(solver = "SCS", eps = 1e-5)
})

test_that("test SCS exp soc 1", {
  StandardTestMixedCPs.test_exp_soc_1(solver = "SCS", eps = 1e-5)
})

test_that("test SCS pcp 1", {
  StandardTestPCPs.test_pcp_1(solver = "SCS")
})

test_that("test SCS pcp 2", {
  StandardTestPCPs.test_pcp_2(solver = "SCS")
})

test_that("test SCS pcp 3", {
  StandardTestPCPs.test_pcp_3(solver = "SCS", eps = 1e-12)
})

#######################
#                     #
#   CLARABEL Tests    #
#                     #
#######################

if("CLARABEL" %in% INSTALLED_SOLVERS) {
  test_clarabel_setup()
  
  test_that("test clarabel lp 0", {
    StandardTestLPs.test_lp_0(solver = "CLARABEL")
  })
  
  test_that("test clarabel lp 1", {
    StandardTestLPs.test_lp_1(solver = "CLARABEL")
  })
  
  test_that("test clarabel lp 2", {
    StandardTestLPs.test_lp_2(solver = "CLARABEL")
  })
  
  test_that("test clarabel lp 3", {
    StandardTestLPs.test_lp_3(solver = "CLARABEL")
  })
  
  test_that("test clarabel lp 4", {
    StandardTestLPs.test_lp_4(solver = "CLARABEL")
  })
  
  test_that("test clarabel lp 5", {
    StandardTestLPs.test_lp_5(solver = "CLARABEL")
  })
  
  test_that("test clarabel qp 0", {
    StandardTestQPs.test_qp_0(solver = "CLARABEL")
  })
  
  test_that("test clarabel qp 0 linear obj", {
    StandardTestQPs.test_qp_0(solver = "CLARABEL", use_quad_obj = FALSE)
  })
  
  test_that("test clarabel socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "CLARABEL")
  })
  
  test_that("test clarabel socp 1", {
    StandardTestSOCPs.test_socp_1(solver = "CLARABEL")
  })
  
  test_that("test clarabel socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "CLARABEL")
  })
  
  test_that("test clarabel socp 3", {
    # Axis 1
    StandardTestSOCPs.test_socp_3ax1(solver = "CLARABEL")
    # Axis 2
    StandardTestSOCPs.test_socp_3ax2(solver = "CLARABEL")
  })
  
  test_that("test clarabel expcone 1", {
    StandardTestECPs.test_expcone_1(solver = "CLARABEL")
  })
  
  test_that("test clarabel exp soc 1", {
    StandardTestMixedCPs.test_exp_soc_1(solver = "CLARABEL")
  })
  
  test_that("test clarabel pcp 0", {
    StandardTestSOCPs.test_socp_0(solver = "CLARABEL")
  })
  
  test_that("test clarabel pcp 1", {
    StandardTestSOCPs.test_socp_1(solver = "CLARABEL")
  })
  
  test_that("test clarabel pcp 2", {
    StandardTestSOCPs.test_socp_2(solver = "CLARABEL")
  })
} else
  warning("CLARABEL is not installed. Skipping tests.")

#######################
#                     #
#     MOSEK Tests     #
#                     #
#######################

if("MOSEK" %in% INSTALLED_SOLVERS) {
  test_that("test mosek lp 0", {
    StandardTestLPs.test_lp_0(solver = "MOSEK")
  })
  
  test_that("test mosek lp 1", {
    # default settings
    StandardTestLPs.test_lp_1(solver = "MOSEK")
    # require a basic feasible solution
    StandardTestLPs.test_lp_1(solver = "MOSEK", tolerance = 1e-6, bfs = TRUE)
  })
  
  test_that("test mosek lp 2", {
    StandardTestLPs.test_lp_2(solver = "MOSEK")
  })
  
  test_that("test mosek lp 3", {
    StandardTestLPs.test_lp_3(solver = "MOSEK")
  })
  
  test_that("test mosek lp 4", {
    StandardTestLPs.test_lp_4(solver = "MOSEK")
  })
  
  test_that("test mosek lp 5", {
    StandardTestLPs.test_lp_5(solver = "MOSEK")
  })
  
  test_that("test mosek socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "MOSEK")
  })
  
  test_that("test mosek socp 1", {
    StandardTestSOCPs.test_socp_1(solver = "MOSEK")
  })
  
  test_that("test mosek socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "MOSEK")
  })
  
  test_that("test mosek socp 3", {
    # Axis 1
    StandardTestSOCPs.test_socp_3ax1(solver = "MOSEK")
    # Axis 2
    StandardTestSOCPs.test_socp_3ax2(solver = "MOSEK")
  })
  
  test_that("test mosek sdp 1", {
    # Minimization
    StandardTestSDPs.test_sdp_1min(solver = "MOSEK")
    # Maximization
    StandardTestSDPs.test_sdp_1max(solver = "MOSEK")
  })
  
  test_that("test mosek sdp 2", {
    StandardTestSDPs.test_sdp_2(solver = "MOSEK")
  })
  
  test_that("test mosek expcone 1", {
    StandardTestECPs.test_expcone_1(solver = "MOSEK")
  })
  
  test_that("test mosek exp soc 1", {
    StandardTestMixedCPs.test_exp_soc_1(solver = "MOSEK")
  })
  
  test_that("test mosek pcp 1", {
    StandardTestPCPs.test_pcp_1(solver = "MOSEK", tolerance = 1e-2)
  })
  
  test_that("test mosek pcp 2", {
    StandardTestPCPs.test_pcp_2(solver = "MOSEK")
  })
  
  test_that("test mosek pcp 3", {
    StandardTestPCPs.test_pcp_3(solver = "MOSEK")
  })
  
  test_that("test mosek mi lp 0", {
    StandardTestLPs.test_mi_lp_0(solver = "MOSEK")
  })
  
  test_that("test mosek mi lp 1", {
    StandardTestLPs.test_mi_lp_1(solver = "MOSEK")
  })
  
  test_that("test mosek mi lp 2", {
    StandardTestLPs.test_mi_lp_2(solver = "MOSEK")
  })
  
  test_that("test mosek mi lp 3", {
    StandardTestLPs.test_mi_lp_3(solver = "MOSEK")
  })
  
  test_that("test mosek mi lp 5", {
    StandardTestLPs.test_mi_lp_5(solver = "MOSEK")
  })
  
  test_that("test mosek mi socp 1", {
    StandardTestSOCPs.test_mi_socp_1(solver = "MOSEK", tolerance = 1e-3)
  })
  
  test_that("test mosek mi socp 2", {
    StandardTestSOCPs.test_mi_socp_2(solver = "MOSEK")
  })
  
  test_that("test mosek mi pcp 0", {
    StandardTestPCPs.test_mi_pcp_0(solver = "MOSEK")
  })
  
  test_that("test mosek params", {
    if("MOSEK" %in% installed_solvers()) {
      require(Rmosek)
      n <- 10
      m <- 4
      set.seed(0)
      A <- matrix(rnorm(m*n), nrow = m, ncol = n)
      x <- matrix(rnorm(n))
      y <- A %*% x
      
      # Solve a simple basis pursuit problem for testing purposes.
      z <- Variable(n)
      objective <- Minimize(norm1(z))
      constraints <- list(A %*% z == y)
      problem <- Problem(objective, constraints)
      
      invalid_mosek_params <- list(dparam.basis_tol_x = 1e-8)
      
      expect_error(do.call("solve", list(a = problem, solver = "MOSEK", mosek_params = invalid_mosek_params)))
      expect_error(do.call("solve", list(solver = "MOSEK", invalid_kwarg = NULL)))
      mosek_params <- list("dparam.basis_tol_x" = 1e-8, "MSK_IPAR_INTPNT_MAX_ITERATIONS" = 20)
      result <- solve(problem, solver = "MOSEK", mosek_params = mosek_params)
    }
  })
} else
  warning("MOSEK is not installed. Skipping tests.")

#######################
#                     #
#    CVXOPT Tests     #
#                     #
#######################

if("CVXOPT" %in% INSTALLED_SOLVERS) {
  test_that("test cvxopt options", {
    # Test that all the CVXOPT solver options works.
    # 'maxiters'
    # maximum number of iterations (default: 100).
    # 'abstol'
    # absolute accuracy (default: 1e-7).
    # 'reltol'
    # relative accuracy (default: 1e-6).
    # 'feastol'
    # tolerance for feasibility conditions (default: 1e-7).
    # 'refinement'
    # number of iterative refinement steps when solving KKT equations
    # (default: 0 if the problem has no second-order cone
    #  or matrix inequality constraints; 1 otherwise).
    EPS <- 1e-7
    prob <- Problem(Minimize(norm1(x) + 1.0), list(x == 0))
    result <- solve(prob, solver = "CVXOPT", feastol = EPS, abstol = EPS, reltol = EPS,
                    max_iters = 20, verbose = TRUE, kktsolver = "chol", refinement = 2)
    expect_equal(result$value, 1.0, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
    
    setup_dummy_factor <- function(c, G, h, dims, A, b) {
      # See CVXR reduction_solvers.R section on KKTSolver for an action implementation of a setup-factor function.
      stop("Unimplemented: This setup-factor function was called")
    }
    
    expect_error(solve(prob, solver = "CVXOPT", kktsolver = setup_dummy_factor),
                 "Unimplemented: This setup-factor function was called", fixed = TRUE)
  })
  
  test_that("test cvxopt lp 0", {
    StandardTestLPs.test_lp_0(solver = "CVXOPT")
  })
  
  test_that("test cvxopt lp 1", {
    # default settings
    StandardTestLPs.test_lp_1(solver = "CVXOPT")
    # require a basic feasible solution
    StandardTestLPs.test_lp_1(solver = "CVXOPT", tolerance = 1e-6, bfs = TRUE)
  })
  
  test_that("test cvxopt lp 2", {
    StandardTestLPs.test_lp_2(solver = "CVXOPT")
  })
  
  test_that("test cvxopt lp 3", {
    StandardTestLPs.test_lp_3(solver = "CVXOPT")
  })
  
  test_that("test cvxopt lp 4", {
    StandardTestLPs.test_lp_4(solver = "CVXOPT")
  })
  
  test_that("test cvxopt lp 5", {
    StandardTestLPs.test_lp_5(solver = "CVXOPT", kktsolver = setup_ldl_factor)
    StandardTestLPs.test_lp_5(solver = "CVXOPT", kktsolver = "chol")
  })
  
  test_that("test cvxopt socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "CVXOPT")
  })
  
  test_that("test cvxopt socp 1", {
    StandardTestSOCPs.test_socp_1(solver = "CVXOPT")
  })
  
  test_that("test cvxopt socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "CVXOPT")
  })
  
  test_that("test cvxopt socp 3", {
    # Axis 1
    StandardTestSOCPs.test_socp_3ax1(solver = "CVXOPT")
    # Axis 2
    StandardTestSOCPs.test_socp_3ax2(solver = "CVXOPT")
  })
  
  test_that("test cvxopt sdp 1", {
    # Minimization
    StandardTestSDPs.test_sdp_1min(solver = "CVXOPT")
    # Maximization
    StandardTestSDPs.test_sdp_1max(solver = "CVXOPT")
  })
  
  test_that("test cvxopt sdp 2", {
    StandardTestSDPs.test_sdp_2(solver = "CVXOPT")
  })
} else
  warning("CVXOPT is not installed. Skipping tests.")

#######################
#                     #
#     SDPA Tests      #
#                     #
#######################

if("SDPA" %in% INSTALLED_SOLVERS) {
  test_that("test sdpa lp 0", {
    StandardTestLPs.test_lp_0(solver = "SDPA")
  })
  
  test_that("test sdpa lp 1", {
    StandardTestLPs.test_lp_1(solver = "SDPA")
  })
  
  test_that("test sdpa lp 2", {
    StandardTestLPs.test_lp_2(solver = "SDPA")
  })
  
  test_that("test sdpa lp 3", {
    StandardTestLPs.test_lp_3(solver = "SDPA")
  })
  
  test_that("test sdpa lp 4", {
    StandardTestLPs.test_lp_4(solver = "SDPA")
  })
  
  test_that("test sdpa lp 5", {
    print("Skipping test. Known limitation of SDPA for degenerate LPs.")
    return()
    
    # This also tests the ability to pass solver options.
    StandardTestLPs.test_lp_5(solver = "SDPA", gammaStar = 0.86, epsilonDash = 8.0e-6, betaStar = 0.18, betaBar = 0.15)
  })
  
  test_that("test sdpa sdp 1", {
    # Minimization
    StandardTestSDPs.test_sdp_1min(solver = "SDPA")
    # Maximization
    StandardTestSDPs.test_sdp_1max(solver = "SDPA")
  })
  
  test_that("test sdpa sdp 2", {
    StandardTestSDPs.test_sdp_2(solver = "SDPA")
  })
} else
  warning("SDPA is not installed. Skipping tests.")

#######################
#                     #
#      CBC Tests      #
#                     #
#######################

test_cbc_setup()

if("CBC" %in% INSTALLED_SOLVERS) {
  test_that("test that all CBC solver options work", {
    prob <- Problem(Minimize(norm2(x)), list(x == Variable(2, boolean = TRUE)))
    if("CBC" %in% INSTALLED_SOLVERS) {
      for(i in 1:2) {
        # Some cut-generators seem to be buggy for now -> set to false
        # prob.solve(solver=cvx.CBC, verbose=True, GomoryCuts=True, MIRCuts=True,
        #            MIRCuts2=True, TwoMIRCuts=True, ResidualCapacityCuts=True,
        #            KnapsackCuts=True, FlowCoverCuts=True, CliqueCuts=True,
        #            LiftProjectCuts=True, AllDifferentCuts=False, OddHoleCuts=True,
        #            RedSplitCuts=False, LandPCuts=False, PreProcessCuts=False,
        #            ProbingCuts=True, SimpleRoundingCuts=True)
        result <- solve(prob, solver = "CBC", verbose = TRUE, maximumSeconds = 100)
      }
      expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
    } else {
      expect_error(solve(prob, solver = "CBC"), "The solver CBC is not installed", fixed = TRUE)
    }
  })
  
  test_that("test cbc lp 0", {
    StandardTestLPs.test_lp_0(solver = "CBC", duals = FALSE)
  })
  
  test_that("test cbc lp 1", {
    StandardTestLPs.test_lp_1(solver = "CBC", duals = FALSE)
  })
  
  test_that("test cbc lp 2", {
    StandardTestLPs.test_lp_2(solver = "CBC", duals = FALSE)
  })
  
  test_that("test cbc lp 3", {
    StandardTestLPs.test_lp_3(solver = "CBC")
  })
  
  test_that("test cbc lp 4", {
    StandardTestLPs.test_lp_4(solver = "CBC")
  })
  
  test_that("test cbc lp 5", {
    StandardTestLPs.test_lp_5(solver = "CBC", duals = FALSE)
  })
  
  test_that("test cbc mi lp 0", {
    StandardTestLPs.test_mi_lp_0(solver = "CBC")
  })
  
  test_that("test cbc mi lp 1", {
    StandardTestLPs.test_mi_lp_1(solver = "CBC")
  })
  
  test_that("test cbc mi lp 2", {
    StandardTestLPs.test_mi_lp_2(solver = "CBC")
  })
  
  test_that("test cbc mi lp 3", {
    StandardTestLPs.test_mi_lp_3(solver = "CBC")
  })
  
  test_that("test cbc mi lp 5", {
    # TODO: CyLP <= 0.92.4 has no working integer infeasibility detection, so
    # if detected, this test should be skipped. See CVXPY test_conic_solvers.py under CBC tests.
    # https://github.com/coin-or/CyLP/pull/150
    StandardTestLPs.test_mi_lp_5(solver = "CBC")
  })
} else
  warning("CBC is not installed. Skipping tests.")

#######################
#                     #
#     GLPK Tests      #
#                     #
#######################

if("GLPK" %in% INSTALLED_SOLVERS) {
  test_that("test glpk lp 0", {
    StandardTestLPs.test_lp_0(solver = "GLPK")
  })
  
  test_that("test glpk lp 1", {
    StandardTestLPs.test_lp_1(solver = "GLPK")
  })
  
  test_that("test glpk lp 2", {
    StandardTestLPs.test_lp_2(solver = "GLPK")
  })
  
  test_that("test glpk lp 3", {
    StandardTestLPs.test_lp_3(solver = "GLPK")
  })
  
  test_that("test glpk lp 4", {
    StandardTestLPs.test_lp_4(solver = "GLPK")
  })
  
  test_that("test glpk lp 5", {
    StandardTestLPs.test_lp_5(solver = "GLPK")
  })
  
  test_that("test glpk lp 6", {
    StandardTestLPs.test_lp_6(solver = "GLPK")
  })
  
  test_that("test glpk mi lp 0", {
    StandardTestLPs.test_mi_lp_0(solver = "GLPK_MI")
  })
  
  test_that("test glpk mi lp 1", {
    StandardTestLPs.test_mi_lp_1(solver = "GLPK_MI")
  })
  
  test_that("test glpk mi lp 2", {
    StandardTestLPs.test_mi_lp_2(solver = "GLPK_MI")
  })
  
  test_that("test glpk mi lp 3", {
    StandardTestLPs.test_mi_lp_3(solver = "GLPK_MI")
  })
  
  test_that("test glpk mi lp 4", {
    StandardTestLPs.test_mi_lp_4(solver = "GLPK_MI")
  })
  
  test_that("test glpk mi lp 5", {
    StandardTestLPs.test_mi_lp_5(solver = "GLPK_MI")
  })
  
  test_that("test glpk options", {
    sth <- lp_1()
    require(cvxopt)
    # TODO: assert "tm_lim" not in cvxopt.glpk.options
    result <- solve(sth, solver = "GLPK", tm_lim = 100)
    # TODO: assert "tm_lim" not in cvxopt.glpk.options
    verify_objective(sth, tolerance = 1e-4)
    check_primal_feasibility(sth, tolerance = 1e-4)
    check_complementarity(sth, tolerance = 1e-4)
    check_dual_domains(sth, tolerance = 1e-4)
  })
  
  test_that("test glpk mi options", {
    sth <- mi_lp_1()
    require(cvxopt)
    # TODO: assert "tm_lim" not in cvxopt.glpk.options
    result <- solve(sth, solver = "GLPK_MI", tm_lim = 100, verbose = TRUE)
    # TODO: assert "tm_lim" not in cvxopt.glpk.options
    verify_objective(sth, tolerance = 1e-4)
    verify_primal_values(sth, tolerance = 1e-4)
  })
} else
  warning("GLPK is not installed. Skipping tests.")

#######################
#                     #
#     GLOP Tests      #
#                     #
#######################

if("GLOP" %in% INSTALLED_SOLVERS) {
  test_that("test glop lp 0", {
    StandardTestLPs.test_lp_0(solver = "GLOP")
  })
  
  test_that("test glop lp 1", {
    StandardTestLPs.test_lp_1(solver = "GLOP")
  })
  
  test_that("test glop lp 2", {
    StandardTestLPs.test_lp_2(solver = "GLOP")
  })
  
  test_that("test glop lp 3 no preprocessing", {
    params <- GlopParameters()
    params$use_preprocessing <- FALSE
    StandardTestLPs.test_lp_3(solver = "GLOP", parameters_proto = params)
  })
  
  # With preprocessing enabled, Glop internally detects
  # INFEASIBLE_OR_UNBOUNDED. This status is translated to
  # MPSOLVER_INFEASIBLE. See
  # https://github.com/google/or-tools/blob/b37d9c786b69128f3505f15beca09e89bf078a89/ortools/linear_solver/glop_utils.cc#L25-L38.
  test_that("test glop lp 3", {
    print("Skipping test. Known limitation of the GLOP interface")
    StandardTestLPs.test_lp_3(solver = "GLOP")
  })
  
  test_that("test glop lp 4", {
    StandardTestLPs.test_lp_4(solver = "GLOP")
  })
  
  test_that("test glop lp 5", {
    StandardTestLPs.test_lp_5(solver = "GLOP")
  })
  
  test_that("test glop lp 6 no preprocessing", {
    params <- GlopParameters()
    params$use_preprocessing <- FALSE
    StandardTestLPs.test_lp_6(solver = "GLOP", parameters_proto = params)
  })
  
  # Same issue as with test_glop_lp_3.
  test_that("test glop lp 6", {
    print("Skipping test. Known limitation of the GLOP interface")
    StandardTestLPs.test_lp_6(solver = "GLOP")
  })
  
  test_that("test glop bad parameters", {
    x <- Variable(1)
    prob <- Problem(Maximize(x), list(x <= 1))
    expect_error(solve(prob, solver = "GLOP", parameters_proto = "not a proto"))
  })
  
  test_that("test glop time limit", {
    sth <- lp_1()
    # Checks that the option doesn't error. A better test would be to solve
    # a large instance and check that the time limit is hit.
    result <- solve(sth, solver = "GLOP", time_limit_sec = 1.0)
  })
} else
  warning("GLOP is not installed. Skipping tests.")

#######################
#                     #
#     PDLP Tests      #
#                     #
#######################

if("PDLP" %in% INSTALLED_SOLVERS) {
  test_that("test pdlp lp 0", {
    StandardTestLPs.test_lp_0(solver = "PDLP")
  })
  
  test_that("test pdlp lp 1", {
    StandardTestLPs.test_lp_1(solver = "PDLP")
  })
  
  test_that("test pdlp lp 2", {
    StandardTestLPs.test_lp_2(solver = "PDLP")
  })
  
  test_that("test pdlp lp 3", {
    sth <- lp_3()
    result <- solve(sth@prob, solver = "PDLP")
    expect_equal(result$status, INFEASIBLE_OR_UNBOUNDED)
  })
  
  # We get the precise status when presolve is disabled.
  test_that("test pdlp lp 3 no presolve", {
    params <- PrimalDualHybridGradientParams()
    params$presolve_options$use_glop <- FALSE
    StandardTestLPs.test_lp_3(solver = "PDLP", parameters_proto = params)
  })
  
  test_that("test pdlp lp 4", {
    sth <- lp_4()
    result <- solve(sth@prob, solver = "PDLP")
    expect_equal(result$status, INFEASIBLE_OR_UNBOUNDED)
  })
  
  test_that("test pdlp lp 4 no presolve", {
    params <- PrimalDualHybridGradientParams()
    params$presolve_options$use_glop <- FALSE
    StandardTestLPs.test_lp_4(solver = "PDLP", parameters_proto = params)
  })
  
  test_that("test pdlp lp 5", {
    StandardTestLPs.test_lp_5(solver = "PDLP")
  })
  
  test_that("test pdlp lp 6", {
    sth <- lp_6()
    result <- solve(sth@prob, solver = "PDLP")
    expect_equal(result$status, INFEASIBLE_OR_UNBOUNDED)
  })
  
  test_that("test pdlp bad parameters", {
    x <- Variable(1)
    prob <- Problem(Maximize(x), list(x <= 1))
    expect_error(solve(prob, solver = "PDLP", parameters_proto = "not a proto"))
  })
  
  test_that("test pdlp time limit", {
    sth <- lp_1()
    # Checks that the option doesn't error. A better test would be to solve a large
    # instance and check that the time limit is hit.
    result <- solve(sth, solver = "PDLP", time_limit_sec = 1.0)
  })
} else
  warning("PDLP is not installed. Skipping tests.")

#######################
#                     #
#    CPLEX Tests      #
#                     #
#######################

test_cplex_setup()

if("CPLEX" %in% INSTALLED_SOLVERS) {
  test_that("test cplex warm start", {
    # Make sure that warm starting CPLEX behaves as expected.
    # Note: This only checks output, not whether or not CPLEX is warm starting internally.
    if("CPLEX" %in% INSTALLED_SOLVERS) {
      Ap <- Parameter(2, 2)
      bp <- Parameter(2)
      hp <- Parameter(2)
      cp <- Parameter(2)
      
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(1,1))
      
      objective <- Maximize(cp[1]*x[1] + cp[2]*x[2])
      constraints <- list(x[1] <= hp[1], x[2] <= hp[2], Ap %*% x == bp)
      prob <- Problem(objective, constraints)
      result <- solve(prob, solver = "CPLEX", warm_start = TRUE)
      expect_equal(result$value, 3)
      expect_equal(result$getValue(x), matrix(c(1,2)), tolerance = TOL)
      
      # Change Ap and bp from original values.
      value(Ap) <- rbind(c(0,0), c(0,1))   # Changed.
      value(bp) <- matrix(c(0,1))          # Changed.
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(1,1))
      
      # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "CPLEX", warm_start = TRUE)
      expect_equal(result$value, 3)
      expect_equal(result$getValue(x), matrix(c(2,1)), tolerance = TOL)
      
      # Change hp from the original values.
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(1,1))   # Changed.
      value(cp) <- matrix(c(1,1))
      
      # Without setting update_ineq_constrs = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "CPLEX", warm_start = TRUE)
      expect_equal(result$value, 2)
      expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
      
      # Change cp from the original values.
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(2,1))   # Changed.
      
      # Without setting update_objective = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "CPLEX", warm_start = TRUE)
      expect_equal(result$value, 4)
      expect_equal(result$getValue(x), matrix(c(1,2)), tolerance = TOL)
    } else {
      prob <- Problem(Minimize(norm1(x)), list(x == 0))
      expect_error(solve(prob, solver = "CPLEX", warm_start = TRUE), 
                   "The solver CPLEX is not installed", fixed = TRUE)
    }
  })
  
  test_that("test cplex params", {
    if("CPLEX" %in% INSTALLED_SOLVERS) {
      n <- 10
      m <- 4
      set.seed(0)
      
      Ac <- matrix(rnorm(m*n), nrow = m, ncol = n)
      xc <- matrix(rnorm(n), nrow = n, ncol = 1)
      yc <- A %*% x
      
      # Solve a simple basis pursuit problem for testing purposes.
      zv <- Variable(n)
      objective <- Minimize(norm1(zv))
      constraints <- list(Ac %*% zv == yc)
      problem <- Problem(objective, constraints)
      
      invalid_cplex_params <- list(bogus = "foo")
      expect_error(solve(problem, solver = "CPLEX", cplex_params = invalid_cplex_params))
      expect_error(solve(problem, solver = "CPLEX", invalid_kwarg = NULL))
      
      cplex_params <- list(advance = 0,   # int param.
                           simplex.limits.iterations = 1000,   # long param.
                           timelimit = 1000.0,   # double param.
                           workdir = "mydir"   # string param.
      )
      result <- solve(problem, solver = "CPLEX", cplex_params = cplex_params)
    }
  })
  
  test_that("test cplex lp 0", {
    StandardTestLPs.test_lp_0(solver = "CPLEX")
  })
  
  test_that("test cplex lp 1", {
    StandardTestLPs.test_lp_1(solver = "CPLEX")
  })
  
  test_that("test cplex lp 2", {
    StandardTestLPs.test_lp_2(solver = "CPLEX")
  })
  
  test_that("test cplex lp 3", {
    # CPLEX initially produces an INFEASIBLE_OR_UNBOUNDED status.
    sth <- lp_3()
    expect_warning(result <- solve(prob, solver = "CPLEX"))
    expect_equal(result$status, INFEASIBLE_OR_UNBOUNDED)
    # Determine the precise status with reoptimize = TRUE.
    StandardTestLPs.test_lp_3(solver = "CPLEX", reoptimize = TRUE)
  })
  
  test_that("test cplex lp 4", {
    StandardTestLPs.test_lp_4(solver = "CPLEX")
  })
  
  test_that("test cplex lp 5", {
    StandardTestLPs.test_lp_5(solver = "CPLEX")
  })
  
  test_that("test cplex socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "CPLEX")
  })
  
  test_that("test cplex socp 1", {
    # Parameters are set due to a minor dual-variable related presolve bug in CPLEX,
    # which will be fixed in the next CPLEX release.
    StandardTestSOCPs.test_socp_1(solver = "CPLEX", tolerance = 1e-2, 
                                  cplex_params = list(preprocessing.presolve = 0, preprocessing.reduce = 2))
  })
  
  test_that("test cplex socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "CPLEX")
  })
  
  test_that("test cplex socp 3", {
    # Axis 1.
    StandardTestSOCPs.test_socp_3ax1(solver = "CPLEX")
    # Axis 2.
    StandardTestSOCPs.test_socp_3ax2(solver = "CPLEX")
  })
  
  test_that("test cplex mi lp 0", {
    StandardTestLPs.test_mi_lp_0(solver = "CPLEX")
  })
  
  test_that("test cplex mi lp 1", {
    StandardTestLPs.test_mi_lp_1(solver = "CPLEX")
  })
  
  test_that("test cplex mi lp 2", {
    StandardTestLPs.test_mi_lp_2(solver = "CPLEX")
  })
  
  test_that("test cplex mi lp 3", {
    StandardTestLPs.test_mi_lp_3(solver = "CPLEX")
  })
  
  test_that("test cplex mi lp 5", {
    StandardTestLPs.test_mi_lp_5(solver = "CPLEX")
  })
  
  test_that("test cplex mi socp 1", {
    StandardTestSOCPs.test_mi_socp_1(solver = "CPLEX", tolerance = 1e-3)
  })
  
  test_that("test cplex mi socp 2", {
    StandardTestSOCPs.test_mi_socp_2(solver = "CPLEX")
  })
} else
  warning("CPLEX is not installed. Skipping tests.")

#######################
#                     #
#    GUROBI Tests     #
#                     #
#######################

test_gurobi_setup()

if("GUROBI" %in% INSTALLED_SOLVERS) {
  test_that("test gurobi warm start", {
    # Make sure that warm starting GUROBI behaves as expected.
    # Note: This only checks output, not whether or not GUROBI is warm starting internally.
    if("GUROBI" %in% installed_solvers()) {
      library(gurobi)
      
      Ap <- Parameter(2, 2)
      bp <- Parameter(2)
      hp <- Parameter(2)
      cp <- Parameter(2)
      
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(1,1))
      
      objective <- Maximize(cp[1]*x[1] + cp[2]*x[2])
      constraints <- list(x[1]^2 <= hp[1]^2, x[2] <= hp[2], Ap %*% x == bp)
      prob <- Problem(objective, constraints)
      result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
      expect_equal(result$value, 3, tolerance = TOL)
      expect_equal(result$getValue(x), matrix(c(1,2)), tolerance = TOL)
      
      # Change Ap and bp from original values.
      value(Ap) <- rbind(c(0,0), c(0,1))   # Changed.
      value(bp) <- matrix(c(0,1))          # Changed.
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(1,1))
      
      # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
      expect_equal(result$value, 3, tolerance = TOL)
      expect_equal(result$getValue(x), matrix(c(2,1)), tolerance = TOL)
      
      # Change hp from the original values.
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(1,1))   # Changed.
      value(cp) <- matrix(c(1,1))
      
      # Without setting update_ineq_constrs = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
      expect_equal(result$value, 2, tolerance = TOL)
      expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
      
      # Change cp from the original values.
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(2,1))   # Changed.
      
      # Without setting update_objective = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
      expect_equal(result$value, 4)
      expect_equal(result$getValue(x), matrix(c(1,2)), tolerance = TOL)
      
      # Try creating a new problem and setting value(x).
      init_value <- matrix(c(2,3))
      value(x) <- init_value
      # objective <- Maximize(cp[1]*x[1] + cp[2]*x[2])
      # constraints <- list(x[1]^2 <= hp[1]^2, x[2] <= hp[2], Ap %*% x == bp)
      prob <- Problem(objective, constraints)
      result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
      expect_equal(result$value, 4)
      expect_equal(result$getValue(x), matrix(c(1,2)), tolerance = TOL)
      
      # Check that "start" value was set appropriately.
      model <- result$solver_stats$extra_stats
      model_x <- getVars(model)
      x_opt <- result$getValue(x)
      for(i in 1:size(x)) {
        expect_equal(init_value[i], model_x[i]$start)
        expect_true(is.allclose(x_opt[i], model_x[i]$x))
      }
      
      # Test with matrix variable.
      z <- Variable()
      Y <- Varabiel(3, 2)
      Y_val <- matrix(0:5, nrow = 3, ncol = 2, byrow = TRUE)
      value(Y) <- Y_val + 1
      objective <- Maximize(z + sum(Y))
      constraints <- list(Y <= Y_val, z <= 2)
      prob <- Problem(objective, constraints)
      result <- solve(prob, solver = "GUROBI", warm_start = tRUE)
      expect_equal(result$value, sum(Y_val) + 2)
      expect_equal(result$getValue(z), 2, tolerance = TOL)
      expect_equal(result$getValue(Y), Y_val, tolerance = TOL)
      
      # Check that "start" value was set appropriately.
      model <- result$solver_stats$extra_stats
      model_x <- getVars(model)
      expect_equal(model_x[1]$start, "UNDEFINED")
      expect_true(is.allclose(model_x[1]$x, 2))
      Y_opt <- result$getValue(Y)
      for(i in 1:size(Y)) {
        row <- (i - 1) %% nrow(Y) + 1
        col <- floor((i - 1) / nrow(Y)) + 1
        expect_equal(Y_val[row, col] + 1 == model_x[i]$start)
        expect_true(is.allclose(Y_opt[row, col], model_x[i]$x))
      }
    } else {
      prob <- Problem(Minimize(norm1(x, 1)), list(x == 0))
      expect_error(solve(prob, solver = "GUROBI", warm_start = TRUE),
                   "The solver GUROBI is not installed", fixed = TRUE)
    }
  })
  
  test_that("test gurobi time limit no solution", {
    # Make sure that if Gurobi terminates due to a time limit before finding a solution:
    #   1) no error is raised,
    #   2) solver stats are returned.
    # The test is skipped if something changes on Gurobi's side so that:
    #   - a solution is found despite a time limit of zero,
    #   - a different termination criteria is hit first.
    if(GUROBI %in% installed_solvers()) {
      require(gurobi)
      objective <- Minimize(x[1])
      constraints <- list(square(x[1]) <= 1)
      prob <- Problem(objective, constraints)
      # TODO.
    } else {
      prob <- Problem(Minimize(norm1(x, 1)), list(x == 0))
      expect_error(solve(prob, solver = "GUROBI", TimeLimit = 0),
                   "The solver GUROBI is not installed", fixed = TRUE)
    }
  })
  
  test_that("test gurobi environment", {
    # Tests that GUROBI environments can be passed to Model.
    # GUROBI environments can include licensing and model parameter data.
    if("GUROBI" %in% installed_solvers()) {
      require(gurobi)
      
      # Set a few parameters to random values close to their defaults.
      params <- list(MIPGap = runif(1),                 # range {0, INFINITY}
                     AggFill = sample(0:9, size = 1),   # range {-1, MAXINT}
                     PerturbValue = runif(1)            # range: {0, INFINITY}
      )
      
      # Create a custom environment and set some parameters.
      custom_env <- gurobi::Env()
      for(k in names(params)) {
        v <- params[[k]]
        setParam(custom_env, k, v)
      }
      
      # Testing Conic Solver Interface.
      tmp <- StandardTestSOCPs.test_socp_0(solver = "GUROBI", env = custom_env)
      sth <- tmp[[1]]
      result <- tmp[[2]]
      model <- result$solver_stats$extra_stats
      for(k in names(params)) {
        v <- params[[k]]
        tmp <- getParamInfo(model, k)
        p_val <- tmp[[3]]
        expect_equal(v, p_val)
      }
    }
  })
  
  test_that("test gurobi lp 0", {
    StandardTestLPs.test_lp_0(solver = "GUROBI")
  })
  
  test_that("test gurobi lp 1", {
    StandardTestLPs.test_lp_1(solver = "GUROBI")
  })
  
  test_that("test gurobi lp 2", {
    StandardTestLPs.test_lp_2(solver = "GUROBI")
  })
  
  test_that("test gurobi lp 3", {
    # GUROBI initially produces an INFEASIBLE_OR_UNBOUNDED status.
    sth <- lp_3()
    expect_warning(result <- solve(prob, solver = "GUROBI"))
    expect_equal(result$status, INFEASIBLE_OR_UNBOUNDED)
    # The user disables presolve and so makes reoptimization unnecessary.
    StandardTestLPs.test_lp_3(solver = "GUROBI", InfUnbdInfo = 1)
    # The user determines the precise status with reoptimize = TRUE
    StandardTestLPs.test_lp_3(solver = "GUROBI", reoptimize = TRUE)
  })
  
  test_that("test gurobi lp 4", {
    StandardTestLPs.test_lp_4(solver = "GUROBI", reoptimize = TRUE)
  })
  
  test_that("test gurobi lp 5", {
    StandardTestLPs.test_lp_5(solver = "GUROBI")
  })
  
  test_that("test gurobi socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "GUROBI")
  })
  
  test_that("test gurobi socp 1", {
    StandardTestSOCPs.test_socp_1(solver = "GUROBI")
  })
  
  test_that("test gurobi socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "GUROBI")
  })
  
  test_that("test gurobi socp 3", {
    # Axis 1.
    StandardTestSOCPs.test_socp_3ax1(solver = "GUROBI")
    # Axis 2.
    StandardTestSOCPs.test_socp_3ax2(solver = "GUROBI")
  })
  
  test_that("test gurobi mi lp 0", {
    StandardTestLPs.test_mi_lp_0(solver = "GUROBI")
  })
  
  test_that("test gurobi mi lp 1", {
    StandardTestLPs.test_mi_lp_1(solver = "GUROBI")
  })
  
  test_that("test gurobi mi lp 2", {
    StandardTestLPs.test_mi_lp_2(solver = "GUROBI")
  })
  
  test_that("test gurobi mi lp 3", {
    StandardTestLPs.test_mi_lp_3(solver = "GUROBI")
  })
  
  test_that("test gurobi mi lp 5", {
    StandardTestLPs.test_mi_lp_5(solver = "GUROBI")
  })
  
  test_that("test gurobi mi socp 1", {
    StandardTestSOCPs.test_mi_socp_1(solver = "GUROBI", tolerance = 1e-2)
  })
  
  test_that("test gurobi mi socp 2", {
    StandardTestSOCPs.test_mi_socp_2(solver = "GUROBI")
  })
} else
  warning("GUROBI is not installed. Skipping tests.")

#######################
#                     #
#    XPRESS Tests     #
#                     #
#######################

if("XPRESS" %in% INSTALLED_SOLVERS) {
  test_that("test xpress warm start", {
    # Make sure that warm starting Xpress behaves as expected.
    # Note: Xpress does not have warm start yuet, it will re-solve problems from scratch.
    if("XPRESS" %in% installed_solvers()) {
      Ap <- Parameter(2, 2)
      bp <- Parameter(2)
      hp <- Parameter(2)
      cp <- Parameter(2)
      
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(1,1))
      
      objective <- Maximize(cp[1]*x[1] + cp[2]*x[2])
      constraints <- list(x[1] <= hp[1], x[2] <= hp[2], Ap %*% x == bp)
      prob <- Problem(objective, constraints)
      result <- solve(prob, solver = "XPRESS", warm_start = TRUE)
      expect_equal(result$value, 3)
      expect_equal(result$getValue(x), matrix(c(1,2)), tolerance = TOL)
      
      # Change Ap and bp from original values.
      value(Ap) <- rbind(c(0,0), c(0,1))   # Changed.
      value(bp) <- matrix(c(0,1))          # Changed.
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(1,1))
      
      # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "XPRESS", warm_start = TRUE)
      expect_equal(result$value, 3)
      expect_equal(result$getValue(x), matrix(c(2,1)), tolerance = TOL)
      
      # Change hp from the original values.
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(1,1))   # Changed.
      value(cp) <- matrix(c(1,1))
      
      # Without setting update_ineq_constrs = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "XPRESS", warm_start = TRUE)
      expect_equal(result$value, 2)
      expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
      
      # Change cp from the original values.
      value(Ap) <- rbind(c(1,0), c(0,0))
      value(bp) <- matrix(c(1,0))
      value(hp) <- matrix(c(2,2))
      value(cp) <- matrix(c(2,1))   # Changed.
      
      # Without setting update_objective = FALSE, the results should change to the correct answer.
      result <- solve(prob, solver = "XPRESS", warm_start = TRUE)
      expect_equal(result$value, 4)
      expect_equal(result$getValue(x), matrix(c(1,2)), tolerance = TOL)
    } else {
      prob <- Problem(Minimize(norm1(x)), list(x == 0))
      expect_error(solve(prob, solver = "XPRESS", warm_start = TRUE),
                   "The solver XPRESS is not installed", fixed = TRUE)
    }
  })
  
  test_that("test xpress params", {
    if("XPRESS" %in% installed_solvers()) {
      n <- 10
      m <- 4
      Ac <- matrix(rnorm(m*n), nrow = m, ncol = n)
      xc <- matrix(rnorm(n))
      yc <- A %*% x
      
      # Solve a simple basis pursuit problem for testing purposes.
      zv <- Variable(n)
      objective <- Minimize(norm1(zv))
      constraints <- list(Ac %*% zv == yc)
      problem <- Problem(objective, constraints)
      
      params <- list(lpiterlimit = 1000,   # maximum number of simplex iterations
                     maxtime = 1000        # time limit
      )
      result <- do.call("solve", c(list(a = problem, solver = "XPRESS"), params))
    }
  })
  
  test_that("test xpress iis none", {
    if("XPRESS" %in% installed_solvers()) {
      Ac <- rbind(c(2,1), c(1,2), c(-3,-3))
      bc <- matrix(c(2,2,-5))
      
      xv <- Variable(2)
      objective <- Minimize(norm2(xv))
      constraints <- list(Ac %*% xv <= bc)
      problem <- Problem(objective, constraints)
      
      params <- list(save_iis = 0)
      result <- do.call("solve", c(list(a = problem, solver = "XPRESS"), params))
    }
  })
  
  test_that("test xpress iis full", {
    if("XPRESS" %in% installed_solvers()) {
      Ac <- rbind(c(2,1), c(1,2), c(-3,-3))
      bc <- matrix(c(2,2,-5))
      
      xv <- Variable(2)
      objective <- Minimize(norm2(xv))
      constraints <- list(Ac %*% xv <= bc)
      problem <- Problem(objective, constraints)
      
      params <- list(save_iis = -1)
      result <- do.call("solve", c(list(a = problem, solver = "XPRESS"), params))
      expect_true("XPRESS_IIS" %in% result$solution$attr)
    }
  })
  
  test_that("test xpress lp 0", {
    StandardTestLPs.test_lp_0(solver = "XPRESS")
  })
  
  test_that("test xpress lp 1", {
    StandardTestLPs.test_lp_1(solver = "XPRESS")
  })
  
  test_that("test xpress lp 2", {
    StandardTestLPs.test_lp_2(solver = "XPRESS")
  })
  
  test_that("test xpress lp 3", {
    StandardTestLPs.test_lp_3(solver = "XPRESS")
  })
  
  test_that("test xpress lp 4", {
    StandardTestLPs.test_lp_4(solver = "XPRESS")
  })
  
  test_that("test xpress socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "XPRESS")
  })
  
  test_that("test xpress socp 1", {
    StandardTestSOCPs.test_socp_1(solver = "XPRESS")
  })
  
  test_that("test xpress socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "XPRESS")
  })
  
  test_that("test xpress mi lp 0", {
    StandardTestLPs.test_mi_lp_0(solver = "XPRESS")
  })
  
  test_that("test xpress mi lp 1", {
    StandardTestLPs.test_mi_lp_1(solver = "XPRESS")
  })
  
  test_that("test xpress mi lp 2", {
    StandardTestLPs.test_mi_lp_2(solver = "XPRESS")
  })
  
  test_that("test xpress mi lp 3", {
    StandardTestLPs.test_mi_lp_3(solver = "XPRESS")
  })
  
  test_that("test xpress mi lp 5", {
    StandardTestLPs.test_mi_lp_5(solver = "XPRESS")
  })
  
  test_that("test xpress mi socp 1", {
    StandardTestSOCPs.test_mi_socp_1(solver = "XPRESS")
  })
  
  test_that("test xpress mi socp 2", {
    StandardTestSOCPs.test_mi_socp_2(solver = "XPRESS")
  })
} else
  warning("XPRESS it not installed. Skipping tests.")


#######################
#                     #
#     NAG Tests       #
#                     #
#######################

if("NAG" %in% INSTALLED_SOLVERS) {
  test_that("test nag lp 0", {
    StandardTestLPs.test_lp_0(solver = "NAG")
  })
  
  test_that("test nag lp 1", {
    StandardTestLPs.test_lp_1(solver = "NAG")
  })
  
  test_that("test nag lp 2", {
    StandardTestLPs.test_lp_2(solver = "NAG")
  })
  
  test_that("test nag lp 3", {
    StandardTestLPs.test_lp_3(solver = "NAG")
  })
  
  test_that("test nag lp 4", {
    StandardTestLPs.test_lp_4(solver = "NAG")
  })
  
  test_that("test nag lp 5", {
    StandardTestLPs.test_lp_5(solver = "NAG")
  })
  
  test_that("test nag socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "NAG")
  })
  
  test_that("test nag socp 1", {
    StandardTestSOCPs.test_socp_1(solver = "NAG")
  })
  
  test_that("test nag socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "NAG")
  })
  
  test_that("test nag socp 3", {
    # Axis 1.
    StandardTestSOCPs.test_socp_3ax1(solver = "NAG")
    # Axis 2.
    StandardTestSOCPs.test_socp_3ax2(solver = "NAG")
  })
} else
  warning("NAG is not installed. Skipping tests.")

#######################
#                     #
#     SCIP Tests      #
#                     #
#######################

if("SCIP" %in% INSTALLED_SOLVERS) {
  test_that("test scip lp 0", {
    StandardTestLPs.test_lp_0(solver = "SCIP")
  })
  
  test_that("test scip lp 1", {
    StandardTestLPs.test_lp_1(solver = "SCIP")
  })
  
  test_that("test scip lp 2", {
    StandardTestLPs.test_lp_2(solver = "SCIP")
  })
  
  test_that("test scip lp 3", {
    StandardTestLPs.test_lp_3(solver = "SCIP")
  })
  
  test_that("test scip lp 4", {
    StandardTestLPs.test_lp_4(solver = "SCIP")
  })
  
  test_that("test scip socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "SCIP")
  })
  
  test_that("test scip socp 1", {
    StandardTestSOCPs.test_socp_1(solver = "SCIP")
  })
  
  test_that("test scip socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "SCIP")
  })
  
  test_that("test scip socp 3", {
    # Axis 1.
    StandardTestSOCPs.test_socp_3ax1(solver = "SCIP")
    # Axis 2.
    StandardTestSOCPs.test_socp_3ax2(solver = "SCIP")
  })
  
  test_that("test scip mi lp 0", {
    StandardTestLPs.test_mi_lp_0(solver = "SCIP")
  })
  
  test_that("test scip mi lp 1", {
    StandardTestLPs.test_mi_lp_1(solver = "SCIP")
  })
  
  test_that("test scip mi lp 2", {
    StandardTestLPs.test_mi_lp_2(solver = "SCIP")
  })
  
  test_that("test scip mi lp 3", {
    StandardTestLPs.test_mi_lp_3(solver = "SCIP")
  })
  
  test_that("test scip mi lp 5", {
    StandardTestLPs.test_mi_lp_5(solver = "SCIP")
  })
  
  test_that("test scip mi socp 1", {
    StandardTestSOCPs.test_mi_socp_1(solver = "SCIP", tolerance = 1e-3)
  })
  
  test_that("test scip mi socp 2", {
    StandardTestSOCPs.test_mi_socp_2(solver = "SCIP")
  })
  
  get_simple_problem <- function() {
    # Example problem that can be used within additional tests.
    x <- Variable()
    y <- Variable()
    constraints <- list(x >= 0, y >= 1, x + y <= 4)
    obj <- Maximize(x)
    prob <- Problem(obj, constraints)
    return(prob)
  }
  
  test_that("test scip test params - no params set", {
    prob <- get_simple_problem()
    result <- solve(prob, solver = "SCIP")
    # Important that passes without raising an error also check obj.
    expect_equal(result$value, 3)
  })
  
  test_that("test scip test params - valid params", {
    prob <- get_simple_problem()
    result <- solve(prob, solver = "SCIP", gp = FALSE)
    # Important that passes without raising an error also check obj.
    expect_equal(result$value, 3)
  })
  
  test_that("test scip test params - valid scip params", {
    prob <- get_simple_problem()
    result <- solve(prob, solver = "SCIP", scip_params = list('lp/fastmip' = 1, 'limits/gap' = 0.1))
    # Important that passes without raising an error also check obj.
    expect_equal(result$value, 3)
  })
  
  test_that("test scip test params - invalid params", {
    prob <- get_simple_problem()
    # Since an invalid NON-scip param is passed, an error is expected t be raised when calling solve.
    expect_error(solve(prob, solver = "SCIP", a = "what?"),
                 "One or more solver params in ['a'] are not valid: 'Not a valid parameter name'")
  })
  
  test_that("test scip test params - invalid scip params", {
    prob <- get_simple_problem()
    # Since an invalid NON-scip param is passed, an error is expected t be raised when calling solve.
    expect_error(solve(prob, solver = "SCIP", scip_params = list(a = "what?")),
                 "One or more solver params in ['a'] are not valid: 'Not a valid parameter name'")
  })
} else
  warning("SCIP is not installed. Skipping tests.")

#######################
#                     #
#    All Solvers      #
#                     #
#######################

test_all_solvers_setup()

test_that("test installed solvers", {
  # Test the list of installed solvers.
  prob <- Problem(Minimize(norm1(x) + 1.0), list(x == 0))
  for(solver in names(SOLVER_MAP_CONIC)) {
    if(solver %in% INSTALLED_SOLVERS) {
      result <- solve(prob, solver = solver)
      expect_equal(result$value, 1.0, tolerance = TOL)
      expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
    } else {
      expect_error(solve(prob, solver = solver), 
                   paste("The solver", solver, "is not installed"), fixed = TRUE)
    }
  }
  
  for(solver in names(SOLVER_MAP_QP)) {
    if(solver %in% INSTALLED_SOLVERS) {
      result <- solve(prob, solver = solver)
      expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
    } else {
      expect_error(solve(prob, solver = solver), 
                   paste("The solver", solver, "is not installed"), fixed = TRUE)
    }
  }
})

test_that("test mixed integer behavior", {
  x <- Variable(2, name = "x", integer = TRUE)
  objective <- Minimize(sum(x))
  prob <- Problem(objective, list(x >= 0))
  if(identical(INSTALLED_MI_SOLVERS, list("ECOS_BB"))) {
    expect_error(solve(prob), "You need a mixed-integer solver for this model")
  } else {
    result <- solve(prob)
    expect_equal(result$getValue(x), matrix(c(0,0)))
  }
})

#######################
#                     #
#   ECOS BB Tests     #
#                     #
#######################

test_that("test ecos bb explicit only", {
  x <- Variable(1, name = "x", integer = TRUE)
  objective <- Minimize(sum(x))
  prob <- Problem(objective, list(x >= 0))
  if(!identical(INSTALLED_MI_SOLVERS, list("ECOS_BB"))) {
    result <- solve(prob)
    expect_false(result$solver_stats$solver_name == "ECOS_BB")
  } else
    expect_error(solve(prob), "You need a mixed-integer solver for this model")
})

test_that("test ecos bb lp 0", {
  StandardTestLPs.test_lp_0(solver = "ECOS_BB")
})

test_that("test ecos bb lp 1", {
  # Default settings.
  StandardTestLPs.test_lp_1(solver = "ECOS_BB")
  # Require a basic feasible solution.
  StandardTestLPs.test_lp_1(solver = "ECOS_BB")
})

test_that("test ecos bb lp 2", {
  StandardTestLPs.test_lp_2(solver = "ECOS_BB")
})

test_that("test ecos bb lp 3", {
  StandardTestLPs.test_lp_3(solver = "ECOS_BB")
})

test_that("test ecos bb lp 4", {
  StandardTestLPs.test_lp_4(solver = "ECOS_BB")
})

test_that("test ecos bb lp 5", {
  StandardTestLPs.test_lp_5(solver = "ECOS_BB")
})

test_that("test ecos bb socp 0", {
  StandardTestSOCPs.test_socp_0(solver = "ECOS_BB")
})

test_that("test ecos bb socp 1", {
  StandardTestSOCPs.test_socp_1(solver = "ECOS_BB")
})

test_that("test ecos bb socp 2", {
  StandardTestSOCPs.test_socp_2(solver = "ECOS_BB")
})

test_that("test ecos bb socp 3", {
  # Axis 1.
  StandardTestSOCPs.test_socp_3ax1(solver = "ECOS_BB")
  # Axis 2.
  StandardTestSOCPs.test_socp_3ax2(solver = "ECOS_BB")
})

test_that("test ecos bb expcone 1", {
  StandardTestECPs.test_expcone_1(solver = "ECOS_BB")
})

test_that("test ecos bb exp soc 1", {
  StandardTestMixedCPs.test_exp_soc_1(solver = "ECOS_BB")
})

test_that("test ecos bb mi lp 0", {
  StandardTestLPs.test_mi_lp_0(solver = "ECOS_BB")
})

test_that("test ecos bb mi lp 2", {
  print("Known bug in ECOS BB. Skipping test.")
  return()
  StandardTestLPs.test_mi_lp_2(solver = "ECOS_BB")
})

test_that("test ecos bb mi lp 3", {
  StandardTestLPs.test_mi_lp_3(solver = "ECOS_BB")
})

test_that("test ecos bb mi lp 5", {
  StandardTestLPs.test_mi_lp_5(solver = "ECOS_BB")
})

test_that("test ecos bb mi socp 1", {
  print("Known bug in ECOS BB. Skipping test.")
  return()
  StandardTestSOCPs.test_mi_socp_1(solver = "ECOS_BB")
})

#######################
#                     #
#     COPT Tests      #
#                     #
#######################

if("COPT" %in% INSTALLED_SOLVERS) {
  test_that("test copt lp 0", {
    StandardTestLPs.test_lp_0(solver = "COPT")
  })
  
  test_that("test copt lp 1", {
    StandardTestLPs.test_lp_1(solver = "COPT")
  })
  
  test_that("test copt lp 2", {
    StandardTestLPs.test_lp_2(solver = "COPT")
  })
  
  test_that("test copt lp 3", {
    StandardTestLPs.test_lp_3(solver = "COPT")
  })
  
  test_that("test copt lp 4", {
    StandardTestLPs.test_lp_4(solver = "COPT")
  })
  
  test_that("test copt lp 5", {
    StandardTestLPs.test_lp_5(solver = "COPT")
  })
  
  test_that("test copt socp 0", {
    StandardTestSOCPs.test_socp_0(solver = "COPT")
  })
  
  test_that("test copt socp 1", {
    StandardTestSOCPs.test_socp_1(solver = "COPT", tolerance = 1e-3)
  })
  
  test_that("test copt socp 2", {
    StandardTestSOCPs.test_socp_2(solver = "COPT")
  })
  
  test_that("test copt socp 3", {
    # Axis 1.
    StandardTestSOCPs.test_socp_3ax1(solver = "COPT")
    # Axis 2.
    StandardTestSOCPs.test_socp_3ax2(solver = "COPT")
  })
  
  test_that("test copt mi lp 0", {
    StandardTestLPs.test_mi_lp_0(solver = "COPT")
  })
  
  test_that("test copt mi lp 1", {
    StandardTestLPs.test_mi_lp_1(solver = "COPT")
  })
  
  test_that("test copt mi lp 2", {
    StandardTestLPs.test_mi_lp_2(solver = "COPT")
  })
  
  test_that("test copt mi lp 3", {
    StandardTestLPs.test_mi_lp_3(solver = "COPT")
  })
  
  test_that("test copt mi lp 5", {
    StandardTestLPs.test_mi_lp_5(solver = "COPT")
  })
  
  test_that("test copt mi socp 1", {
    # COPT does not support MISOCP.
    expect_error(StandardTestSOCPs.test_mi_socp_1(solver = "COPT"), "do not support")
  })
  
  test_that("test copt sdp 1min", {
    StandardTestSDPs.test_sdp_1min(solver = "COPT")
  })
  
  test_that("test copt sdp 1max", {
    StandardTestSDPs.test_sdp_1max(solver = "COPT")
  })
  
  test_that("test copt sdp 2", {
    StandardTestSDPs.test_sdp_2(solver = "COPT")
  })
  
  test_that("test copt params", {
    n <- 10
    m <- 4
    set.seed(0)
    A <- matrix(rnorm(m*n), nrow = m, ncol = n)
    x <- matrix(rnorm(n))
    y <- A %*% x
    
    # Solve a simple basis pursuit problem for testing purposes.
    z <- Variable(n)
    objective <- Minimize(norm1(z))
    constraints <- list(A %*% z == y)
    problem <- Problem(objective, constraints)
    
    expect_error(solve(problem, solver = "COPT", invalid_kwarg = NULL))
    
    # Valid arg.
    result <- solve(problem, solver = "COPT", feastol = 1e-9)
  })
} else
  warning("COPT is not installed. Skipping tests.")
