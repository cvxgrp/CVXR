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
#   Clarabel Tests    #
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

