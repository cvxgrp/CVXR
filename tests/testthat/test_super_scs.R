TOL <- 0.01

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2,2, name = "A")
B <- Variable(2,2, name = "B")
C <- Variable(3,2, name = "C")

test_that("test log problem", {
  if("SUPER_SCS" %in% installed_solvers()) {
    # Log in objective.
    obj <- Maximize(sum(log(x)))
    constr <- list(x <= c(1, exp(1)))
    p <- Problem(obj, constr)
    result <- solve(p, solver = "SUPER_SCS")
    expect_equal(result$value, 1, tolerance = TOL)
    expect_equal(result$getValue(x), list(1, exp(1)), tolerance = TOL)
    
    # Log in constraint.
    obj <- Minimize(sum(x))
    constr <- list(log(x) >= 0, x <= c(1,1))
    p <- Problem(obj, constr)
    result <- solve(p, solver = "SUPER_SCS")
    expect_equal(result$value, 2, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(1, nrow = 2, ncol = 1), tolerance = TOL)
  }
})

test_that("test sigma_max", {
  if("SUPER_SCS" %in% installed_solvers()) {
    const <- Constant(rbind(1:3, 4:6))
    constr <- list(C == const)
    prob <- Problem(Minimize(cvxr_norm(C, 2)), constr)
    result <- solve(prob, solver = "SUPER_SCS")
    expect_equal(result$value, value(cvxr_norm(const, 2)), tolerance = TOL)
    expect_equal(result$getValue(C), value(const), tolerance = TOL)
  }
})

test_that("test SDP variable", {
  if("SUPER_SCS" %in% installed_solvers()) {
    const <- Constant(rbind(1:3, 4:6, 7:9))
    X <- Variable(3,3, PSD = TRUE)
    prob <- Problem(Minimize(0), list(X == const))
    result <- solve(prob, verbose = TRUE, solver = "SUPER_SCS")
    expect_equal(result$status, INFEASIBLE)
  }
})

test_that("test complex matrices", {
  if("SUPER_SCS" %in% installed_solvers()) {
    # Complex-valued matrix.
    K <- matrix(rnorm(4), nrow = 2, ncol = 2) + 1i*matrix(rnorm(4), nrow = 2, ncol = 2)   # Example matrix.
    n1 <- sum(svd(K)$d)   # Trace norm of K.
    
    # Dual problem.
    X <- Variable(2,2, complex = TRUE)
    Y <- Variable(2,2, complex = TRUE)
    Z <- Variable(2,2)
    
    # X, Y >= 0 so trace is real.
    objective <- Minimize(Re(0.5*matrix_trace(X) + 0.5*matrix_trace(Y)))
    constraints <- list(bmat(list(list(X, -Conj(t(K))), list(-K, Y))) %>>% 0, X %>>% 0, Y %>>% 0)
    problem <- Problem(objective, constraints)
    sol_scs <- solve(problem, solver = "SUPER_SCS")
    expect_equal(dim(dual_value(constraints[[1]])), c(4,4))
    expect_equal(dim(dual_value(constraints[[2]])), c(2,2))
    expect_equal(dim(dual_value(constraints[[3]])), c(2,2))
    expect_equal(sol_scs, n1, tolerance = TOL)
  }
})

test_that("test a problem with entr", {
  if("SUPER_SCS" %in% installed_solvers()) {
    for(n in c(5,10,25)) {
      print(n)
      x <- Variable(n)
      obj <- Maximize(sum(entr(x)))
      p <- Problem(obj, list(sum(x) == 1))
      result <- solve(p, solver = "SUPER_SCS", verbose = TRUE)
      expect_equal(result$getValue(x), n*(1.0/n), tolerance = TOL)
    }
  }
})

test_that("test a problem with exp", {
  if("SUPER_SCS" %in% installed_solvers()) {
    for(n in c(5,10,25)) {
      print(n)
      x <- Variable(n)
      obj <- Minimize(sum(exp(x)))
      p <- Problem(obj, list(sum(x) == 1))
      result <- solve(p, solver = "SUPER_SCS", verbose = TRUE)
      expect_equal(result$getValue(x), n*(1.0/n), tolerance = TOL)
    }
  }
})

test_that("test a problem with log", {
  if("SUPER_SCS" %in% installed_solvers()) {
    for(n in c(5,10,25)) {
      print(n)
      x <- Variable(n)
      obj <- Maximize(sum(log(x)))
      p <- Problem(obj, list(sum(x) == 1))
      result <- solve(p, solver = "SUPER_SCS", verbose = TRUE)
      expect_equal(result$getValue(x), n*(1.0/n), tolerance = TOL)
    }
  }
})

test_that("test warm starting", {
  if("SUPER_SCS" %in% installed_solvers()) {
    x <- Variable(10)
    obj <- Minimize(sum(exp(x)))
    prob <- Problem(obj, list(sum(x) == 1))
    result <- solve(prob, solver = "SUPER_SCS", eps = 1e-4)
    time <- result$solver_stats$solve_time
    result2 <- solve(prob, solver = "SUPER_SCS", warm_start = TRUE, eps = 1e-4)
    time2 <- result2$solver_stats$solve_time
    expect_equal(result2$value, result$value, tolerance = TOL)
  }
})
