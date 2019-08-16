context("test-g04-cbc")
TOL <- 1e-6

a <- Variable(name='a')
b <- Variable(name='b')
c <- Variable(name='c')

x <- Variable(2, name='x')
y <- Variable(3, name='y')
z <- Variable(2, name='z')

A <- Variable(2, 2, name='A')
B <- Variable(2, 2, name='B')
C <- Variable(3, 2, name='C')

test_that("Test basic LPs", {
  if("CBC" %in% installed_solvers()) {
    prob <- Problem(Minimize(0), list(x == 2))
    result <- solve(prob, verbose = FALSE, solver = "CBC")
    expect_equal(result$value, 0, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(2,2)), tolerance = TOL)
    
    prob <- Problem(Minimize(-a), list(a <= 1))
    result <- solve(prob, verbose = FALSE, solver = "CBC")
    expect_equal(result$value, -1, tolerance = TOL)
    expect_equal(result$getValue(a), 1, tolerance = TOL)
  }
})

test_that("Test a basic LP", {
  if("CBC" %in% installed_solvers()) {
    prob <- Problem(Minimize(p_norm(x, 1)), list(x == 0))
    result <- solve(prob, verbose = FALSE, solver = "CBC")
    expect_equal(result$value, 0, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] +   x[2] <= 3,
                          x[1] + 2*x[2] <= 3,
                          x[1] >= 0,
                          x[2] >= 0)
    prob <- Problem(objective, constraints)
    result <- solve(prob, verbose = FALSE, solver = "CBC")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
  } else {
    prob <- Problem(Minimize(p_norm(x, 1)), list(x == 0))
    expect_error(solve(prob, verbose = FALSE, solver = "CBC"))
  }
})

test_that("Test a basic MILP", {
  if("CBC" %in% installed_solvers()) {
    bool_var <- Variable(boolean = TRUE)
    int_var <- Variable(integer = TRUE)
    prob <- Problem(Minimize(p_norm(x, 1)),
                    list(x == bool_var, bool_var == 0))
    result <- solve(prob, verbose = FALSE, solver = "CBC")
    expect_equal(result$value, 0, tolerance = TOL)
    expect_equal(result$getValue(bool_var), 0, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] +   x[2] <= int_var,
                          x[1] + 2*x[2] <= 3*bool_var,
                          x[1] >= 0,
                          x[2] >= 0,
                        int_var == 3*bool_var,
                        int_var == 3)
    prob <- Problem(objective, constraints)
    result <- solve(prob, verbose = FALSE, solver = "CBC")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(int_var), 3, tolerance = TOL)
    expect_equal(result$getValue(bool_var), 1, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
  } else {
    prob <- Problem(Minimize(p_norm(x, 1)), list(x == 0))
    expect_error(solve(prob, verbose = FALSE, solver = "CBC"))
  }
})
