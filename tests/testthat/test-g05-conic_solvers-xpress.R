context("test-g05-conic_solvers-xpress")

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

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
