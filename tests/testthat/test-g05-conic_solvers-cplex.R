context("test-g05-conic_solvers-cplex")
TOL <- 1e-4

CPLEX_AVAILABLE <- "CPLEX" %in% installed_solvers()
if(CPLEX_AVAILABLE)
  require(Rcplex)

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test cplex warm start", {
  # Make sure that warm starting CPLEX behaves as expected.
  # Note: This only checks output, not whether or not CPLEX is warm starting internally.
  if("CPLEX" %in% installed_solvers()) {
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
  if("CPLEX" %in% installed_solvers()) {
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
  result <- solve(prob, solver = "CPLEX")
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
