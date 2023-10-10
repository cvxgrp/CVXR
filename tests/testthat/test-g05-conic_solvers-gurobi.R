context("test-g05-conic_solvers-gurobi")
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
