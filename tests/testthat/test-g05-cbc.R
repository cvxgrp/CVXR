context("test-g05-cbc")
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

CBC_AVAILABLE <- "CBC" %in% installed_solvers()

test_that("Test basic LPs", {
    skip_on_cran()
    skip_if_not(CBC_AVAILABLE, "Skipping CBC test as it is not available.!")
    ## TODO: This _used to be_ a bug in the rcbc library. NOW FIXED!
    prob <- Problem(Minimize(0), list(x == 2))
    result <- solve(prob, verbose = FALSE, solver = "CBC")
    expect_equal(result$value, 0, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(2,2)), tolerance = TOL)

    prob <- Problem(Minimize(-a), list(a <= 1))
    result <- solve(prob, verbose = FALSE, solver = "CBC")
    expect_equal(result$value, -1, tolerance = TOL)
    expect_equal(result$getValue(a), 1, tolerance = TOL)
})

test_that("Test a basic LP", {
    skip_on_cran()
    skip_if_not(CBC_AVAILABLE, "Skipping CBC test as it is not available.!")
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
})

test_that("Test a basic MILP", {
  skip_on_cran()
  skip_if_not(CBC_AVAILABLE, "Skipping CBC test as it is not available.!")
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
})

test_that("Test a hard knapsack problem", {
  skip_on_cran()
  skip_if_not(CBC_AVAILABLE, "Skipping CBC test as it is not available.!")
    ## Instance "knapPI_1_50_1000_1" from "http://www.diku.dk/~pisinger/genhard.c"
  n <- 50
  c <- 995
  z <- 8373
  coeffs <- rbind(c(1,   94, 485, 0), c(2,  506, 326, 0), c(3,  416, 248, 0),
                  c(4,  992, 421, 0), c(5,  649, 322, 0), c(6,  237, 795, 0),
                  c(7,  457,  43, 1), c(8,  815, 845, 0), c(9,  446, 955, 0),
                  c(10, 422, 252, 0), c(11, 791,   9, 1), c(12, 359, 901, 0),
                  c(13, 667, 122, 1), c(14, 598,  94, 1), c(15,   7, 738, 0),
                  c(16, 544, 574, 0), c(17, 334, 715, 0), c(18, 766, 882, 0),
                  c(19, 994, 367, 0), c(20, 893, 984, 0), c(21, 633, 299, 0),
                  c(22, 131, 433, 0), c(23, 428, 682, 0), c(24, 700,  72, 1),
                  c(25, 617, 874, 0), c(26, 874, 138, 1), c(27, 720, 856, 0),
                  c(28, 419, 145, 0), c(29, 794, 995, 0), c(30, 196, 529, 0),
                  c(31, 997, 199, 1), c(32, 116, 277, 0), c(33, 908,  97, 1),
                  c(34, 539, 719, 0), c(35, 707, 242, 0), c(36, 569, 107, 0),
                  c(37, 537, 122, 0), c(38, 931,  70, 1), c(39, 726,  98, 1),
                  c(40, 487, 600, 0), c(41, 772, 645, 0), c(42, 513, 267, 0),
                  c(43,  81, 972, 0), c(44, 943, 895, 0), c(45,  58, 213, 0),
                  c(46, 303, 748, 0), c(47, 764, 487, 0), c(48, 536, 923, 0),
                  c(49, 724,  29, 1), c(50, 789, 674, 0))  # index, p / w / x

  X <- Variable(n, boolean = TRUE)
  prob <- Problem(Maximize(sum(multiply(coeffs[,2], X))),
                  list(sum(multiply(coeffs[,3], X)) <= c))
  result <- solve(prob, verbose = FALSE, solver = "CBC")
  expect_equal(result$value, z, tolerance = TOL)   # objective
})

test_that("Test that all CBC solver options work", {
  skip_on_cran()
  skip_if_not(CBC_AVAILABLE, "Skipping CBC test as it is not available.!")
  prob <- Problem(Minimize(p_norm(x, 1)),
                  list(x == Variable(2, boolean = TRUE)))
  for(i in 1:2) {
      # Some cut-generators seem to be buggy for now -> set to false
      # result <- solve(prob, verbose = TRUE, solver = "CBC", GomoryCuts = TRUE, MIRCuts = TRUE,
      #            MIRCuts2 = TRUE, TwoMIRCuts = TRUE, ResidualCapacityCuts = TRUE,
      #            KnapsackCuts = TRUE, FlowCoverCuts = TRUE, CliqueCuts = TRUE,
      #            LiftProjectCuts = TRUE, AllDifferentCuts = FALSE, OddHoleCuts = TRUE,
      #            RedSplitCuts = FALSE, LandPCuts = FALSE, PreProcessCuts = FALSE,
      #            ProbingCuts = TRUE, SimpleRoundingCuts = TRUE)
      result <- solve(prob, verbose = TRUE, solver = "CBC", maximumSeconds = 100)
      expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
  }
})
