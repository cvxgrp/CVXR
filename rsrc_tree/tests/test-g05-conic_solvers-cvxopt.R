context("test-g05-conic_solvers-cvxopt")
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
