context("test-g05-conic_solvers-mosek")
TOL <- 1e-6

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
