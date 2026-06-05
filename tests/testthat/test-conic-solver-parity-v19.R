## CVXPY 1.9 conic-solver parity additions.

infeasible_status <- function(prob) {
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
  invisible(prob)
}

## @cvxpy test_conic_solvers.py::TestClarabel::test_infeasible_lp_ineq_constraints test_conic_solvers.py::TestSCS::test_infeasible_lp_ineq_constraints test_conic_solvers.py::TestClarabel::test_infeasible_lp_eq_constraints test_conic_solvers.py::TestSCS::test_infeasible_lp_eq_constraints
test_that("conic solvers report infeasible LP certificates", {
  for (solver in c("CLARABEL", "SCS")) {
    x <- Variable(2)
    A <- matrix(c(1, 0, 0, 1, -1, -1), 3, 2, byrow = TRUE)
    ineq <- Problem(Minimize(0), list(A %*% x <= c(0, 0, -1)))
    psolve(ineq, solver = solver)
    infeasible_status(ineq)
    expect_null(value(x))

    y <- Variable(2)
    B <- matrix(c(1, 0, 0, 1, 1, 1), 3, 2, byrow = TRUE)
    eq <- Problem(Minimize(0), list(B %*% y == c(0, 0, 1), y >= 0))
    psolve(eq, solver = solver)
    infeasible_status(eq)
    expect_null(value(y))
  }
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_infeasible_soc test_conic_solvers.py::TestSCS::test_infeasible_soc test_conic_solvers.py::TestClarabel::test_infeasible_exp_cone test_conic_solvers.py::TestSCS::test_infeasible_exp_cone test_conic_solvers.py::TestClarabel::test_infeasible_psd_cone test_conic_solvers.py::TestSCS::test_infeasible_psd_cone test_conic_solvers.py::TestClarabel::test_infeasible_soc_exp_mixed test_conic_solvers.py::TestSCS::test_infeasible_soc_exp_mixed
test_that("conic solvers report infeasible SOC, exponential, and PSD cones", {
  for (solver in c("CLARABEL", "SCS")) {
    x <- Variable(2)
    soc <- Problem(Minimize(0), list(SOC(-1, x)))
    psolve(soc, solver = solver)
    infeasible_status(soc)
    expect_null(value(x))

    a <- Variable()
    b <- Variable()
    c <- Variable()
    exp_prob <- Problem(Minimize(0), list(ExpCone(a, b, c), a == 1, b == 1, c == 1))
    psolve(exp_prob, solver = solver)
    infeasible_status(exp_prob)
    expect_true(all(vapply(variables(exp_prob), function(v) is.null(value(v)), logical(1))))

    X <- Variable(c(2, 2))
    psd_prob <- Problem(Minimize(0), list(X %>>% 0, X[1, 1] == -1))
    psolve(psd_prob, solver = solver)
    infeasible_status(psd_prob)
    expect_null(value(X))

    u <- Variable()
    v <- Variable()
    mixed <- Problem(Maximize(0), list(u == v, SOC(u, Constant(1)), ExpCone(u, v, 1)))
    psolve(mixed, solver = solver)
    infeasible_status(mixed)
    expect_true(all(vapply(variables(mixed), function(var) is.null(value(var)), logical(1))))
  }
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_infeasible_power_cone_3d test_conic_solvers.py::TestSCS::test_infeasible_power_cone_3d test_conic_solvers.py::TestClarabel::test_infeasible_power_cone_nd
test_that("conic solvers report infeasible power cones where supported", {
  for (solver in c("CLARABEL", "SCS")) {
    z <- Variable()
    p3d <- Problem(Minimize(0), list(PowCone3D(-1, 1, z, 0.5)))
    psolve(p3d, solver = solver)
    infeasible_status(p3d)
    expect_null(value(z))
  }

  W <- Variable(3)
  z <- Variable()
  pnd <- Problem(Minimize(0), list(PowConeND(W, z, Constant(matrix(c(0.2, 0.3, 0.5), 3, 1))),
                                   W == 1, z == 2))
  psolve(pnd, solver = "CLARABEL")
  infeasible_status(pnd)
  expect_true(all(vapply(variables(pnd), function(v) is.null(value(v)), logical(1))))
})

## @cvxpy test_conic_solvers.py::test_offset_in_opt_val
test_that("partial_optimize preserves constant objective offsets", {
  x <- Variable()
  t <- Variable()
  inner <- Problem(Minimize(t + 1000), list(t >= x))
  f <- partial_optimize(inner, opt_vars = list(t))
  value(x) <- 0
  expect_equal(as.numeric(value(f)), 1000, tolerance = 1e-6)
})
