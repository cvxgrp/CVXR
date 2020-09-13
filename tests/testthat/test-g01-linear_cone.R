context("test-g01-linear_cone")
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

accepts <- CVXR:::accepts
perform <- CVXR:::perform
invert <- CVXR:::invert
reduction_solve <- CVXR:::reduction_solve
ConeMatrixStuffing <- CVXR:::ConeMatrixStuffing
CvxAttr2Constr <- CVXR:::CvxAttr2Constr
FlipObjective <- CVXR:::FlipObjective
ExpCone <- CVXR:::ExpCone
SOC <- CVXR:::SOC

# solvers <- c(CVXR:::ECOS(), CVXR:::GUROBI(), CVXR:::MOSEK(), CVXR:::SCS(), CVXR:::CVXOPT(), CVXR:::GLPK())
solvers <- c(CVXR:::ECOS())

test_that("test scalar LP problems", {
  skip_on_cran()
  for(solver in solvers) {
    p <- Problem(Minimize(3*a), list(a >= 2))
    expect_true(accepts(ConeMatrixStuffing(), p))
    result <- solve(p, solver = name(solver))
    p_new <- perform(ConeMatrixStuffing(), p)
    result_new <- solve(p_new[[2]], solver = name(solver))
    expect_equal(result$value, result_new$value, tolerance = TOL)
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value, tolerance = TOL)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
    expect_equal(inv_sltn@opt_val, result$value, tolerance = TOL)
    expect_equal(inv_sltn@primal_vars[[as.character(id(a))]][1], result$getValue(a), tolerance = TOL)

    # TODO: Maximize.
    p <- Problem(Minimize(-3*a + b), list(a <= 2, b == a, b <= 5))
    result <- solve(p, solver = name(solver))
    expect_true(accepts(ConeMatrixStuffing(), p))
    p_new <- perform(ConeMatrixStuffing(), p)
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value, tolerance = TOL)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
    expect_equal(inv_sltn@opt_val, result$value, tolerance = TOL)
    expect_equal(inv_sltn@primal_vars[[as.character(id(a))]][1], result$getValue(a), tolerance = TOL)
    expect_equal(inv_sltn@primal_vars[[as.character(id(b))]][1], result$getValue(b), tolerance = TOL)

    # With a constant in the objective.
    p <- Problem(Minimize(3*a - b + 100), list(a >= 2, b + 5*c - 2 == a, b <= 5 + c))
    expect_true(accepts(ConeMatrixStuffing(), p))
    result <- solve(p, solver = name(solver))
    p_new <- perform(ConeMatrixStuffing(), p)
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value - 100, tolerance = TOL)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
    expect_equal(inv_sltn@opt_val, result$value, tolerance = TOL)
    expect_equal(inv_sltn@primal_vars[[as.character(id(a))]][1], result$getValue(a))
    expect_equal(inv_sltn@primal_vars[[as.character(id(b))]][1], result$getValue(b))

    # Unbounded problems.
    # TODO: Maximize.
    p <- Problem(Minimize(-a), list(a >= 2))
    expect_true(accepts(ConeMatrixStuffing(), p))
    if(name(solver) != "GUROBI") {
      result <- solve(p, solver = name(solver))
      p_new <- perform(ConeMatrixStuffing(), p)
      expect_true(accepts(solver, p_new[[2]]))
      sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
      expect_equal(sltn@opt_val, result$value, tolerance = TOL)
      inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
      expect_equal(inv_sltn@opt_val, result$value, tolerance = TOL)
    }

    # Infeasible problems.
    p <- Problem(Maximize(a), list(a >= 2, a <= 1))
    result <- solve(p, solver = name(solver))
    expect_true(accepts(FlipObjective(), p))
    p_min <- perform(FlipObjective(), p)
    expect_true(accepts(ConeMatrixStuffing(), p_min[[2]]))
    p_new <- perform(ConeMatrixStuffing(), p_min[[2]])
    result_new <- solve(p_new[[2]], solver = name(solver))
    expect_equal(result$value, -result_new$value, tolerance = TOL)
    expect_true(accepts(solver, p_new[[2]]))
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, -result$value, tolerance = TOL)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
    expect_equal(inv_sltn@opt_val, -result$value, tolerance = TOL)
    inv_flipped_sltn <- invert(FlipObjective(), inv_sltn, p_min[[3]])
    expect_equal(inv_flipped_sltn@opt_val, result$value, tolerance = TOL)
  }
})

test_that("test vector LP problems", {
  skip_on_cran()
  for(solver in solvers) {
    c <- Constant(matrix(1:2))
    p <- Problem(Minimize(t(c) %*% x), list(x >= c))
    result <- solve(p, solver = name(solver))
    expect_true(accepts(ConeMatrixStuffing(), p))
    p_new <- perform(ConeMatrixStuffing(), p)
    # result_new <- solve(p_new[[2]], solver = name(solver))
    # expect_equal(result$value, result_new$value)
    expect_true(accepts(solver, p_new[[2]]))
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value, tolerance = TOL)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])

    p_new1 <- perform(ConeMatrixStuffing(), p)
    expect_true(accepts(solver, p_new1[[2]]))
    sltn <- reduction_solve(solver, p_new1[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value, tolerance = TOL)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new1[[3]])

    expect_equal(inv_sltn@opt_val, result$value, tolerance = TOL)
    expect_equal(inv_sltn@primal_vars[[as.character(id(x))]], result$getValue(x))

    A <- value(Constant(cbind(c(3,5), c(1,2))))
    I <- Constant(diag(2))
    p <- Problem(Minimize(t(c) %*% x + a), list(A %*% x >= c(-1, 1), 4*I %*% z == x, z >= c(2, 2), a >= 2))
    expect_true(accepts(ConeMatrixStuffing(), p))
    result <- solve(p, solver = name(solver))
    p_new <- perform(ConeMatrixStuffing(), p)
    result_new <- solve(p_new[[2]], solver = name(solver))
    expect_equal(result$value, result_new$value, tolerance = TOL)
    expect_true(accepts(solver, p_new[[2]]))
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value, tolerance = 0.1)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
    expect_equal(inv_sltn@opt_val, result$value, tolerance = 0.1)
    for(var in variables(p))
      expect_equal(inv_sltn@primal_vars[[as.character(id(var))]], as.matrix(result$getValue(var)), tolerance = 0.1)
  }
})

test_that("test matrix LP problems", {
  skip_on_cran()
  for(solver in solvers) {
    Tmat <- value(Constant(matrix(1, nrow = 2, ncol = 2)))
    p <- Problem(Minimize(1 + a), list(A == Tmat + a, a >= 0))
    expect_true(accepts(ConeMatrixStuffing(), p))
    result <- solve(p, solver = name(solver))
    p_new <- perform(ConeMatrixStuffing(), p)
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value - 1, tolerance = TOL)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
    expect_equal(inv_sltn@opt_val, result$value, tolerance = TOL)
    for(var in variables(p))
      expect_equal(inv_sltn@primal_vars[[as.character(id(var))]], as.matrix(result$getValue(var)), tolerance = TOL)

    Tmat <- value(Constant(matrix(2, nrow = 2, ncol = 3)))
    p <- Problem(Minimize(1), list(A >= Tmat %*% C, A == B, C == t(Tmat)))
    expect_true(accepts(ConeMatrixStuffing(), p))
    result <- solve(p, solver = name(solver))
    p_new <- perform(ConeMatrixStuffing(), p)
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value - 1, tolerance = TOL)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
    expect_equal(inv_sltn@opt_val, result$value, tolerance = TOL)
    for(var in variables(p))
      expect_equal(inv_sltn@primal_vars[[as.character(id(var))]], as.matrix(result$getValue(var)), tolerance = TOL)
  }
})

test_that("test SOCP problems", {
  skip_on_cran()
  for(solver in solvers) {
    # Basic.
    p <- Problem(Minimize(b), list(p_norm(x, p = 2) <= b))
    pmod <- Problem(Minimize(b), list(SOC(b, x)))
    expect_true(accepts(ConeMatrixStuffing(), pmod))
    p_new <- perform(ConeMatrixStuffing(), pmod)
    if(accepts(solver, p_new[[2]])) {
      result <- solve(p, solver = name(solver))
      sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
      expect_equal(sltn@opt_val, result$value, tolerance = TOL)
      inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
      expect_equal(inv_sltn@opt_val, result$value, tolerance = TOL)
      for(var in variables(p))
        expect_equal(inv_sltn@primal_vars[[as.character(id(var))]], as.matrix(result$getValue(var)), tolerance = TOL)
    }

    # More complex.
    p <- Problem(Minimize(b), list(p_norm(x/2 + y[1:2], p = 2) <= b + 5, x >= 1, y == 5))
    pmod <- Problem(Minimize(b), list(SOC(b + 5, x/2 + y[1:2]), x >= 1, y == 5))
    expect_true(accepts(ConeMatrixStuffing(), pmod))
    result <- solve(p, solver = name(solver))
    p_new <- perform(ConeMatrixStuffing(), pmod)
    sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
    expect_equal(sltn@opt_val, result$value, tolerance = 1e-2)
    inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
    expect_equal(inv_sltn@opt_val, result$value, tolerance = 1e-2)
    for(var in variables(p))
      expect_equal(inv_sltn@primal_vars[[as.character(id(var))]], as.matrix(result$getValue(var)), tolerance = 1e-2)
  }
})

test_that("test exponential cone problems", {
  skip_on_cran()
  for(solver in solvers) {
    # Basic.
    p <- Problem(Minimize(b), list(exp(a) <= b, a >= 1))
    pmod <- Problem(Minimize(b), list(ExpCone(a, Constant(1), b), a >= 1))
    expect_true(accepts(ConeMatrixStuffing(), pmod))
    p_new <- perform(ConeMatrixStuffing(), pmod)
    if(accepts(solver, p_new[[2]])) {
      result <- solve(p, solver = name(solver))
      sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
      expect_equal(sltn@opt_val, result$value, tolerance = 0.1)
      inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
      expect_equal(inv_sltn@opt_val, result$value, tolerance = 0.1)
      for(var in variables(p))
        expect_equal(inv_sltn@primal_vars[[as.character(id(var))]], as.matrix(result$getValue(var)), tolerance = 0.1)
    }

    # More complex.
    # TODO: CVXOPT fails here.
    if(name(solver) != "CVXOPT") {
      p <- Problem(Minimize(b), list(exp(a/2 + c) <= b+5, a >= 1, c >= 5))
      pmod <- Problem(Minimize(b), list(ExpCone(a/2 + c, Constant(1), b+5), a >= 1, c >= 5))
      expect_true(accepts(ConeMatrixStuffing(), pmod))
      result <- solve(p, solver = name(solver))
      p_new <- perform(ConeMatrixStuffing(), pmod)
      sltn <- reduction_solve(solver, p_new[[2]], FALSE, FALSE, 1e-8, 1e-8, 1e-8, 10000, list())
      expect_equal(sltn@opt_val, result$value, tolerance = 1)
      inv_sltn <- invert(ConeMatrixStuffing(), sltn, p_new[[3]])
      expect_equal(inv_sltn@opt_val, result$value, tolerance = 1)
      for(var in variables(pmod))
        expect_equal(inv_sltn@primal_vars[[as.character(id(var))]], as.matrix(result$getValue(var)), tolerance = 1)
    }
  }
})

test_that("test positive semidefinite constraints", {
  skip_on_cran()
  C <- Variable(3,3)
  obj <- Maximize(C[1,3])
  constraints <- list(diag(C) == 1, C[1,2] == 0.6, C[2,3] == -0.3, C == t(C), C %>>% 0)
  prob <- Problem(obj, constraints)
  expect_true(accepts(FlipObjective(), prob))
  p_min <- perform(FlipObjective(), prob)
  expect_true(accepts(ConeMatrixStuffing(), p_min[[2]]))

  C <- Variable(2,2)
  obj <- Maximize(C[1,2])
  constraints <- list(C == 1, C %>>% diag(2,2))
  prob <- Problem(obj, constraints)
  expect_true(accepts(FlipObjective(), prob))
  p_min <- perform(FlipObjective(), prob)
  expect_true(accepts(ConeMatrixStuffing(), p_min[[2]]))

  C <- Variable(2,2, symmetric = TRUE)
  obj <- Minimize(C[1,1])
  constraints <- list(C %<<% diag(2,2))
  tmp <- perform(CvxAttr2Constr(), Problem(obj, constraints))
  expect_true(accepts(ConeMatrixStuffing(), tmp[[2]]))
})
