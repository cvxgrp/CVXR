## Tests for CPLEX solver (QP path via Rcplex)

# ── Constant / infrastructure tests (no solver needed) ──────────────

## @cvxpy NONE
test_that("CPLEX_SOLVER constant is exported", {
  expect_equal(CPLEX_SOLVER, "CPLEX")
})

## @cvxpy NONE
test_that("CPLEX appears in installed_solvers() when Rcplex available", {
  skip_if_not_installed("Rcplex")
  expect_true("CPLEX" %in% installed_solvers())
})

# ── LP tests ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("CPLEX solves basic LP", {
  skip_if_not_installed("Rcplex")
  ## min x1 + 2*x2  s.t. x1 + x2 >= 3, x >= 0
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x >= 0))
  val <- psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 3.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(3, 0), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("CPLEX solves LP with equality constraint", {
  skip_if_not_installed("Rcplex")
  ## min x1 + 2*x2  s.t. x1 + x2 == 4, x1 <= 3, x >= 0
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] == 4, x[1] <= 3, x >= 0))
  val <- psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 5.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(3, 1), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("CPLEX solves Maximize LP", {
  skip_if_not_installed("Rcplex")
  ## max x1 + 2*x2  s.t. x1 + x2 <= 3, x1 <= 2, x2 <= 2, x >= 0
  x <- Variable(2)
  prob <- Problem(Maximize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] <= 3, x[1] <= 2, x[2] <= 2, x >= 0))
  val <- psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 5.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-5)
})

# ── QP tests ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("CPLEX solves basic QP", {
  skip_if_not_installed("Rcplex")
  ## min quad_form(y, 2*I) + sum(y)  s.t. y1 + y2 == 1, y >= 0
  ## quad_form(y, 2I) = y'(2I)y = 2*(y1^2 + y2^2). At y=(0.5,0.5): 1+1 = 2.
  y <- Variable(2)
  prob <- Problem(Minimize(quad_form(y, 2 * diag(2)) + sum(y)),
                  list(sum(y) == 1, y >= 0))
  val <- psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 2.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(y)), c(0.5, 0.5), tolerance = 1e-5)
})

# ── MIP tests ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("CPLEX solves boolean MIP", {
  skip_if_not_installed("Rcplex")
  ## max 3*b1 + 2*b2 + b3  s.t. b1 + b2 + b3 <= 2
  b <- Variable(3, boolean = TRUE)
  prob <- Problem(Maximize(3 * b[1] + 2 * b[2] + b[3]),
                  list(b[1] + b[2] + b[3] <= 2))
  val <- psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 5.0, tolerance = 1e-5)
  bval <- as.numeric(value(b))
  expect_equal(bval[1], 1, tolerance = 1e-5)
  expect_equal(bval[2], 1, tolerance = 1e-5)
  expect_equal(bval[3], 0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("CPLEX solves integer MIP", {
  skip_if_not_installed("Rcplex")
  ## min w1 + 2*w2 + 3*w3  s.t. w1 + w2 + w3 >= 5, w >= 0
  w <- Variable(3, integer = TRUE)
  prob <- Problem(Minimize(w[1] + 2 * w[2] + 3 * w[3]),
                  list(w[1] + w[2] + w[3] >= 5, w >= 0))
  val <- psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 5.0, tolerance = 1e-5)
  wval <- as.numeric(value(w))
  expect_equal(wval[1], 5, tolerance = 1e-5)
  expect_equal(wval[2], 0, tolerance = 1e-5)
  expect_equal(wval[3], 0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("CPLEX solves MIQP", {
  skip_if_not_installed("Rcplex")
  ## min quad_form(z, I) + sum(z)  s.t. sum(z) >= 3, z integer, z >= 0
  z <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(quad_form(z, diag(2)) + sum(z)),
                  list(sum(z) >= 3, z >= 0))
  val <- psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  ## z=(1,2) or z=(2,1): obj = 1+4+3 = 8 or 4+1+3 = 8
  expect_equal(val, 8.0, tolerance = 1e-5)
})

# ── Infeasible ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("CPLEX detects infeasible problem", {
  skip_if_not_installed("Rcplex")
  z <- Variable(1)
  prob <- Problem(Minimize(z), list(z >= 5, z <= 3))
  val <- psolve(prob, solver = "CPLEX")
  expect_true(status(prob) %in%
                c(INFEASIBLE, INFEASIBLE_OR_UNBOUNDED))
})

# ── Dual values ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("CPLEX returns correct LP duals", {
  skip_if_not_installed("Rcplex")
  ## min x1 + 2*x2  s.t. x1 + x2 == 4 (c1), x1 <= 3 (c2), x >= 0
  x <- Variable(2)
  c1 <- x[1] + x[2] == 4
  c2 <- x[1] <= 3
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(c1, c2, x >= 0))
  psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  ## CVXPY verified: eq_dual = -2.0, ineq_dual = 1.0
  expect_equal(as.numeric(dual_value(c1)), -2.0, tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(c2)), 1.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("CPLEX returns correct QP duals", {
  skip_if_not_installed("Rcplex")
  ## min 0.5*quad_form(y, [[2,0],[0,1]]) + y1
  ##     s.t. y1 + y2 == 1 (c1), y1 >= 0.1 (c2), y >= 0
  y <- Variable(2)
  Q <- matrix(c(2, 0, 0, 1), 2, 2)
  c1 <- y[1] + y[2] == 1
  c2 <- y[1] >= 0.1
  prob <- Problem(Minimize(0.5 * quad_form(y, Q) + y[1]),
                  list(c1, c2, y >= 0))
  psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(as.numeric(value(y)), c(0.1, 0.9), tolerance = 1e-3)
  ## CVXPY verified: eq_dual = -0.9, ineq_dual = 0.3
  expect_equal(as.numeric(dual_value(c1)), -0.9, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2)), 0.3, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("CPLEX Maximize LP duals", {
  skip_if_not_installed("Rcplex")
  ## max x1 + 2*x2  s.t. x1+x2<=3 (c1), x1<=2 (c2), x2<=2 (c3), x>=0
  x <- Variable(2)
  c1 <- x[1] + x[2] <= 3
  c2 <- x[1] <= 2
  c3 <- x[2] <= 2
  prob <- Problem(Maximize(x[1] + 2 * x[2]),
                  list(c1, c2, c3, x >= 0))
  psolve(prob, solver = "CPLEX")
  expect_equal(status(prob), OPTIMAL)
  ## CVXPY verified: c1 dual=1.0, c2 dual=0.0, c3 dual=1.0
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(c2)), 0.0, tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(c3)), 1.0, tolerance = 1e-4)
})

# ── Cross-solver parity ─────────────────────────────────────────────

## @cvxpy NONE
test_that("CPLEX matches Clarabel on LP values and duals", {
  skip_if_not_installed("Rcplex")
  skip_if_not_installed("clarabel")
  ## min 3*x1 + x2  s.t. x1 + x2 == 5 (c1), x1 >= 1 (c2), x2 >= 1 (c3)
  x <- Variable(2)
  c1_cplex <- x[1] + x[2] == 5
  c2_cplex <- x[1] >= 1
  c3_cplex <- x[2] >= 1
  prob_cplex <- Problem(Minimize(3 * x[1] + x[2]),
                        list(c1_cplex, c2_cplex, c3_cplex))
  psolve(prob_cplex, solver = "CPLEX")

  x2 <- Variable(2)
  c1_clar <- x2[1] + x2[2] == 5
  c2_clar <- x2[1] >= 1
  c3_clar <- x2[2] >= 1
  prob_clar <- Problem(Minimize(3 * x2[1] + x2[2]),
                       list(c1_clar, c2_clar, c3_clar))
  psolve(prob_clar, solver = "CLARABEL")

  expect_equal(value(prob_cplex), value(prob_clar), tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), as.numeric(value(x2)), tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(c1_cplex)),
               as.numeric(dual_value(c1_clar)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2_cplex)),
               as.numeric(dual_value(c2_clar)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("CPLEX matches OSQP on QP values and duals", {
  skip_if_not_installed("Rcplex")
  skip_if_not_installed("osqp")
  ## min quad_form(y, I) + sum(y)  s.t. y1+y2==2 (c1), y>=0
  y <- Variable(2)
  c1_cplex <- y[1] + y[2] == 2
  prob_cplex <- Problem(Minimize(quad_form(y, diag(2)) + sum(y)),
                        list(c1_cplex, y >= 0))
  psolve(prob_cplex, solver = "CPLEX")

  y2 <- Variable(2)
  c1_osqp <- y2[1] + y2[2] == 2
  prob_osqp <- Problem(Minimize(quad_form(y2, diag(2)) + sum(y2)),
                       list(c1_osqp, y2 >= 0))
  psolve(prob_osqp, solver = "OSQP")

  expect_equal(value(prob_cplex), value(prob_osqp), tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), as.numeric(value(y2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c1_cplex)),
               as.numeric(dual_value(c1_osqp)), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("CPLEX matches HiGHS on MIP values", {
  skip_if_not_installed("Rcplex")
  skip_if_not_installed("highs")
  ## min w1 + 2*w2  s.t. w1 + w2 >= 4, w integer, w >= 0
  w <- Variable(2, integer = TRUE)
  prob_cplex <- Problem(Minimize(w[1] + 2 * w[2]),
                        list(w[1] + w[2] >= 4, w >= 0))
  val_cplex <- psolve(prob_cplex, solver = "CPLEX")

  w2 <- Variable(2, integer = TRUE)
  prob_highs <- Problem(Minimize(w2[1] + 2 * w2[2]),
                        list(w2[1] + w2[2] >= 4, w2 >= 0))
  val_highs <- psolve(prob_highs, solver = "HIGHS")

  expect_equal(val_cplex, val_highs, tolerance = 1e-4)
})
