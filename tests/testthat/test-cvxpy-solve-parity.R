## CVXPY solve parity tests
## Mirrors tests from CVXPY test_problem.py, test_conic_solvers.py
## Verifies optimal values, variable values, expression values, and status.
##
## CVXPY source: /Users/naras/GitHub/cvxpy branch claude

# ── LP: Scalar problems ─────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py test_scalar_lp()

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("scalar LP: min 3a, a >= 2", {
  a <- Variable(1)
  p <- Problem(Minimize(3 * a), list(a >= 2))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 6, tolerance = 1e-4)
  expect_equal(as.numeric(value(a)), 2, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("scalar LP: max 3a - b, with equalities", {
  a <- Variable(1)
  b <- Variable(1)
  p <- Problem(Maximize(3 * a - b),
               list(a <= 2, b == a, b <= 5))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 4.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(a)), 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(b)), 2, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("scalar LP: min with constant offset", {
  a <- Variable(1)
  b <- Variable(1)
  c_var <- Variable(1)
  p <- Problem(Minimize(3 * a - b + 100),
               list(a >= 2,
                    b + 5 * c_var - 2 == a,
                    b <= 5 + c_var))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 101 + 1.0 / 6, tolerance = 1e-3)
  expect_equal(as.numeric(value(a)), 2, tolerance = 1e-3)
  expect_equal(as.numeric(value(b)), 5 - 1.0 / 6, tolerance = 1e-3)
  expect_equal(as.numeric(value(c_var)), -1.0 / 6, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("scalar LP: status and dual value populated", {
  a <- Variable(1)
  constr <- list(a <= 2)
  p <- Problem(Maximize(a), constr)
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, value(p))
  expect_equal(status(p), OPTIMAL)
  expect_true(!is.null(value(a)))
  expect_true(!is.null(dual_value(constr[[1]])))
})

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("scalar LP: unbounded", {
  a <- Variable(1)
  p <- Problem(Maximize(a), list(a >= 2))
  psolve(p, verbose = FALSE)
  expect_equal(status(p), UNBOUNDED)
  expect_true(is.infinite(value(p)))
  expect_true(value(p) > 0)  # Maximize unbounded → +Inf (can maximize without bound)
})

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("scalar LP: infeasible maximize", {
  a <- Variable(1)
  p <- Problem(Maximize(a), list(a >= 2, a <= 1))
  psolve(p, verbose = FALSE)
  expect_equal(status(p), INFEASIBLE)
  expect_true(is.infinite(value(p)))
  expect_true(value(p) < 0)
})

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("scalar LP: infeasible minimize", {
  a <- Variable(1)
  p <- Problem(Minimize(-a), list(a >= 2, a <= 1))
  psolve(p, verbose = FALSE)
  expect_equal(status(p), INFEASIBLE)
  expect_true(is.infinite(value(p)))
  expect_true(value(p) > 0)
})

# ── LP: norm1 ────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py test_norm1()

## @cvxpy test_problem.py::TestProblem::test_norm1
test_that("norm1: scalar min", {
  a <- Variable(1)
  p <- Problem(Minimize(norm1(a)), list(a <= -2))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(a)), -2, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_norm1
test_that("norm1: vector with offset", {
  x <- Variable(2)
  z <- Variable(2)
  p <- Problem(Minimize(norm1(x - z) + 5),
               list(x >= c(2, 3), z <= c(-1, -4)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 15, tolerance = 1e-3)
  xv <- as.numeric(value(x))
  zv <- as.numeric(value(z))
  expect_equal(xv[2] - zv[2], 7, tolerance = 1e-3)
})

# ── LP: norm_inf ─────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py test_norm_inf()

## @cvxpy test_problem.py::TestProblem::test_norm_inf
test_that("norm_inf: scalar min", {
  a <- Variable(1)
  p <- Problem(Minimize(norm_inf(a)), list(a >= 2))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(a)), 2, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_norm_inf
test_that("norm_inf: combined scalars", {
  a <- Variable(1)
  b <- Variable(1)
  c_var <- Variable(1)
  p <- Problem(Minimize(3 * norm_inf(a + 2 * b) + c_var),
               list(a >= 2, b <= -1, c_var == 3))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 3, tolerance = 1e-3)
  expect_equal(as.numeric(value(a)) + 2 * as.numeric(value(b)), 0, tolerance = 1e-3)
  expect_equal(as.numeric(value(c_var)), 3, tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_norm_inf
test_that("norm_inf: vector", {
  x <- Variable(2)
  z <- Variable(2)
  p <- Problem(Minimize(norm_inf(x - z) + 5),
               list(x >= c(2, 3), z <= c(-1, -4)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 12, tolerance = 1e-3)
})

# ── LP: abs ──────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py test_abs()

## @cvxpy test_problem.py::TestProblem::test_abs
test_that("abs: sum of matrix", {
  A <- Variable(c(2, 2))
  p <- Problem(Minimize(sum(abs(A))), list(-2 >= A))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 8, tolerance = 1e-4)
  expect_equal(as.numeric(value(A)), rep(-2, 4), tolerance = 1e-4)
})

# ── LP: multiply ─────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py test_multiply()

## @cvxpy test_problem.py::TestProblem::test_multiply
test_that("multiply: elementwise with norm_inf", {
  A <- Variable(c(2, 2))
  c_mat <- matrix(c(1, 2, -1, -2), 2, 2)
  expr <- c_mat * A
  p <- Problem(Minimize(norm_inf(expr)), list(A == 5))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 10, tolerance = 1e-4)
})

# ── SOCP ─────────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py, test_conic_solvers.py

## @cvxpy NONE
test_that("SOCP: p_norm minimization", {
  x <- Variable(3)
  p <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, sqrt(3), tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("SOCP: sum_squares", {
  x <- Variable(2)
  p <- Problem(Minimize(sum_squares(x)), list(x >= 1))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("SOCP: quad_over_lin", {
  x <- Variable(2)
  a <- Variable(1)
  p <- Problem(Minimize(quad_over_lin(x, a)),
               list(x >= 1, a == 2))
  result <- psolve(p, verbose = FALSE)
  ## x = (1,1), a = 2 → (1+1)/2 = 1
  expect_equal(result, 1, tolerance = 1e-4)
})

# ── ExpCone ──────────────────────────────────────────────────────────
## CVXPY SOURCE: test_conic_solvers.py test_exp(), test_log(), test_entr()

## @cvxpy test_conic_solvers.py::TestSCS::test_exp
test_that("exp: minimize sum(exp(x)), sum(x) == 1, n = 5", {
  x <- Variable(5)
  p <- Problem(Minimize(sum(exp(x))), list(sum(x) == 1))
  psolve(p, verbose = FALSE)
  expect_equal(as.numeric(value(x)), rep(1 / 5, 5), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_exp
test_that("exp: minimize sum(exp(x)), sum(x) == 1, n = 10", {
  x <- Variable(10)
  p <- Problem(Minimize(sum(exp(x))), list(sum(x) == 1))
  psolve(p, verbose = FALSE)
  expect_equal(as.numeric(value(x)), rep(1 / 10, 10), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_log
test_that("log: maximize sum(log(x)), sum(x) == 1, n = 5", {
  x <- Variable(5)
  p <- Problem(Maximize(sum(log(x))), list(sum(x) == 1))
  psolve(p, verbose = FALSE)
  expect_equal(as.numeric(value(x)), rep(1 / 5, 5), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_log
test_that("log: maximize sum(log(x)), sum(x) == 1, n = 10", {
  x <- Variable(10)
  p <- Problem(Maximize(sum(log(x))), list(sum(x) == 1))
  psolve(p, verbose = FALSE)
  expect_equal(as.numeric(value(x)), rep(1 / 10, 10), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_entr
test_that("entr: maximize sum(entr(x)), sum(x) == 1, n = 5", {
  x <- Variable(5)
  p <- Problem(Maximize(sum(entr(x))), list(sum(x) == 1))
  psolve(p, verbose = FALSE)
  expect_equal(as.numeric(value(x)), rep(1 / 5, 5), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_entr
test_that("entr: maximize sum(entr(x)), sum(x) == 1, n = 10", {
  x <- Variable(10)
  p <- Problem(Maximize(sum(entr(x))), list(sum(x) == 1))
  psolve(p, verbose = FALSE)
  expect_equal(as.numeric(value(x)), rep(1 / 10, 10), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_log_problem
test_that("log: objective and bounds", {
  ## CVXPY SOURCE: test_conic_solvers.py test_log_problem()
  x <- Variable(2)
  p <- Problem(Maximize(sum(log(x))), list(x <= c(1, exp(1))))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 1, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, exp(1)), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_log_problem
test_that("log: in constraint", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(log(x) >= 0, x <= c(1, 1)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 2, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-3)
})

# ── Expression values after solve ────────────────────────────────────
## CVXPY SOURCE: test_problem.py test_expression_values()

## @cvxpy test_problem.py::TestProblem::test_expression_values
test_that("expression values populated after solve", {
  x <- Variable(2)
  z <- Variable(2)
  diff_exp <- x - z
  inf_exp <- norm_inf(diff_exp)
  sum_exp <- 5 + norm1(z) + norm1(x) + inf_exp
  constr_exp <- p_norm(x + z, 2)
  obj <- p_norm(sum_exp, 2)
  p <- Problem(Minimize(obj),
               list(x >= c(2, 3), z <= c(-1, -4), constr_exp <= 2))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 22, tolerance = 1e-2)
  xv <- as.numeric(value(x))
  zv <- as.numeric(value(z))
  expect_equal(xv, c(2, 3), tolerance = 1e-2)
  expect_equal(zv, c(-1, -4), tolerance = 1e-2)
  ## Check intermediate expression values
  expect_equal(as.numeric(value(diff_exp)), xv - zv, tolerance = 1e-3)
  expect_equal(as.numeric(value(inf_exp)),
               max(abs(xv - zv)), tolerance = 1e-3)
  expect_equal(as.numeric(value(constr_exp)),
               sqrt(sum((xv + zv)^2)), tolerance = 1e-3)
})

# ── VStack / HStack ─────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py test_vstack(), test_hstack()

## @cvxpy test_problem.py::TestProblem::test_vstack
test_that("vstack: column assembly", {
  x <- Variable(c(2, 1))
  y <- Variable(c(3, 1))
  cc <- matrix(1, 1, 5)
  p <- Problem(Minimize(cc %*% vstack(x, y)),
               list(x == matrix(c(1, 2), 2, 1),
                    y == matrix(c(3, 4, 5), 3, 1)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 15, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_vstack
test_that("vstack: matrix assembly", {
  A <- Variable(c(2, 2))
  C <- Variable(c(2, 2))
  cc <- matrix(1, 2, 2)
  ## A >= cc (all 1's), C == -2 → sum = 4 + (-8) = -4
  p <- Problem(Minimize(sum(vstack(A, C))),
               list(A >= cc, C == -2))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, -4, tolerance = 1e-3)
})

# ── Trace ────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("trace: in objective", {
  A <- Variable(c(2, 2))
  p <- Problem(Minimize(matrix_trace(A)), list(A >= diag(2)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 2, tolerance = 1e-4)
})

# ── Reshape ──────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py test_reshape()

## @cvxpy test_problem.py::TestProblem::test_reshape
test_that("reshape: matrix to vector", {
  A <- Variable(c(2, 2))
  cc <- c(1, 2, 3, 4)
  expr <- reshape_expr(A, c(4, 1))
  p <- Problem(Minimize(t(expr) %*% cc),
               list(A == matrix(c(-1, -2, 3, 4), 2, 2)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 20, tolerance = 1e-3)
})

# ── DiagVec / DiagMat ────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_diag_prob
test_that("diag: extract and constrain diagonal", {
  A <- Variable(c(3, 3))
  ## DiagMat extracts diagonal from matrix; DiagVec creates diagonal matrix from vector
  p <- Problem(Minimize(sum(A)),
               list(DiagMat(A) == matrix(c(1, 2, 3), 3, 1), A >= 0))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 6, tolerance = 1e-4)
  Av <- value(A)
  expect_equal(diag(Av), c(1, 2, 3), tolerance = 1e-4)
})

# ── Nonneg / Nonpos variable attributes ─────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_pos
test_that("nonneg variable: bounded by 0", {
  x <- Variable(2, nonneg = TRUE)
  p <- Problem(Minimize(sum(x)), list(x <= 3))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-5)
})

## @cvxpy test_problem.py::TestProblem::test_pos
test_that("nonpos variable: bounded by 0", {
  x <- Variable(2, nonpos = TRUE)
  p <- Problem(Maximize(sum(x)), list(x >= -3))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-5)
})

# ── SCS solver variants ─────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("SCS: scalar LP", {
  a <- Variable(1)
  p <- Problem(Minimize(3 * a), list(a >= 2))
  result <- psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(result, 6, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_exp
test_that("SCS: exp minimize", {
  x <- Variable(5)
  p <- Problem(Minimize(sum(exp(x))), list(sum(x) == 1))
  psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(as.numeric(value(x)), rep(1 / 5, 5), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_log
test_that("SCS: log maximize", {
  x <- Variable(5)
  p <- Problem(Maximize(sum(log(x))), list(sum(x) == 1))
  psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(as.numeric(value(x)), rep(1 / 5, 5), tolerance = 1e-3)
})

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("SCS: infeasible", {
  a <- Variable(1)
  p <- Problem(Minimize(a), list(a >= 2, a <= 1))
  psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(status(p), INFEASIBLE)
})

# ── Solve twice: consistency ─────────────────────────────────────────

## @cvxpy test_conic_solvers.py::TestSCS::test_solve_problem_twice
test_that("solve twice gives consistent results", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 1))
  v1 <- psolve(p, verbose = FALSE)
  v2 <- psolve(p, verbose = FALSE)
  expect_equal(v1, v2, tolerance = 1e-6)
})

# ── Dual values: complementary slackness ─────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_dual_variables
test_that("dual values: active vs inactive constraints", {
  x <- Variable(1)
  c1 <- x >= 1     # active at optimal
  c2 <- x >= -10   # inactive
  p <- Problem(Minimize(x), list(c1, c2))
  psolve(p, verbose = FALSE)
  ## Active constraint should have nonzero dual
  dv1 <- as.numeric(dual_value(c1))
  expect_true(abs(dv1) > 0.1)
  ## Inactive constraint should have ~zero dual
  dv2 <- as.numeric(dual_value(c2))
  expect_equal(dv2, 0, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_dual_variables
test_that("dual values: equality constraint", {
  x <- Variable(2)
  eq <- x[1] + x[2] == 5
  p <- Problem(Minimize(x[1] - x[2]), list(eq, x >= 0))
  psolve(p, verbose = FALSE)
  expect_equal(as.numeric(value(x[1])), 0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x[2])), 5, tolerance = 1e-4)
  ## Dual for equality should be non-null
  expect_true(!is.null(dual_value(eq)))
})

# ── Constraint violation ─────────────────────────────────────────────

## @cvxpy NONE
test_that("constraint violation near zero after solve", {
  x <- Variable(2)
  c1 <- x >= c(1, 2)
  c2 <- sum(x) <= 10
  p <- Problem(Minimize(sum(x)), list(c1, c2))
  psolve(p, verbose = FALSE)
  expect_true(all(as.numeric(violation(c1)) <= 1e-5))
  expect_equal(max(as.numeric(violation(c2))), 0, tolerance = 1e-5)
})

# ── Multiple variable problems ───────────────────────────────────────

## @cvxpy NONE
test_that("three variables with mixed constraints", {
  x <- Variable(1)
  y <- Variable(1)
  z <- Variable(1)
  p <- Problem(Minimize(x + 2 * y + 3 * z),
               list(x >= 1, y >= 2, z >= 3, x + y + z <= 10))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 1 + 4 + 9, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 1, tolerance = 1e-4)
  expect_equal(as.numeric(value(y)), 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(z)), 3, tolerance = 1e-4)
})

# ── Matrix variable LP ──────────────────────────────────────────────

## @cvxpy NONE
test_that("matrix variable with element constraints", {
  X <- Variable(c(2, 2))
  p <- Problem(Minimize(sum(X)),
               list(X >= 0, X[1, 1] == 3, X[2, 2] == 5))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 8, tolerance = 1e-4)
  Xv <- value(X)
  expect_equal(Xv[1, 1], 3, tolerance = 1e-4)
  expect_equal(Xv[2, 2], 5, tolerance = 1e-4)
  expect_equal(Xv[1, 2], 0, tolerance = 1e-4)
  expect_equal(Xv[2, 1], 0, tolerance = 1e-4)
})

# ── Cumsum in constraints ────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_cumsum
test_that("cumsum in constraint", {
  x <- Variable(3)
  ## cumsum(x) >= c(1, 3, 6) means x1>=1, x1+x2>=3, x1+x2+x3>=6
  p <- Problem(Minimize(sum(x)), list(cumsum(x) >= c(1, 3, 6)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 6, tolerance = 1e-3)
})

# ── Power atom ───────────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_power
test_that("power: square in objective", {
  x <- Variable(1)
  p <- Problem(Minimize(power(x, 2)), list(x >= 3))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 9, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 3, tolerance = 1e-3)
})

# ── sum_largest ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("sum_largest: top-k elements", {
  x <- Variable(4)
  p <- Problem(Minimize(sum_largest(x, 2)),
               list(x >= c(1, 2, 3, 4)))
  result <- psolve(p, verbose = FALSE)
  ## Minimize top-2 sum → x = (1, 2, 3, 4), top-2 = 3+4 = 7
  expect_equal(result, 7, tolerance = 1e-3)
})

# ── inv_pos ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("inv_pos: minimize 1/x with equality", {
  x <- Variable(1)
  ## Minimize 1/x with x == 2 → result = 0.5
  ## (with only x >= 2, the infimum is 0 as x→∞, not a useful test)
  p <- Problem(Minimize(inv_pos(x)), list(x == 2))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 0.5, tolerance = 1e-3)
})

# ── pos / neg atoms ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("pos atom in objective", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(pos(x))), list(x == c(-1, 2)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 2, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("neg atom in objective", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(neg(x))), list(x == c(-1, 2)))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 1, tolerance = 1e-4)
})

# ── max_elemwise / min_elemwise ──────────────────────────────────────

## @cvxpy NONE
test_that("max_elemwise in objective", {
  x <- Variable(1)
  p <- Problem(Minimize(max_elemwise(x, 0)), list(x >= -5))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 0, tolerance = 1e-5)
  expect_true(as.numeric(value(x)) <= 0 + 1e-4)
})

## @cvxpy NONE
test_that("min_elemwise in objective", {
  x <- Variable(1)
  p <- Problem(Maximize(min_elemwise(x, 5)), list(x <= 10))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 5, tolerance = 1e-4)
})

# ── Kron ─────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("kron: Kronecker product constraint", {
  x <- Variable(c(2, 1))
  I2 <- diag(2)
  ## kron(I2, x) = [x; 0; 0; x] (4x1)
  p <- Problem(Minimize(sum(kron(I2, x))),
               list(x >= 1))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 4, tolerance = 1e-4)
})

# ── upper_tri ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("upper_tri: extract and minimize", {
  X <- Variable(c(2, 2))
  p <- Problem(Minimize(sum(upper_tri(X))),
               list(X == matrix(c(1, 2, 3, 4), 2, 2)))
  result <- psolve(p, verbose = FALSE)
  ## upper_tri extracts off-diagonal upper triangle: X[1,2] = 3
  expect_equal(result, 3, tolerance = 1e-4)
})
