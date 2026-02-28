## CVXPY Parity Tests: MEDIUM-Priority Miscellaneous Gaps
## =====================================================
## Tests covering:
##   - Constraints: duals (equality, inequality, PSD) + violation()
##   - Cone2Cone: SOC, ExpCone, PSD, PowCone preserved
##   - Quad Form: constant x, sparse P, 1x1, quad_over_lin equivalence
##   - Quadratic: sum_squares, matrix_frac, quad_over_lin solving, is_quadratic
##   - Perspective: exp, square, parameter
##   - Power Atom: power(x,2), power(x,0.5), power(x,-1), fractional, curvature, mono, numeric
##
## Reference values verified via `uv run python` against CVXPY 1.8.1.

library(testthat)
library(CVXR)

# =====================================================================
# CONSTRAINT TESTS (4 tests)
# =====================================================================
# CVXPY SOURCE: test_problem.py::test_dual_variables, test_constraints.py

# -- 1. equality_constraint_dual -----------------------------------------------
# CVXPY: min sum(x) s.t. x == [3, 4]
# Dual: [-1, -1] (each unit increase in RHS decreases obj by 1)

## @cvxpy test_problem.py::TestProblem::test_dual_variables
test_that("equality_constraint_dual: equality constraint dual values (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  x <- Variable(2)
  constr <- (x == c(3, 4))
  prob <- Problem(Minimize(sum_entries(x)), list(constr))
  val <- psolve(prob, solver = CLARABEL_SOLVER)

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 7, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(3, 4), tolerance = 1e-5)

  ## Dual values
  dv <- dual_value(constr)
  expect_equal(as.numeric(dv), c(-1, -1), tolerance = 1e-4)
})

# -- 2. inequality_constraint_dual ---------------------------------------------
# CVXPY: min sum(x) s.t. x >= [3, 4]
# Dual: [1, 1] (each unit relaxation of RHS decreases obj by 1)

## @cvxpy test_problem.py::TestProblem::test_dual_variables
test_that("inequality_constraint_dual: inequality constraint dual values (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  x <- Variable(2)
  constr <- (x >= c(3, 4))
  prob <- Problem(Minimize(sum_entries(x)), list(constr))
  val <- psolve(prob, solver = CLARABEL_SOLVER)

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 7, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(3, 4), tolerance = 1e-5)

  ## Dual values
  dv <- dual_value(constr)
  expect_equal(as.numeric(dv), c(1, 1), tolerance = 1e-4)
})

# -- 3. psd_constraint_dual ---------------------------------------------------
# CVXPY: max C[0,0] s.t. C << [[2,0],[0,2]]
# test_problem.py::test_psd_duals
# opt value = 2, dual = [[1,0],[0,0]]

## @cvxpy test_problem.py::TestProblem::test_psd_duals
test_that("psd_constraint_dual: PSD constraint dual values (CVXPY parity)", {
  skip_if_not_installed("scs")

  C <- Variable(c(2, 2), symmetric = TRUE)
  bound <- matrix(c(2, 0, 0, 2), 2, 2)
  constr <- (C %<<% bound)  ## C << [[2,0],[0,2]]
  prob <- Problem(Maximize(C[1, 1]), list(constr))
  val <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 2, tolerance = 1e-2)

  ## Dual value: may be flat vector (n^2) or matrix; reshape if needed
  dv <- dual_value(constr)
  if (is.null(dim(dv))) {
    dv <- matrix(as.numeric(dv), 2, 2)
  }
  expect_equal(dim(dv), c(2, 2))
  ## Diagonal: [1, ~0], off-diagonal: ~0
  expect_equal(dv[1, 1], 1, tolerance = 0.1)
  expect_equal(dv[2, 2], 0, tolerance = 0.1)
  expect_equal(dv[1, 2], 0, tolerance = 0.1)
})

# -- 4. constraint_violation ---------------------------------------------------
# CVXPY: violation() returns elementwise residual for Equality/Inequality

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("constraint_violation: violation() on satisfied and violated constraints", {
  x <- Variable(2)

  ## Equality violation
  constr_eq <- (x == c(1, 2))
  value(x) <- c(1.5, 2.0)
  viol_eq <- violation(constr_eq)
  ## |1.5 - 1| = 0.5, |2.0 - 2| = 0
  expect_equal(as.numeric(viol_eq), c(0.5, 0), tolerance = 1e-10)

  ## Equality satisfied
  value(x) <- c(1, 2)
  viol_eq_ok <- violation(constr_eq)
  expect_equal(as.numeric(viol_eq_ok), c(0, 0), tolerance = 1e-10)

  ## Inequality violated: x=[1.5, 2.5] <= [1, 2]
  constr_ineq <- (x <= c(1, 2))
  value(x) <- c(1.5, 2.5)
  viol_ineq <- violation(constr_ineq)
  ## max(1.5 - 1, 0) = 0.5, max(2.5 - 2, 0) = 0.5
  expect_equal(as.numeric(viol_ineq), c(0.5, 0.5), tolerance = 1e-10)

  ## Inequality satisfied: x=[0.5, 1.5] <= [1, 2]
  value(x) <- c(0.5, 1.5)
  viol_ineq_ok <- violation(constr_ineq)
  expect_equal(as.numeric(viol_ineq_ok), c(0, 0), tolerance = 1e-10)

  ## PSD violation: X has negative eigenvalue
  X <- Variable(c(2, 2), symmetric = TRUE)
  constr_psd <- PSD(X)
  value(X) <- matrix(c(1, 0, 0, -0.5), 2, 2)
  ## violation = max(0, -min_eig) = max(0, 0.5) = 0.5
  viol_psd <- violation(constr_psd)
  expect_equal(as.numeric(viol_psd), 0.5, tolerance = 1e-10)

  ## PSD satisfied: X is PSD
  value(X) <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  viol_psd_ok <- violation(constr_psd)
  expect_equal(as.numeric(viol_psd_ok), 0, tolerance = 1e-10)
})


# =====================================================================
# CONE2CONE TESTS (4 tests)
# =====================================================================
# Verify cone constraints are preserved through the solving pipeline.

# -- 1. soc_to_soc -------------------------------------------------------------
# SOC constraint: norm(x) <= t preserved in solve

## @cvxpy test_constraints.py::TestConstraints::test_soc_constraint
test_that("soc_to_soc: SOC constraint preserved through solving pipeline", {
  skip_if_not_installed("scs")

  x <- Variable(2)
  t_var <- Variable()
  ## min t s.t. SOC(t, x), x == [3, 4]
  ## Expect t = norm(c(3,4)) = 5
  soc_con <- SOC(t_var, x, axis = 2L)
  prob <- Problem(Minimize(t_var), list(soc_con, x == c(3, 4)))
  val <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 5, tolerance = 1e-3)
  expect_equal(as.numeric(value(t_var)), 5, tolerance = 1e-3)
})

# -- 2. expcone_to_expcone ----------------------------------------------------
# ExpCone constraint: y * exp(x/y) <= z

## @cvxpy NONE
test_that("expcone_to_expcone: ExpCone constraint preserved through solving pipeline", {
  skip_if_not_installed("scs")

  xec <- Variable()
  yec <- Variable()
  zec <- Variable()
  ## min z s.t. ExpCone(x, y, z), x==1, y==1
  ## y * exp(x/y) <= z => 1*exp(1) <= z => z >= e
  ec_con <- ExpCone(xec, yec, zec)
  prob <- Problem(Minimize(zec), list(ec_con, xec == 1, yec == 1))
  val <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, exp(1), tolerance = 1e-3)
  expect_equal(as.numeric(value(zec)), exp(1), tolerance = 1e-3)
})

# -- 3. psd_to_psd ------------------------------------------------------------
# PSD constraint: X >> 0 preserved

## @cvxpy test_problem.py::TestProblem::test_psd_constraints
test_that("psd_to_psd: PSD constraint preserved through solving pipeline", {
  skip_if_not_installed("scs")

  X <- Variable(c(2, 2), symmetric = TRUE)
  ## min trace(X) s.t. X >> 0, X[1,1] >= 1, X[2,2] >= 1
  ## Optimal: X = I, trace = 2
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(X %>>% 0, X[1, 1] >= 1, X[2, 2] >= 1))
  val <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 2, tolerance = 1e-2)
  Xval <- value(X)
  expect_equal(Xval[1, 1], 1, tolerance = 1e-2)
  expect_equal(Xval[2, 2], 1, tolerance = 1e-2)
  expect_equal(Xval[1, 2], 0, tolerance = 1e-2)
})

# -- 4. powcone_handling -------------------------------------------------------
# Power cone constraint through pipeline via power atom

## @cvxpy NONE
test_that("powcone_handling: Power cone through solving pipeline (via power atom)", {
  skip_if_not_installed("clarabel")

  ## min power(x, 1.5) s.t. x >= 4
  ## Optimal: x = 4, value = 4^1.5 = 8
  x <- Variable()
  prob <- Problem(Minimize(power(x, 1.5)), list(x >= 4))
  val <- psolve(prob, solver = CLARABEL_SOLVER)

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 8, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 4, tolerance = 1e-3)
})


# =====================================================================
# QUAD FORM TESTS (4 tests)
# =====================================================================
# CVXPY SOURCE: test_quad_form.py, test_problem.py

# -- 1. quad_form_with_constant_x ----------------------------------------------
# CVXPY: quad_form(constant, P) evaluates to scalar
# c = [1, 2], P = [[4, 0], [0, 9]] => c'Pc = 4 + 36 = 40

## @cvxpy test_atoms.py::TestAtoms::test_quad_form
test_that("quad_form_with_constant_x: quad_form(constant, P) evaluates to scalar", {
  c_vec <- c(1, 2)
  P <- matrix(c(4, 0, 0, 9), 2, 2)

  ## quad_form with two constants => Constant expression
  expr <- quad_form(c_vec, P)
  expect_true(is_constant(expr))
  expect_equal(as.numeric(value(expr)), 40, tolerance = 1e-10)

  ## Solve: min quad_form([1,2], P) => 40 (no variables)
  prob <- Problem(Minimize(quad_form(c_vec, P)))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 40, tolerance = 1e-2)
})

# -- 2. quad_form_sparse_P ----------------------------------------------------
# CVXPY: quad_form with sparse Matrix P

## @cvxpy test_quad_form.py::TestNonOptimal::test_sparse_quad_form
test_that("quad_form_sparse_P: quad_form with sparse Matrix P", {
  skip_if_not_installed("osqp")

  P_sparse <- Matrix::Diagonal(2)  # 2x2 identity
  x <- Variable(2)
  cost <- quad_form(x, P_sparse)
  ## min x'Ix s.t. x == [1, 2] => 1 + 4 = 5
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  val <- psolve(prob, solver = "OSQP")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 5, tolerance = 1e-3)
})

# -- 3. quad_form_1x1 ---------------------------------------------------------
# quad_form with 1x1 matrix

## @cvxpy NONE
test_that("quad_form_1x1: quad_form with 1x1 matrix/scalar", {
  skip_if_not_installed("osqp")

  x <- Variable(1)
  P <- matrix(3, 1, 1)
  ## min 3*x^2 s.t. x == 2 => 3*4 = 12
  prob <- Problem(Minimize(quad_form(x, P)), list(x == 2))
  val <- psolve(prob, solver = "OSQP")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 12, tolerance = 1e-3)
})

# -- 4. quad_over_lin_equivalence -----------------------------------------------
# quad_over_lin(x, 1) == sum_squares(x)

## @cvxpy NONE
test_that("quad_over_lin_equivalence: quad_over_lin(x, 1) == sum_squares(x)", {
  skip_if_not_installed("clarabel")

  x1 <- Variable(3)
  prob_qol <- Problem(Minimize(quad_over_lin(x1, 1)), list(x1 == c(1, 2, 3)))
  val_qol <- psolve(prob_qol, solver = CLARABEL_SOLVER)

  x2 <- Variable(3)
  prob_ss <- Problem(Minimize(sum_squares(x2)), list(x2 == c(1, 2, 3)))
  val_ss <- psolve(prob_ss, solver = CLARABEL_SOLVER)

  ## Both should equal 1^2 + 2^2 + 3^2 = 14
  expect_equal(val_qol, 14, tolerance = 1e-3)
  expect_equal(val_ss, 14, tolerance = 1e-3)
  expect_equal(val_qol, val_ss, tolerance = 1e-3)
})


# =====================================================================
# QUADRATIC TESTS (4 tests)
# =====================================================================
# CVXPY SOURCE: test_quadratic.py

# -- 1. sum_squares_solving ----------------------------------------------------
# min sum_squares(x - [1,2,3]) => x = [1,2,3], value = 0

## @cvxpy test_quadratic.py::TestExpressions::test_sum_squares
test_that("sum_squares_solving: sum_squares(x) minimization (CVXPY parity)", {
  skip_if_not_installed("osqp")

  x <- Variable(3)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2, 3))))
  val <- psolve(prob, solver = "OSQP")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 2, 3), tolerance = 1e-3)
})

# -- 2. matrix_frac_solving ----------------------------------------------------
# min matrix_frac(x, P) s.t. x == [1, 1]
# P = [[2,1],[1,3]], P^{-1} = (1/5)*[[3,-1],[-1,2]]
# x'P^{-1}x = (1/5)*(3 - 1 - 1 + 2) = 3/5 = 0.6

## @cvxpy test_quadratic.py::TestExpressions::test_matrix_frac
test_that("matrix_frac_solving: matrix_frac(x, P) with PSD P (CVXPY parity)", {
  skip_if_not_installed("scs")

  x <- Variable(2)
  P <- matrix(c(2, 1, 1, 3), 2, 2)  # col-major
  prob <- Problem(Minimize(matrix_frac(x, P)), list(x == c(1, 1)))
  val <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 0.6, tolerance = 1e-2)

  ## Also test actual optimization:
  ## min x'P^{-1}x s.t. sum(x) >= 2, P = 2*I
  ## P^{-1} = 0.5*I => min 0.5*||x||^2 s.t. 1'x >= 2
  ## Optimal: x = [1,1], val = 1
  x2 <- Variable(2)
  P2 <- 2 * diag(2)
  prob2 <- Problem(Minimize(matrix_frac(x2, P2)), list(sum_entries(x2) >= 2))
  val2 <- psolve(prob2, solver = "SCS")

  expect_equal(status(prob2), OPTIMAL)
  expect_equal(val2, 1, tolerance = 1e-2)
  expect_equal(as.numeric(value(x2)), c(1, 1), tolerance = 1e-2)
})

# -- 3. quad_over_lin_solving --------------------------------------------------
# min quad_over_lin(x, y) s.t. sum(x) >= 4, y <= 2
# Optimal: x = [2, 2], y = 2, value = 8/2 = 4

## @cvxpy test_quadratic.py::TestExpressions::test_quad_over_lin
test_that("quad_over_lin_solving: quad_over_lin(x, y) minimization (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  x <- Variable(2)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(quad_over_lin(x, y)),
                  list(sum_entries(x) >= 4, y <= 2))
  val <- psolve(prob, solver = CLARABEL_SOLVER)

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 4, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(2, 2), tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 2, tolerance = 1e-2)
})

# -- 4. is_quadratic_detection ------------------------------------------------
# CVXPY: is_quadratic() for various expressions

## @cvxpy test_quadratic.py::TestExpressions::test_quadratic_form
test_that("is_quadratic_detection: is_quadratic() returns TRUE for QP objectives", {
  x <- Variable(2)

  ## quad_form(x, I) is quadratic
  q <- quad_form(x, diag(2))
  expect_true(is_quadratic(q))
  expect_true(is_dcp(Minimize(q)))

  ## sum_squares(x) is quadratic
  ss <- sum_squares(x)
  expect_true(is_quadratic(ss))

  ## quad_over_lin(x, 1) is quadratic (constant denominator)
  qol <- quad_over_lin(x, 1)
  expect_true(is_quadratic(qol))

  ## matrix_frac(x, I) is quadratic
  mf <- matrix_frac(x, diag(2))
  expect_true(is_quadratic(mf))

  ## Variable is quadratic (and affine)
  expect_true(is_quadratic(x))
  expect_true(is_affine(x))

  ## power(x, 2) is quadratic
  x_scalar <- Variable()
  pw2 <- power(x_scalar, 2)
  expect_true(is_quadratic(pw2))

  ## Affine expression is quadratic
  aff <- 2 * x + c(3, 4)
  expect_true(is_quadratic(aff))
  expect_true(is_affine(aff))
})


# =====================================================================
# PERSPECTIVE TESTS (3 tests)
# =====================================================================
# CVXPY SOURCE: test_perspective.py

# -- 1. perspective_exp --------------------------------------------------------
# CVXPY: perspective(exp(x), s) with s >= 1, x >= 1
# Optimal: x=1, s=1, value = s*exp(x/s) = exp(1) = 2.71828...

## @cvxpy test_perspective.py::test_exp
test_that("perspective_exp: perspective(exp(x), s) with known solution (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- exp(x)
  obj <- Minimize(perspective(f, s))
  prob <- Problem(obj, list(s >= 1, x >= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER)

  expect_equal(status(prob), OPTIMAL)
  ## Optimal value = exp(1)
  expect_equal(val, exp(1), tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1, tolerance = 1e-2)
  expect_equal(as.numeric(value(s)), 1, tolerance = 1e-2)
})

# -- 2. perspective_square -----------------------------------------------------
# CVXPY: test_quad_atom with r=2
# f = x^2 + 2*x - 4, perspective(f, s) with s <= 0.5, x >= 2
# Reference: quad_over_lin(x, s) + 2*x - 4*s
# CVXPY value: 10

## @cvxpy test_perspective.py::test_quad_atom
test_that("perspective_square: perspective(square(x)+2x-4, s) (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  ## Reference problem: quad_over_lin + affine terms
  ref_x <- Variable()
  ref_s <- Variable(nonneg = TRUE)
  ref_obj <- quad_over_lin(ref_x, ref_s) + 2 * ref_x - 4 * ref_s
  ref_prob <- Problem(Minimize(ref_obj), list(ref_x >= 2, ref_s <= 0.5))
  ref_val <- psolve(ref_prob, solver = CLARABEL_SOLVER)

  ## Perspective problem: perspective(x^2 + 2x - 4, s)
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- square(x) + 2 * x - 4
  obj <- Minimize(perspective(f, s))
  prob <- Problem(obj, list(s <= 0.5, x >= 2))
  val <- psolve(prob, solver = CLARABEL_SOLVER)

  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 10, tolerance = 0.1)
  expect_equal(val, ref_val, tolerance = 0.1)
  expect_equal(as.numeric(value(x)), as.numeric(value(ref_x)), tolerance = 0.1)
  expect_equal(as.numeric(value(s)), as.numeric(value(ref_s)), tolerance = 0.1)
})

# -- 3. perspective_parameter ---------------------------------------------------
# CVXPY: test_parameter
# f = p * x^2, perspective(f, s), s <= 1, x >= 2, p = 99
# Value: 4 * p = 396

## @cvxpy test_perspective.py::test_parameter
test_that("perspective_parameter: perspective with Parameter in f (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  p <- Parameter(nonneg = TRUE)
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- p * square(x)
  obj <- Minimize(perspective(f, s))
  prob <- Problem(obj, list(s <= 1, x >= 2))
  value(p) <- 99

  val <- psolve(prob, solver = CLARABEL_SOLVER)

  expect_equal(status(prob), OPTIMAL)
  ## Value = x^2/s * p at optimum: x=2, s=1 => 4*99 = 396
  expect_equal(val, 396, tolerance = 1)
})


# =====================================================================
# POWER ATOM TESTS (7 tests)
# =====================================================================
# CVXPY SOURCE: test_atoms.py::test_power, test_expressions.py::test_powers,
#               test_problem.py::test_power, test_quadratic.py::test_power

# -- 1. power_2 ---------------------------------------------------------------
# power(x, 2) == square(x)

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power_2: power(x, 2) == square(x) (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  x1 <- Variable()
  prob1 <- Problem(Minimize(power(x1, 2)), list(x1 >= 3))
  val1 <- psolve(prob1, solver = CLARABEL_SOLVER)

  x2 <- Variable()
  prob2 <- Problem(Minimize(square(x2)), list(x2 >= 3))
  val2 <- psolve(prob2, solver = CLARABEL_SOLVER)

  ## Both should give 9
  expect_equal(val1, 9, tolerance = 1e-3)
  expect_equal(val2, 9, tolerance = 1e-3)
  expect_equal(val1, val2, tolerance = 1e-3)
})

# -- 2. power_half -------------------------------------------------------------
# power(x, 0.5) == sqrt(x)

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power_half: power(x, 0.5) == sqrt(x) (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  x1 <- Variable()
  prob1 <- Problem(Maximize(power(x1, 0.5)), list(x1 <= 4, x1 >= 0))
  val1 <- psolve(prob1, solver = CLARABEL_SOLVER)

  x2 <- Variable()
  prob2 <- Problem(Maximize(sqrt(x2)), list(x2 <= 4, x2 >= 0))
  val2 <- psolve(prob2, solver = CLARABEL_SOLVER)

  ## Both should give sqrt(4) = 2
  expect_equal(val1, 2, tolerance = 1e-3)
  expect_equal(val2, 2, tolerance = 1e-3)
  expect_equal(val1, val2, tolerance = 1e-3)
})

# -- 3. power_neg1 -------------------------------------------------------------
# power(x, -1) == inv_pos(x)

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power_neg1: power(x, -1) == inv_pos(x) (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  x1 <- Variable()
  prob1 <- Problem(Minimize(power(x1, -1)), list(x1 >= 2, x1 <= 5))
  val1 <- psolve(prob1, solver = CLARABEL_SOLVER)

  x2 <- Variable()
  prob2 <- Problem(Minimize(inv_pos(x2)), list(x2 >= 2, x2 <= 5))
  val2 <- psolve(prob2, solver = CLARABEL_SOLVER)

  ## Both should give 1/5 = 0.2
  expect_equal(val1, 0.2, tolerance = 1e-3)
  expect_equal(val2, 0.2, tolerance = 1e-3)
  expect_equal(val1, val2, tolerance = 1e-3)
})

# -- 4. power_frac -------------------------------------------------------------
# power(x, p/q) for fractional exponents

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power_frac: power(x, p/q) for fractional exponents (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  ## power(x, 2/3): concave => max power(x, 2/3) s.t. x <= 8
  ## Optimal: x = 8, value = 8^(2/3) = 4
  x1 <- Variable()
  prob1 <- Problem(Maximize(power(x1, 2 / 3)), list(x1 <= 8, x1 >= 0))
  val1 <- psolve(prob1, solver = CLARABEL_SOLVER)

  expect_equal(status(prob1), OPTIMAL)
  expect_equal(val1, 4, tolerance = 1e-2)
  expect_equal(as.numeric(value(x1)), 8, tolerance = 1e-2)

  ## power(x, 3/2): convex => min power(x, 3/2) s.t. x >= 4
  ## Optimal: x = 4, value = 4^(3/2) = 8
  x2 <- Variable()
  prob2 <- Problem(Minimize(power(x2, 3 / 2)), list(x2 >= 4))
  val2 <- psolve(prob2, solver = CLARABEL_SOLVER)

  expect_equal(status(prob2), OPTIMAL)
  expect_equal(val2, 8, tolerance = 1e-2)
  expect_equal(as.numeric(value(x2)), 4, tolerance = 1e-2)
})

# -- 5. power_curvature -------------------------------------------------------
# DCP curvature of power for various exponents

## @cvxpy test_expressions.py::TestExpressions::test_powers
test_that("power_curvature: DCP curvature of power for various exponents (CVXPY parity)", {
  x <- Variable(c(2, 1))
  y <- Variable(c(2, 1))
  expr <- x + y

  ## p > 1 => CONVEX
  expect_equal(curvature(power(expr, 2)), "CONVEX")
  expect_equal(curvature(power(expr, 3)), "CONVEX")
  expect_equal(curvature(power(expr, 2.7)), "CONVEX")

  ## p < 0 => CONVEX
  expect_equal(curvature(power(expr, -1)), "CONVEX")
  expect_equal(curvature(power(expr, -2.3)), "CONVEX")

  ## p == 1 => AFFINE
  expect_equal(curvature(power(expr, 1)), "AFFINE")

  ## p == 0 => CONSTANT
  expect_true(is_constant(power(expr, 0)))

  ## 0 < p < 1 => CONCAVE
  expect_equal(curvature(power(expr, 0.67)), "CONCAVE")
  expect_equal(curvature(power(expr, 0.5)), "CONCAVE")

  ## Nonneg sign for all p != 1
  for (p in c(0, 2, 3, 2.7, 0.67, 0.5, -1, -2.3)) {
    atom <- power(expr, p)
    expect_true(is_nonneg(atom), info = paste("p =", p))
  }
})

# -- 6. power_mono -------------------------------------------------------------
# Monotonicity of power for various exponents

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power_mono: Monotonicity of power for various exponents (CVXPY parity)", {
  ## is_incr/is_decr are internal generics with 0-based idx
  is_incr <- CVXR:::is_incr
  is_decr <- CVXR:::is_decr

  x <- Variable(nonneg = TRUE)

  ## p == 0 => constant, increasing and decreasing both TRUE (trivially)
  atom0 <- power(x, 0)
  expect_true(is_incr(atom0, 1L))
  expect_true(is_decr(atom0, 1L))

  ## 0 < p < 1 => increasing
  atom_half <- power(x, 0.5)
  expect_true(is_incr(atom_half, 1L))
  expect_false(is_decr(atom_half, 1L))

  ## p == 1 => increasing
  atom1 <- power(x, 1)
  expect_true(is_incr(atom1, 1L))
  expect_false(is_decr(atom1, 1L))

  ## p > 1 => increasing (for nonneg domain)
  atom2 <- power(x, 2)
  expect_true(is_incr(atom2, 1L))
  expect_false(is_decr(atom2, 1L))

  atom3 <- power(x, 3)
  expect_true(is_incr(atom3, 1L))
  expect_false(is_decr(atom3, 1L))

  ## p < 0 => decreasing
  atom_neg1 <- power(x, -1)
  expect_false(is_incr(atom_neg1, 1L))
  expect_true(is_decr(atom_neg1, 1L))
})

# -- 7. power_numeric ---------------------------------------------------------
# Numeric evaluation of power atom on constants

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power_numeric: Numeric evaluation of power atom (CVXPY parity)", {
  ## power(-1, 2) = 1 (CVXPY: assert cp.power(-1, 2).value == 1)
  expect_equal(as.numeric(value(power(Constant(-1), 2))), 1, tolerance = 1e-10)

  ## power(4, 0.5) = 2
  expect_equal(as.numeric(value(power(Constant(4), 0.5))), 2, tolerance = 1e-10)

  ## power(2, -1) = 0.5
  expect_equal(as.numeric(value(power(Constant(2), -1))), 0.5, tolerance = 1e-10)

  ## power(8, 1/3) = 2
  expect_equal(as.numeric(value(power(Constant(8), 1 / 3))), 2, tolerance = 1e-4)

  ## power(3, 0) = 1
  expect_equal(as.numeric(value(power(Constant(3), 0))), 1, tolerance = 1e-10)

  ## power(27, 2/3) = 9
  expect_equal(as.numeric(value(power(Constant(27), 2 / 3))), 9, tolerance = 1e-4)

  ## Vector evaluation: power([1,4,9], 0.5) = [1,2,3]
  v <- Constant(matrix(c(1, 4, 9), 3, 1))
  result <- as.numeric(value(power(v, 0.5)))
  expect_equal(result, c(1, 2, 3), tolerance = 1e-10)
})


# =====================================================================
# SEMIDEFINITE VARIABLE TESTS (2 tests)
# =====================================================================
# CVXPY SOURCE: test_semidefinite_vars.py

# -- 1. test_symm -----------------------------------------------------------
# CVXPY: M PSD, x1 PSD, x2 PSD; M + C1 == x1, M + C2 == x2;
# minimize trace(M). Verify M is symmetric.

## @cvxpy test_semidefinite_vars.py::TestSemidefiniteVariable::test_symm
test_that("test_symm: PSD variable M is symmetric after solve (CVXPY parity)", {
  skip_if_not_installed("scs")

  M <- Variable(c(3, 3), PSD = TRUE)
  C1 <- matrix(c(0, 0, 0.5, 0, 0, 0, 0.5, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  C2 <- matrix(c(0, 0, 0, 0, 0, 0.5, 0, 0.5, 1), nrow = 3, ncol = 3, byrow = TRUE)
  x1 <- Variable(c(3, 3), PSD = TRUE)
  x2 <- Variable(c(3, 3), PSD = TRUE)
  constraints <- list(M + C1 == x1, M + C2 == x2)
  objective <- Minimize(matrix_trace(M))
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  Mval <- value(M)
  ## Verify M is symmetric: M == t(M)
  expect_equal(Mval, t(Mval), tolerance = 1e-4)
})

# -- 2. test_sdp_problem --------------------------------------------------------
# CVXPY: Three sub-problems testing PSD variable behavior

## @cvxpy test_semidefinite_vars.py::TestSemidefiniteVariable::test_sdp_problem
test_that("test_sdp_problem: PSD in objective — min sum(square(X - F)) (CVXPY parity)", {
  skip_if_not_installed("scs")

  ## Part 1: PSD in objective
  ## X is PSD, F = [[1,0],[0,-1]]. Nearest PSD to F is [[1,0],[0,0]].
  ## sum(square(X - F)) = (0-(-1))^2 = 1
  X <- Variable(c(2, 2), PSD = TRUE)
  F_mat <- matrix(c(1, 0, 0, -1), 2, 2)
  obj <- Minimize(sum_entries(square(X - F_mat)))
  prob <- Problem(obj, list())
  result <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(result, 1, tolerance = 1e-2)
  ## X should be approximately [[1,0],[0,0]]
  Xval <- value(X)
  expect_equal(Xval[1, 1], 1, tolerance = 0.1)
  expect_equal(Xval[2, 2], 0, tolerance = 0.1)
})

## @cvxpy test_semidefinite_vars.py::TestSemidefiniteVariable::test_sdp_problem
test_that("test_sdp_problem: PSD in constraint (CVXPY parity)", {
  skip_if_not_installed("scs")

  ## Part 2: PSD in constraint
  ## Y unconstrained, but Y == Z where Z is PSD
  ## F = [[1,0],[0,-1]]. Closest PSD to F is [[1,0],[0,0]].
  Y <- Variable(c(2, 2))
  Z <- Variable(c(2, 2), PSD = TRUE)
  F_mat <- matrix(c(1, 0, 0, -1), 2, 2)
  obj <- Minimize(sum_entries(square(Y - F_mat)))
  prob <- Problem(obj, list(Y == Z))
  result <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(result, 1, tolerance = 0.1)
})

## @cvxpy test_semidefinite_vars.py::TestSemidefiniteVariable::test_sdp_problem
test_that("test_sdp_problem: index into PSD variable (CVXPY parity)", {
  skip_if_not_installed("scs")

  ## Part 3: Index into semidef
  ## min (X[1,1]-1)^2 + (X[2,1]-2)^2 + (X[2,2]-4)^2
  ## X is PSD, so X[1,2] == X[2,1]. Optimal: X = [[1,2],[2,4]] (rank 1, PSD)
  ## Note: CVXPY uses 0-based: X[0,0], X[1,0], X[1,1]
  ## R uses 1-based: X[1,1], X[2,1], X[2,2]
  X <- Variable(c(2, 2), PSD = TRUE)
  obj <- Minimize(
    square(X[1, 1] - 1) + square(X[2, 1] - 2) + square(X[2, 2] - 4)
  )
  prob <- Problem(obj, list())
  result <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(result, 0, tolerance = 1e-2)
  ## X should be [[1,2],[2,4]]
  Xval <- value(X)
  expect_equal(Xval[1, 1], 1, tolerance = 0.1)
  expect_equal(Xval[2, 1], 2, tolerance = 0.1)
  expect_equal(Xval[1, 2], 2, tolerance = 0.1)
  expect_equal(Xval[2, 2], 4, tolerance = 0.1)
})
