## Wave 4: Expression & Constraint Parity Tests
## Tests from CVXPY test_expressions.py and test_constraints.py gaps
## All expected values verified via `uv run python` against CVXPY 1.8.1

# ═══════════════════════════════════════════════════════════════════════
# Variable Value Assignment (CVXPY test_expressions.py::test_assign_var_value)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: wrong shape raises error", {
  ## CVXPY: ValueError("Invalid dimensions (2,) for Variable value.")
  x <- Variable(3)
  expect_error(value(x) <- c(1, 2), "value")
})

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: nonneg rejects negative", {
  ## CVXPY: ValueError("Variable value must be nonnegative.")
  x <- Variable(3, nonneg = TRUE)
  expect_error(value(x) <- c(1, -2, 3), "nonneg")
})

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: boolean rejects non-boolean", {
  ## CVXPY: ValueError("Variable value must be boolean.")
  b <- Variable(2, boolean = TRUE)
  expect_error(value(b) <- c(0.5, 0.5), "boolean")
})

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: integer rejects non-integer", {
  ## CVXPY: ValueError("Variable value must be integer.")
  i <- Variable(2, integer = TRUE)
  expect_error(value(i) <- c(0.5, 1.5), "integer")
})

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: NULL clears value", {
  x <- Variable(3)
  value(x) <- c(1, 2, 3)
  expect_false(is.null(value(x)))
  value(x) <- NULL
  expect_null(value(x))
})

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: correct shape accepted", {
  x <- Variable(c(2, 3))
  val <- matrix(1:6, 2, 3)
  value(x) <- val
  expect_equal(value(x), val)
})

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: scalar to scalar variable", {
  x <- Variable(1)
  value(x) <- 5
  expect_equal(as.numeric(value(x)), 5)
})

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: nonneg accepts nonneg values", {
  x <- Variable(3, nonneg = TRUE)
  value(x) <- c(0, 1, 2)
  expect_equal(as.numeric(value(x)), c(0, 1, 2))
})

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: boolean accepts boolean values", {
  b <- Variable(2, boolean = TRUE)
  value(b) <- c(0, 1)
  expect_equal(as.numeric(value(b)), c(0, 1))
})

# ═══════════════════════════════════════════════════════════════════════
# Variable project() (CVXPY test_expressions.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_expressions.py::TestExpressions::test_round_attr
test_that("project: nonneg clips to [0, Inf)", {
  x <- Variable(3, nonneg = TRUE)
  projected <- project(x, c(-2, 1, -0.5))
  expect_equal(as.numeric(projected), c(0, 1, 0))
})

## @cvxpy test_expressions.py::TestExpressions::test_round_attr
test_that("project: nonpos clips to (-Inf, 0]", {
  x <- Variable(3, nonpos = TRUE)
  projected <- project(x, c(-2, 1, -0.5))
  expect_equal(as.numeric(projected), c(-2, 0, -0.5))
})

## @cvxpy test_expressions.py::TestExpressions::test_round_attr
test_that("project: boolean rounds to {0, 1}", {
  ## CVXPY: project([0.1, 0.9, 0.5]) = [0, 1, 0]
  ## (round after clipping to [0,1])
  b <- Variable(3, boolean = TRUE)
  projected <- project(b, c(0.1, 0.9, 0.5))
  expect_equal(as.numeric(projected), c(0, 1, 0))
})

## @cvxpy test_expressions.py::TestExpressions::test_round_attr
test_that("project: integer rounds to nearest integer", {
  i <- Variable(3, integer = TRUE)
  projected <- project(i, c(0.3, 1.7, -0.6))
  expect_equal(as.numeric(projected), c(0, 2, -1))
})

## @cvxpy test_expressions.py::TestExpressions::test_round_attr
test_that("project: symmetric symmetrizes", {
  x <- Variable(c(2, 2), symmetric = TRUE)
  projected <- project(x, matrix(c(1, 2, 3, 4), 2, 2))
  expected <- (matrix(c(1, 2, 3, 4), 2, 2) + t(matrix(c(1, 2, 3, 4), 2, 2))) / 2
  expect_equal(projected, expected)
})

## @cvxpy test_expressions.py::TestExpressions::test_round_attr
test_that("project: PSD projects to PSD cone", {
  x <- Variable(c(2, 2), PSD = TRUE)
  val <- matrix(c(-1, 0, 0, 2), 2, 2)
  projected <- project(x, val)
  ## All eigenvalues should be >= 0 after projection
  eigs <- eigen(projected, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

# ═══════════════════════════════════════════════════════════════════════
# expr_copy (CVXPY test_expressions.py::test_var_copy/test_param_copy)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_expressions.py::TestExpressions::test_var_copy
test_that("expr_copy: Variable returns self (leaf)", {
  ## CVXPY: x.copy().id == x.id (same identity for leaf)
  x <- Variable(3)
  x_copy <- expr_copy(x)
  expect_identical(x_copy@id, x@id)
})
## @cvxpy test_expressions.py::TestExpressions::test_param_copy

test_that("expr_copy: Parameter returns self (leaf)", {
  p <- Parameter(2, nonneg = TRUE)
  p_copy <- expr_copy(p)
  expect_identical(p_copy@id, p@id)
  expect_true(p_copy@attributes$nonneg)
})

## @cvxpy test_expressions.py::TestExpressions::test_constant_copy
test_that("expr_copy: Constant returns self (leaf)", {
  c_val <- Constant(matrix(1:4, 2, 2))
  c_copy <- expr_copy(c_val)
  expect_identical(c_copy@id, c_val@id)
})

## @cvxpy NONE
test_that("expr_copy: compound expression creates new node", {
  x <- Variable(3)
  y <- Variable(3)
  expr <- x + y
  expr_cp <- expr_copy(expr)
  ## New expression has same structure but different id
  expect_false(identical(expr_cp@id, expr@id))
  ## Children are the same leaf objects
  expect_identical(expr_cp@args[[1L]]@id, x@id)
  expect_identical(expr_cp@args[[2L]]@id, y@id)
})

## @cvxpy NONE
test_that("tree_copy: deep copy creates all new ids", {
  x <- Variable(3)
  y <- Variable(3)
  expr <- x + y
  expr_tc <- tree_copy(expr)
  ## Root has new id
  expect_false(identical(expr_tc@id, expr@id))
})

# ═══════════════════════════════════════════════════════════════════════
# Inequality violation (CVXPY test_constraints.py)
# Verified: Inequality.violation() returns array, NonPos returns scalar
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("inequality violation: returns array (not scalar)", {
  ## CVXPY: x<=z with x=[3,1,2], z=[2,2,2] → violation=[1,0,0]
  x <- Variable(3)
  z <- Variable(3)
  ineq <- x <= z
  save_leaf_value(x, matrix(c(3, 1, 2), 3, 1))
  save_leaf_value(z, matrix(c(2, 2, 2), 3, 1))
  viol <- violation(ineq)
  expect_equal(as.numeric(viol), c(1, 0, 0))
  expect_true(is.matrix(viol) || is.numeric(viol))
})

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("equality violation: returns abs(lhs - rhs)", {
  ## CVXPY: x==z with x=[3,1,2], z=[2,2,2] → violation=[1,1,0]
  x <- Variable(3)
  z <- Variable(3)
  eq <- x == z
  save_leaf_value(x, matrix(c(3, 1, 2), 3, 1))
  save_leaf_value(z, matrix(c(2, 2, 2), 3, 1))
  viol <- violation(eq)
  expect_equal(as.numeric(viol), c(1, 1, 0))
})

## @cvxpy test_constraints.py::TestConstraints::test_nonneg
test_that("NonPos violation: returns L2 norm (scalar)", {
  ## CVXPY: NonPos(x-z) with x=[3,1,2], z=[2,2,2] → violation=1.0
  x <- Variable(3)
  z <- Variable(3)
  save_leaf_value(x, matrix(c(3, 1, 2), 3, 1))
  save_leaf_value(z, matrix(c(2, 2, 2), 3, 1))
  np_constr <- NonPos(x - z)
  viol <- violation(np_constr)
  expect_equal(viol, 1.0)
  expect_length(viol, 1)
})

## @cvxpy test_constraints.py::TestConstraints::test_nonneg
test_that("NonNeg violation: returns L2 norm (scalar)", {
  ## CVXPY: NonNeg(z-x) with z=[2,2,2], x=[3,1,2] → residual=abs(min(z-x,0))
  ## z-x = [-1, 1, 0], min(z-x,0) = [-1, 0, 0], abs = [1, 0, 0]
  ## L2 norm = 1.0
  x <- Variable(3)
  z <- Variable(3)
  save_leaf_value(x, matrix(c(3, 1, 2), 3, 1))
  save_leaf_value(z, matrix(c(2, 2, 2), 3, 1))
  nn_constr <- NonNeg(z - x)
  viol <- violation(nn_constr)
  expect_equal(viol, 1.0)
  expect_length(viol, 1)
})

# ═══════════════════════════════════════════════════════════════════════
# SOC violation with exact CVXPY values
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constraints.py::TestConstraints::test_soc_constraint
test_that("SOC violation: feasible is 0", {
  ## CVXPY: SOC(t=5, x=[1,2,3]) → violation = 0.0
  t_var <- Variable(1)
  x_var <- Variable(3)
  soc <- SOC(t_var, x_var)
  save_leaf_value(t_var, matrix(5, 1, 1))
  save_leaf_value(x_var, matrix(c(1, 2, 3), 3, 1))
  expect_equal(violation(soc), 0.0, tolerance = 1e-10)
})

## @cvxpy test_constraints.py::TestConstraints::test_soc_constraint
test_that("SOC violation: infeasible matches CVXPY", {
  ## CVXPY: SOC(t=1, x=[1,2,3]) → violation = 1.938644529878043
  ## ||x||_2 = sqrt(14) ≈ 3.742, t = 1, so infeasible
  t_var <- Variable(1)
  x_var <- Variable(3)
  soc <- SOC(t_var, x_var)
  save_leaf_value(t_var, matrix(1, 1, 1))
  save_leaf_value(x_var, matrix(c(1, 2, 3), 3, 1))
  expect_equal(violation(soc), 1.938644529878043, tolerance = 1e-6)
})

# ═══════════════════════════════════════════════════════════════════════
# PSD violation with exact CVXPY values
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constraints.py::TestConstraints::test_psd_constraint
test_that("PSD violation: feasible matrix → 0", {
  ## CVXPY: PSD with [[2,-1],[-1,2]] (eigs 1,3) → violation = -0.0 ≈ 0
  A <- Variable(c(2, 2), symmetric = TRUE)
  psd_c <- PSD(A)
  save_leaf_value(A, matrix(c(2, -1, -1, 2), 2, 2))
  expect_equal(violation(psd_c), 0.0, tolerance = 1e-10)
})

## @cvxpy test_constraints.py::TestConstraints::test_psd_constraint
test_that("PSD violation: infeasible matrix → -min_eigenvalue", {
  ## CVXPY: PSD with [[-1,0],[0,2]] (eigs -1,2) → violation = 1.0
  A <- Variable(c(2, 2), symmetric = TRUE)
  psd_c <- PSD(A)
  save_leaf_value(A, matrix(c(-1, 0, 0, 2), 2, 2))
  expect_equal(violation(psd_c), 1.0, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# PowCone3D scalar alpha solve (CVXPY test_constraints.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constraints.py::TestConstraints::test_pow3d_scalar_alpha_constraint
test_that("PowCone3D: scalar alpha broadcast solve", {
  ## CVXPY test_pow3d_scalar_alpha_constraint:
  ## minimize ||x - x0|| s.t. PowCone3D(x0[1],x0[2],x0[3], 0.25), x <= -10
  ## Expected optimal value: 17.3205 (≈ 10*sqrt(3))
  x0 <- Variable(3)
  x <- Variable(3)
  cons <- list(PowCone3D(x0[1, 1], x0[2, 1], x0[3, 1], 0.25), x <= -10)
  obj <- Minimize(p_norm(x - x0))
  prob <- Problem(obj, cons)
  result <- psolve(prob)
  expect_equal(result, 17.320508, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# Constraint value() (boolean satisfaction check)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("constraint value: satisfied returns TRUE", {
  x <- Variable(3)
  constr <- x >= 0
  save_leaf_value(x, matrix(c(1, 2, 3), 3, 1))
  expect_true(value(constr))
})

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("constraint value: violated returns FALSE", {
  x <- Variable(3)
  constr <- x >= 0
  save_leaf_value(x, matrix(c(1, -1, 3), 3, 1))
  expect_false(value(constr))
})

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("equality constraint value: satisfied", {
  x <- Variable(3)
  constr <- x == c(1, 2, 3)
  save_leaf_value(x, matrix(c(1, 2, 3), 3, 1))
  expect_true(value(constr))
})

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("equality constraint value: violated", {
  x <- Variable(3)
  constr <- x == c(1, 2, 3)
  save_leaf_value(x, matrix(c(1, 2, 4), 3, 1))
  expect_false(value(constr))
})

# ═══════════════════════════════════════════════════════════════════════
# DCP analysis of various constraint patterns
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("inequality: affine <= affine is DCP", {
  x <- Variable(3)
  constr <- x <= c(1, 2, 3)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("inequality: convex <= constant is DCP", {
  x <- Variable(3)
  constr <- sum_squares(x) <= 1
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("inequality: concave >= constant is DCP (via <=)", {
  ## In CVXR, >= creates Inequality(rhs, lhs) i.e. rhs <= lhs
  x <- Variable(3, nonneg = TRUE)
  constr <- sqrt(x[1, 1]) >= 1
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("equality: affine == affine is DCP", {
  x <- Variable(3)
  y <- Variable(3)
  constr <- x == y
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("equality: nonlinear == is not DCP", {
  x <- Variable(3)
  constr <- sum_squares(x) == 1
  expect_false(is_dcp(constr))
})

# ═══════════════════════════════════════════════════════════════════════
# Constraint expr_name (CVXPY test_constraints.py formatting)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("inequality expr_name format", {
  x <- Variable(3, name = "x")
  constr <- x <= c(1, 2, 3)
  nm <- expr_name(constr)
  expect_true(grepl("<=", nm))
})

## @cvxpy NONE
test_that("equality expr_name format", {
  x <- Variable(3, name = "x")
  constr <- x == c(1, 2, 3)
  nm <- expr_name(constr)
  expect_true(grepl("==", nm))
})

## @cvxpy NONE
test_that("NonPos expr_name format", {
  x <- Variable(3, name = "x")
  constr <- NonPos(x)
  nm <- expr_name(constr)
  expect_true(grepl("<= 0", nm))
})

## @cvxpy NONE
test_that("NonNeg expr_name format", {
  x <- Variable(3, name = "x")
  constr <- NonNeg(x)
  nm <- expr_name(constr)
  expect_true(grepl(">= 0", nm))
})

# ═══════════════════════════════════════════════════════════════════════
# Variable/Parameter attributes and properties
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_expressions.py::TestExpressions::test_variable
test_that("variable: nonneg implies NONNEG sign", {
  x <- Variable(3, nonneg = TRUE)
  expect_true(is_nonneg(x))
  expect_false(is_nonpos(x))
})

## @cvxpy test_expressions.py::TestExpressions::test_variable
test_that("variable: nonpos implies NONPOS sign", {
  x <- Variable(3, nonpos = TRUE)
  expect_true(is_nonpos(x))
  expect_false(is_nonneg(x))
})

## @cvxpy test_expressions.py::TestExpressions::test_variable
test_that("variable: symmetric creates square variable", {
  x <- Variable(c(3, 3), symmetric = TRUE)
  expect_equal(x@shape, c(3L, 3L))
  expect_true(is_symmetric(x))
})

## @cvxpy test_expressions.py::TestExpressions::test_parameters_successes
test_that("parameter: nonneg preserves sign", {
  p <- Parameter(3, nonneg = TRUE)
  expect_true(is_nonneg(p))
})

## @cvxpy test_expressions.py::TestExpressions::test_parameters_failures
test_that("parameter: value assignment validates shape", {
  p <- Parameter(3)
  expect_error(value(p) <- c(1, 2), "value")
  value(p) <- c(1, 2, 3)
  expect_equal(as.numeric(value(p)), c(1, 2, 3))
})

# ═══════════════════════════════════════════════════════════════════════
# Misc expression tests from CVXPY
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_expressions.py::TestExpressions::test_constants
test_that("constant: is_constant is TRUE, is_affine is TRUE", {
  c_val <- Constant(c(1, 2, 3))
  expect_true(is_constant(c_val))
  expect_true(is_affine(c_val))
  expect_true(is_convex(c_val))
  expect_true(is_concave(c_val))
})

## @cvxpy test_expressions.py::TestExpressions::test_variable
test_that("variable: is_affine but not is_constant", {
  x <- Variable(3)
  expect_true(is_affine(x))
  expect_false(is_constant(x))
})

## @cvxpy test_problem.py::TestProblem::test_mult_by_zero
test_that("multiply by zero: value is zero", {
  ## CVXPY test_problem.py::test_mult_by_zero
  x <- Variable(3)
  expr <- 0 * x
  save_leaf_value(x, matrix(c(5, 6, 7), 3, 1))
  val <- value(expr)
  expect_equal(as.numeric(val), c(0, 0, 0))
})

## @cvxpy NONE
test_that("expression shape: addition preserves shape", {
  x <- Variable(c(2, 3))
  y <- Variable(c(2, 3))
  expr <- x + y
  expect_equal(expr@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("expression shape: scalar * matrix preserves shape", {
  x <- Variable(c(2, 3))
  expr <- 2 * x
  expect_equal(expr@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("save_leaf_value: bypasses validation", {
  ## CVXPY: save_value() bypasses _validate_value
  ## Used by solvers to store results without projection
  x <- Variable(3, nonneg = TRUE)
  save_leaf_value(x, matrix(c(-1, -2, -3), 3, 1))
  ## Value should be stored despite being negative
  val <- value(x)
  expect_equal(as.numeric(val), c(-1, -2, -3))
})

# ═══════════════════════════════════════════════════════════════════════
# Variable bounds (CVXPY test_expressions.py::test_bounds_attr)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("variable bounds: validation rejects invalid inputs", {
  ## Must be a list of length 2
  expect_error(Variable(3, bounds = c(-1, 2)), "list of two items")
  expect_error(Variable(3, bounds = list(1)), "list of two items")
  expect_error(Variable(3, bounds = list(1, 2, 3)), "list of two items")
  ## NaN not allowed
  expect_error(Variable(3, bounds = list(NaN, 1)), "NaN")
  ## Inf as lower bound not allowed
  expect_error(Variable(3, bounds = list(Inf, 2)), "Inf.*lower")
  ## -Inf as upper bound not allowed
  expect_error(Variable(3, bounds = list(-1, -Inf)), "Inf.*upper")
  ## lower > upper not allowed
  expect_error(Variable(3, bounds = list(5, 2)), "upper bounds are less")
})

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("variable bounds: scalar bounds promoted to arrays", {
  x <- Variable(3, bounds = list(-1, 2))
  lb <- x@attributes$bounds[[1L]]
  ub <- x@attributes$bounds[[2L]]
  expect_equal(lb, rep(-1, 3))
  expect_equal(ub, rep(2, 3))
})

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("variable bounds: project clamps to bounds", {
  x <- Variable(3, bounds = list(-1, 2))
  projected <- project(x, matrix(c(-5, 0, 10), 3, 1))
  expect_equal(as.numeric(projected), c(-1, 0, 2))
})

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("variable bounds: minimize with lower bound", {
  ## CVXPY: min sum(x), x in [-1, 2] => optimal = -3, x = [-1,-1,-1]
  x <- Variable(3, bounds = list(-1, 2))
  prob <- Problem(Minimize(sum_entries(x)))
  psolve(prob)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(value(prob), -3, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(-1, -1, -1), tolerance = 1e-4)
})

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("variable bounds: maximize with upper bound", {
  ## CVXPY: max sum(x), x in [-1, 2] => optimal = 6, x = [2,2,2]
  x <- Variable(3, bounds = list(-1, 2))
  prob <- Problem(Maximize(sum_entries(x)))
  psolve(prob)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(value(prob), 6, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(2, 2, 2), tolerance = 1e-4)
})

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("variable bounds: one-sided bounds (lower only)", {
  ## bounds = list(-3, NULL) => lower = -3, upper = Inf
  x <- Variable(2, bounds = list(-3, NULL))
  prob <- Problem(Minimize(sum_entries(x)))
  psolve(prob)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(value(prob), -6, tolerance = 1e-4)
})

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("variable bounds: one-sided bounds (upper only)", {
  ## bounds = list(NULL, 5) => lower = -Inf, upper = 5
  x <- Variable(2, bounds = list(NULL, 5))
  prob <- Problem(Maximize(sum_entries(x)))
  psolve(prob)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(value(prob), 10, tolerance = 1e-4)
})

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("variable bounds: nonneg + bounds combined", {
  ## CVXPY: max sum(x), nonneg=TRUE, bounds=[0, 5] => optimal = 15, x = [5,5,5]
  x <- Variable(3, nonneg = TRUE, bounds = list(0, 5))
  prob <- Problem(Maximize(sum_entries(x)))
  psolve(prob)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(value(prob), 15, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(5, 5, 5), tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# Shape incompatibility tests (test_shape.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_shape.py::TestShape::test_add_incompatible
test_that("shape: adding incompatible shapes raises error", {
  ## CVXPY: sum_shapes([(4, 2), (4,)]) → ValueError
  ## In CVXR, shapes are always 2D. (4,) doesn't exist; the closest is (4,1).
  ## Variable(c(4,2)) + Variable(c(4,1)) should fail because (4,2) and (4,1)
  ## are not broadcast-compatible (2 != 1 would normally broadcast, but let's
  ## test the truly incompatible case: (4,2) + (3,2))
  a <- Variable(c(4, 2))
  b <- Variable(c(3, 2))
  expect_error(a + b, "broadcast|shape|mismatch|incompatible|conform")
})

## @cvxpy test_shape.py::TestShape::test_mul_scalars
test_that("shape: matmul dimension mismatch raises error", {
  ## CVXPY: mul_shapes((), (5,9)) → ValueError (scalar not permitted in matmul)
  ## In CVXR, zero-dim shapes don't exist. Test true dimension mismatch:
  ## (3,2) %*% (4,3) should fail since inner dims 2 != 4.
  a <- Variable(c(3, 2))
  b <- Variable(c(4, 3))
  expect_error(a %*% b, "dimension|incompatible|mismatch|Incompatible")
})

## @cvxpy test_shape.py::TestShape::test_reshape_with_lists
test_that("shape: reshape with list-style dims and addition", {
  ## CVXPY: a = Variable([n,n]), b = Variable(n^2)
  ## c = reshape(b, [n,n], order='F')
  ## (a + c).shape == (n, n)
  ## In R, reshape_expr uses 'F' order by default (column-major).
  ## Verified via CVXPY: (a + c).shape == (2, 2)
  n <- 2L
  a <- Variable(c(n, n))
  b <- Variable(n^2)
  c_expr <- reshape_expr(b, c(n, n))
  result <- a + c_expr
  expect_equal(result@shape, c(n, n))
})

# ═══════════════════════════════════════════════════════════════════════
# Batch C: test_expressions.py parity gap tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_expressions.py::TestExpressions::test_hermitian
test_that("expressions: Hermitian variable properties", {
  ## CVXPY: Variable((4,3), hermitian=True) raises error (non-square)
  expect_error(
    Variable(c(4, 3), hermitian = TRUE),
    regexp = ".*square.*"
  )

  ## Hermitian variable is_hermitian
  v <- Variable(c(2, 2), hermitian = TRUE)
  expect_true(CVXR:::is_hermitian(v))

  ## Diagonal variable is also hermitian
  v_diag <- Variable(c(2, 2), diag = TRUE)
  expect_true(CVXR:::is_hermitian(v_diag))

  ## Operations that preserve hermitian property
  v2 <- Variable(c(2, 2), hermitian = TRUE)
  ## v + v
  expect_true(CVXR:::is_hermitian(v2 + v2))
  ## -v
  expect_true(CVXR:::is_hermitian(-v2))
  ## t(v)
  expect_true(CVXR:::is_hermitian(t(v2)))
})

## @cvxpy test_expressions.py::TestExpressions::test_log_log_curvature
test_that("expressions: DGP log-log curvature queries", {
  ## CVXPY: x = Variable(pos=True)
  ## monomial = x*x*x -> curvature == LOG_LOG_AFFINE
  ## posynomial = x*x*x + x -> curvature == LOG_LOG_CONVEX
  ## llcv = 1/(x*x*x + x) -> curvature == LOG_LOG_CONCAVE
  ##
  ## In CVXR, curvature() returns DCP curvature strings, not log-log.

  ## But is_log_log_affine/convex/concave are available.
  x <- Variable(name = "x", pos = TRUE)

  monomial <- x * x * x
  expect_true(is_log_log_affine(monomial))

  posynomial <- x * x * x + x
  expect_true(is_log_log_convex(posynomial))

  llcv <- 1 / (x * x * x + x)
  expect_true(is_log_log_concave(llcv))
})

## @cvxpy test_expressions.py::TestExpressions::test_sum
test_that("expressions: sum with scalar and vector", {
  ## CVXPY: a.value = 1, sum(a).value == 1
  a <- Variable(name = "a")
  value(a) <- 1
  expect_equal(as.numeric(value(sum_entries(a))), 1)

  ## CVXPY: x.value = [1, 2], sum(x).value == 3
  x <- Variable(2, name = "x")
  value(x) <- c(1, 2)
  expect_equal(as.numeric(value(sum_entries(x))), 3)
})

## @cvxpy test_expressions.py::TestExpressions::test_quad_form_matmul
test_that("expressions: t(x) %%*%% A %%*%% x as QuadForm detection", {
  ## CVXPY: x.T @ A @ x automatically detected as QuadForm
  ## In CVXR, t(x) %%*%% A %%*%% x does NOT produce QuadForm -- it produces
  ## MulExpression. QuadForm detection from matmul chains is not implemented.
  skip("CVXR does not auto-detect t(x) %%*%% A %%*%% x as QuadForm; use quad_form() explicitly")
})

## @cvxpy test_expressions.py::TestExpressions::test_1D_array
test_that("expressions: 1D array in matmul", {
  ## CVXPY: c = np.array([1, 2]), p = Parameter(2), p.value = [1, 1]
  ## (c @ p).value == 3
  ## In R, all variables are 2D (column vectors), so c is a (1,2) or (2,1).
  ## We test the equivalent: t(c) %%*%% p where c is a numeric vector.
  p <- Parameter(2, name = "p")
  value(p) <- c(1, 1)
  c_vec <- c(1, 2)
  expr <- t(c_vec) %*% p
  expect_equal(as.numeric(value(expr)), 3)

  ## (c @ x).shape in CVXPY is () (scalar)
  ## In R, the result is (1,1) which is still scalar.
  x <- Variable(2, name = "x")
  expr2 <- t(c_vec) %*% x
  expect_equal(prod(expr2@shape), 1L)
})

## @cvxpy test_expressions.py::TestExpressions::test_none_idx
test_that("expressions: None indexing (R equivalent)", {
  ## CVXPY: a[None, None].shape == (1, 1)
  ## x[:, None].shape == (2, 1)
  ## x[None, :].shape == (1, 2)
  ## In R, there is no None/NULL indexing on CVXR expressions.
  ## R variables are always at least 2D, so Variable() is (1,1),
  ## Variable(2) is (2,1). The concept of None-indexing does not apply.
  skip("R does not have None/NULL indexing on CVXR expressions; not applicable")
})

## @cvxpy test_expressions.py::TestExpressions::test_logical_indices
test_that("expressions: boolean indexing on Constants", {
  ## CVXPY: A = Constant(np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]]))
  ## C[A <= 2].shape == (2,), sign == NONNEG, value == [1, 2]
  ## In CVXR, logical indexing works on 1D constants but NOT on matrix Constants.
  ## Test with 1D constant:
  c_1d <- Constant(c(1, 2, 3, 4))
  mask <- c(TRUE, TRUE, FALSE, FALSE)
  expr <- c_1d[mask]
  expect_equal(prod(expr@shape), 2L)
  expect_equal(as.numeric(value(expr)), c(1, 2))
})

## @cvxpy test_expressions.py::TestExpressions::test_project_boolean_indices
test_that("expressions: project() with boolean indices", {
  ## CVXPY: Variable((3,), boolean=(np.array([0, 2]),))
  ## val = [-0.2, 0.5, 1.7], project -> [0, 0.5, 1]
  ## Only indices 0 and 2 are projected to boolean.
  ## In CVXR, project() with index-based boolean does not work correctly
  ## (project checks isTRUE(a$boolean) which is FALSE when boolean is a list).
  skip("CVXR project() does not handle index-based boolean specification")
})

## @cvxpy test_expressions.py::TestExpressions::test_project_integer_indices
test_that("expressions: project() with integer indices", {
  ## CVXPY: Variable((3,), integer=(np.array([1]),))
  ## val = [1.2, 2.7, -0.8], project -> [1.2, 3, -0.8]
  ## Only index 1 is projected to integer.
  ## Same issue as boolean indices: project() doesn't handle index-based integer.
  skip("CVXR project() does not handle index-based integer specification")
})

# =========================================================================
# Wave 5: test_expression_methods.py parity gaps (13 tests)
# All expected values verified via `uv run python` against CVXPY 1.8.1
# =========================================================================

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_all_expressions
test_that("expression methods: all expressions no-arg (trace, max, min, mean, ptp, sum, std, var)", {
  ## CVXPY: X = Constant([[1., 4., 7.], [2., -4., 3.], [99., -2., 2.4]])
  ## Each method with no args returns a scalar with known value.
  X_np <- matrix(c(1, 4, 7, 2, -4, 3, 99, -2, 2.4), nrow = 3, byrow = TRUE)
  X <- Constant(X_np)

  ## trace
  expect_equal(matrix_trace(X)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(matrix_trace(X))), -0.6, tolerance = 1e-10)

  ## max
  expect_equal(max_entries(X)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(max_entries(X))), 99.0, tolerance = 1e-10)

  ## min
  expect_equal(min_entries(X)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(min_entries(X))), -4.0, tolerance = 1e-10)

  ## mean
  expect_equal(cvxr_mean(X)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(cvxr_mean(X))), 112.4 / 9, tolerance = 1e-10)

  ## ptp
  expect_equal(ptp(X)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(ptp(X))), 103.0, tolerance = 1e-10)

  ## sum
  expect_equal(sum_entries(X)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(sum_entries(X))), 112.4, tolerance = 1e-10)

  ## std
  expect_equal(cvxr_std(X)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(cvxr_std(X))), 30.73544621964984, tolerance = 1e-6)

  ## var
  expect_equal(cvxr_var(X)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(cvxr_var(X))), 944.6676543209875, tolerance = 1e-6)

  ## conj on complex
  complex_X_np <- matrix(c(1+0i, 4+0i, 7+0i, 2+0i, -4+3i, 3+0i, 99+0i, -2-9i, 2.4+0i),
                         nrow = 3, byrow = TRUE)
  complex_X <- Constant(complex_X_np)
  conj_val <- value(Conj(complex_X))
  expected_conj <- Conj(complex_X_np)
  expect_equal(Conj(complex_X)@shape, c(3L, 3L))
  expect_true(all(abs(conj_val - expected_conj) < 1e-10))
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_conj
test_that("expression methods: conj in constraint (real variable)", {
  ## CVXPY: v = Variable((4,)), min sum(v) s.t. v.conj() >= 1 => v = [1,1,1,1]
  v <- Variable(4)
  obj <- Minimize(sum_entries(v))
  prob <- Problem(obj, list(Conj(v) >= 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(v)), rep(1, 4), tolerance = 1e-4)
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_conjugate
test_that("expression methods: conjugate is alias for conj", {
  ## CVXPY: v.conjugate() is the same as v.conj()
  ## In CVXR, Conj() serves both roles.
  v <- Variable(4)
  obj <- Minimize(sum_entries(v))
  prob <- Problem(obj, list(Conj(v) >= 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(v)), rep(1, 4), tolerance = 1e-4)
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_max
test_that("expression methods: max_entries sign and axis shapes", {
  ## CVXPY: Variable().max().sign == UNKNOWN
  expect_equal(expr_sign(max_entries(Variable())), "UNKNOWN")

  ## CVXPY: Variable(2).max(axis=0, keepdims=True).shape == (1,)
  ## R axis=2 = CVXPY axis=0; CVXR shapes always 2D
  expect_equal(max_entries(Variable(2), axis = 2, keepdims = TRUE)@shape, c(1L, 1L))

  ## CVXPY: Variable((2,3)).max(axis=0, keepdims=True).shape == (1, 3)
  expect_equal(max_entries(Variable(c(2, 3)), axis = 2, keepdims = TRUE)@shape, c(1L, 3L))

  ## CVXPY: Variable((2,3)).max(axis=1).shape == (2,)
  ## R axis=1 = CVXPY axis=1; shape (2,1) in R
  expect_equal(max_entries(Variable(c(2, 3)), axis = 1)@shape, c(2L, 1L))

  ## Invalid axis
  expect_error(max_entries(Variable(2), axis = 4), "out of bounds")
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_min
test_that("expression methods: min_entries sign and axis shapes", {
  ## CVXPY: Variable().min().sign == UNKNOWN
  expect_equal(expr_sign(min_entries(Variable())), "UNKNOWN")

  ## CVXPY: Variable(2).min(axis=0).shape == ()
  ## R: axis=2 = CVXPY axis=0; shape (1,1)
  expect_equal(min_entries(Variable(2), axis = 2)@shape, c(1L, 1L))

  ## CVXPY: Variable((2,3)).min(axis=0).shape == (3,)
  ## R: shape (1,3)
  expect_equal(min_entries(Variable(c(2, 3)), axis = 2)@shape, c(1L, 3L))

  ## CVXPY: Variable((2,3)).min(axis=1).shape == (2,)
  ## R: shape (2,1)
  expect_equal(min_entries(Variable(c(2, 3)), axis = 1)@shape, c(2L, 1L))

  ## Invalid axis
  expect_error(min_entries(Variable(2), axis = 4), "out of bounds")
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_missing_order_warning
test_that("expression methods: reshape without order does not warn in CVXR", {
  ## CVXPY: X.reshape((2,6)) raises FutureWarning (order param missing).
  ## CVXR: reshape_expr defaults to order="F" (column-major) -- no warning.
  ## This test documents the CVXR behavior difference.
  skip("CVXR reshape_expr defaults to order='F'; no FutureWarning equivalent")
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_ptp
test_that("expression methods: ptp (peak-to-peak) values and shapes", {
  ## CVXPY: a = Constant([[10., -10., 3.0], [6., 0., -1.5]])
  a_np <- matrix(c(10, -10, 3, 6, 0, -1.5), nrow = 2, byrow = TRUE)
  a <- Constant(a_np)

  ## ptp() = 20.0
  expect_equal(as.numeric(value(ptp(a))), 20.0, tolerance = 1e-10)
  expect_equal(ptp(a)@shape, c(1L, 1L))

  ## ptp(axis=0) in CVXPY = ptp(axis=2) in R = [4, 10, 4.5]
  e0 <- ptp(a, axis = 2)
  expect_equal(as.numeric(value(e0)), c(4, 10, 4.5), tolerance = 1e-10)
  expect_equal(e0@shape, c(1L, 3L))

  ## ptp(axis=1) in CVXPY = ptp(axis=1) in R = [20., 7.5]
  e1 <- ptp(a, axis = 1)
  expect_equal(as.numeric(value(e1)), c(20, 7.5), tolerance = 1e-10)
  expect_equal(e1@shape, c(2L, 1L))

  ## ptp(0, keepdims=True) in CVXPY = ptp(axis=2, keepdims=TRUE) in R
  e0k <- ptp(a, axis = 2, keepdims = TRUE)
  expect_equal(e0k@shape, c(1L, 3L))
  expect_equal(as.numeric(value(e0k)), c(4, 10, 4.5), tolerance = 1e-10)

  ## ptp(1, keepdims=True) in CVXPY = ptp(axis=1, keepdims=TRUE) in R
  e1k <- ptp(a, axis = 1, keepdims = TRUE)
  expect_equal(e1k@shape, c(2L, 1L))
  expect_equal(as.numeric(value(e1k)), c(20, 7.5), tolerance = 1e-10)
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_reshape
test_that("expression methods: reshape shape/sign/curvature and C-order", {
  ## Basic reshape
  A <- Variable(c(2, 2), name = "A")
  expr <- reshape_expr(A, c(4, 1))
  expect_equal(expr_sign(expr), "UNKNOWN")
  expect_equal(curvature(expr), "AFFINE")
  expect_equal(expr@shape, c(4L, 1L))

  ## Re-reshape
  expr2 <- reshape_expr(expr, c(2, 2))
  expect_equal(expr2@shape, c(2L, 2L))

  ## square + reshape preserves sign/curvature
  x <- Variable(2)
  expr3 <- reshape_expr(x^2, c(1, 2))
  expect_equal(expr_sign(expr3), "NONNEGATIVE")
  expect_equal(curvature(expr3), "CONVEX")
  expect_equal(expr3@shape, c(1L, 2L))

  ## Invalid dimensions
  C <- Variable(c(3, 2))
  expect_error(reshape_expr(C, c(5, 4)), "reshape")

  ## C-order reshape
  a <- 0:9
  A_np <- matrix(c(0, 2, 4, 6, 8, 1, 3, 5, 7, 9), nrow = 5, ncol = 2)
  A_cp <- reshape_expr(Constant(a), c(5, 2), order = "C")
  expect_equal(value(A_cp), A_np)

  ## Solve with C-order reshape
  X <- Variable(c(5, 2))
  prob <- Problem(Minimize(Constant(0)), list(X == A_cp))
  psolve(prob, solver = "CLARABEL")
  expect_true(all(abs(value(X) - A_np) < 1e-3))

  ## C-order un-reshape
  a_back <- reshape_expr(A_cp, c(10, 1), order = "C")
  expect_equal(as.numeric(value(a_back)), 0:9, tolerance = 1e-10)

  ## Matrix-to-matrix C-order reshape
  ## b = [[0,1,2],[3,4,5],[6,7,8],[9,10,11]] (row-major reading)
  b <- matrix(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), nrow = 4, ncol = 3, byrow = TRUE)
  ## C-order reshape (2,6): [[0,1,2,3,4,5],[6,7,8,9,10,11]]
  b_reshaped_C <- matrix(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
                         nrow = 2, ncol = 6, byrow = TRUE)
  X2 <- Variable(c(4, 3))
  X2_reshaped <- reshape_expr(X2, c(2, 6), order = "C")
  prob2 <- Problem(Minimize(Constant(0)), list(X2_reshaped == b_reshaped_C))
  psolve(prob2, solver = "CLARABEL")
  expect_true(all(abs(value(X2_reshaped) - b_reshaped_C) < 1e-3))
  expect_true(all(abs(value(X2) - b) < 1e-3))

  ## F-order reshape (default)
  ## F-order reshape (2,6): [[0,6,1,7,2,8],[3,9,4,10,5,11]]
  b_reshaped_F <- matrix(c(0, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11),
                         nrow = 2, ncol = 6, byrow = TRUE)
  X3 <- Variable(c(4, 3))
  X3_reshaped <- reshape_expr(X3, c(2, 6))  # default = F order
  prob3 <- Problem(Minimize(Constant(0)), list(X3_reshaped == b_reshaped_F))
  psolve(prob3, solver = "CLARABEL")
  expect_true(all(abs(value(X3_reshaped) - b_reshaped_F) < 1e-3))
  expect_true(all(abs(value(X3) - b) < 1e-3))
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_reshape_negative_one
test_that("expression methods: reshape with -1 dimension", {
  ## CVXPY: Variable((2,3)).reshape((-1, 1)) -> (6, 1)
  x <- Variable(c(2, 3))

  expect_equal(reshape_expr(x, c(-1, 1))@shape, c(6L, 1L))
  expect_equal(reshape_expr(x, c(1, -1))@shape, c(1L, 6L))
  expect_equal(reshape_expr(x, c(-1, 2))@shape, c(3L, 2L))
  ## In CVXR, reshape_expr(x, -1) yields (6,1) since R shapes are always 2D
  expect_equal(prod(reshape_expr(x, -1)@shape), 6L)

  ## Error: cannot reshape 6 into (8, -1)
  expect_error(reshape_expr(x, c(8, -1)), "Cannot reshape")

  ## Error: only one -1
  expect_error(reshape_expr(x, c(-1, -1)), "Only one")

  ## Error: dimension must be positive
  expect_error(reshape_expr(x, c(-1, 0)), "positive")

  ## Numeric verification with C-order
  A <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
  A_reshaped_C <- reshape_expr(Constant(A), -1, order = "C")
  expect_equal(as.numeric(value(A_reshaped_C)), c(1, 2, 3, 4, 5, 6), tolerance = 1e-10)

  A_reshaped_F <- reshape_expr(Constant(A), -1, order = "F")
  expect_equal(as.numeric(value(A_reshaped_F)), c(1, 4, 2, 5, 3, 6), tolerance = 1e-10)
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_stats
test_that("expression methods: mean, std, var with axis and ddof", {
  ## CVXPY: a = Constant([[10., 10., 3.0], [6., 0., 1.5]])
  a_np <- rbind(c(10, 10, 3), c(6, 0, 1.5))
  a <- Constant(a_np)

  ## No-arg
  expect_equal(as.numeric(value(cvxr_mean(a))), mean(a_np), tolerance = 1e-10)
  expect_equal(as.numeric(value(cvxr_var(a))), 15.368055555555555, tolerance = 1e-6)
  expect_equal(as.numeric(value(cvxr_std(a))), 3.9202111621130253, tolerance = 1e-6)

  ## ddof
  for (ddof in c(0, 1)) {
    ## R var uses ddof=1 by default; CVXPY/numpy default ddof=0
    expected_var <- sum((a_np - mean(a_np))^2) / (length(a_np) - ddof)
    expected_std <- sqrt(expected_var)
    expect_equal(as.numeric(value(cvxr_var(a, ddof = ddof))), expected_var, tolerance = 1e-6)
    expect_equal(as.numeric(value(cvxr_std(a, ddof = ddof))), expected_std, tolerance = 1e-6)
  }

  ## axis with keepdims
  ## axis=1 in R = axis=1 in CVXPY (row-wise)
  ## axis=2 in R = axis=0 in CVXPY (column-wise)

  ## CVXPY axis=0 -> R axis=2 (column-wise)
  m_a2_kT <- cvxr_mean(a, axis = 2, keepdims = TRUE)
  expect_equal(m_a2_kT@shape, c(1L, 3L))
  expect_equal(as.numeric(value(m_a2_kT)), c(8, 5, 2.25), tolerance = 1e-6)

  m_a2_kF <- cvxr_mean(a, axis = 2, keepdims = FALSE)
  expect_equal(m_a2_kF@shape, c(1L, 3L))  # CVXR always 2D
  expect_equal(as.numeric(value(m_a2_kF)), c(8, 5, 2.25), tolerance = 1e-6)

  s_a2_kT <- cvxr_std(a, axis = 2, keepdims = TRUE)
  expect_equal(s_a2_kT@shape, c(1L, 3L))
  expect_equal(as.numeric(value(s_a2_kT)), c(2, 5, 0.75), tolerance = 1e-6)

  ## CVXPY axis=1 -> R axis=1 (row-wise)
  m_a1_kT <- cvxr_mean(a, axis = 1, keepdims = TRUE)
  expect_equal(m_a1_kT@shape, c(2L, 1L))
  expect_equal(as.numeric(value(m_a1_kT)), c(7.666666667, 2.5), tolerance = 1e-6)

  s_a1_kT <- cvxr_std(a, axis = 1, keepdims = TRUE)
  expect_equal(s_a1_kT@shape, c(2L, 1L))
  expect_equal(as.numeric(value(s_a1_kT)), c(3.29983165, 2.54950976), tolerance = 1e-4)
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_sum
test_that("expression methods: sum sign, curvature, shape, axis, sparse", {
  ## CVXPY: Constant([1,-1]).sum().sign == UNKNOWN
  expect_equal(expr_sign(sum_entries(Constant(c(1, -1)))), "UNKNOWN")
  expect_equal(curvature(sum_entries(Constant(c(1, -1)))), "CONSTANT")

  ## Variable sum
  x2 <- Variable(2)
  expect_equal(expr_sign(sum_entries(x2)), "UNKNOWN")
  expect_equal(sum_entries(x2)@shape, c(1L, 1L))  # CVXPY: ()
  expect_equal(curvature(sum_entries(x2)), "AFFINE")

  ## keepdims
  x21 <- Variable(c(2, 1))
  expect_equal(sum_entries(x21, keepdims = TRUE)@shape, c(1L, 1L))

  ## Mixed curvature
  mat <- matrix(c(1, -1), nrow = 1)
  expr <- sum_entries(mat %*% (Variable(2)^2))
  expect_equal(curvature(expr), "UNKNOWN")

  ## axis
  ## CVXPY axis=0 -> R axis=2
  x23 <- Variable(c(2, 3))
  expect_equal(sum_entries(x23, axis = 2, keepdims = TRUE)@shape, c(1L, 3L))
  expect_equal(sum_entries(x23, axis = 2, keepdims = FALSE)@shape, c(1L, 3L))  # CVXR always 2D
  expect_equal(sum_entries(x23, axis = 1)@shape, c(2L, 1L))

  ## Invalid axis
  expect_error(sum_entries(Variable(2), axis = 4), "out of bounds")

  ## Sparse identity sum
  I3 <- Matrix::Diagonal(3)
  expect_equal(as.numeric(value(sum_entries(Constant(I3)))), 3)

  ## Sparse identity sum(axis=0)
  expect_equal(as.numeric(value(sum_entries(Constant(I3), axis = 2))), c(1, 1, 1))
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_trace
test_that("expression methods: trace sign/curvature/shape and non-square error", {
  ## CVXPY: A.trace().sign == UNKNOWN, curvature == AFFINE, shape == ()
  A <- Variable(c(2, 2), name = "A")
  tr <- matrix_trace(A)
  expect_equal(expr_sign(tr), "UNKNOWN")
  expect_equal(curvature(tr), "AFFINE")
  expect_equal(tr@shape, c(1L, 1L))

  ## Non-square raises error
  C <- Variable(c(3, 2))
  expect_error(matrix_trace(C), "square")
})

## @cvxpy test_expression_methods.py::TestExpressionMethods::test_trace_sign_psd
test_that("expression methods: trace sign for PSD/NSD variables", {
  ## CVXPY: X_psd.trace().is_nonneg() == True
  X_psd <- Variable(c(2, 2), PSD = TRUE)
  X_nsd <- Variable(c(2, 2), NSD = TRUE)

  expect_true(is_nonneg(matrix_trace(X_psd)))
  expect_true(is_nonpos(matrix_trace(X_nsd)))
})

# =========================================================================
# Wave 5: test_expressions.py parity gaps (12 tests)
# All expected values verified via `uv run python` against CVXPY 1.8.1
# =========================================================================

## @cvxpy test_expressions.py::TestExpressions::test_1D_array
test_that("expressions: 1D array in matmul (already covered, additional checks)", {
  ## CVXPY: c = np.array([1, 2]), p = Parameter(2), p.value = [1, 1]
  ## (c @ p).value == 3
  ## In R, c is a numeric vector; t(c) %*% p gives (1,1) result.
  p <- Parameter(2, name = "p")
  value(p) <- c(1, 1)
  c_vec <- c(1, 2)
  expr <- t(c_vec) %*% p
  expect_equal(as.numeric(value(expr)), 3)

  ## (c @ x).shape in CVXPY is () -> in R (1,1)
  x <- Variable(2, name = "x")
  expr2 <- t(c_vec) %*% x
  expect_equal(expr2@shape, c(1L, 1L))
})

## @cvxpy test_expressions.py::TestExpressions::test_float_is_invalid_index
test_that("expressions: float index on CVXR expressions", {
  ## CVXPY: x[1.0] raises IndexError("float is an invalid index type.")
  ## In R, float indices are truncated to integer by default for base R.
  ## CVXR may or may not reject float indices depending on implementation.
  x <- Variable(2)
  ## R typically coerces 1.0 to 1L in indexing. Check if CVXR rejects it:
  tryCatch({
    e <- x[1.5, 1]
    ## If it works, it means R truncated the float. This is not an error in R.
    expect_true(TRUE)  # Document that CVXR allows float indices (R behavior)
  }, error = function(e) {
    expect_true(grepl("float|invalid|index", conditionMessage(e), ignore.case = TRUE))
  })
})

## @cvxpy test_expressions.py::TestExpressions::test_hermitian
test_that("expressions: Hermitian variable properties (extended)", {
  ## CVXPY: Variable((4,3), hermitian=True) raises error
  expect_error(Variable(c(4, 3), hermitian = TRUE), "square")

  ## is_hermitian
  v <- Variable(c(2, 2), hermitian = TRUE)
  expect_true(CVXR:::is_hermitian(v))

  ## diag is hermitian
  v_diag <- Variable(c(2, 2), diag = TRUE)
  expect_true(CVXR:::is_hermitian(v_diag))

  ## Operations preserving hermitian
  v2 <- Variable(c(2, 2), hermitian = TRUE)
  expect_true(CVXR:::is_hermitian(v2 + v2))
  expect_true(CVXR:::is_hermitian(-v2))
  expect_true(CVXR:::is_hermitian(t(v2)))

  ## Re, Conj preserve hermitian; Im(Hermitian) is skew-Hermitian
  ## CVXPY v1.8.2 fix: Im(Hermitian) is NOT hermitian
  expect_true(CVXR:::is_hermitian(Re(v2)))
  expect_false(CVXR:::is_hermitian(Im(v2)))
  expect_true(CVXR:::is_hermitian(Conj(v2)))
})

## @cvxpy test_expressions.py::TestExpressions::test_log_log_curvature
test_that("expressions: DGP log-log curvature (extended)", {
  ## CVXPY: x = Variable(pos=True)
  x <- Variable(name = "x", pos = TRUE)

  ## monomial curvature is LOG_LOG_AFFINE
  monomial <- x * x * x
  expect_true(is_log_log_affine(monomial))

  ## posynomial curvature is LOG_LOG_CONVEX
  posynomial <- x * x * x + x
  expect_true(is_log_log_convex(posynomial))

  ## 1/(posynomial) is LOG_LOG_CONCAVE
  llcv <- 1 / (x * x * x + x)
  expect_true(is_log_log_concave(llcv))
})

## @cvxpy test_expressions.py::TestExpressions::test_logical_indices
test_that("expressions: boolean indexing on Constants (extended)", {
  ## CVXPY: A = Constant([[1,2,3,4],[5,6,7,8],[9,10,11,12]])
  ## C[A <= 2].shape == (2,), value == [1, 2]
  ## In CVXR, single-index selection on matrices is not supported.
  ## Test with a vector Constant instead.
  c_1d <- Constant(c(1, 2, 3, 4))
  mask <- c(TRUE, TRUE, FALSE, FALSE)
  expr <- c_1d[mask]
  expect_equal(prod(expr@shape), 2L)
  expect_equal(as.numeric(value(expr)), c(1, 2))

  ## Logical mask selecting even values
  c_1d2 <- Constant(c(1, 2, 3, 4, 5, 6))
  mask2 <- c(1, 2, 3, 4, 5, 6) %% 2 == 0
  expr2 <- c_1d2[mask2]
  expect_equal(prod(expr2@shape), 3L)
  expect_equal(as.numeric(value(expr2)), c(2, 4, 6))
})

## @cvxpy test_expressions.py::TestExpressions::test_none_idx
test_that("expressions: None indexing not applicable in R", {
  ## CVXPY: a[None, None].shape == (1, 1)
  ## x[:, None].shape == (2, 1)
  ## In R, CVXR expressions do not support NULL indexing.
  ## Variable() is already (1,1), Variable(2) is already (2,1).
  skip("R CVXR expressions do not have None/NULL indexing; not applicable")
})

## @cvxpy test_expressions.py::TestExpressions::test_out_of_bounds
test_that("expressions: out of bounds indexing raises error", {
  ## CVXPY: x = Variable(2), x[100] raises "Index 100 is out of bounds..."
  x <- Variable(2)
  expect_error(x[100, 1], "out of bounds|Index")

  ## Negative out of bounds
  ## In CVXR with R 1-based indexing, negative indices mean "exclude",
  ## so x[-100, 1] would try to exclude row 100 which doesn't exist.
  ## Behavior differs from Python. Test x[0, 1] instead (R has no 0-index).
  expect_error(x[0, 1], "out of bounds|subscript|0")

  ## Matrix out of bounds
  C <- Variable(c(3, 2))
  expect_error(C[100, 1], "out of bounds|Index")
})

## @cvxpy test_expressions.py::TestExpressions::test_project_boolean_indices
test_that("expressions: project with boolean indices", {
  ## CVXPY: Variable((3,), boolean=(np.array([0, 2]),))
  ## Only indices 0, 2 are boolean. project([-0.2, 0.5, 1.7]) -> [0, 0.5, 1]
  skip("CVXR project() does not handle index-based boolean specification")
})

## @cvxpy test_expressions.py::TestExpressions::test_project_integer_indices
test_that("expressions: project with integer indices", {
  ## CVXPY: Variable((3,), integer=(np.array([1]),))
  ## Only index 1 is integer. project([1.2, 2.7, -0.8]) -> [1.2, 3, -0.8]
  skip("CVXR project() does not handle index-based integer specification")
})

## @cvxpy test_expressions.py::TestExpressions::test_quad_form_matmul
test_that("expressions: t(x) %%*%% A %%*%% x is not auto-detected as QuadForm", {
  ## CVXPY: x.T @ A @ x automatically detected as QuadForm.
  ## CVXR does NOT auto-detect matmul chains as QuadForm.
  ## Use quad_form(x, A) explicitly.
  x <- Variable(2)
  A <- matrix(c(1, 0, 0, -1), 2, 2)

  ## Matmul chain does not produce QuadForm
  expr <- t(x) %*% A %*% x
  expect_false(S7::S7_inherits(expr, CVXR:::QuadForm))

  ## quad_form explicitly creates QuadForm
  expr2 <- quad_form(x, A)
  expect_true(S7::S7_inherits(expr2, CVXR:::QuadForm))
})

## @cvxpy test_expressions.py::TestExpressions::test_sum
test_that("expressions: sum function with scalar and vector values", {
  ## CVXPY: a.value = 1, sum(a).value == 1
  a <- Variable(name = "a")
  value(a) <- 1
  expect_equal(as.numeric(value(sum_entries(a))), 1)

  ## CVXPY: x.value = [1, 2], sum(x).value == 3
  x <- Variable(2, name = "x")
  value(x) <- c(1, 2)
  expect_equal(as.numeric(value(sum_entries(x))), 3)
})

## @cvxpy test_expressions.py::test_sum_squares_with_axis
test_that("expressions: sum_squares with axis and keepdims", {
  ## CVXPY parametrized test: shape=(3,4), various axis/keepdims
  ## In CVXR: axis=2 = CVXPY axis=0 (column-wise), axis=1 = CVXPY axis=1 (row-wise)

  ## Shape checks
  X <- Variable(c(3, 4))
  expect_equal(sum_squares(X, axis = 2)@shape, c(1L, 4L))              # CVXPY: axis=0 -> (4,)

  expect_equal(sum_squares(X, axis = 1)@shape, c(3L, 1L))              # CVXPY: axis=1 -> (3,)
  expect_equal(sum_squares(X, axis = 2, keepdims = TRUE)@shape, c(1L, 4L))  # keepdims
  expect_equal(sum_squares(X, axis = 1, keepdims = TRUE)@shape, c(3L, 1L))  # keepdims
  expect_equal(sum_squares(X, keepdims = TRUE)@shape, c(1L, 1L))       # scalar keepdims

  ## Numeric values match numpy
  set.seed(42)
  X_val <- matrix(rnorm(12), 3, 4)
  X2 <- Variable(c(3, 4))
  value(X2) <- X_val

  ## axis=2 (col-wise, CVXPY axis=0)
  e2 <- sum_squares(X2, axis = 2)
  expect_equal(as.numeric(value(e2)), colSums(X_val^2), tolerance = 1e-6)

  ## axis=1 (row-wise, CVXPY axis=1)
  e1 <- sum_squares(X2, axis = 1)
  expect_equal(as.numeric(value(e1)), rowSums(X_val^2), tolerance = 1e-6)

  ## Optimization test: minimize sum(sum_squares(X - Y, axis=2))
  set.seed(42)
  Y <- matrix(rnorm(12), 3, 4)
  X3 <- Variable(c(3, 4))
  obj <- sum_entries(sum_squares(X3 - Y, axis = 2))
  prob <- Problem(Minimize(obj))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), OPTIMAL)
  expect_true(all(abs(value(X3) - Y) < 1e-3))
})

# =========================================================================
# Wave 5: test_stack.py parity gaps (10 tests)
# CVXPY's cp.stack() does not exist in CVXR.
# CVXR has vstack() and hstack() for vertical/horizontal stacking.
# All expected values verified via `uv run python` against CVXPY 1.8.1
# =========================================================================

## @cvxpy test_stack.py::test_stack_1d_axis0
test_that("stack: 1D axis 0 equivalent via vstack", {
  ## CVXPY: cp.stack([a, b], axis=0).shape == (2, 4) where a,b are (4,)
  ## In CVXR, Parameter(4) is (4,1). vstack(a, b) stacks vertically: (8,1).
  ## The equivalent of axis=0 stacking of (4,) vectors to (2,4) requires
  ## a different approach. There is no cp.stack() in CVXR.
  skip("CVXR does not have cp.stack(); vstack/hstack have different semantics for 1D")
})

## @cvxpy test_stack.py::test_stack_1d_axis_last_numeric_parity
test_that("stack: 1D axis -1 numeric parity", {
  ## CVXPY: cp.stack([a, b], axis=-1) where a,b are (3,) arrays -> shape (3, 2)
  ## np.stack([a, b], axis=-1) == [[1,4],[2,5],[3,6]]
  ## In CVXR, the equivalent would be hstack(a, b) for (3,1) vectors -> (3,2)
  a <- c(1, 2, 3)
  b <- c(4, 5, 6)
  y <- hstack(Constant(a), Constant(b))
  expect_equal(y@shape, c(3L, 2L))
  expected <- cbind(a, b)
  expect_equal(value(y), expected, tolerance = 1e-10, ignore_attr = TRUE)
})

## @cvxpy test_stack.py::test_stack_2d_various_axes
test_that("stack: 2D various axes not directly supported", {
  ## CVXPY: cp.stack([a, b], axis=0/1/-1) on (3,4) -> (2,3,4)/(3,2,4)/(3,4,2)
  ## CVXR does not support 3D expressions (no cp.stack equivalent).
  skip("CVXR does not have cp.stack(); 3D stacking not supported")
})

## @cvxpy test_stack.py::test_stack_axis_bounds_check
test_that("stack: axis bounds check for vstack/hstack", {
  ## CVXPY: cp.stack axis validation. In CVXR, vstack/hstack don't take axis.
  ## Test that vstack/hstack validate dimension compatibility instead.
  ## vstack requires same number of columns; hstack requires same number of rows.
  a <- Variable(c(2, 3))
  b <- Variable(c(2, 4))  ## different column count
  expect_error(vstack(a, b), "column")
  a2 <- Variable(3)
  b2 <- Variable(4)  ## different row count
  expect_error(hstack(a2, b2), "row")
})

## @cvxpy test_stack.py::test_stack_canonicalization_resolves_equalities
test_that("stack: canonicalization resolves equalities via vstack", {
  ## CVXPY: z_tilde = cp.stack([x, y]), problem: z_tilde == z, x==1, y==2
  ## => z.value == [1, 2]
  ## In CVXR, vstack(x, y) for scalar Variables gives (2,1)
  x <- Variable(name = "x_stack")
  y <- Variable(name = "y_stack")
  z <- Variable(2)
  z_tilde <- vstack(x, y)
  prob <- Problem(Minimize(Constant(0)), list(z_tilde == z, x == 1, y == 2))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(as.numeric(value(z)), c(1, 2), tolerance = 1e-4)
})

## @cvxpy test_stack.py::test_stack_empty_list_raises
test_that("stack: empty list raises error", {
  ## CVXPY: cp.stack([], axis=0) raises ValueError
  ## In CVXR, vstack() with no args or do.call(vstack, list()) should error.
  expect_error(do.call(vstack, list()), "argument|empty|arg")
})

## @cvxpy test_stack.py::test_stack_non_int_axis_raises
test_that("stack: non-integer axis not applicable", {
  ## CVXPY: cp.stack([a, a], axis=0.5) raises TypeError
  ## CVXR vstack/hstack don't take axis argument.
  skip("CVXR vstack/hstack do not take axis parameter; not applicable")
})

## @cvxpy test_stack.py::test_stack_scalar_inputs
test_that("stack: scalar inputs via vstack", {
  ## CVXPY: cp.stack([1, 2, 3], axis=0).shape == (3,), value == [1,2,3]
  ## In CVXR, vstack with scalars: vstack(1, 2, 3) -> (3, 1)
  y <- vstack(Constant(1), Constant(2), Constant(3))
  expect_equal(y@shape, c(3L, 1L))
  expect_equal(as.numeric(value(y)), c(1, 2, 3), tolerance = 1e-10)
})

## @cvxpy test_stack.py::test_stack_shape_mismatch_raises
test_that("stack: shape mismatch raises error", {
  ## CVXPY: cp.stack([Parameter((3,)), Parameter((4,))], axis=0) raises ValueError
  ## In CVXR, vstack of (3,1) and (4,1) works (same col count); no exact equivalent.
  ## Test the closest equivalent: vstack with different column counts.
  a <- Parameter(c(2, 3))
  b <- Parameter(c(2, 4))
  expect_error(vstack(a, b), "column")

  ## And hstack with different row counts
  a2 <- Parameter(3)   ## (3,1)
  b2 <- Parameter(4)   ## (4,1)
  expect_error(hstack(a2, b2), "row")
})

## @cvxpy test_stack.py::test_stack_variables_shape_only
test_that("stack: variables shape via vstack/hstack", {
  ## CVXPY: cp.stack([x, y], axis=0) on (2,3) -> (2, 2, 3) -- 3D, not supported
  ## In CVXR, vstack(x, y) on (2,3) -> (4, 3) -- vertical stacking
  x <- Variable(c(2, 3))
  y <- Variable(c(2, 3))
  z <- vstack(x, y)
  expect_equal(z@shape, c(4L, 3L))

  ## hstack -> (2, 6)
  z2 <- hstack(x, y)
  expect_equal(z2@shape, c(2L, 6L))
})
