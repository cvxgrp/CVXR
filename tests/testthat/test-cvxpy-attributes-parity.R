## Tests for Variable/Parameter attributes parity with CVXPY
##
## 12 CRITICAL gaps from the gap analysis:
##  1. test_variable_bounds
##  2. test_boolean_var_value
##  3. test_integer_var_value
##  4. test_multiple_attributes
##  5. test_nonneg_PSD
##  6. test_nonpos_NSD
##  7. test_integer_bounds
##  8. test_nonpos_nonneg_variable
##  9. test_parameter_multiple_attributes
## 10. test_parameter_bounds_and_attributes
## 11. test_parameter_psd_and_attributes
## 12. test_bool_int_variable
##
## All expected values verified against CVXPY 1.8.1 using `uv run python`.
##
## Both BUG-A (boolean/integer stripped) and BUG-B (nonneg/nonpos else-if)
## in CvxAttr2Constr were fixed at v1.8.0-9040. All tests now pass.

# ====================================================================
# 1. test_variable_bounds
# ====================================================================

## @cvxpy test_attributes.py::TestAttributes::test_variable_bounds
test_that("Variable bounds create correct constraints (CVXPY parity)", {
  ## Lower bound only: minimize x with x >= 2 -> optimal at 2
  x <- Variable(bounds = list(2, NULL))
  prob <- Problem(Minimize(x))
  psolve(prob, solver = "SCS")
  expect_equal(as.numeric(value(x)), 2, tolerance = 1e-3)

  ## Upper bound only: maximize x with x <= 5 -> optimal at 5
  x2 <- Variable(bounds = list(NULL, 5))
  prob2 <- Problem(Maximize(x2))
  psolve(prob2, solver = "SCS")
  expect_equal(as.numeric(value(x2)), 5, tolerance = 1e-3)

  ## Both bounds: minimize -> lower, maximize -> upper
  x3 <- Variable(bounds = list(1, 5))
  prob3 <- Problem(Minimize(x3))
  psolve(prob3, solver = "SCS")
  expect_equal(as.numeric(value(x3)), 1, tolerance = 1e-3)

  x4 <- Variable(bounds = list(1, 5))
  prob4 <- Problem(Maximize(x4))
  psolve(prob4, solver = "SCS")
  expect_equal(as.numeric(value(x4)), 5, tolerance = 1e-3)

  ## Vector bounds with per-element lower bounds
  x5 <- Variable(3, bounds = list(c(1, 2, 3), c(4, 5, 6)))
  prob5 <- Problem(Minimize(sum(x5)))
  psolve(prob5, solver = "SCS")
  expect_equal(as.numeric(value(x5)), c(1, 2, 3), tolerance = 1e-3)

  ## Maximize -> hits upper bounds
  x6 <- Variable(3, bounds = list(c(1, 2, 3), c(4, 5, 6)))
  prob6 <- Problem(Maximize(sum(x6)))
  psolve(prob6, solver = "SCS")
  expect_equal(as.numeric(value(x6)), c(4, 5, 6), tolerance = 1e-3)

  ## Scalar bounds broadcast to vector
  x7 <- Variable(2, bounds = list(1, 2))
  prob7 <- Problem(Minimize(sum(x7)))
  psolve(prob7, solver = "SCS")
  expect_equal(as.numeric(value(x7)), c(1, 1), tolerance = 1e-3)

  ## Mixed: array lower, scalar upper
  x8 <- Variable(2, bounds = list(c(4, 5), 6))
  prob8 <- Problem(Maximize(sum(x8)))
  psolve(prob8, solver = "SCS")
  expect_equal(as.numeric(value(x8)), c(6, 6), tolerance = 1e-3)

  ## Mixed: scalar lower, array upper
  x9 <- Variable(3, bounds = list(1, c(2, 3, 4)))
  prob9 <- Problem(Maximize(sum(x9)))
  psolve(prob9, solver = "SCS")
  expect_equal(as.numeric(value(x9)), c(2, 3, 4), tolerance = 1e-3)
})

## @cvxpy test_attributes.py::TestAttributes::test_variable_bounds
test_that("Variable bounds validation errors (CVXPY parity)", {
  ## Bounds must be list of 2
  expect_error(Variable(2, bounds = list(c(0, 1, 2))),
               "two items")

  ## Upper < lower is invalid
  expect_error(Variable(2, bounds = list(c(2, 3), c(1, 4))),
               "upper bounds.*less than")

  ## -Inf as upper bound is invalid
  expect_error(Variable(2, bounds = list(NULL, -Inf)),
               "Inf.*upper bound")

  ## Inf as lower bound is invalid
  expect_error(Variable(2, bounds = list(Inf, Inf)),
               "Inf.*lower bound")

  ## NaN is invalid
  expect_error(Variable(2, bounds = list(NaN, NaN)),
               "NaN.*bound")
})

## @cvxpy test_attributes.py::TestAttributes::test_variable_bounds
test_that("Variable bounds with -Inf/Inf are effectively unconstrained", {
  ## -Inf lower, finite upper
  x <- Variable(2, bounds = list(-Inf, c(1, 2)))
  prob <- Problem(Maximize(sum(x)))
  psolve(prob, solver = "SCS")
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-3)

  ## Finite lower, Inf upper
  x2 <- Variable(3, bounds = list(3, Inf))
  prob2 <- Problem(Minimize(sum(x2)))
  psolve(prob2, solver = "SCS")
  expect_equal(as.numeric(value(x2)), c(3, 3, 3), tolerance = 1e-3)
})

# ====================================================================
# 2. test_boolean_var_value
# ====================================================================

## @cvxpy test_attributes.py::TestAttributes::test_boolean_var_value
test_that("Boolean variable solving yields 0 or 1 (CVXPY parity)", {
  ## Minimize boolean variable -> 0
  x <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(x))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), 0)

  ## Maximize boolean variable -> 1
  x2 <- Variable(boolean = TRUE)
  prob2 <- Problem(Maximize(x2))
  psolve(prob2, solver = "GLPK_MI")
  expect_equal(status(prob2), "optimal")
  expect_equal(round(as.numeric(value(x2))), 1)

  ## Boolean with constraint x >= 0.5 -> x = 1
  x3 <- Variable(boolean = TRUE)
  prob3 <- Problem(Minimize(x3), list(x3 >= 0.5))
  psolve(prob3, solver = "GLPK_MI")
  expect_equal(status(prob3), "optimal")
  expect_equal(round(as.numeric(value(x3))), 1)

  ## Boolean vector: minimize sum with x >= 0.5 -> all 1
  x4 <- Variable(3, boolean = TRUE)
  prob4 <- Problem(Minimize(sum(x4)), list(x4 >= 0.5))
  psolve(prob4, solver = "GLPK_MI")
  expect_equal(status(prob4), "optimal")
  expect_equal(round(as.numeric(value(x4))), c(1, 1, 1))

  ## Boolean vector: maximize sum -> all 1
  x5 <- Variable(3, boolean = TRUE)
  prob5 <- Problem(Maximize(sum(x5)))
  psolve(prob5, solver = "GLPK_MI")
  expect_equal(status(prob5), "optimal")
  expect_equal(round(as.numeric(value(x5))), c(1, 1, 1))
})

# ====================================================================
# 3. test_integer_var_value
# ====================================================================

## @cvxpy test_attributes.py::TestAttributes::test_integer_var_value
test_that("Integer variable solving yields integer values (CVXPY parity)", {
  ## Minimize integer variable with explicit constraints -> rounds up
  ## Using explicit constraints (not bounds attribute) to avoid BUG-A
  x <- Variable(integer = TRUE)
  prob <- Problem(Minimize(x), list(x >= 2.7, x <= 10))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), 3)

  ## Maximize integer variable with explicit constraints -> rounds down
  x2 <- Variable(integer = TRUE)
  prob2 <- Problem(Maximize(x2), list(x2 >= 0, x2 <= 7.3))
  psolve(prob2, solver = "GLPK_MI")
  expect_equal(status(prob2), "optimal")
  expect_equal(round(as.numeric(value(x2))), 7)

  ## Integer vector: minimize sum with lower bounds
  y <- Variable(3, integer = TRUE)
  prob3 <- Problem(Minimize(sum(y)), list(y >= 0.5, y <= 5))
  psolve(prob3, solver = "GLPK_MI")
  expect_equal(status(prob3), "optimal")
  expect_equal(round(as.numeric(value(y))), c(1, 1, 1))
})

# ====================================================================
# 4. test_multiple_attributes
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_multiple_attributes
test_that("Variable with multiple attributes: attribute storage (CVXPY parity)", {
  ## Attributes are stored correctly even in combinations
  x1 <- Variable(boolean = TRUE, nonneg = TRUE)
  expect_true(x1@attributes$boolean)
  expect_true(x1@attributes$nonneg)

  x2 <- Variable(integer = TRUE, nonneg = TRUE)
  expect_true(x2@attributes$integer)
  expect_true(x2@attributes$nonneg)

  x3 <- Variable(integer = TRUE, nonpos = TRUE)
  expect_true(x3@attributes$integer)
  expect_true(x3@attributes$nonpos)
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_multiple_attributes
test_that("Variable with boolean+nonneg solving (CVXPY parity)", {
  ## CvxAttr2Constr now preserves boolean/integer when convex attributes trigger replacement.
  x <- Variable(boolean = TRUE, nonneg = TRUE)
  prob <- Problem(Minimize(x), list(x >= 0.5))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), 1)
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_multiple_attributes
test_that("Variable with integer+nonneg solving (CVXPY parity)", {
  ## CvxAttr2Constr now preserves integer when nonneg triggers replacement.
  x <- Variable(integer = TRUE, nonneg = TRUE)
  prob <- Problem(Minimize(x), list(x >= 1.5))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), 2)
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_multiple_attributes
test_that("Variable with integer+nonpos solving (CVXPY parity)", {
  ## CvxAttr2Constr now preserves integer when nonpos triggers replacement.
  x <- Variable(integer = TRUE, nonpos = TRUE)
  prob <- Problem(Maximize(x), list(x <= -1.5))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), -2)
})

# ====================================================================
# 5. test_nonneg_PSD
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_nonneg_PSD
test_that("PSD variable that is also nonneg (CVXPY parity)", {
  ## PSD + nonneg: CVXPY enforces both PSD constraint and nonneg constraint.
  ## In CVXR, only PSD is enforced (BUG-A: nonneg branch is skipped because
  ## PSD branch takes priority in if-else). However, for minimize-trace with
  ## diagonal >= 1, the optimal PSD solution (identity matrix) happens to
  ## have nonneg entries, so the test passes coincidentally.
  X <- Variable(c(2, 2), PSD = TRUE, nonneg = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(X[1, 1] >= 1, X[2, 2] >= 1))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(prob)), 2, tolerance = 1e-3)

  ## Check solution: diagonal = 1, off-diagonal ~0
  Xval <- value(X)
  expect_equal(Xval[1, 1], 1, tolerance = 1e-3)
  expect_equal(Xval[2, 2], 1, tolerance = 1e-3)
  expect_equal(Xval[1, 2], 0, tolerance = 1e-2)
  expect_equal(Xval[2, 1], 0, tolerance = 1e-2)

  ## Check attributes are stored
  expect_true(X@attributes$PSD)
  expect_true(X@attributes$nonneg)
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_nonneg_PSD
test_that("PSD+nonneg: nonneg constraint is actually enforced (CVXPY parity)", {
  ## nonneg/nonpos constraints are now added independently of PSD/NSD.
  X <- Variable(c(2, 2), PSD = TRUE, nonneg = TRUE)
  prob <- Problem(Minimize(matrix_trace(X) + 2 * X[1, 2]),
                  list(X[1, 1] >= 1, X[2, 2] >= 1,
                       X[1, 1] <= 2, X[2, 2] <= 2))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  Xval <- value(X)
  ## With nonneg: off-diagonal must be >= 0
  expect_true(Xval[1, 2] >= -1e-4)
})

# ====================================================================
# 6. test_nonpos_NSD
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_nonpos_NSD
test_that("NSD variable that is also nonpos (CVXPY parity)", {
  ## NSD + nonpos: CVXPY enforces both NSD constraint and nonpos constraint.
  ## Same as test 5: NSD branch takes priority, nonpos is coincidentally
  ## satisfied for maximize-trace with diagonal <= -1.
  X <- Variable(c(2, 2), NSD = TRUE, nonpos = TRUE)
  prob <- Problem(Maximize(matrix_trace(X)),
                  list(X[1, 1] <= -1, X[2, 2] <= -1))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(prob)), -2, tolerance = 1e-3)

  ## Check solution: diagonal = -1, off-diagonal ~0
  Xval <- value(X)
  expect_equal(Xval[1, 1], -1, tolerance = 1e-3)
  expect_equal(Xval[2, 2], -1, tolerance = 1e-3)
  expect_equal(Xval[1, 2], 0, tolerance = 1e-2)
  expect_equal(Xval[2, 1], 0, tolerance = 1e-2)

  ## Check attributes
  expect_true(X@attributes$NSD)
  expect_true(X@attributes$nonpos)
})

# ====================================================================
# 7. test_integer_bounds
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_integer_bounds
test_that("Integer variable with bounds, integer endpoints (CVXPY parity)", {
  ## Integer + bounds [2, 7]: minimize -> 2, maximize -> 7
  ## Works because bounds endpoints are integers, so even without integer
  ## attribute being preserved, the continuous solution is integer.
  x <- Variable(integer = TRUE, bounds = list(2, 7))
  prob <- Problem(Minimize(x))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), 2)

  x2 <- Variable(integer = TRUE, bounds = list(2, 7))
  prob2 <- Problem(Maximize(x2))
  psolve(prob2, solver = "GLPK_MI")
  expect_equal(status(prob2), "optimal")
  expect_equal(round(as.numeric(value(x2))), 7)

  ## Check attributes are stored
  expect_true(x@attributes$integer)
  expect_true(is.list(x@attributes$bounds))
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_integer_bounds
test_that("Integer variable with bounds, non-integer endpoints (CVXPY parity)", {
  ## CvxAttr2Constr now preserves integer when bounds triggers replacement.
  x <- Variable(integer = TRUE, bounds = list(2.3, 6.7))
  prob <- Problem(Minimize(x))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), 3)

  x2 <- Variable(integer = TRUE, bounds = list(2.3, 6.7))
  prob2 <- Problem(Maximize(x2))
  psolve(prob2, solver = "GLPK_MI")
  expect_equal(status(prob2), "optimal")
  expect_equal(round(as.numeric(value(x2))), 6)
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_integer_bounds
test_that("Integer variable with explicit constraint bounds (CVXPY parity workaround)", {
  ## Workaround for BUG-A: use explicit constraints instead of bounds attribute.
  ## This correctly preserves the integer attribute.
  x <- Variable(integer = TRUE)
  prob <- Problem(Minimize(x), list(x >= 2.3, x <= 6.7))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), 3)

  x2 <- Variable(integer = TRUE)
  prob2 <- Problem(Maximize(x2), list(x2 >= 2.3, x2 <= 6.7))
  psolve(prob2, solver = "GLPK_MI")
  expect_equal(status(prob2), "optimal")
  expect_equal(round(as.numeric(value(x2))), 6)
})

# ====================================================================
# 8. test_nonpos_nonneg_variable
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_nonpos_nonneg_variable
test_that("Variable with both nonpos and nonneg forces zero (CVXPY parity)", {
  ## nonneg/nonpos constraints are now added independently (both applied).
  x <- Variable(nonneg = TRUE, nonpos = TRUE)
  prob <- Problem(Maximize(x))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 0, tolerance = 1e-4)
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_nonpos_nonneg_variable
test_that("Variable with both nonpos and nonneg: attribute storage (CVXPY parity)", {
  ## Verify that both attributes are stored, even if enforcement is buggy
  x <- Variable(nonneg = TRUE, nonpos = TRUE)
  expect_true(x@attributes$nonneg)
  expect_true(x@attributes$nonpos)
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_nonpos_nonneg_variable
test_that("Variable with both nonpos and nonneg: workaround (CVXPY parity)", {
  ## Workaround: use explicit constraints to enforce x >= 0 and x <= 0
  x <- Variable()
  prob <- Problem(Maximize(x), list(x >= 0, x <= 0))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 0, tolerance = 1e-4)

  x2 <- Variable()
  prob2 <- Problem(Minimize(x2), list(x2 >= 0, x2 <= 0))
  psolve(prob2, solver = "SCS")
  expect_equal(status(prob2), "optimal")
  expect_equal(as.numeric(value(x2)), 0, tolerance = 1e-4)
})

# ====================================================================
# 9. test_parameter_multiple_attributes
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_parameter_multiple_attributes
test_that("Parameter with nonneg validates values correctly (CVXPY parity)", {
  ## Nonneg parameter: accepts positive values, rejects negative
  p <- Parameter(nonneg = TRUE)
  value(p) <- 5
  expect_equal(as.numeric(value(p)), 5)

  expect_error(value(p) <- -1, "nonneg")

  ## Nonpos parameter: accepts negative values, rejects positive
  p2 <- Parameter(nonpos = TRUE)
  value(p2) <- -3
  expect_equal(as.numeric(value(p2)), -3)

  expect_error(value(p2) <- 1, "nonpos")

  ## Nonneg + nonpos parameter: only accepts zero
  ## Note: Parameter value validation uses project() which for multiple
  ## attributes (n_attr > 1) returns val as-is, then checks tolerance.
  ## With both nonneg+nonpos, n_attr=2, projection skipped, val=0 passes.
  p3 <- Parameter(nonneg = TRUE, nonpos = TRUE)
  value(p3) <- 0
  expect_equal(as.numeric(value(p3)), 0)
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_parameter_multiple_attributes
test_that("Parameter with PSD validates matrix values (CVXPY parity)", {
  ## PSD parameter: accepts PSD matrices
  P <- Parameter(c(2, 2), PSD = TRUE)
  value(P) <- diag(2)
  expect_true(!is.null(value(P)))

  ## PSD parameter: rejects NSD matrices
  expect_error(value(P) <- -diag(2), "positive semidefinite")
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_parameter_multiple_attributes
test_that("Parameter with NSD validates matrix values (CVXPY parity)", {
  ## NSD parameter: accepts NSD matrices
  P <- Parameter(c(2, 2), NSD = TRUE)
  value(P) <- -diag(2)
  expect_true(!is.null(value(P)))

  ## NSD parameter: rejects PSD matrices
  expect_error(value(P) <- diag(2), "negative semidefinite")
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_parameter_multiple_attributes
test_that("Parameter with symmetric validates matrix symmetry (CVXPY parity)", {
  ## Symmetric parameter: accepts symmetric matrices
  P <- Parameter(c(3, 3), symmetric = TRUE)
  val <- matrix(c(1, 2, 3, 2, 5, 6, 3, 6, 9), 3, 3)
  value(P) <- val
  expect_true(!is.null(value(P)))

  ## Symmetric parameter: rejects non-symmetric matrices
  nonsym <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, 3)
  expect_error(value(P) <- nonsym, "symmetric")
})

# ====================================================================
# 10. test_parameter_bounds_and_attributes
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_parameter_bounds_and_attributes
test_that("Parameter with bounds validates values (CVXPY parity)", {
  ## Parameter with bounds [0, 10]: accepts in-range, rejects out-of-range
  p <- Parameter(bounds = list(0, 10))
  value(p) <- 5
  expect_equal(as.numeric(value(p)), 5)

  expect_error(value(p) <- -1)
  expect_error(value(p) <- 11)

  ## Parameter with bounds + nonneg: attributes stored correctly
  p2 <- Parameter(nonneg = TRUE, bounds = list(0, 10))
  value(p2) <- 5
  expect_equal(as.numeric(value(p2)), 5)
  expect_true(p2@attributes$nonneg)
  expect_true(is.list(p2@attributes$bounds))
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_parameter_bounds_and_attributes
test_that("Parameter with bounds used in optimization (CVXPY parity)", {
  ## Parameter value affects the problem solution
  p <- Parameter(bounds = list(0, 10))
  value(p) <- 3
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= p))
  psolve(prob, solver = "SCS")
  expect_equal(as.numeric(value(x)), 3, tolerance = 1e-3)

  ## Change parameter value and re-solve with fresh problem
  value(p) <- 7
  x2 <- Variable()
  prob2 <- Problem(Minimize(x2), list(x2 >= p))
  psolve(prob2, solver = "SCS")
  expect_equal(as.numeric(value(x2)), 7, tolerance = 1e-3)
})

# ====================================================================
# 11. test_parameter_psd_and_attributes
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_parameter_psd_and_attributes
test_that("Parameter with PSD + nonneg accepts valid matrices (CVXPY parity)", {
  ## PSD + nonneg: matrix must be PSD AND all entries nonneg
  P <- Parameter(c(2, 2), PSD = TRUE, nonneg = TRUE)
  expect_true(P@attributes$PSD)
  expect_true(P@attributes$nonneg)

  ## Valid: PSD and all entries nonneg
  val <- matrix(c(2, 1, 1, 2), 2, 2)
  value(P) <- val
  expect_true(!is.null(value(P)))
})

## @cvxpy test_attributes.py::TestMultipleAttributes::test_parameter_psd_and_attributes
test_that("Parameter with PSD used in quad_form (CVXPY parity)", {
  ## PSD parameter in quadratic form: x^T P x with x >= 1
  P <- Parameter(c(2, 2), PSD = TRUE)
  value(P) <- diag(2)
  x <- Variable(2)
  cost <- quad_form(x, P)
  prob <- Problem(Minimize(cost), list(x >= 1))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  ## x = [1, 1], cost = 1^2 + 1^2 = 2
  expect_equal(as.numeric(value(prob)), 2, tolerance = 1e-3)

  ## Change P to 2*I -> cost doubles
  value(P) <- 2 * diag(2)
  x2 <- Variable(2)
  cost2 <- quad_form(x2, P)
  prob2 <- Problem(Minimize(cost2), list(x2 >= 1))
  psolve(prob2, solver = "SCS")
  expect_equal(status(prob2), "optimal")
  ## x = [1, 1], cost = 2*(1^2 + 1^2) = 4
  expect_equal(as.numeric(value(prob2)), 4, tolerance = 1e-3)
})

# ====================================================================
# 12. test_bool_int_variable
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_bool_int_variable
test_that("Variable with both boolean and integer (CVXPY parity)", {
  ## Boolean is a subset of integer, so both should work.
  ## Neither boolean nor integer is in CONVEX_ATTRIBUTES, so
  ## CvxAttr2Constr does NOT strip them. Both are preserved.
  x <- Variable(boolean = TRUE, integer = TRUE)
  expect_true(x@attributes$boolean)
  expect_true(x@attributes$integer)

  ## Minimize -> 0 (boolean forces {0, 1})
  prob <- Problem(Minimize(x))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(round(as.numeric(value(x))), 0)

  ## Maximize -> 1
  x2 <- Variable(boolean = TRUE, integer = TRUE)
  prob2 <- Problem(Maximize(x2))
  psolve(prob2, solver = "GLPK_MI")
  expect_equal(status(prob2), "optimal")
  expect_equal(round(as.numeric(value(x2))), 1)

  ## With constraint x >= 0.5 -> x = 1
  x3 <- Variable(boolean = TRUE, integer = TRUE)
  prob3 <- Problem(Minimize(x3), list(x3 >= 0.5))
  psolve(prob3, solver = "GLPK_MI")
  expect_equal(status(prob3), "optimal")
  expect_equal(round(as.numeric(value(x3))), 1)
})

# ====================================================================
# Additional integration tests
# ====================================================================

## @cvxpy test_attributes.py::TestMultipleAttributes::test_integer_bounds test_attributes.py::TestMultipleAttributes::test_nonneg_PSD
test_that("Nonneg variable with bounds: both constraints apply (CVXPY parity)", {
  ## nonneg + bounds: both constraints should apply.
  ## In this case nonneg is redundant since lower bound = 2 > 0.
  ## BUG-A technically applies but is invisible because nonneg is redundant.
  x <- Variable(nonneg = TRUE, bounds = list(2, 5))
  prob <- Problem(Minimize(x))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 2, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Boolean variable in mixed-integer LP (CVXPY parity)", {
  ## Boolean + integer variables as separate variables in same problem
  ## (not combining attributes on a single variable)
  x <- Variable(2)
  boolvar <- Variable(boolean = TRUE)
  intvar <- Variable(integer = TRUE)
  prob <- Problem(
    Minimize(-4 * x[1] - 5 * x[2]),
    list(
      2 * x[1] + x[2] <= intvar,
      x[1] + 2 * x[2] <= 3 * boolvar,
      x >= 0,
      intvar == 3 * boolvar,
      intvar == 3
    )
  )
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(prob)), -9, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-3)
  expect_equal(round(as.numeric(value(boolvar))), 1)
  expect_equal(round(as.numeric(value(intvar))), 3)
})

## @cvxpy NONE
test_that("Variable attributes are stored correctly in @attributes (CVXPY parity)", {
  ## Test that all attributes are properly stored
  x1 <- Variable(nonneg = TRUE)
  expect_true(x1@attributes$nonneg)
  expect_false(x1@attributes$nonpos)
  expect_false(x1@attributes$boolean)
  expect_false(x1@attributes$integer)

  x2 <- Variable(nonpos = TRUE)
  expect_true(x2@attributes$nonpos)
  expect_false(x2@attributes$nonneg)

  x3 <- Variable(boolean = TRUE)
  expect_true(x3@attributes$boolean)

  x4 <- Variable(integer = TRUE)
  expect_true(x4@attributes$integer)

  X5 <- Variable(c(3, 3), PSD = TRUE)
  expect_true(X5@attributes$PSD)
  expect_false(X5@attributes$NSD)

  X6 <- Variable(c(3, 3), NSD = TRUE)
  expect_true(X6@attributes$NSD)
  expect_false(X6@attributes$PSD)

  X7 <- Variable(c(3, 3), symmetric = TRUE)
  expect_true(X7@attributes$symmetric)

  ## PSD/NSD/symmetric require square matrix
  expect_error(Variable(c(2, 3), PSD = TRUE), "square matrix")
  expect_error(Variable(c(2, 3), NSD = TRUE), "square matrix")
  expect_error(Variable(c(2, 3), symmetric = TRUE), "square matrix")
})

## @cvxpy NONE
test_that("Variable sign queries reflect attributes (CVXPY parity)", {
  ## nonneg -> is_nonneg TRUE
  x_nn <- Variable(nonneg = TRUE)
  expect_true(is_nonneg(x_nn))
  expect_false(is_nonpos(x_nn))

  ## nonpos -> is_nonpos TRUE
  x_np <- Variable(nonpos = TRUE)
  expect_true(is_nonpos(x_np))
  expect_false(is_nonneg(x_np))

  ## boolean implies nonneg (CVXPY leaf.py lines 262-269)
  x_bool <- Variable(boolean = TRUE)
  expect_true(is_nonneg(x_bool))

  ## plain variable: neither nonneg nor nonpos
  x_plain <- Variable()
  expect_false(is_nonneg(x_plain))
  expect_false(is_nonpos(x_plain))
})

## @cvxpy NONE
test_that("Variable PSD/NSD queries reflect attributes (CVXPY parity)", {
  X_psd <- Variable(c(2, 2), PSD = TRUE)
  expect_true(is_psd(X_psd))
  expect_false(is_nsd(X_psd))
  expect_true(is_symmetric(X_psd))

  X_nsd <- Variable(c(2, 2), NSD = TRUE)
  expect_true(is_nsd(X_nsd))
  expect_false(is_psd(X_nsd))
  expect_true(is_symmetric(X_nsd))
})

## @cvxpy NONE
test_that("Parameter attributes stored correctly (CVXPY parity)", {
  p_nn <- Parameter(nonneg = TRUE)
  expect_true(p_nn@attributes$nonneg)

  p_np <- Parameter(nonpos = TRUE)
  expect_true(p_np@attributes$nonpos)

  P_psd <- Parameter(c(2, 2), PSD = TRUE)
  expect_true(P_psd@attributes$PSD)

  P_nsd <- Parameter(c(2, 2), NSD = TRUE)
  expect_true(P_nsd@attributes$NSD)

  P_sym <- Parameter(c(3, 3), symmetric = TRUE)
  expect_true(P_sym@attributes$symmetric)
})

## @cvxpy test_attributes.py::TestAttributes::test_variable_bounds
test_that("Bounds with QP solver OSQP (CVXPY parity)", {
  ## Bounds should work correctly with QP solvers
  x <- Variable(3, bounds = list(c(1, 2, 3), c(4, 5, 6)))
  Q <- diag(c(1, 2, 3))
  c_vec <- c(1, 1, 1)
  prob <- Problem(Minimize(quad_form(x, Q) + t(c_vec) %*% x))
  psolve(prob, solver = "OSQP")
  expect_equal(status(prob), "optimal")
  ## The optimal is at x_i = max(lower_i, -c_i/(2*Q_ii))
  ## x1: max(1, -1/2) = 1, x2: max(2, -1/4) = 2, x3: max(3, -1/6) = 3
  ## All at lower bounds since the unconstrained minimum is negative
  xval <- as.numeric(value(x))
  expect_true(all(xval >= c(1, 2, 3) - 1e-3))
  expect_true(all(xval <= c(4, 5, 6) + 1e-3))
})

## @cvxpy NONE
test_that("Nonneg variable solving (CVXPY parity)", {
  ## Basic nonneg: minimize x subject to x >= 0 (implicit) and x >= 5
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(x), list(x >= 5))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 5, tolerance = 1e-3)

  ## Nonneg prevents going below 0 even without explicit constraint
  x2 <- Variable(nonneg = TRUE)
  prob2 <- Problem(Minimize(x2))
  psolve(prob2, solver = "SCS")
  expect_equal(status(prob2), "optimal")
  expect_equal(as.numeric(value(x2)), 0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Nonpos variable solving (CVXPY parity)", {
  ## Basic nonpos: maximize x subject to x <= 0 (implicit) and x <= -5
  x <- Variable(nonpos = TRUE)
  prob <- Problem(Maximize(x), list(x <= -5))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), -5, tolerance = 1e-3)

  ## Nonpos prevents going above 0 even without explicit constraint
  x2 <- Variable(nonpos = TRUE)
  prob2 <- Problem(Maximize(x2))
  psolve(prob2, solver = "SCS")
  expect_equal(status(prob2), "optimal")
  expect_equal(as.numeric(value(x2)), 0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("PSD variable solving (CVXPY parity)", {
  ## PSD variable: minimize trace with diagonal constrained >= 1
  X <- Variable(c(2, 2), PSD = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(X[1, 1] >= 1, X[2, 2] >= 1))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(prob)), 2, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("NSD variable solving (CVXPY parity)", {
  ## NSD variable: maximize trace with diagonal constrained <= -1
  X <- Variable(c(2, 2), NSD = TRUE)
  prob <- Problem(Maximize(matrix_trace(X)),
                  list(X[1, 1] <= -1, X[2, 2] <= -1))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(prob)), -2, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Symmetric variable solving (CVXPY parity)", {
  ## Symmetric variable: X = X^T is enforced
  X <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(X[1, 1] >= 1, X[2, 2] >= 1, X >= -10))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  Xval <- value(X)
  ## Check symmetry
  expect_equal(Xval[1, 2], Xval[2, 1], tolerance = 1e-6)
})

# ═══════════════════════════════════════════════════════════════════════
# Scalar boolean variable with nonpos variable (test_attributes.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_attributes.py::TestAttributes::test_scalar_bool
test_that("Scalar boolean + nonpos variable optimization", {
  ## CVXPY: x = Variable(nonpos=True), n = Variable(boolean=True)
  ## Maximize x, subject to n == 1
  ## x is nonpos -> max x = 0; n is constrained to 1.
  ## Verified via CVXPY: status=optimal, x=0.0, n=1.0
  x <- Variable(nonpos = TRUE)
  n <- Variable(boolean = TRUE)
  prob <- Problem(Maximize(x), list(n == 1))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 0, tolerance = 1e-4)
  expect_equal(as.numeric(value(n)), 1, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# Scalar integer variable with nonpos variable (test_attributes.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_attributes.py::TestAttributes::test_scalar_int
test_that("Scalar integer + nonpos variable optimization", {
  ## CVXPY: x = Variable(nonpos=True), n = Variable(integer=True)
  ## Maximize x, subject to n == 1
  ## x is nonpos -> max x = 0; n is constrained to 1.
  ## Verified via CVXPY: status=optimal, x=0.0, n=1.0
  x <- Variable(nonpos = TRUE)
  n <- Variable(integer = TRUE)
  prob <- Problem(Maximize(x), list(n == 1))
  psolve(prob, solver = "GLPK_MI")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 0, tolerance = 1e-4)
  expect_equal(as.numeric(value(n)), 1, tolerance = 1e-4)
})
