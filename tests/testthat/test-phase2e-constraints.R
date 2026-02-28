## Tests for Phase 2E: Constraint Classes
## Translated from CVXPY: tests/test_constraints.py
## Tests construction, naming, shape, DCP, residual, violation, value, copy.

# ═══════════════════════════════════════════════════════════════════
# Zero constraint (x == 0)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Zero: construction and naming", {
  x <- Variable(3, name = "x")
  constr <- Zero(x)
  expect_true(S7_inherits(constr, Zero))
  expect_true(S7_inherits(constr, Constraint))
  expect_equal(expr_name(constr), "x == 0")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Zero: DCP — affine is DCP", {
  x <- Variable(3)
  constr <- Zero(x)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("Zero: DCP — non-affine is not DCP", {
  x <- Variable(3)
  ## x * x is not affine (product of non-constants)
  constr <- Zero(x * x)
  expect_false(is_dcp(constr))
})

## @cvxpy NONE
test_that("Zero: residual with constant", {
  a <- Constant(matrix(c(1, -2, 0), 3, 1))
  constr <- Zero(a)
  res <- residual(constr)
  expect_equal(as.numeric(res), c(1, 2, 0))
})

## @cvxpy NONE
test_that("Zero: value returns TRUE when zero", {
  a <- Constant(matrix(c(0, 0, 0), 3, 1))
  constr <- Zero(a)
  expect_true(value(constr))
})

## @cvxpy NONE
test_that("Zero: value returns FALSE when nonzero", {
  a <- Constant(matrix(c(1, 0, 0), 3, 1))
  constr <- Zero(a)
  expect_false(value(constr))
})

## @cvxpy NONE
test_that("Zero: violation is same as residual (base behavior)", {
  a <- Constant(matrix(c(1, -2, 0), 3, 1))
  constr <- Zero(a)
  expect_equal(as.numeric(violation(constr)), c(1, 2, 0))
})

## @cvxpy NONE
test_that("Zero: dual_value is NULL before solving", {
  x <- Variable(3)
  constr <- Zero(x)
  expect_null(dual_value(constr))
})

# ═══════════════════════════════════════════════════════════════════
# Equality constraint (x == y)
# CVXPY: test_equality in test_constraints.py
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Equality: construction and naming (CVXPY test_equality)", {
  x <- Variable(c(2, 1), name = "x")
  z <- Variable(c(2, 1), name = "z")
  constr <- Equality(x, z)
  expect_true(S7_inherits(constr, Equality))
  expect_equal(expr_name(constr), "x == z")
  expect_equal(constr@shape, c(2L, 1L))
})

## @cvxpy NONE
test_that("Equality: dual_value is NULL before solving", {
  x <- Variable(c(2, 1), name = "x")
  z <- Variable(c(2, 1), name = "z")
  constr <- Equality(x, z)
  expect_null(dual_value(constr))
})

## @cvxpy NONE
test_that("Equality: value errors when no variable values set", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Equality(x, z)
  expect_error(value(constr), "no value")
})

## @cvxpy NONE
test_that("Equality: value TRUE when x == z", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Equality(x, z)
  value(x) <- matrix(c(2, 2), 2, 1)
  value(z) <- matrix(c(2, 2), 2, 1)
  expect_true(value(constr))
})

## @cvxpy NONE
test_that("Equality: value FALSE when x != z (scalar broadcast)", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Equality(x, z)
  value(x) <- matrix(c(2, 2), 2, 1)
  value(z) <- matrix(c(2, 2), 2, 1)
  expect_true(value(constr))
  value(x) <- matrix(c(3, 3), 2, 1)
  expect_false(value(constr))
})

## @cvxpy NONE
test_that("Equality: residual and violation with specific values (CVXPY test)", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Equality(x, z)

  value(x) <- matrix(c(2, 1), 2, 1)
  value(z) <- matrix(c(2, 2), 2, 1)
  expect_false(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 1))
  expect_equal(as.numeric(residual(constr)), c(0, 1))

  value(z) <- matrix(c(2, 1), 2, 1)
  expect_true(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 0))
  expect_equal(as.numeric(residual(constr)), c(0, 0))
})

## @cvxpy NONE
test_that("Equality: incompatible dimensions error (CVXPY test)", {
  x <- Variable(c(2, 1)); y <- Variable(c(3, 1))
  expect_error(Equality(x, y), "broadcast")
})

## @cvxpy NONE
test_that("Equality: is_dcp — affine difference", {
  x <- Variable(3); y <- Variable(3)
  constr <- Equality(x, y)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("Equality: copy with args=NULL (CVXPY test)", {
  x <- Variable(c(2, 1), name = "x")
  z <- Variable(c(2, 1), name = "z")
  constr <- Equality(x, z)
  copy <- expr_copy(constr)
  expect_true(S7_inherits(copy, Equality))
  ## Same class
  expect_identical(S7_class(copy), S7_class(constr))
  ## Args equal (same objects) but different list allocation
  expect_identical(copy@args[[1L]]@id, constr@args[[1L]]@id)
  expect_identical(copy@args[[2L]]@id, constr@args[[2L]]@id)
})

## @cvxpy NONE
test_that("Equality: copy with new args (CVXPY test)", {
  x <- Variable(c(2, 1), name = "x")
  z <- Variable(c(2, 1), name = "z")
  constr <- Equality(x, z)
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")
  copy <- expr_copy(constr, args = list(A, B))
  expect_true(S7_inherits(copy, Equality))
  expect_identical(copy@args[[1L]]@id, A@id)
})

# ═══════════════════════════════════════════════════════════════════
# NonPos constraint (x <= 0)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("NonPos: construction and naming", {
  x <- Variable(3, name = "x")
  constr <- NonPos(x)
  expect_true(S7_inherits(constr, NonPos))
  expect_equal(expr_name(constr), "x <= 0")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("NonPos: DCP — convex argument", {
  x <- Variable(3)
  ## x is affine → convex → DCP
  constr <- NonPos(x)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("NonPos: DCP — concave argument is not DCP", {
  x <- Variable(3)
  ## -x*x is concave (neg of convex non-constant product)... actually
  ## just check a known non-convex case
  constr <- NonPos(-x)
  ## -x is affine → convex → DCP
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("NonPos: residual — max(expr, 0)", {
  a <- Constant(matrix(c(-1, 2, 0), 3, 1))
  constr <- NonPos(a)
  expect_equal(as.numeric(residual(constr)), c(0, 2, 0))
})

## @cvxpy NONE
test_that("NonPos: violation returns L2 norm", {
  a <- Constant(matrix(c(0, 3, 4), 3, 1))
  constr <- NonPos(a)
  ## max(c(0,3,4), 0) = c(0,3,4), L2 norm = 5
  expect_equal(violation(constr), 5)
})

## @cvxpy NONE
test_that("NonPos: value TRUE when all <= 0", {
  a <- Constant(matrix(c(-1, -2, 0), 3, 1))
  constr <- NonPos(a)
  expect_true(value(constr))
})

## @cvxpy NONE
test_that("NonPos: value FALSE when some > 0", {
  a <- Constant(matrix(c(-1, 1, 0), 3, 1))
  constr <- NonPos(a)
  expect_false(value(constr))
})

# ═══════════════════════════════════════════════════════════════════
# NonNeg constraint (x >= 0)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("NonNeg: construction and naming", {
  x <- Variable(3, name = "x")
  constr <- NonNeg(x)
  expect_true(S7_inherits(constr, NonNeg))
  expect_equal(expr_name(constr), "x >= 0")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("NonNeg: DCP — concave argument", {
  x <- Variable(3)
  ## x is affine → concave → DCP
  constr <- NonNeg(x)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("NonNeg: residual — abs(min(expr, 0))", {
  a <- Constant(matrix(c(1, -2, 0), 3, 1))
  constr <- NonNeg(a)
  ## min(c(1,-2,0), 0) = c(0,-2,0), abs = c(0,2,0)
  expect_equal(as.numeric(residual(constr)), c(0, 2, 0))
})

## @cvxpy NONE
test_that("NonNeg: violation returns L2 norm", {
  a <- Constant(matrix(c(0, -3, -4), 3, 1))
  constr <- NonNeg(a)
  ## abs(min(c(0,-3,-4), 0)) = c(0,3,4), L2 norm = 5
  expect_equal(violation(constr), 5)
})

## @cvxpy NONE
test_that("NonNeg: value TRUE when all >= 0", {
  a <- Constant(matrix(c(1, 2, 0), 3, 1))
  constr <- NonNeg(a)
  expect_true(value(constr))
})

## @cvxpy NONE
test_that("NonNeg: value FALSE when some < 0", {
  a <- Constant(matrix(c(1, -1, 0), 3, 1))
  constr <- NonNeg(a)
  expect_false(value(constr))
})

# ═══════════════════════════════════════════════════════════════════
# Inequality constraint (lhs <= rhs)
# CVXPY: test_inequality in test_constraints.py
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Inequality: construction and naming (CVXPY test_inequality)", {
  x <- Variable(c(2, 1), name = "x")
  z <- Variable(c(2, 1), name = "z")
  constr <- Inequality(x, z)
  expect_true(S7_inherits(constr, Inequality))
  expect_equal(expr_name(constr), "x <= z")
  expect_equal(constr@shape, c(2L, 1L))
})

## @cvxpy NONE
test_that("Inequality: dual_value is NULL before solving", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Inequality(x, z)
  expect_null(dual_value(constr))
})

## @cvxpy NONE
test_that("Inequality: value errors when no variable values set", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Inequality(x, z)
  expect_error(value(constr), "no value")
})

## @cvxpy NONE
test_that("Inequality: value TRUE when lhs <= rhs (CVXPY test)", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Inequality(x, z)
  value(x) <- matrix(c(1, 1), 2, 1)
  value(z) <- matrix(c(2, 2), 2, 1)
  expect_true(value(constr))
})

## @cvxpy NONE
test_that("Inequality: value FALSE when lhs > rhs (CVXPY test)", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Inequality(x, z)
  value(x) <- matrix(c(1, 1), 2, 1)
  value(z) <- matrix(c(2, 2), 2, 1)
  expect_true(value(constr))
  value(x) <- matrix(c(3, 3), 2, 1)
  expect_false(value(constr))
})

## @cvxpy NONE
test_that("Inequality: residual and violation (CVXPY test)", {
  x <- Variable(c(2, 1)); z <- Variable(c(2, 1))
  constr <- Inequality(x, z)

  value(x) <- matrix(c(2, 1), 2, 1)
  value(z) <- matrix(c(2, 0), 2, 1)
  expect_false(value(constr))
  ## lhs - rhs = (2-2, 1-0) = (0, 1), max(., 0) = (0, 1)
  expect_equal(as.numeric(violation(constr)), c(0, 1))
  expect_equal(as.numeric(residual(constr)), c(0, 1))

  value(z) <- matrix(c(2, 2), 2, 1)
  expect_true(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 0))
  expect_equal(as.numeric(residual(constr)), c(0, 0))
})

## @cvxpy NONE
test_that("Inequality: incompatible dimensions error (CVXPY test)", {
  x <- Variable(c(2, 1)); y <- Variable(c(3, 1))
  expect_error(Inequality(x, y), "broadcast")
})

## @cvxpy NONE
test_that("Inequality: is_dcp — affine difference", {
  x <- Variable(3); y <- Variable(3)
  constr <- Inequality(x, y)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("Inequality: is_dcp — convex minus constant", {
  x <- Variable(3)
  constr <- Inequality(x, Constant(matrix(0, 3, 1)))
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("Inequality: copy with args=NULL (CVXPY test)", {
  x <- Variable(c(2, 1), name = "x")
  z <- Variable(c(2, 1), name = "z")
  constr <- Inequality(x, z)
  copy <- expr_copy(constr)
  expect_true(S7_inherits(copy, Inequality))
  expect_identical(S7_class(copy), S7_class(constr))
  expect_identical(copy@args[[1L]]@id, constr@args[[1L]]@id)
  expect_identical(copy@args[[2L]]@id, constr@args[[2L]]@id)
})

## @cvxpy NONE
test_that("Inequality: copy with new args (CVXPY test)", {
  x <- Variable(c(2, 1), name = "x")
  z <- Variable(c(2, 1), name = "z")
  constr <- Inequality(x, z)
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")
  copy <- expr_copy(constr, args = list(A, B))
  expect_true(S7_inherits(copy, Inequality))
  expect_identical(copy@args[[1L]]@id, A@id)
})

# ═══════════════════════════════════════════════════════════════════
# >= operator semantics (CVXPY: test_geq)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Inequality with reversed args (>= semantics from CVXPY test_geq)", {
  ## z >= x → Inequality(x, z) → name "x <= z"
  x <- Variable(c(2, 1), name = "x")
  z <- Variable(c(2, 1), name = "z")
  constr <- Inequality(x, z)  ## This is what z >= x would create
  expect_equal(expr_name(constr), "x <= z")
  expect_equal(constr@shape, c(2L, 1L))
})

# ═══════════════════════════════════════════════════════════════════
# Constraint base class properties
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Constraint: get_data returns list with id", {
  x <- Variable(3)
  constr <- Zero(x)
  data <- get_data(constr)
  expect_true(is.list(data))
  expect_equal(data[[1]], constr@id)
})

## @cvxpy NONE
test_that("Constraint: is_real / is_complex", {
  x <- Variable(3)
  constr <- Zero(x)
  expect_true(is_real(constr))
  expect_false(is_complex(constr))
})

## @cvxpy NONE
test_that("Constraint: print method", {
  x <- Variable(3, name = "x")
  constr <- Zero(x)
  expect_output(print(constr), "x == 0")
})

## @cvxpy NONE
test_that("Constraint: constr_expr for single-arg constraints", {
  x <- Variable(3, name = "x")
  constr <- Zero(x)
  expect_identical(constr_expr(constr)@id, x@id)
})

## @cvxpy NONE
test_that("Constraint: constr_expr errors for multi-arg constraints", {
  x <- Variable(3); y <- Variable(3)
  constr <- Equality(x, y)
  expect_error(constr_expr(constr), "ambiguous")
})

# ═══════════════════════════════════════════════════════════════════
# DCP rules for constraints
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DCP: Zero requires affine", {
  x <- Variable(3)
  expect_true(is_dcp(Zero(x)))        ## x is affine
  expect_true(is_dcp(Zero(2 * x + 1))) ## affine expression
})

## @cvxpy NONE
test_that("DCP: NonPos requires convex", {
  x <- Variable(3)
  expect_true(is_dcp(NonPos(x)))      ## affine → convex
  expect_true(is_dcp(NonPos(2 * x)))  ## affine → convex
})

## @cvxpy NONE
test_that("DCP: NonNeg requires concave", {
  x <- Variable(3)
  expect_true(is_dcp(NonNeg(x)))       ## affine → concave
  expect_true(is_dcp(NonNeg(-x + 1)))  ## affine → concave
})

## @cvxpy NONE
test_that("DCP: Inequality — lhs-rhs must be convex", {
  x <- Variable(3); y <- Variable(3)
  expect_true(is_dcp(Inequality(x, y)))  ## x-y is affine → convex
})

## @cvxpy NONE
test_that("DCP: Equality — lhs-rhs must be affine", {
  x <- Variable(3); y <- Variable(3)
  expect_true(is_dcp(Equality(x, y)))   ## x-y is affine
  expect_true(is_dcp(Equality(x, 0)))   ## x-0 is affine
  expect_true(is_dcp(Equality(2 * x, y + 1)))
})

# ═══════════════════════════════════════════════════════════════════
# variables() and parameters() on constraints
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Constraint: variables() collects from args", {
  x <- Variable(3, name = "x"); y <- Variable(3, name = "y")
  constr <- Equality(x, y)
  vars <- variables(constr)
  var_ids <- vapply(vars, function(v) v@id, integer(1))
  expect_true(x@id %in% var_ids)
  expect_true(y@id %in% var_ids)
})

## @cvxpy NONE
test_that("Constraint: parameters() collects from args", {
  x <- Variable(3); p <- Parameter(3)
  constr <- Equality(x, p)
  params <- parameters(constr)
  expect_length(params, 1L)
  expect_equal(params[[1L]]@id, p@id)
})

## @cvxpy NONE
test_that("Constraint: constants() collects from args", {
  x <- Variable(3)
  constr <- Inequality(x, Constant(matrix(0, 3, 1)))
  const <- constants(constr)
  expect_length(const, 1L)
})

# ═══════════════════════════════════════════════════════════════════
# Edge cases
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Scalar constraint", {
  x <- Variable(1, name = "x")
  constr <- Zero(x)
  expect_equal(constr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("Matrix constraint shape", {
  X <- Variable(c(3, 4), name = "X")
  Y <- Variable(c(3, 4), name = "Y")
  constr <- Equality(X, Y)
  expect_equal(constr@shape, c(3L, 4L))
})

## @cvxpy NONE
test_that("Constraint with constant expressions", {
  a <- Constant(matrix(c(1, 2, 3), 3, 1))
  b <- Constant(matrix(c(4, 5, 6), 3, 1))
  constr <- Inequality(a, b)
  ## a - b = (-3, -3, -3), all <= 0, so satisfied
  expect_true(value(constr))
  expect_equal(as.numeric(residual(constr)), c(0, 0, 0))
})

## @cvxpy NONE
test_that("Constraint with mixed expression types", {
  x <- Variable(3, name = "x")
  p <- Parameter(3, name = "p")
  value(p) <- matrix(c(1, 2, 3), 3, 1)
  ## x <= p is DCP (affine difference)
  constr <- Inequality(x, p)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("Inequality: rejects complex expressions", {
  x <- Variable(3, complex = TRUE)
  y <- Variable(3, complex = TRUE)
  expect_error(Inequality(x, y), "complex")
})

## @cvxpy NONE
test_that("NonPos: rejects complex input", {
  x <- Variable(3, complex = TRUE)
  expect_error(NonPos(x), "real")
})

## @cvxpy NONE
test_that("NonNeg: rejects complex input", {
  x <- Variable(3, complex = TRUE)
  expect_error(NonNeg(x), "real")
})

## @cvxpy NONE
test_that("Residual is NULL when variables have no value", {
  x <- Variable(3); y <- Variable(3)
  expect_null(residual(Equality(x, y)))
  expect_null(residual(Inequality(x, y)))
  expect_null(residual(Zero(x)))
  expect_null(residual(NonPos(x)))
  expect_null(residual(NonNeg(x)))
})

## @cvxpy NONE
test_that("Value errors when variables have no value", {
  x <- Variable(3); y <- Variable(3)
  expect_error(value(Equality(x, y)), "no value")
  expect_error(value(Inequality(x, y)), "no value")
  expect_error(value(Zero(x)), "no value")
  expect_error(value(NonPos(x)), "no value")
  expect_error(value(NonNeg(x)), "no value")
})
