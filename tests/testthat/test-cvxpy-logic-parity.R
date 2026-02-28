## CVXPY Logic Atoms Parity Tests
## ================================
## These tests mirror CVXPY's test_logic.py to verify that CVXR logic atoms
## (Not, And, Or, Xor, implies, iff) match CVXPY behavior exactly.
##
## Source: /Users/naras/GitHub/cvxpy/cvxpy/tests/test_logic.py
## R implementation: rsrc_tree/atoms/elementwise/logic.R
## R canonicalizers: rsrc_tree/reductions/dcp2cone/canonicalizers/logic_canon.R
##
## Solver: HIGHS (MIP-capable, required for boolean variables; matches CVXPY)
##
## Key differences from Python:
##   - R uses `!` for Not (Python uses `~`)
##   - R uses `Xor()` functional only (Python `^` is power in R)
##   - R names: "!x" (Python: "~x"), " XOR " separator (Python: " ^ ")
##   - R Variable(boolean=TRUE) shape is (1,1); Python is ()

library(testthat)
library(CVXR)

# ── Helper: solve feasibility and return value of expr ──────────────
.eval_logic <- function(expr, constraints = list()) {
  y <- Variable(expr@shape)
  prob <- Problem(Minimize(Constant(0)), c(list(y == expr), constraints))
  psolve(prob, solver = "HIGHS")
  value(y)
}

# ═══════════════════════════════════════════════════════════════════════
# 1. Construction tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicProperties::test_expression_types
test_that("Not: construction with boolean variable", {
  x <- Variable(boolean = TRUE)
  expr <- Not(x)
  expect_true(S7_inherits(expr, Not))
  expect_true(S7_inherits(expr, LogicExpression))
  expect_equal(length(expr@args), 1L)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy test_logic.py::TestLogicProperties::test_expression_types
test_that("And: construction with two boolean variables", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- And(x, y)
  expect_true(S7_inherits(expr, And))
  expect_true(S7_inherits(expr, NaryLogicExpression))
  expect_equal(length(expr@args), 2L)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy test_logic.py::TestLogicProperties::test_expression_types
test_that("Or: construction with two boolean variables", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- Or(x, y)
  expect_true(S7_inherits(expr, Or))
  expect_true(S7_inherits(expr, NaryLogicExpression))
  expect_equal(length(expr@args), 2L)
})

## @cvxpy test_logic.py::TestLogicProperties::test_expression_types
test_that("Xor: construction with two boolean variables", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- Xor(x, y)
  expect_true(S7_inherits(expr, Xor))
  expect_true(S7_inherits(expr, NaryLogicExpression))
  expect_equal(length(expr@args), 2L)
})

## @cvxpy NONE
test_that("LogicExpression on vector boolean variables", {
  x <- Variable(3, boolean = TRUE)
  y <- Variable(3, boolean = TRUE)
  expect_equal(And(x, y)@shape, c(3L, 1L))
  expect_equal(Or(x, y)@shape, c(3L, 1L))
  expect_equal(Not(x)@shape, c(3L, 1L))
  expect_equal(Xor(x, y)@shape, c(3L, 1L))
})

# ═══════════════════════════════════════════════════════════════════════
# 2. Validation tests: error on non-boolean arguments
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicValidation::test_non_boolean_raises
test_that("And rejects non-boolean variable", {
  x <- Variable()
  y <- Variable(boolean = TRUE)
  expect_error(And(x, y), "boolean")
})

## @cvxpy test_logic.py::TestLogicValidation::test_non_boolean_raises test_logic.py::TestLogicOperators::test_non_boolean_raises
test_that("Or rejects non-boolean variable", {
  x <- Variable()
  y <- Variable(boolean = TRUE)
  expect_error(Or(x, y), "boolean")
})

## @cvxpy test_logic.py::TestLogicOperators::test_non_boolean_raises
test_that("Xor rejects non-boolean variable", {
  x <- Variable()
  y <- Variable(boolean = TRUE)
  expect_error(Xor(x, y), "boolean")
})

## @cvxpy test_logic.py::TestLogicValidation::test_non_boolean_raises test_logic.py::TestLogicOperators::test_non_boolean_raises
test_that("Not rejects non-boolean variable", {
  x <- Variable()
  expect_error(Not(x), "boolean")
})

## @cvxpy test_logic.py::TestLogicValidation::test_too_few_args
test_that("And requires at least 2 arguments", {
  x <- Variable(boolean = TRUE)
  expect_error(And(x), "2 arguments")
})

## @cvxpy test_logic.py::TestLogicValidation::test_too_few_args
test_that("Or requires at least 2 arguments", {
  x <- Variable(boolean = TRUE)
  expect_error(Or(x), "2 arguments")
})

## @cvxpy test_logic.py::TestLogicValidation::test_too_few_args
test_that("Xor requires at least 2 arguments", {
  x <- Variable(boolean = TRUE)
  expect_error(Xor(x), "2 arguments")
})

## @cvxpy NONE
test_that("Numeric constant 0/1 is accepted as boolean-valued in R", {
  ## In R's CVXR, Constant(1) has value 1 (which is 0 or 1),
  ## so .is_boolean_arg() accepts it. This differs from CVXPY which checks
  ## is_boolean_valued on the Constant class. We test the R behavior.
  x <- Variable(boolean = TRUE)
  expect_no_error(And(x, Constant(TRUE)))
  expect_no_error(Or(x, Constant(FALSE)))
})

## @cvxpy test_logic.py::TestLogicValidation::test_numeric_constant_raises
test_that("Non-boolean numeric constant (e.g. 2) is rejected", {
  x <- Variable(boolean = TRUE)
  expect_error(And(x, Constant(2)), "boolean")
  expect_error(Or(x, Constant(0.5)), "boolean")
})

# ═══════════════════════════════════════════════════════════════════════
# 3. Operator syntax: !, &, | dispatch
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicOperators::test_operators
test_that("! operator creates Not expression", {
  x <- Variable(boolean = TRUE)
  expr <- !x
  expect_true(S7_inherits(expr, Not))
})

## @cvxpy test_logic.py::TestLogicOperators::test_operators
test_that("& operator creates And expression", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- x & y
  expect_true(S7_inherits(expr, And))
})

## @cvxpy test_logic.py::TestLogicOperators::test_operators
test_that("| operator creates Or expression", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- x | y
  expect_true(S7_inherits(expr, Or))
})

## @cvxpy test_logic.py::TestLogicOperators::test_non_boolean_raises
test_that("! on non-boolean raises error", {
  x <- Variable()
  expect_error(!x, "boolean")
})

## @cvxpy test_logic.py::TestLogicOperators::test_non_boolean_raises
test_that("& on non-boolean raises error", {
  x <- Variable()
  y <- Variable(boolean = TRUE)
  expect_error(x & y, "boolean")
})

## @cvxpy test_logic.py::TestLogicOperators::test_non_boolean_raises
test_that("| on non-boolean raises error", {
  x <- Variable()
  y <- Variable(boolean = TRUE)
  expect_error(x | y, "boolean")
})

# ═══════════════════════════════════════════════════════════════════════
# 4. DCP properties: curvature, sign, monotonicity
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicProperties::test_dcp_compliance
test_that("All logic atoms are DCP-compliant", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expect_true(is_dcp(Not(x)))
  expect_true(is_dcp(And(x, y)))
  expect_true(is_dcp(Or(x, y)))
  expect_true(is_dcp(Xor(x, y)))
})

## @cvxpy NONE
test_that("Logic atoms are both convex and concave (affine-like)", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  for (expr in list(Not(x), And(x, y), Or(x, y), Xor(x, y))) {
    expect_true(is_atom_convex(expr))
    expect_true(is_atom_concave(expr))
  }
})

## @cvxpy NONE
test_that("Logic atom signs: all nonneg, not nonpos", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  for (expr in list(Not(x), And(x, y), Or(x, y), Xor(x, y))) {
    sf <- sign_from_args(expr)
    expect_true(sf$is_nonneg)
    expect_false(sf$is_nonpos)
  }
})

## @cvxpy test_logic.py::TestLogicProperties::test_monotonicity
test_that("Not is decreasing, not increasing", {
  x <- Variable(boolean = TRUE)
  expr <- Not(x)
  expect_true(is_decr(expr, 1L))
  expect_false(is_incr(expr, 1L))
})

## @cvxpy test_logic.py::TestLogicProperties::test_monotonicity
test_that("And is increasing in all arguments", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- And(x, y)
  expect_true(is_incr(expr, 1L))
  expect_true(is_incr(expr, 2L))
  expect_false(is_decr(expr, 1L))
  expect_false(is_decr(expr, 2L))
})

## @cvxpy test_logic.py::TestLogicProperties::test_monotonicity
test_that("Or is increasing in all arguments", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- Or(x, y)
  expect_true(is_incr(expr, 1L))
  expect_true(is_incr(expr, 2L))
  expect_false(is_decr(expr, 1L))
  expect_false(is_decr(expr, 2L))
})

## @cvxpy test_logic.py::TestLogicProperties::test_monotonicity
test_that("Xor is neither increasing nor decreasing", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- Xor(x, y)
  expect_false(is_incr(expr, 1L))
  expect_false(is_incr(expr, 2L))
  expect_false(is_decr(expr, 1L))
  expect_false(is_decr(expr, 2L))
})

# ═══════════════════════════════════════════════════════════════════════
# 5. Name display: expr_name() output
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicName::test_basic_names
test_that("Basic expr_name for logic atoms", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)
  ## R uses ! instead of ~, XOR instead of ^
  expect_equal(expr_name(Not(x)), "!x")
  expect_equal(expr_name(And(x, y)), "x & y")
  expect_equal(expr_name(Or(x, y)), "x | y")
  expect_match(expr_name(Xor(x, y)), "XOR")
})

## @cvxpy test_logic.py::TestLogicName::test_nary_names
test_that("N-ary names with 3 arguments", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)
  z <- Variable(name = "z", boolean = TRUE)
  expect_equal(expr_name(And(x, y, z)), "x & y & z")
  expect_equal(expr_name(Or(x, y, z)), "x | y | z")
})

## @cvxpy test_logic.py::TestLogicName::test_not_parenthesizes_binary
test_that("Not parenthesizes binary sub-expressions", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)
  expect_equal(expr_name(Not(And(x, y))), "!(x & y)")
  expect_equal(expr_name(Not(Or(x, y))), "!(x | y)")
})

## @cvxpy test_logic.py::TestLogicName::test_not_no_parens_for_leaf
test_that("Not does not parenthesize leaf variables", {
  x <- Variable(name = "x", boolean = TRUE)
  expect_equal(expr_name(Not(x)), "!x")
  expect_equal(expr_name(Not(Not(x))), "!!x")
})

## @cvxpy test_logic.py::TestLogicName::test_precedence
test_that("And parenthesizes lower-precedence Or and Xor children", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)
  z <- Variable(name = "z", boolean = TRUE)
  ## & is higher precedence than | and XOR
  expect_equal(expr_name(And(Or(x, y), z)), "(x | y) & z")
})

## @cvxpy test_logic.py::TestLogicName::test_precedence
test_that("Or has lowest precedence - never parenthesizes children", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)
  z <- Variable(name = "z", boolean = TRUE)
  expect_equal(expr_name(Or(And(x, y), z)), "x & y | z")
})

# ═══════════════════════════════════════════════════════════════════════
# 6. Numeric evaluation: truth tables via solver
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicSolve::test_not_truth_table
test_that("Not truth table: Not(0)=1, Not(1)=0", {
  for (val in c(0, 1)) {
    x <- Variable(boolean = TRUE)
    result <- .eval_logic(Not(x), list(x == val))
    expect_equal(round(drop(result)), 1 - val,
                 info = paste("Not(", val, ")"))
  }
})

## @cvxpy test_logic.py::TestLogicSolve::test_and_truth_table
test_that("And truth table: full 2-input table", {
  truth <- list(
    c(0, 0, 0), c(0, 1, 0), c(1, 0, 0), c(1, 1, 1)
  )
  for (row in truth) {
    x <- Variable(boolean = TRUE)
    y <- Variable(boolean = TRUE)
    result <- .eval_logic(And(x, y), list(x == row[1], y == row[2]))
    expect_equal(round(drop(result)), row[3],
                 info = paste("And(", row[1], ",", row[2], ")"))
  }
})

## @cvxpy test_logic.py::TestLogicSolve::test_or_truth_table
test_that("Or truth table: full 2-input table", {
  truth <- list(
    c(0, 0, 0), c(0, 1, 1), c(1, 0, 1), c(1, 1, 1)
  )
  for (row in truth) {
    x <- Variable(boolean = TRUE)
    y <- Variable(boolean = TRUE)
    result <- .eval_logic(Or(x, y), list(x == row[1], y == row[2]))
    expect_equal(round(drop(result)), row[3],
                 info = paste("Or(", row[1], ",", row[2], ")"))
  }
})

## @cvxpy test_logic.py::TestLogicSolve::test_xor_truth_table
test_that("Xor truth table: full 2-input table", {
  truth <- list(
    c(0, 0, 0), c(0, 1, 1), c(1, 0, 1), c(1, 1, 0)
  )
  for (row in truth) {
    x <- Variable(boolean = TRUE)
    y <- Variable(boolean = TRUE)
    result <- .eval_logic(Xor(x, y), list(x == row[1], y == row[2]))
    expect_equal(round(drop(result)), row[3],
                 info = paste("Xor(", row[1], ",", row[2], ")"))
  }
})

## @cvxpy test_logic.py::TestLogicSolve::test_implies_truth_table
test_that("implies truth table: full 2-input table", {
  truth <- list(
    c(0, 0, 1), c(0, 1, 1), c(1, 0, 0), c(1, 1, 1)
  )
  for (row in truth) {
    x <- Variable(boolean = TRUE)
    y <- Variable(boolean = TRUE)
    result <- .eval_logic(implies(x, y), list(x == row[1], y == row[2]))
    expect_equal(round(drop(result)), row[3],
                 info = paste("implies(", row[1], ",", row[2], ")"))
  }
})

## @cvxpy test_logic.py::TestLogicSolve::test_iff_truth_table
test_that("iff truth table: full 2-input table", {
  truth <- list(
    c(0, 0, 1), c(0, 1, 0), c(1, 0, 0), c(1, 1, 1)
  )
  for (row in truth) {
    x <- Variable(boolean = TRUE)
    y <- Variable(boolean = TRUE)
    result <- .eval_logic(iff(x, y), list(x == row[1], y == row[2]))
    expect_equal(round(drop(result)), row[3],
                 info = paste("iff(", row[1], ",", row[2], ")"))
  }
})

# ═══════════════════════════════════════════════════════════════════════
# 7. Solving tests: actual optimization with each logic atom
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicInConstraint::test_and_maximize
test_that("Maximize And(x, y) -> x=1, y=1", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Maximize(And(x, y)))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(x))), 1)
  expect_equal(round(drop(value(y))), 1)
  expect_equal(round(drop(value(prob))), 1)
})

## @cvxpy test_logic.py::TestLogicInConstraint::test_or_minimize
test_that("Minimize Or(x, y) -> x=0, y=0", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(Or(x, y)))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(x))), 0)
  expect_equal(round(drop(value(y))), 0)
  expect_equal(round(drop(value(prob))), 0)
})

## @cvxpy test_logic.py::TestLogicInConstraint::test_xor_constraint
test_that("Xor constraint: Xor(x, y)==1, x==1 -> y=0", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(Constant(0)), list(Xor(x, y) == 1, x == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(y))), 0)
})

## @cvxpy NONE
test_that("Not constraint: Not(x)==1 -> x=0", {
  x <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(Constant(0)), list(Not(x) == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(x))), 0)
})

## @cvxpy test_logic.py::TestLogicInConstraint::test_and_maximize
test_that("Maximize x with And constraint", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Maximize(x), list(And(x, y) == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(x))), 1)
  expect_equal(round(drop(value(y))), 1)
})

## @cvxpy test_logic.py::TestLogicInConstraint::test_or_minimize
test_that("Minimize Or as objective", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(Or(x, y)), list(x == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  ## Or(1, y) is minimized when y=0, but Or(1,0) = 1
  expect_equal(round(drop(value(prob))), 1)
})

# ═══════════════════════════════════════════════════════════════════════
# 8. implies/iff convenience functions
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicProperties::test_namespace
test_that("implies returns Or expression", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- implies(x, y)
  expect_true(S7_inherits(expr, Or))
})

## @cvxpy test_logic.py::TestLogicProperties::test_namespace
test_that("iff returns Not expression", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- iff(x, y)
  expect_true(S7_inherits(expr, Not))
})

## @cvxpy NONE
test_that("Maximize x+y with implies(x,y)==1", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Maximize(x + y), list(implies(x, y) == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  ## implies(x, y) = Or(Not(x), y)
  ## Maximize x+y subject to implies(x,y)==1 -> x=1, y=1 (both 1 satisfies implication)
  expect_equal(round(drop(value(x))), 1)
  expect_equal(round(drop(value(y))), 1)
})

## @cvxpy NONE
test_that("Maximize x with implies(x,y)==1, y==0 -> x=0", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Maximize(x), list(implies(x, y) == 1, y == 0))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  ## implies(x, 0) = Or(Not(x), 0) = Not(x); for this to be 1, x must be 0
  expect_equal(round(drop(value(x))), 0)
})

## @cvxpy NONE
test_that("Minimize x+y with iff(x,y)==1 -> both 0", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(x + y), list(iff(x, y) == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  ## iff(x,y) = Not(Xor(x,y)); equals 1 when x==y
  ## Minimize x+y with x==y -> x=y=0
  expect_equal(round(drop(value(x))), 0)
  expect_equal(round(drop(value(y))), 0)
  expect_equal(round(drop(value(prob))), 0)
})

## @cvxpy NONE
test_that("Maximize x+y with iff(x,y)==1 -> both 1", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Maximize(x + y), list(iff(x, y) == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  ## Maximize x+y with x==y -> x=y=1
  expect_equal(round(drop(value(x))), 1)
  expect_equal(round(drop(value(y))), 1)
  expect_equal(round(drop(value(prob))), 2)
})

# ═══════════════════════════════════════════════════════════════════════
# 9. N-ary operations: And/Or/Xor with 3+ arguments
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicSolve::test_nary
test_that("3-ary And: all-true -> 1, one-false -> 0", {
  x1 <- Variable(boolean = TRUE)
  x2 <- Variable(boolean = TRUE)
  x3 <- Variable(boolean = TRUE)

  ## And(1, 1, 1) = 1
  result <- .eval_logic(And(x1, x2, x3), list(x1 == 1, x2 == 1, x3 == 1))
  expect_equal(round(drop(result)), 1)

  ## And(1, 0, 1) = 0
  result <- .eval_logic(And(x1, x2, x3), list(x1 == 1, x2 == 0, x3 == 1))
  expect_equal(round(drop(result)), 0)
})

## @cvxpy test_logic.py::TestLogicSolve::test_nary
test_that("3-ary Or: all-false -> 0, one-true -> 1", {
  x1 <- Variable(boolean = TRUE)
  x2 <- Variable(boolean = TRUE)
  x3 <- Variable(boolean = TRUE)

  ## Or(0, 0, 0) = 0
  result <- .eval_logic(Or(x1, x2, x3), list(x1 == 0, x2 == 0, x3 == 0))
  expect_equal(round(drop(result)), 0)

  ## Or(0, 1, 0) = 1
  result <- .eval_logic(Or(x1, x2, x3), list(x1 == 0, x2 == 1, x3 == 0))
  expect_equal(round(drop(result)), 1)
})

## @cvxpy test_logic.py::TestLogicSolve::test_nary
test_that("3-ary Xor parity: 0->0, 1->1, 2->0, 3->1", {
  x1 <- Variable(boolean = TRUE)
  x2 <- Variable(boolean = TRUE)
  x3 <- Variable(boolean = TRUE)

  ## Xor(0, 0, 0) = 0 (0 ones = even parity)
  result <- .eval_logic(Xor(x1, x2, x3), list(x1 == 0, x2 == 0, x3 == 0))
  expect_equal(round(drop(result)), 0)

  ## Xor(1, 0, 0) = 1 (1 one = odd parity)
  result <- .eval_logic(Xor(x1, x2, x3), list(x1 == 1, x2 == 0, x3 == 0))
  expect_equal(round(drop(result)), 1)

  ## Xor(1, 1, 0) = 0 (2 ones = even parity)
  result <- .eval_logic(Xor(x1, x2, x3), list(x1 == 1, x2 == 1, x3 == 0))
  expect_equal(round(drop(result)), 0)

  ## Xor(1, 1, 1) = 1 (3 ones = odd parity)
  result <- .eval_logic(Xor(x1, x2, x3), list(x1 == 1, x2 == 1, x3 == 1))
  expect_equal(round(drop(result)), 1)
})

# ═══════════════════════════════════════════════════════════════════════
# 10. Composition: nested logic expressions
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicSolve::test_composition
test_that("Or(And(x1, x2), Not(x3)) truth table", {
  for (a in c(0, 1)) {
    for (b in c(0, 1)) {
      for (cc in c(0, 1)) {
        x1 <- Variable(boolean = TRUE)
        x2 <- Variable(boolean = TRUE)
        x3 <- Variable(boolean = TRUE)
        expected <- as.integer((a & b) | (!cc))
        result <- .eval_logic(
          Or(And(x1, x2), Not(x3)),
          list(x1 == a, x2 == b, x3 == cc)
        )
        expect_equal(round(drop(result)), expected,
                     info = paste("Or(And(", a, ",", b, "), Not(", cc, "))"))
      }
    }
  }
})

## @cvxpy test_logic.py::TestLogicOperators::test_invert_on_logic_expr
test_that("!(x & y) via operator (NAND)", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  expr <- !(x & y)
  expect_true(S7_inherits(expr, Not))

  ## NAND(1,1) = 0
  result <- .eval_logic(expr, list(x == 1, y == 1))
  expect_equal(round(drop(result)), 0)

  ## NAND(1,0) = 1
  result2 <- .eval_logic(!(x & y), list(x == 1, y == 0))
  expect_equal(round(drop(result2)), 1)
})

## @cvxpy test_logic.py::TestLogicOperators::test_composed_operators
test_that("XNOR via operators: (x & y) | (!x & !y)", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  for (a in c(0, 1)) {
    for (b in c(0, 1)) {
      expr <- (x & y) | (!x & !y)
      result <- .eval_logic(expr, list(x == a, y == b))
      expected <- as.integer(a == b)
      expect_equal(round(drop(result)), expected,
                   info = paste("XNOR(", a, ",", b, ")"))
    }
  }
})

## @cvxpy NONE
test_that("De Morgan's law: !(x & y) == (!x | !y)", {
  for (a in c(0, 1)) {
    for (b in c(0, 1)) {
      x <- Variable(boolean = TRUE)
      y <- Variable(boolean = TRUE)
      lhs <- .eval_logic(Not(And(x, y)), list(x == a, y == b))
      rhs <- .eval_logic(Or(Not(x), Not(y)), list(x == a, y == b))
      expect_equal(round(drop(lhs)), round(drop(rhs)),
                   info = paste("De Morgan (", a, ",", b, ")"))
    }
  }
})

## @cvxpy NONE
test_that("De Morgan's law: !(x | y) == (!x & !y)", {
  for (a in c(0, 1)) {
    for (b in c(0, 1)) {
      x <- Variable(boolean = TRUE)
      y <- Variable(boolean = TRUE)
      lhs <- .eval_logic(Not(Or(x, y)), list(x == a, y == b))
      rhs <- .eval_logic(And(Not(x), Not(y)), list(x == a, y == b))
      expect_equal(round(drop(lhs)), round(drop(rhs)),
                   info = paste("De Morgan (", a, ",", b, ")"))
    }
  }
})

# ═══════════════════════════════════════════════════════════════════════
# 11. Vectorized: logic on vector boolean variables
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicSolve::test_vector
test_that("Vectorized And: elementwise truth table", {
  x <- Variable(3, boolean = TRUE)
  y <- Variable(3, boolean = TRUE)
  fix <- list(x == c(1, 0, 1), y == c(1, 1, 0))
  result <- .eval_logic(And(x, y), fix)
  expect_equal(round(drop(result)), c(1, 0, 0))
})

## @cvxpy test_logic.py::TestLogicSolve::test_vector
test_that("Vectorized Or: elementwise truth table", {
  x <- Variable(3, boolean = TRUE)
  y <- Variable(3, boolean = TRUE)
  fix <- list(x == c(1, 0, 1), y == c(1, 1, 0))
  result <- .eval_logic(Or(x, y), fix)
  expect_equal(round(drop(result)), c(1, 1, 1))
})

## @cvxpy test_logic.py::TestLogicSolve::test_vector
test_that("Vectorized Not: elementwise truth table", {
  x <- Variable(3, boolean = TRUE)
  fix <- list(x == c(1, 0, 1))
  result <- .eval_logic(Not(x), fix)
  expect_equal(round(drop(result)), c(0, 1, 0))
})

## @cvxpy test_logic.py::TestLogicOperators::test_vector_operators
test_that("Vectorized Xor: elementwise truth table", {
  x <- Variable(3, boolean = TRUE)
  y <- Variable(3, boolean = TRUE)
  fix <- list(x == c(1, 0, 1), y == c(1, 1, 0))
  result <- .eval_logic(Xor(x, y), fix)
  expect_equal(round(drop(result)), c(0, 1, 1))
})

## @cvxpy test_logic.py::TestLogicOperators::test_vector_operators
test_that("Vectorized operators: &, |, !", {
  x <- Variable(3, boolean = TRUE)
  y <- Variable(3, boolean = TRUE)
  fix <- list(x == c(1, 0, 1), y == c(1, 1, 0))

  result_and <- .eval_logic(x & y, fix)
  expect_equal(round(drop(result_and)), c(1, 0, 0))

  result_or <- .eval_logic(x | y, fix)
  expect_equal(round(drop(result_or)), c(1, 1, 1))

  result_not <- .eval_logic(!x, fix)
  expect_equal(round(drop(result_not)), c(0, 1, 0))
})

## @cvxpy NONE
test_that("Vectorized optimization: maximize sum(And(x, y))", {
  x <- Variable(3, boolean = TRUE)
  y <- Variable(3, boolean = TRUE)
  prob <- Problem(Maximize(sum_entries(And(x, y))))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(prob))), 3)
  expect_equal(round(drop(value(x))), c(1, 1, 1))
  expect_equal(round(drop(value(y))), c(1, 1, 1))
})

## @cvxpy NONE
test_that("Vectorized optimization: minimize sum(Or(x, y))", {
  x <- Variable(3, boolean = TRUE)
  y <- Variable(3, boolean = TRUE)
  prob <- Problem(Minimize(sum_entries(Or(x, y))))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(prob))), 0)
  expect_equal(round(drop(value(x))), c(0, 0, 0))
  expect_equal(round(drop(value(y))), c(0, 0, 0))
})

# ═══════════════════════════════════════════════════════════════════════
# Additional: Boolean constant arguments
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicBoolConstant::test_scalar_bool_constants
test_that("Logic atoms accept R logical TRUE/FALSE constants", {
  x <- Variable(boolean = TRUE)

  ## And(x, TRUE): when x=1 -> 1; when x=0 -> 0
  result <- .eval_logic(And(x, Constant(TRUE)), list(x == 1))
  expect_equal(round(drop(result)), 1)

  result2 <- .eval_logic(And(x, Constant(FALSE)), list(x == 1))
  expect_equal(round(drop(result2)), 0)

  ## Or(x, TRUE) -> always 1
  result3 <- .eval_logic(Or(x, Constant(TRUE)), list(x == 0))
  expect_equal(round(drop(result3)), 1)

  ## Not(TRUE) -> 0
  result4 <- .eval_logic(Not(Constant(TRUE)), list())
  expect_equal(round(drop(result4)), 0)

  ## Xor(x, TRUE): when x=1 -> 0; when x=0 -> 1
  result5 <- .eval_logic(Xor(x, Constant(TRUE)), list(x == 1))
  expect_equal(round(drop(result5)), 0)

  result6 <- .eval_logic(Xor(x, Constant(TRUE)), list(x == 0))
  expect_equal(round(drop(result6)), 1)
})

# ═══════════════════════════════════════════════════════════════════════
# Additional: Logic in constraints
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicInConstraint::test_and_maximize
test_that("And constraint: And(x, y) == 1 forces both to 1", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(Constant(0)), list(And(x, y) == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(x))), 1)
  expect_equal(round(drop(value(y))), 1)
})

## @cvxpy test_logic.py::TestLogicInConstraint::test_or_minimize
test_that("Or constraint: Or(x, y) == 0 forces both to 0", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(Constant(0)), list(Or(x, y) == 0))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(x))), 0)
  expect_equal(round(drop(value(y))), 0)
})

## @cvxpy test_logic.py::TestLogicInConstraint::test_xor_constraint
test_that("Xor constraint in optimization context", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  ## Maximize x subject to Xor(x, y) == 1
  ## Xor==1 means exactly one is 1. Maximize x -> x=1, y=0
  prob <- Problem(Maximize(x), list(Xor(x, y) == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(x))), 1)
  expect_equal(round(drop(value(y))), 0)
})

## @cvxpy NONE
test_that("implies constraint restricts feasible region", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)
  ## Maximize x - y subject to implies(x, y) == 1
  ## implies(x, y) = !x | y; satisfied when x=0 or y=1
  ## Max x-y: try x=1, y=0 -> implies fails. Best is x=0, y=0 (val=0) or x=1, y=1 (val=0)
  prob <- Problem(Maximize(x - y), list(implies(x, y) == 1))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
  expect_equal(round(drop(value(prob))), 0)
})

# ═══════════════════════════════════════════════════════════════════════
# Boolean array constant as mask in And (test_logic.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicBoolConstant::test_bool_array_constant
test_that("And/& with boolean array constant as mask", {
  ## CVXPY: x = Variable(3, boolean=True), mask = [True, False, True]
  ## fix = [x == [1, 1, 0]]
  ## And(x, mask) → [1, 0, 0] (elementwise: 1&T=1, 1&F=0, 0&T=0)
  ## x & mask → same result
  ## Verified via CVXPY: [1, 0, 0]
  x <- Variable(3, boolean = TRUE)
  mask <- Constant(c(1, 0, 1))  # boolean mask as 0/1
  fix <- list(x == c(1, 1, 0))

  ## Functional And(x, mask)
  result1 <- .eval_logic(And(x, mask), fix)
  expect_equal(as.numeric(result1), c(1, 0, 0), tolerance = 1e-4)

  ## Operator dispatch: x & mask
  result2 <- .eval_logic(x & mask, fix)
  expect_equal(as.numeric(result2), c(1, 0, 0), tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# Operators with scalar boolean constant (test_logic.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_logic.py::TestLogicBoolConstant::test_operator_with_bool_constant
test_that("Logic operators with scalar boolean constants", {
  ## CVXPY: x = Variable(boolean=True), fix = [x == 1]
  ## x & True → 1, x | False → 1, x ^ True → 0
  ## Verified via CVXPY: 1.0, 1.0, 0.0
  x <- Variable(boolean = TRUE)
  fix <- list(x == 1)

  ## x & TRUE → And(x, TRUE) → 1 & 1 = 1
  result_and <- .eval_logic(x & Constant(1), fix)
  expect_equal(as.numeric(result_and), 1, tolerance = 1e-4)

  ## x | FALSE → Or(x, FALSE) → 1 | 0 = 1
  result_or <- .eval_logic(x | Constant(0), fix)
  expect_equal(as.numeric(result_or), 1, tolerance = 1e-4)

  ## Xor(x, TRUE) → 1 XOR 1 = 0
  result_xor <- .eval_logic(Xor(x, Constant(1)), fix)
  expect_equal(as.numeric(result_xor), 0, tolerance = 1e-4)
})
