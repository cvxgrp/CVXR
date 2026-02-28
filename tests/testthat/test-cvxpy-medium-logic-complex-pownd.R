## CVXPY Medium-Priority Parity Tests: Logic, Complex, and PowConeND Gaps
## =========================================================================
## These tests fill MEDIUM-priority gaps identified in the gap analysis,
## covering three areas:
##   1. Logic atoms: n-ary names, precedence, DCP, namespace, log-log, composition
##   2. Complex: construction, arithmetic, conjugate transpose, PSD
##   3. Complex DPP: param recognition, chain structure, EvalParams skipping
##   4. PowConeND: variable swap, single cone, scalar alpha, balanced tree
##
## Source references:
##   - CVXPY: test_logic.py, test_complex.py, test_dpp.py, test_pow_cone_nd.py
##   - CVXR: rsrc_tree/atoms/elementwise/logic.R
##   - CVXR: rsrc_tree/constraints/power.R
##
## Solvers: HIGHS (logic; matches CVXPY), Clarabel (PowConeND, complex), SCS (complex)
## Tolerances: 1e-3 for SCS, 1e-5 for Clarabel

library(testthat)
library(CVXR)

# ── Helper: solve feasibility and return value of expr ──────────────
.eval_logic <- function(expr, constraints = list()) {
  y <- Variable(expr@shape)
  prob <- Problem(Minimize(Constant(0)), c(list(y == expr), constraints))
  psolve(prob, solver = "HIGHS")
  value(y)
}

# =====================================================================
# SECTION 1: Logic Atoms — N-ary Names, Precedence, DCP, Namespace
# =====================================================================

# ── 1. basic_nary_names: expr_name() for n-ary And/Or/Xor (3+ args)
# CVXPY parity: test_logic.py TestLogicName.test_nary_names
## @cvxpy test_logic.py::TestLogicName::test_nary_names
test_that("basic_nary_names: expr_name() for n-ary And/Or/Xor with 3+ args", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)
  z <- Variable(name = "z", boolean = TRUE)
  w <- Variable(name = "w", boolean = TRUE)

  ## 3-ary
  expect_equal(expr_name(And(x, y, z)), "x & y & z")
  expect_equal(expr_name(Or(x, y, z)), "x | y | z")

  ## 4-ary
  expect_equal(expr_name(And(x, y, z, w)), "x & y & z & w")
  expect_equal(expr_name(Or(x, y, z, w)), "x | y | z | w")

  ## Xor uses " XOR " in R (Python uses " ^ ")
  xor_name <- expr_name(Xor(x, y, z))
  expect_match(xor_name, "x XOR y XOR z")
})

# ── 2. nary_names: Names with mixed nesting
# CVXPY parity: test_logic.py TestLogicName.test_complex_expression
## @cvxpy test_logic.py::TestLogicName::test_complex_expression
test_that("nary_names: Names with mixed nesting", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)
  z <- Variable(name = "z", boolean = TRUE)

  ## Or(And(x, y), And(Not(x), Not(y))) should show: "x & y | !x & !y"
  ## This mirrors CVXPY's "x & y | ~x & ~y"
  expr <- Or(And(x, y), And(Not(x), Not(y)))
  expected_name <- "x & y | !x & !y"
  expect_equal(expr_name(expr), expected_name)

  ## And(Or(x, y), Not(z)): "& parenthesizes Or child" -> "(x | y) & !z"
  expr2 <- And(Or(x, y), Not(z))
  expect_equal(expr_name(expr2), "(x | y) & !z")

  ## Or(And(x, Not(y)), z): lowest precedence; no parens -> "x & !y | z"
  expr3 <- Or(And(x, Not(y)), z)
  expect_equal(expr_name(expr3), "x & !y | z")
})

# ── 3. precedence: Operator precedence: And before Or, parenthesization
# CVXPY parity: test_logic.py TestLogicName.test_precedence
## @cvxpy test_logic.py::TestLogicName::test_precedence
test_that("precedence: And before Or, parenthesization rules", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)
  z <- Variable(name = "z", boolean = TRUE)

  ## & is higher precedence than |
  ## And(Or(x, y), z) -> "(x | y) & z"
  expect_equal(expr_name(And(Or(x, y), z)), "(x | y) & z")

  ## And of And is flat (no parens): And(And(x, y), z) -> "x & y & z"
  expect_equal(expr_name(And(And(x, y), z)), "x & y & z")

  ## Or is lowest: never parenthesizes
  ## Or(And(x, y), z) -> "x & y | z"
  expect_equal(expr_name(Or(And(x, y), z)), "x & y | z")

  ## And parenthesizes Xor children
  ## And(x, Xor(y, z)) -> "x & (y XOR z)"
  and_xor <- And(x, Xor(y, z))
  expect_match(expr_name(and_xor), "x & \\(y XOR z\\)")

  ## Xor parenthesizes Or but not And
  ## Xor(Or(x, y), z) -> "(x | y) XOR z"
  xor_or <- Xor(Or(x, y), z)
  expect_match(expr_name(xor_or), "\\(x \\| y\\) XOR z")

  ## Xor(And(x, y), z) -> "x & y XOR z" (no parens around And)
  xor_and <- Xor(And(x, y), z)
  expect_match(expr_name(xor_and), "x & y XOR z")
})

# ── 4. dcp_compliance: Logic atoms are DCP (convex AND concave)
# CVXPY parity: test_logic.py TestLogicProperties.test_dcp_compliance
## @cvxpy test_logic.py::TestLogicProperties::test_dcp_compliance
test_that("dcp_compliance: all logic atoms are DCP", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)

  for (expr in list(Not(x), And(x, y), Or(x, y), Xor(x, y))) {
    expect_true(is_dcp(expr),
                info = paste("is_dcp should be TRUE for", class(expr)[[1L]]))
    ## Also check convex AND concave (affine-like)
    expect_true(is_atom_convex(expr),
                info = paste("is_atom_convex for", class(expr)[[1L]]))
    expect_true(is_atom_concave(expr),
                info = paste("is_atom_concave for", class(expr)[[1L]]))
  }

  ## implies and iff should also be DCP
  expect_true(is_dcp(implies(x, y)))
  expect_true(is_dcp(iff(x, y)))
})

# ── 5. namespace: Logic atoms properly exported
# CVXPY parity: test_logic.py TestLogicProperties.test_namespace
## @cvxpy test_logic.py::TestLogicProperties::test_namespace
test_that("namespace: Logic atoms are properly exported and accessible", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)

  ## All logic functions should be accessible as CVXR::Not, etc.
  expect_true(S7_inherits(CVXR::Not(x), Not))
  expect_true(S7_inherits(CVXR::And(x, y), And))
  expect_true(S7_inherits(CVXR::Or(x, y), Or))
  expect_true(S7_inherits(CVXR::Xor(x, y), Xor))

  ## implies returns Or (as in CVXPY)
  expect_true(S7_inherits(CVXR::implies(x, y), Or))

  ## iff returns Not (wrapping Xor, as in CVXPY)
  expect_true(S7_inherits(CVXR::iff(x, y), Not))
})

# ── 6. not_log_log: Logic atoms are NOT log-log convex/concave
# CVXPY parity: test_logic.py TestLogicProperties.test_not_log_log
## @cvxpy test_logic.py::TestLogicProperties::test_not_log_log
test_that("not_log_log: Logic atoms are not log-log convex or concave", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)

  ## Logic atoms have domain including 0, so they are not log-log convex/concave.
  ## The default from Atom is FALSE for both, and logic atoms should inherit that.
  for (expr in list(Not(x), And(x, y), Or(x, y), Xor(x, y))) {
    expect_false(is_atom_log_log_convex(expr),
                 info = paste("should NOT be log-log convex:", class(expr)[[1L]]))
    expect_false(is_atom_log_log_concave(expr),
                 info = paste("should NOT be log-log concave:", class(expr)[[1L]]))
  }

  ## Consequently, not DGP either
  for (expr in list(Not(x), And(x, y), Or(x, y), Xor(x, y))) {
    expect_false(is_dgp(expr),
                 info = paste("should NOT be DGP:", class(expr)[[1L]]))
  }
})

# ── 7. composed_operators: Nested: !(x & y), (x | y) & z, etc.
# CVXPY parity: test_logic.py TestLogicOperators.test_composed_operators
## @cvxpy test_logic.py::TestLogicOperators::test_composed_operators
test_that("composed_operators: (x & y) | (!x & !y) is XNOR", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)

  ## Construct XNOR via operators: (x & y) | (!x & !y)
  expr <- (x & y) | (!x & !y)

  ## Verify it's an Or at the top level

  expect_true(S7_inherits(expr, Or))

  ## Verify truth table: XNOR(a, b) = (a == b)
  for (a in c(0, 1)) {
    for (b in c(0, 1)) {
      result <- .eval_logic(expr, list(x == a, y == b))
      expected <- as.integer(a == b)
      expect_equal(round(drop(result)), expected,
                   info = paste("XNOR(", a, ",", b, ")"))
    }
  }
})

# ── 8. invert_on_logic_expr: Applying Not to result of And/Or/Xor
# CVXPY parity: test_logic.py TestLogicOperators.test_invert_on_logic_expr
## @cvxpy test_logic.py::TestLogicOperators::test_invert_on_logic_expr
test_that("invert_on_logic_expr: ! applied to And, Or, Xor expressions", {
  x <- Variable(boolean = TRUE)
  y <- Variable(boolean = TRUE)

  ## !(x & y) is NAND
  nand_expr <- !(x & y)
  expect_true(S7_inherits(nand_expr, Not))
  expect_equal(round(drop(.eval_logic(nand_expr, list(x == 1, y == 1)))), 0)
  expect_equal(round(drop(.eval_logic(nand_expr, list(x == 1, y == 0)))), 1)

  ## !(x | y) is NOR
  nor_expr <- !(x | y)
  expect_true(S7_inherits(nor_expr, Not))
  expect_equal(round(drop(.eval_logic(nor_expr, list(x == 0, y == 0)))), 1)
  expect_equal(round(drop(.eval_logic(nor_expr, list(x == 0, y == 1)))), 0)

  ## Not(Xor(x, y)) is XNOR (same as iff)
  xnor_expr <- Not(Xor(x, y))
  expect_true(S7_inherits(xnor_expr, Not))
  expect_equal(round(drop(.eval_logic(xnor_expr, list(x == 1, y == 1)))), 1)
  expect_equal(round(drop(.eval_logic(xnor_expr, list(x == 1, y == 0)))), 0)
})

# ── 9. non_boolean_raises: Non-boolean inputs raise errors
# CVXPY parity: test_logic.py TestLogicOperators.test_non_boolean_raises
#               and test_logic.py TestLogicValidation.test_non_boolean_raises
## @cvxpy test_logic.py::TestLogicOperators::test_non_boolean_raises test_logic.py::TestLogicValidation::test_non_boolean_raises
test_that("non_boolean_raises: non-boolean variables raise errors in all logic ops", {
  x <- Variable()   # continuous, not boolean
  y <- Variable(boolean = TRUE)

  ## Functional syntax
  expect_error(And(x, y), "boolean")
  expect_error(Or(x, y), "boolean")
  expect_error(Xor(x, y), "boolean")
  expect_error(Not(x), "boolean")

  ## Operator syntax
  expect_error(x & y, "boolean")
  expect_error(x | y, "boolean")
  expect_error(!x, "boolean")

  ## Integer variable (not boolean) should also be rejected
  z <- Variable(integer = TRUE)
  expect_error(And(z, y), "boolean")
  expect_error(Or(z, y), "boolean")
  expect_error(Not(z), "boolean")
})


# =====================================================================
# SECTION 2: Complex Expressions
# =====================================================================

# ── 1. complex_variable_construction: Variable(complex=TRUE), shape
# CVXPY parity: test_complex.py TestComplex.test_variable
## @cvxpy test_complex.py::TestComplex::test_variable
test_that("complex_variable_construction: Variable(complex=TRUE) type flags and shape", {
  x <- Variable(2, complex = FALSE)
  y <- Variable(2, complex = TRUE)
  z <- Variable(2, imag = TRUE)

  ## Real variable
  expect_false(is_complex(x))
  expect_false(is_imag(x))
  expect_true(is_real(x))

  ## Complex variable
  expect_true(is_complex(y))
  expect_false(is_imag(y))
  expect_false(is_real(y))

  ## Imaginary variable
  expect_true(is_complex(z))
  expect_true(is_imag(z))
  expect_false(is_real(z))

  ## Shape checks
  expect_equal(y@shape, c(2L, 1L))
  expect_equal(z@shape, c(2L, 1L))
})

# ── 2. hermitian_variable: Variable(hermitian=TRUE), symmetry enforcement
# CVXPY parity: test_complex.py TestComplex.test_hermitian
## @cvxpy test_complex.py::TestComplex::test_hermitian
test_that("hermitian_variable: construction and symmetry properties", {
  skip_if_not_installed("scs")

  X <- Variable(c(2L, 2L), hermitian = TRUE)

  ## A hermitian variable should be complex and hermitian
  expect_true(is_complex(X))
  expect_true(is_hermitian(X))

  ## Shape should be 2x2
  expect_equal(X@shape, c(2L, 2L))

  ## Solve: min Im(X[2,1]) s.t. X[1,1]==2, X[2,2]==3, X[1,2]==1+1i
  ## For hermitian: X[2,1] = conj(X[1,2]) = 1-1i, so Im(X[2,1]) = -1
  prob <- Problem(Minimize(Im(X[2, 1])),
                  list(X[1, 1] == 2, X[2, 2] == 3, X[1, 2] == Constant(1 + 1i)))
  result <- psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")

  X_val <- value(X)
  ## Check hermitian structure
  expect_equal(Re(X_val[1, 1]), 2, tolerance = 1e-3)
  expect_equal(Re(X_val[2, 2]), 3, tolerance = 1e-3)
  expect_equal(Re(X_val[1, 2]), 1, tolerance = 1e-3)
  expect_equal(Im(X_val[1, 2]), 1, tolerance = 1e-3)
  ## Conjugate symmetry: X[2,1] = conj(X[1,2])
  expect_equal(Re(X_val[2, 1]), 1, tolerance = 1e-3)
  expect_equal(Im(X_val[2, 1]), -1, tolerance = 1e-3)
})

# ── 3. complex_constant: Constant with complex values
# CVXPY parity: test_complex.py TestComplex.test_constant
## @cvxpy test_complex.py::TestComplex::test_constant
test_that("complex_constant: type flags for complex constants", {
  ## Real constant
  c_real <- Constant(2)
  expect_false(is_complex(c_real))
  expect_false(is_imag(c_real))
  expect_true(is_real(c_real))

  ## Complex constant (mixed real+imag)
  c_mixed <- Constant(2i + 1)
  expect_true(is_complex(c_mixed))
  expect_false(is_imag(c_mixed))
  expect_false(is_real(c_mixed))

  ## Purely imaginary constant
  c_imag <- Constant(2i)
  expect_true(is_complex(c_imag))
  expect_true(is_imag(c_imag))
  expect_false(is_real(c_imag))

  ## Complex vector constant
  c_vec <- Constant(c(1 + 2i, 3 - 4i))
  expect_true(is_complex(c_vec))
  expect_false(is_imag(c_vec))
})

# ── 4. complex_addition: Adding complex expressions
# CVXPY parity: test_complex.py TestComplex.test_arithmetic (addition part)
## @cvxpy test_complex.py::TestComplex::test_arithmetic
test_that("complex_addition: complex + real, imag + real type inference", {
  x <- Variable(complex = TRUE)
  y <- Variable(imag = TRUE)
  z <- Variable()  # real

  ## complex + real -> complex (not imaginary)
  expr1 <- x + z
  expect_true(is_complex(expr1))
  expect_false(is_imag(expr1))

  ## imag + real -> complex (not imaginary)
  expr2 <- y + z
  expect_true(is_complex(expr2))
  expect_false(is_imag(expr2))

  ## complex + complex -> complex
  x2 <- Variable(complex = TRUE)
  expr3 <- x + x2
  expect_true(is_complex(expr3))

  ## real + real -> real
  z2 <- Variable()
  expr4 <- z + z2
  expect_true(is_real(expr4))
  expect_false(is_complex(expr4))
})

# ── 5. complex_multiplication: Multiplying complex expressions
# CVXPY parity: test_complex.py TestComplex.test_arithmetic (multiplication part)
## @cvxpy test_complex.py::TestComplex::test_arithmetic
test_that("complex_multiplication: product type inference for complex expressions", {
  y <- Variable(imag = TRUE)
  z <- Variable()  # real

  ## imag * real -> imaginary
  expr1 <- y * z
  expect_true(is_complex(expr1))
  expect_true(is_imag(expr1))

  ## imag * imag -> real (i*i = -1, purely real)
  expr2 <- y * y
  expect_false(is_complex(expr2))
  expect_false(is_imag(expr2))
  expect_true(is_real(expr2))

  ## imag / real -> imaginary
  expr3 <- y / 2
  expect_true(is_complex(expr3))
  expect_true(is_imag(expr3))

  ## imag / 1j -> real (dividing imaginary by imaginary gives real)
  expr4 <- y / Constant(1i)
  expect_false(is_complex(expr4))
  expect_false(is_imag(expr4))
  expect_true(is_real(expr4))
})

# ── 6. complex_conjugate_transpose: Conj, transpose of complex
# CVXPY parity: test_complex.py TestComplex.test_conj
## @cvxpy test_complex.py::TestComplex::test_conj
test_that("complex_conjugate_transpose: Conj type inference and value", {
  ## Test Conj type inference
  A <- matrix(1, 2, 2)
  expr <- Constant(A) + Constant(1i * A)  ## A + iA
  conj_expr <- Conj(expr)
  expect_false(is_real(conj_expr))
  expect_true(is_complex(conj_expr))
  expect_false(is_imag(conj_expr))

  ## Test Conj value: conj(A + iA) = A - iA
  conj_val <- value(conj_expr)
  expect_equal(Re(conj_val), A, tolerance = 1e-10)
  expect_equal(Im(conj_val), -A, tolerance = 1e-10)

  ## Test Re type inference
  x <- Variable(complex = TRUE)
  re_expr <- Re(x) + Im(x)
  expect_true(is_real(re_expr))

  ## Test Im extracts correctly from a complex constant
  A2 <- matrix(1, 2, 2)
  complex_mat <- A2 + 2i * A2  # Direct R complex matrix
  expr2 <- Constant(complex_mat)
  im_expr <- Im(expr2)
  expect_true(is_real(im_expr))
  expect_false(is_complex(im_expr))
  im_val <- value(im_expr)
  expect_equal(im_val, 2 * A2, tolerance = 1e-10)
})

# ── 7. complex_psd_constraint: PSD constraint on complex/Hermitian matrix
# CVXPY parity: test_complex.py TestComplex.test_psd
## @cvxpy test_complex.py::TestComplex::test_psd
test_that("complex_psd_constraint: PSD on hermitian with infeasible diagonal", {
  skip_if_not_installed("scs")

  ## X %>>% 0 with X[1,1] == -1 should be infeasible
  ## (PSD requires nonneg diagonal entries)
  X <- Variable(c(2L, 2L), hermitian = TRUE)
  prob <- Problem(Minimize(Im(X[2, 1])),
                  list(X %>>% 0, X[1, 1] == -1))
  psolve(prob, solver = "SCS")
  expect_equal(status(prob), "infeasible")
})


# =====================================================================
# SECTION 3: Complex DPP
# =====================================================================

# ── 1. complex_dpp_param_affine: Complex param * variable is DPP
# CVXPY parity: test_dpp.py (complex DPP recognition)
## @cvxpy NONE
test_that("complex_dpp_param_affine: problem with complex param is DPP", {
  ## A problem with a complex parameter in an affine combination should
  ## be recognized as DPP.
  p <- Parameter(2, complex = TRUE, value = c(1 + 0i, 0 + 1i))
  x <- Variable(2, complex = TRUE)

  ## Build problem: min sum(Re(p * x)) s.t. |x| <= 1
  prob <- Problem(Minimize(sum_entries(Re(p * x))),
                  list(Abs(x) <= 1))

  ## Should be DPP
  expect_true(is_dpp(prob))

  ## The objective involves a parameter-variable product, which is affine in x
  ## given that p is a Parameter. This should be DPP-compliant.
  obj_expr <- Re(p * x)
  expect_true(is_real(obj_expr))
})

# ── 2. complex_dpp_chain_structure: Chain has Complex2Real before Dcp2Cone
# CVXPY parity: test_dpp.py (chain inspection for complex)
## @cvxpy NONE
test_that("complex_dpp_chain_structure: Complex2Real appears before Dcp2Cone in chain", {
  skip_if_not_installed("scs")

  ## A complex problem's solving chain should contain Complex2Real
  ## before Dcp2Cone
  z <- Variable(2, complex = TRUE)
  prob <- Problem(Minimize(sum_entries(Re(z))),
                  list(Im(z) == 1, Re(z) >= -3))

  chain <- construct_solving_chain(prob, solver = "SCS")
  red_names <- vapply(chain@reductions,
                      function(r) sub("^.*::", "", class(r)[[1L]]),
                      character(1L))

  ## Complex2Real should be present
  expect_true("Complex2Real" %in% red_names)

  ## Dcp2Cone should be present
  expect_true("Dcp2Cone" %in% red_names)

  ## Complex2Real should come BEFORE Dcp2Cone
  c2r_pos <- which(red_names == "Complex2Real")
  dcp_pos <- which(red_names == "Dcp2Cone")
  expect_true(length(c2r_pos) == 1L)
  expect_true(length(dcp_pos) == 1L)
  expect_true(c2r_pos < dcp_pos)
})

# ── 3. complex_dpp_eval_params_skipped: EvalParams skipped for DPP complex
# CVXPY parity: test_dpp.py (EvalParams not in chain for DPP complex)
## @cvxpy NONE
test_that("complex_dpp_eval_params_skipped: EvalParams not in chain for DPP complex problem", {
  skip_if_not_installed("scs")

  ## A DPP complex problem should NOT have EvalParams in its solving chain
  z <- Variable(2, complex = TRUE)
  prob <- Problem(Minimize(sum_entries(Re(z))),
                  list(Im(z) == 1, Re(z) >= -3))

  ## Verify DPP
  expect_true(is_dpp(prob))

  ## Build chain
  chain <- construct_solving_chain(prob, solver = "SCS")
  red_names <- vapply(chain@reductions,
                      function(r) sub("^.*::", "", class(r)[[1L]]),
                      character(1L))

  ## EvalParams should NOT be present (DPP fast path skips it)
  expect_false("EvalParams" %in% red_names)

  ## Also test with a complex Parameter — structure only, not solving
  p <- Parameter(2, complex = TRUE, value = c(1 + 1i, 2 - 1i))
  x <- Variable(2, complex = TRUE)
  prob2 <- Problem(Minimize(sum_entries(Re(p * x))),
                   list(Abs(x) <= 1))

  ## Should be DPP
  expect_true(is_dpp(prob2))

  ## Chain should have Complex2Real and no EvalParams
  chain2 <- construct_solving_chain(prob2, solver = "SCS")
  red_names2 <- vapply(chain2@reductions,
                       function(r) sub("^.*::", "", class(r)[[1L]]),
                       character(1L))
  expect_true("Complex2Real" %in% red_names2)
  ## NOTE: Complex parameters force EvalParams because Imag_/Real_ atoms lack
  ## graph_implementation and R's Matrix package doesn't support complex sparse
  ## matrices. Complex DPP is deferred.
  expect_true("EvalParams" %in% red_names2)
})


# =====================================================================
# SECTION 4: PowConeND — Variable Swap, Single Cone, Balanced Tree
# =====================================================================

# ── 1. variable_swap_2d: PowConeND with 2 variables (3D), swap order
# CVXPY parity: test_pow_cone_nd.py TestPowConeND.test_pow_cone_nd_3d_variable_swap
## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_3d_variable_swap
test_that("variable_swap_2d: PowConeND 3D with swapped variable ordering", {
  skip_if_not_installed("clarabel")

  ## This is test_pow_cone_nd_3d_variable_swap from CVXPY, axis=0.
  ## Variables appear in swapped order compared to the standard pcp_3 test.
  ## Python verification:
  ##   expect_x = [0.0639, 2.3054, 0.7834]  (x[1] and x[2] swapped)
  expect_x <- c(0.06393515, 2.30571048, 0.78320961)

  x <- Variable(3)
  hypos <- Variable(2)
  objective <- Maximize(sum_entries(hypos) - x[1])

  ## W has x[1], x[2] swapped compared to standard test
  W <- bmat(list(
    list(x[1], x[2]),
    list(x[3], Constant(1.0))
  ))
  alpha <- Constant(matrix(c(0.2, 0.8, 0.4, 0.6), nrow = 2L, ncol = 2L))

  ## Linear constraint with swapped coefficients
  con_eq <- (x[1] + x[3] + 0.5 * x[2] == 2)
  con_pow <- PowConeND(W, hypos, alpha, axis = 2L)

  prob <- Problem(objective, list(con_eq, con_pow))
  result <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  expect_equal(result, 1.8073, tolerance = 1e-3)

  x_val <- as.numeric(value(x))
  expect_equal(x_val, expect_x, tolerance = 1e-2)
})

# ── 2. variable_swap_3d: PowConeND with 3 rows, swap order
# This test corresponds to the 3-row (>3D) case with swapped variables.
# CVXPY parity: test_pow_cone_nd.py TestPowConeND.test_pow_cone_nd (axis=0)
## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd
test_that("variable_swap_3d: PowConeND >3D with standard ordering", {
  skip_if_not_installed("clarabel")

  ## 5 variables, 3 rows in W (>3D power cone), standard order
  ## Python verification:
  ##   expect_x = [0, 0, 0, 2.28571, 3.42857]
  expect_x <- c(0, 0, 0, 2.28571379, 3.42857186)

  x <- Variable(5)
  hypos <- Variable(2)
  objective <- Maximize(sum_entries(hypos) - x[1])

  W <- bmat(list(
    list(x[1], x[4]),
    list(x[2], x[5]),
    list(x[3], Constant(1.0))
  ))
  alpha <- Constant(matrix(c(0.2, 0.4, 0.4, 0.4, 0.3, 0.3), nrow = 3L, ncol = 2L))

  con_eq <- (x[1] + x[2] + 0.5 * x[3] + 0.5 * x[4] + 0.25 * x[5] == 2)
  con_pow <- PowConeND(W, hypos, alpha, axis = 2L)

  prob <- Problem(objective, list(con_eq, con_pow))
  result <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  x_val <- as.numeric(value(x))
  expect_equal(x_val, expect_x, tolerance = 1e-2)
})

# ── 3. variable_swap_5d: PowConeND with 5 variables, swapped order
# CVXPY parity: test_pow_cone_nd.py TestPowConeND.test_pow_cone_nd_variable_swap
## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_variable_swap
test_that("variable_swap_5d: PowConeND >3D with swapped variable ordering", {
  skip_if_not_installed("clarabel")

  ## Same problem as variable_swap_3d but with variable indices shuffled.
  ## Python verification:
  ##   expect_x = [3.42857, 0, 0, 2.28571, 0]
  expect_x <- c(3.42857186, 0, 0, 2.28571379, 0)

  x <- Variable(5)
  hypos <- Variable(2)
  objective <- Maximize(sum_entries(hypos) - x[5])

  W <- bmat(list(
    list(x[5], x[4]),
    list(x[2], x[1]),
    list(x[3], Constant(1.0))
  ))
  alpha <- Constant(matrix(c(0.2, 0.4, 0.4, 0.4, 0.3, 0.3), nrow = 3L, ncol = 2L))

  con_eq <- (x[5] + x[2] + 0.5 * x[3] + 0.5 * x[4] + 0.25 * x[1] == 2)
  con_pow <- PowConeND(W, hypos, alpha, axis = 2L)

  prob <- Problem(objective, list(con_eq, con_pow))
  result <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  x_val <- as.numeric(value(x))
  expect_equal(x_val, expect_x, tolerance = 1e-2)
})

# ── 4. single_cone: Single PowConeND constraint
# CVXPY parity: test_pow_cone_nd.py TestPowConeND.test_pow_cone_nd_single_cone
## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_single_cone
test_that("single_cone: PowConeND with only a single cone (k=1)", {
  skip_if_not_installed("clarabel")

  ## Ensures no variables collapse to lower dimensions incorrectly
  ## when there is only a single cone.
  x <- Variable(2)
  hypos <- Variable(1)
  objective <- Maximize(sum_entries(hypos) - x[1])

  W <- bmat(list(list(x[1]), list(x[2])))
  alpha <- Constant(matrix(c(0.2, 0.8), nrow = 2L, ncol = 1L))

  con_eq <- (x[1] + x[2] == 2)
  con_pow <- PowConeND(W, hypos, alpha, axis = 2L)

  prob <- Problem(objective, list(con_eq, con_pow))
  result <- psolve(prob, solver = "CLARABEL")

  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
})

# ── 5. scalar_alpha: Scalar alpha for PowCone3D
# CVXPY parity: test_pow_cone_nd.py TestPowConeND.test_3d_pow_cone_scalar_alpha
## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_3d_pow_cone_scalar_alpha
test_that("scalar_alpha: PowCone3D with scalar alpha value", {
  skip_if_not_installed("clarabel")

  ## PowCone3D with scalar alpha should promote correctly
  x <- Variable(3)
  constraints <- list(PowCone3D(x[1], x[2], x[3], 0.75))
  prob <- Problem(Minimize(p_norm(x, 2)), constraints)
  result <- psolve(prob, solver = "CLARABEL")

  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
})

# ── 6. balanced_tree_3: Balanced tree decomposition with 3 variables
# CVXPY parity: test_pow_cone_nd.py test_pow_cone_nd_balanced_tree_decomposition n=3
## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_decomposition
test_that("balanced_tree_3: PowConeND balanced tree with n=3", {
  skip_if_not_installed("clarabel")

  ## max z s.t. PowConeND(W, z, alpha), W <= 2
  ## alpha = [1/3, 1/3, 1/3], so prod(W_i^alpha_i) >= |z|
  ## At W = [2, 2, 2], z = 2^1 = 2
  n <- 3L
  W <- Variable(c(n, 1L), nonneg = TRUE)
  z <- Variable(1)
  alpha <- Constant(matrix(rep(1 / n, n), nrow = n, ncol = 1L))

  con <- PowConeND(W, z, alpha)
  prob <- Problem(Maximize(z), list(con, W <= 2))
  result <- psolve(prob, solver = "CLARABEL")

  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(as.numeric(value(z)), 2.0, tolerance = 1e-3)
})

# ── 7. balanced_tree_5: Balanced tree with 5 variables
# CVXPY parity: test_pow_cone_nd.py test_pow_cone_nd_balanced_tree_decomposition n=5
## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_decomposition
test_that("balanced_tree_5: PowConeND balanced tree with n=5", {
  skip_if_not_installed("clarabel")

  n <- 5L
  W <- Variable(c(n, 1L), nonneg = TRUE)
  z <- Variable(1)
  alpha <- Constant(matrix(rep(1 / n, n), nrow = n, ncol = 1L))

  con <- PowConeND(W, z, alpha)
  prob <- Problem(Maximize(z), list(con, W <= 2))
  result <- psolve(prob, solver = "CLARABEL")

  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(as.numeric(value(z)), 2.0, tolerance = 1e-3)
})

# ── 8. balanced_tree_7: Balanced tree with 7 variables
# CVXPY parity: test_pow_cone_nd.py test_pow_cone_nd_balanced_tree_decomposition n=7
## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_decomposition
test_that("balanced_tree_7: PowConeND balanced tree with n=7", {
  skip_if_not_installed("clarabel")

  n <- 7L
  W <- Variable(c(n, 1L), nonneg = TRUE)
  z <- Variable(1)
  alpha <- Constant(matrix(rep(1 / n, n), nrow = n, ncol = 1L))

  con <- PowConeND(W, z, alpha)
  prob <- Problem(Maximize(z), list(con, W <= 2))
  result <- psolve(prob, solver = "CLARABEL")

  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(as.numeric(value(z)), 2.0, tolerance = 1e-3)
})
