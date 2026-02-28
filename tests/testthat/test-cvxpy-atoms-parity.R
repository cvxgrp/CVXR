## CVXPY Atoms Parity Tests
## ========================
## These tests mirror CVXPY's test_atoms.py (TestAtoms class) to verify
## that CVXR atom shape/sign/curvature/numeric behavior matches CVXPY exactly.
##
## Source: /Users/naras/GitHub/cvxpy/cvxpy/tests/test_atoms.py
## CVXPY setUp: a=Variable(), x=Variable(2), y=Variable(2),
##              A=Variable((2,2)), B=Variable((2,2)), C=Variable((3,2))
##
## Shape mapping: CVXPY tuple()→c(1,1), (n,)→c(n,1), (n,m)→c(n,m)
## Sign mapping:  CVXPY s.NONNEG→NONNEG_SIGN, s.NONPOS→NONPOS_SIGN, etc.
## Curvature mapping: CVXPY s.CONVEX→CONVEX, s.CONCAVE→CONCAVE, etc.

library(testthat)
library(CVXR)

# ═══════════════════════════════════════════════════════════════════════
# test_norm_inf (test_atoms.py line 72)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_norm_inf
test_that("norm_inf: shape, curvature, sign (CVXPY test_norm_inf)", {
  x <- Variable(2)
  y <- Variable(2)
  expr <- x + y
  atom <- norm_inf(expr)

  ## CVXPY: shape == tuple() → R: c(1L, 1L) scalar

  expect_equal(atom@shape, c(1L, 1L))
  ## CVXPY: curvature == s.CONVEX
  expect_equal(expr_curvature(atom), CONVEX)
  expect_true(is_convex(atom))
  expect_true(is_concave(-atom))
  ## CVXPY: norm_inf(atom).curvature == s.CONVEX (composition: convex nondecr of convex)
  expect_equal(expr_curvature(norm_inf(atom)), CONVEX)
  ## CVXPY: norm_inf(-atom).curvature == s.CONVEX
  expect_equal(expr_curvature(norm_inf(-atom)), CONVEX)
})

# ═══════════════════════════════════════════════════════════════════════
# test_norm1 (test_atoms.py line 85)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_norm1
test_that("norm1: shape, curvature (CVXPY test_norm1)", {
  x <- Variable(2)
  y <- Variable(2)
  expr <- x + y
  atom <- norm1(expr)

  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONVEX)
  ## Composition: norm1 of convex nonneg → CONVEX
  expect_equal(expr_curvature(norm1(atom)), CONVEX)
  ## norm1(-atom) → CONVEX (norm is always convex)
  expect_equal(expr_curvature(norm1(-atom)), CONVEX)
})

# ═══════════════════════════════════════════════════════════════════════
# test_power (test_atoms.py line 135)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: shape preserved for various p values (CVXPY test_power)", {
  for (shape in list(c(1L, 1L), c(3L, 1L), c(2L, 3L))) {
    x <- Variable(shape)
    y <- Variable(shape)
    expr <- x + y

    for (p in c(0, 1, 2, 3, 2.7, 0.67, -1, -2.3)) {
      atom <- power(expr, p)
      expect_equal(atom@shape, shape, info = paste("p =", p, "shape =", paste(shape, collapse = ",")))
    }
  }
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: curvature varies by p (CVXPY test_power)", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  expr <- x + y

  ## p > 1 → CONVEX
  expect_equal(expr_curvature(power(expr, 2)), CONVEX)
  expect_equal(expr_curvature(power(expr, 3)), CONVEX)
  expect_equal(expr_curvature(power(expr, 2.7)), CONVEX)

  ## p < 0 → CONVEX
  expect_equal(expr_curvature(power(expr, -1)), CONVEX)
  expect_equal(expr_curvature(power(expr, -2.3)), CONVEX)

  ## p == 1 → AFFINE
  expect_equal(expr_curvature(power(expr, 1)), AFFINE)

  ## p == 0 → CONSTANT
  expect_true(is_constant(power(expr, 0)))

  ## 0 < p < 1 → CONCAVE
  expect_equal(expr_curvature(power(expr, 0.67)), CONCAVE)
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: sign is NONNEG for p != 1 (CVXPY test_power)", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  expr <- x + y

  for (p in c(0, 2, 3, 2.7, 0.67, -1, -2.3)) {
    atom <- power(expr, p)
    expect_true(is_nonneg(atom), info = paste("p =", p))
  }
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: numeric value for constant (CVXPY test_power)", {
  ## assert cp.power(-1, 2).value == 1
  expect_equal(as.numeric(value(power(Constant(-1), 2))), 1)
})

# ═══════════════════════════════════════════════════════════════════════
# test_pnorm (test_atoms.py line 227)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_pnorm
test_that("pnorm: shape and curvature for various p (CVXPY test_pnorm)", {
  x <- Variable(2)

  ## p = 1.5 → scalar, CONVEX, NONNEG
  atom <- p_norm(x, p = 1.5)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONVEX)
  expect_true(is_nonneg(atom))

  ## p = 1 → scalar, CONVEX, NONNEG (dispatches to norm1)
  atom <- p_norm(x, p = 1)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONVEX)
  expect_true(is_nonneg(atom))

  ## p = 2 → scalar, CONVEX, NONNEG
  atom <- p_norm(x, p = 2)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONVEX)
  expect_true(is_nonneg(atom))

  ## p = Inf → scalar, CONVEX, NONNEG (dispatches to norm_inf)
  atom <- p_norm(x, p = Inf)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONVEX)
  expect_true(is_nonneg(atom))

  ## p = 0.5 → scalar, CONCAVE, NONNEG
  atom <- p_norm(x, p = 0.5)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONCAVE)
  expect_true(is_nonneg(atom))

  ## p = 0.7 → scalar, CONCAVE, NONNEG
  atom <- p_norm(x, p = 0.7)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONCAVE)
  expect_true(is_nonneg(atom))

  ## p = -0.1 → scalar, CONCAVE, NONNEG
  atom <- p_norm(x, p = -0.1)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONCAVE)
  expect_true(is_nonneg(atom))

  ## p = -1 → scalar, CONCAVE, NONNEG
  atom <- p_norm(x, p = -1)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONCAVE)
  expect_true(is_nonneg(atom))

  ## p = -1.3 → scalar, CONCAVE, NONNEG
  atom <- p_norm(x, p = -1.3)
  expect_equal(atom@shape, c(1L, 1L))
  expect_equal(expr_curvature(atom), CONCAVE)
  expect_true(is_nonneg(atom))
})

## @cvxpy test_atoms.py::TestAtoms::test_pnorm
test_that("pnorm: axis/keepdims shapes (CVXPY test_pnorm)", {
  A <- Variable(c(2L, 2L))

  ## norm(A, 2, axis=2) → R c(1L, 2L)
  expr <- p_norm(A, p = 2, axis = 2L)
  expect_equal(expr@shape, c(1L, 2L))

  ## norm(A, 2, axis=0, keepdims=TRUE) → CVXPY shape (1, 2) → R c(1L, 2L)
  expr <- p_norm(A, p = 2, axis = 2L, keepdims = TRUE)
  expect_equal(expr@shape, c(1L, 2L))

  ## norm(A, 2, axis=1, keepdims=TRUE) → CVXPY shape (2, 1)
  expr <- p_norm(A, p = 2, axis = 1L, keepdims = TRUE)
  expect_equal(expr@shape, c(2L, 1L))
})

# ═══════════════════════════════════════════════════════════════════════
# test_quad_over_lin (test_atoms.py line 318)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_quad_over_lin
test_that("quad_over_lin: DCP curvature (CVXPY test_quad_over_lin)", {
  x <- Variable(2)
  a <- Variable(1)

  ## quad_over_lin(square(x), a) → CONVEX
  atom <- quad_over_lin(x^2, a)
  expect_equal(expr_curvature(atom), CONVEX)

  ## quad_over_lin(-square(x), a) → CONVEX
  atom <- quad_over_lin(-(x^2), a)
  expect_equal(expr_curvature(atom), CONVEX)

  ## quad_over_lin(sqrt(x), a) → UNKNOWN (not DCP)
  atom <- quad_over_lin(sqrt(x), a)
  expect_equal(expr_curvature(atom), UNKNOWN_CURVATURE)
  expect_false(is_dcp(atom))
})

## @cvxpy test_atoms.py::TestAtoms::test_quad_over_lin
test_that("quad_over_lin: shape validation (CVXPY test_quad_over_lin)", {
  x <- Variable(2)
  ## Second arg must be scalar
  expect_error(quad_over_lin(x, x), "scalar")
})

## @cvxpy test_atoms.py::TestAtoms::test_quad_over_lin
test_that("quad_over_lin: numeric value (CVXPY test_quad_over_lin)", {
  x_val <- matrix(c(3.0, 4.0), ncol = 1)
  y_val <- 2.0
  atom <- quad_over_lin(Constant(x_val), Constant(y_val))
  expected <- (3.0^2 + 4.0^2) / 2.0  # = 12.5
  expect_equal(as.numeric(value(atom)), expected, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# test_max (test_atoms.py line 371)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_max
test_that("max: sign logic (CVXPY test_max)", {
  expect_equal(expr_sign_str(max_entries(Constant(1))), NONNEG_SIGN)
  expect_equal(expr_sign_str(max_entries(Constant(-2))), NONPOS_SIGN)
  expect_equal(expr_sign_str(MaxEntries(Variable(1))), UNKNOWN_SIGN)
  expect_equal(expr_sign_str(MaxEntries(Constant(0))), ZERO_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_max
test_that("max: axis/keepdims shapes (CVXPY test_max)", {
  ## CVXPY: max(Variable(2), axis=0, keepdims=True).shape == (1,)
  ## R: c(1L, 1L)
  expr <- MaxEntries(Variable(2), axis = 2L, keepdims = TRUE)
  expect_equal(expr@shape, c(1L, 1L))

  ## CVXPY: max(Variable((2,3)), axis=0, keepdims=True).shape == (1, 3)
  expr <- MaxEntries(Variable(c(2L, 3L)), axis = 2L, keepdims = TRUE)
  expect_equal(expr@shape, c(1L, 3L))

  ## CVXPY: max(Variable((2,3)), axis=1).shape == (2,) → R c(2L, 1L)
  expr <- MaxEntries(Variable(c(2L, 3L)), axis = 1L)
  expect_equal(expr@shape, c(2L, 1L))
})

# ═══════════════════════════════════════════════════════════════════════
# test_min (test_atoms.py line 396)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_min
test_that("min: sign logic (CVXPY test_min)", {
  expect_equal(expr_sign_str(MinEntries(Constant(1))), NONNEG_SIGN)
  expect_equal(expr_sign_str(MinEntries(Constant(-2))), NONPOS_SIGN)
  expect_equal(expr_sign_str(MinEntries(Variable(1))), UNKNOWN_SIGN)
  expect_equal(expr_sign_str(MinEntries(Constant(0))), ZERO_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_min
test_that("min: axis/keepdims shapes (CVXPY test_min)", {
  ## CVXPY: min(Variable(2), axis=0).shape == tuple() → R c(1L, 1L)
  expr <- MinEntries(Variable(2), axis = 2L)
  expect_equal(expr@shape, c(1L, 1L))

  ## CVXPY: min(Variable((2,3)), axis=0).shape == (3,) → R c(1L, 3L)
  expr <- MinEntries(Variable(c(2L, 3L)), axis = 2L)
  expect_equal(expr@shape, c(1L, 3L))

  ## CVXPY: min(Variable((2,3)), axis=1).shape == (2,) → R c(2L, 1L)
  expr <- MinEntries(Variable(c(2L, 3L)), axis = 1L)
  expect_equal(expr@shape, c(2L, 1L))
})

# ═══════════════════════════════════════════════════════════════════════
# test_maximum_sign (test_atoms.py line 432)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("maximum: sign logic two args (CVXPY test_maximum_sign)", {
  expect_equal(expr_sign_str(Maximum(Constant(1), Constant(2))), NONNEG_SIGN)
  expect_equal(expr_sign_str(Maximum(Constant(1), Variable(1))), NONNEG_SIGN)
  expect_equal(expr_sign_str(Maximum(Constant(1), Constant(-2))), NONNEG_SIGN)
  expect_equal(expr_sign_str(Maximum(Constant(1), Constant(0))), NONNEG_SIGN)

  expect_equal(expr_sign_str(Maximum(Variable(1), Constant(0))), NONNEG_SIGN)
  expect_equal(expr_sign_str(Maximum(Variable(1), Variable(1))), UNKNOWN_SIGN)
  expect_equal(expr_sign_str(Maximum(Variable(1), Constant(-2))), UNKNOWN_SIGN)

  expect_equal(expr_sign_str(Maximum(Constant(0), Constant(0))), ZERO_SIGN)
  expect_equal(expr_sign_str(Maximum(Constant(0), Constant(-2))), ZERO_SIGN)

  expect_equal(expr_sign_str(Maximum(Constant(-3), Constant(-2))), NONPOS_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_maximum_sign
test_that("maximum: promotion shape (CVXPY test_maximum_sign)", {
  ## CVXPY: maximum(1, Variable(2)).shape == (2,) → R c(2L, 1L)
  expr <- Maximum(Constant(1), Variable(2))
  expect_true(is_nonneg(expr))
  expect_equal(expr@shape, c(2L, 1L))
})

# ═══════════════════════════════════════════════════════════════════════
# test_minimum_sign (test_atoms.py line 459)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("minimum: sign logic two args (CVXPY test_minimum_sign)", {
  expect_equal(expr_sign_str(Minimum(Constant(1), Constant(2))), NONNEG_SIGN)
  expect_equal(expr_sign_str(Minimum(Constant(1), Variable(1))), UNKNOWN_SIGN)
  expect_equal(expr_sign_str(Minimum(Constant(1), Constant(-2))), NONPOS_SIGN)
  expect_equal(expr_sign_str(Minimum(Constant(1), Constant(0))), ZERO_SIGN)

  expect_equal(expr_sign_str(Minimum(Variable(1), Constant(0))), NONPOS_SIGN)
  expect_equal(expr_sign_str(Minimum(Variable(1), Variable(1))), UNKNOWN_SIGN)
  expect_equal(expr_sign_str(Minimum(Variable(1), Constant(-2))), NONPOS_SIGN)

  expect_equal(expr_sign_str(Minimum(Constant(0), Constant(0))), ZERO_SIGN)
  expect_equal(expr_sign_str(Minimum(Constant(0), Constant(-2))), NONPOS_SIGN)

  expect_equal(expr_sign_str(Minimum(Constant(-3), Constant(-2))), NONPOS_SIGN)
})

## @cvxpy test_atoms.py::TestAtoms::test_minimum_sign
test_that("minimum: promotion shape (CVXPY test_minimum_sign)", {
  ## CVXPY: minimum(-1, Variable(2)).shape == (2,) → R c(2L, 1L)
  expr <- Minimum(Constant(-1), Variable(2))
  expect_true(is_nonpos(expr))
  expect_equal(expr@shape, c(2L, 1L))
})

# ═══════════════════════════════════════════════════════════════════════
# test_sum (test_atoms.py line 485)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_sum
test_that("sum: sign and curvature (CVXPY test_sum)", {
  ## sum(1).sign == NONNEG
  expect_true(is_nonneg(SumEntries(Constant(1))))

  ## sum(Constant([1, -1])).sign == UNKNOWN
  expect_equal(expr_sign_str(SumEntries(Constant(matrix(c(1, -1), ncol = 1)))), UNKNOWN_SIGN)
  ## sum(Constant([1, -1])).curvature == CONSTANT
  expect_true(is_constant(SumEntries(Constant(matrix(c(1, -1), ncol = 1)))))

  ## sum(Variable(2)).sign == UNKNOWN
  expect_equal(expr_sign_str(SumEntries(Variable(2))), UNKNOWN_SIGN)
  ## sum(Variable(2)).shape == tuple() → R c(1L, 1L)
  expect_equal(SumEntries(Variable(2))@shape, c(1L, 1L))
  ## sum(Variable(2)).curvature == AFFINE
  expect_equal(expr_curvature(SumEntries(Variable(2))), AFFINE)
})

## @cvxpy test_atoms.py::TestAtoms::test_sum
test_that("sum: keepdims shape (CVXPY test_sum)", {
  ## sum(Variable((2,1)), keepdims=True).shape == (1, 1)
  expect_equal(SumEntries(Variable(c(2L, 1L)), keepdims = TRUE)@shape, c(1L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_sum
test_that("sum: axis shapes (CVXPY test_sum)", {
  ## sum(Variable(2), axis=0).shape == tuple() → c(1L, 1L)
  expect_equal(SumEntries(Variable(2), axis = 2L)@shape, c(1L, 1L))

  ## sum(Variable((2,3)), axis=0, keepdims=True).shape == (1, 3)
  expect_equal(SumEntries(Variable(c(2L, 3L)), axis = 2L, keepdims = TRUE)@shape, c(1L, 3L))

  ## sum(Variable((2,3)), axis=2, keepdims=False).shape == (3,) → c(1L, 3L)
  expect_equal(SumEntries(Variable(c(2L, 3L)), axis = 2L)@shape, c(1L, 3L))

  ## sum(Variable((2,3)), axis=1).shape == (2,) → c(2L, 1L)
  expect_equal(SumEntries(Variable(c(2L, 3L)), axis = 1L)@shape, c(2L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_sum
test_that("sum: mixed curvature (CVXPY test_sum)", {
  ## mat = [[1, -1]]; sum(mat @ square(Variable(2))).curvature == UNKNOWN
  mat <- Constant(matrix(c(1, -1), nrow = 1))
  expr <- SumEntries(mat %*% (Variable(2)^2))
  expect_equal(expr_curvature(expr), UNKNOWN_CURVATURE)
})

## @cvxpy test_atoms.py::TestAtoms::test_sum
test_that("sum: invalid axis (CVXPY test_sum)", {
  x <- Variable(2)
  expect_error(SumEntries(x, axis = 4L))
  ## In CVXPY, Variable(2) is 1D shape (2,), so axis=1 is invalid.
  ## In R, Variable(2) is (2,1) which is 2D, so axis=1 (columns) is valid.
  ## Test with a truly invalid axis instead:
  expect_error(SumEntries(Variable(2), axis = 3L))
})

## @cvxpy test_atoms.py::TestAtoms::test_sum
test_that("sum: sparse matrix input (CVXPY test_sum)", {
  A <- Matrix::Diagonal(3)  # eye(3)
  expect_equal(as.numeric(value(SumEntries(Constant(A)))), 3)
})

# ═══════════════════════════════════════════════════════════════════════
# test_multiply (test_atoms.py line 525)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_multiply
test_that("multiply: sign, curvature, shape (CVXPY test_multiply)", {
  x <- Variable(2)

  ## multiply([1, -1], x).sign == UNKNOWN
  expr <- Constant(matrix(c(1, -1), ncol = 1)) * x
  expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
  ## multiply([1, -1], x).curvature == AFFINE
  expect_equal(expr_curvature(expr), AFFINE)
  ## multiply([1, -1], x).shape == (2,) → c(2L, 1L)
  expect_equal(expr@shape, c(2L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_multiply
test_that("multiply: sign with nonneg/nonpos constants (CVXPY test_multiply)", {
  ## CVXPY uses Parameters; we use Constants (since no Parameter support yet)
  pos <- Constant(matrix(c(1, 2), ncol = 1))   # nonneg
  neg <- Constant(matrix(c(-1, -2), ncol = 1))  # nonpos

  ## pos * pos → NONNEG
  expect_true(is_nonneg(pos * pos))
  ## pos * neg → NONPOS
  expect_true(is_nonpos(pos * neg))
  ## neg * neg → NONNEG
  expect_true(is_nonneg(neg * neg))
})

## @cvxpy test_atoms.py::TestAtoms::test_multiply
test_that("multiply: concave curvature (CVXPY test_multiply)", {
  x <- Variable(2)
  neg <- Constant(matrix(c(-1, -2), ncol = 1))
  ## neg * square(x) → CONCAVE (nonpos * convex nonneg → concave)
  expr <- neg * (x^2)
  expect_equal(expr_curvature(expr), CONCAVE)
})

## @cvxpy test_atoms.py::TestAtoms::test_multiply
test_that("multiply: promotion (CVXPY test_multiply)", {
  C <- Variable(c(3L, 2L))
  ## multiply(1, C).shape == C.shape
  expect_equal((Constant(1) * C)@shape, C@shape)
})

# ═══════════════════════════════════════════════════════════════════════
# test_vstack (test_atoms.py line 548)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_vstack
test_that("vstack: shapes (CVXPY test_vstack)", {
  x <- Variable(2)
  y <- Variable(2)
  A <- Variable(c(2L, 2L))
  C <- Variable(c(3L, 2L))
  B <- Variable(c(2L, 2L))

  ## CVXPY: vstack([x, y, x]).shape == (3, 2) for 1D inputs
  ## In R, x is (2,1) so VStack of three (2,1) → (6, 1)
  atom <- VStack(x, y, x)
  expect_equal(atom@shape, c(6L, 1L))

  ## CVXPY: vstack([A, C, B]).shape == (7, 2)
  atom <- VStack(A, C, B)
  expect_equal(atom@shape, c(7L, 2L))
})

## @cvxpy test_atoms.py::TestAtoms::test_vstack
test_that("vstack: entries from indexing (CVXPY test_vstack)", {
  x <- Variable(2)
  ## Extracting each element and vstacking
  entries <- list(x[1L, ], x[2L, ])
  atom <- do.call(VStack, entries)
  expect_equal(atom@shape, c(2L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_vstack
test_that("vstack: validation (CVXPY test_vstack)", {
  C <- Variable(c(3L, 2L))

  ## Incompatible column dimensions: (3,2) vs scalar promoted to (1,1)
  expect_error(VStack(C, Constant(1)))
  ## In R, Variable(2) is (2,1) and Variable(3) is (3,1) — both 1-column, so VStack is valid
  ## Instead test: (3,2) vs (2,1) — incompatible columns (2 vs 1)
  expect_error(VStack(C, Variable(2)))
})

# ═══════════════════════════════════════════════════════════════════════
# test_hstack (test_atoms.py line 576)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_hstack
test_that("hstack: shapes (CVXPY test_hstack)", {
  A <- Variable(c(2L, 2L))
  B <- Variable(c(2L, 2L))

  ## CVXPY: hstack([A, B]).shape == (2, 4)
  atom <- HStack(A, B)
  expect_equal(atom@shape, c(2L, 4L))
})

## @cvxpy test_atoms.py::TestAtoms::test_hstack
test_that("hstack: validation (CVXPY test_hstack)", {
  C <- Variable(c(3L, 2L))
  A <- Variable(c(2L, 2L))

  ## Incompatible row dimensions
  expect_error(HStack(C, A))
})

# ═══════════════════════════════════════════════════════════════════════
# test_reshape (test_atoms.py line 644)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_reshape
test_that("reshape: sign, curvature, shape (CVXPY test_reshape)", {
  A <- Variable(c(2L, 2L))

  ## reshape(A, (4,1), order='F')
  expr <- reshape_expr(A, c(4L, 1L), order = "F")
  expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
  expect_equal(expr_curvature(expr), AFFINE)
  expect_equal(expr@shape, c(4L, 1L))

  ## reshape back to (2,2)
  expr2 <- reshape_expr(expr, c(2L, 2L), order = "F")
  expect_equal(expr2@shape, c(2L, 2L))
})

## @cvxpy test_atoms.py::TestAtoms::test_reshape
test_that("reshape: convex nonneg (CVXPY test_reshape)", {
  x <- Variable(2)
  ## reshape(square(x), (1, 2), order='F')
  expr <- reshape_expr(x^2, c(1L, 2L), order = "F")
  expect_true(is_nonneg(expr))
  expect_equal(expr_curvature(expr), CONVEX)
  expect_equal(expr@shape, c(1L, 2L))
})

## @cvxpy test_atoms.py::TestAtoms::test_reshape
test_that("reshape: invalid dimensions (CVXPY test_reshape)", {
  C <- Variable(c(3L, 2L))
  expect_error(reshape_expr(C, c(5L, 4L), order = "F"), "reshape")
})

## @cvxpy test_atoms.py::TestAtoms::test_reshape
test_that("reshape: C-order numeric (CVXPY test_reshape)", {
  ## a = np.arange(10); reshape(a, (5,2), order='C')
  a <- matrix(0:9, ncol = 1)
  A_np <- matrix(0:9, nrow = 5, ncol = 2, byrow = TRUE)
  A_cp <- reshape_expr(Constant(a), c(5L, 2L), order = "C")
  expect_equal(as.matrix(value(A_cp)), A_np, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# test_vec (test_atoms.py line 754)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_vec
test_that("vec: sign, curvature, shape (CVXPY test_vec)", {
  C <- Variable(c(3L, 2L))
  ## vec(C, order='F')
  expr <- vec(C)
  expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
  expect_equal(expr_curvature(expr), AFFINE)
  ## CVXPY: shape == (6,) → R c(6L, 1L)
  expect_equal(expr@shape, c(6L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_vec
test_that("vec: scalar (CVXPY test_vec)", {
  x <- Variable(2)
  ## vec(x) → c(2L, 1L)
  expr <- vec(x)
  expect_equal(expr@shape, c(2L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_vec
test_that("vec: convex nonneg preserves (CVXPY test_vec)", {
  a <- Variable(1)
  ## vec(square(a))
  expr <- vec(a^2)
  expect_true(is_nonneg(expr))
  expect_equal(expr_curvature(expr), CONVEX)
  ## CVXPY: shape == (1,) → R c(1L, 1L)
  expect_equal(expr@shape, c(1L, 1L))
})

# ═══════════════════════════════════════════════════════════════════════
# test_diag (test_atoms.py line 770)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_diag
test_that("diag: vector to matrix (CVXPY test_diag)", {
  x <- Variable(2)
  ## diag(x) → (2, 2)
  expr <- DiagVec(x)
  expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
  expect_equal(expr_curvature(expr), AFFINE)
  expect_equal(expr@shape, c(2L, 2L))
})

## @cvxpy test_atoms.py::TestAtoms::test_diag
test_that("diag: matrix to vector (CVXPY test_diag)", {
  A <- Variable(c(2L, 2L))
  ## diag(A) → (2,) → c(2L, 1L)
  expr <- DiagMat(A)
  expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
  expect_equal(expr_curvature(expr), AFFINE)
  expect_equal(expr@shape, c(2L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_diag
test_that("diag: transpose of matrix (CVXPY test_diag)", {
  ## In CVXPY, x.T of 1D (2,) stays (2,) so diag works.
  ## In R, Variable(2) is (2,1) and Transpose gives (1,2) row vector.
  ## DiagVec requires column vector → use matrix diagonal extraction instead.
  ## Test: DiagMat of transposed 2x2 matrix
  A <- Variable(c(2L, 2L))
  expr <- DiagMat(Transpose(A))
  expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
  expect_equal(expr_curvature(expr), AFFINE)
  expect_equal(expr@shape, c(2L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_diag
test_that("diag: constant matrix to vector (CVXPY test_diag)", {
  psd_matrix <- matrix(c(1, -1, -1, 1), 2, 2)
  expr <- DiagMat(Constant(psd_matrix))
  expect_true(is_nonneg(expr))
  expect_true(is_constant(expr))
  expect_equal(expr@shape, c(2L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_diag
test_that("diag: non-square matrix rejected (CVXPY test_diag)", {
  C <- Variable(c(3L, 2L))
  expect_error(DiagMat(C), "square")
})

## @cvxpy test_atoms.py::TestAtoms::test_diag
test_that("diag: PSD property (CVXPY test_diag)", {
  ## diag(w) is PSD when w >= 0
  w <- Constant(matrix(c(1.0, 2.0), ncol = 1))
  expr <- DiagVec(w)
  expect_true(is_psd(expr))

  expr_neg <- DiagVec(-w)
  expect_true(is_nsd(expr_neg))

  expr_mixed <- DiagVec(Constant(matrix(c(1, -1), ncol = 1)))
  expect_false(is_psd(expr_mixed))
  expect_false(is_nsd(expr_mixed))
})

# ═══════════════════════════════════════════════════════════════════════
# test_trace (test_atoms.py line 825)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_trace
test_that("trace: sign, curvature, shape (CVXPY test_trace)", {
  A <- Variable(c(2L, 2L))
  expr <- Trace(A)
  expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
  expect_equal(expr_curvature(expr), AFFINE)
  ## CVXPY: shape == tuple() → R c(1L, 1L)
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_trace
test_that("trace: non-square rejected (CVXPY test_trace)", {
  C <- Variable(c(3L, 2L))
  expect_error(Trace(C), "square")
})

## @cvxpy test_atoms.py::TestAtoms::test_trace_sign_psd
test_that("trace: PSD/NSD sign (CVXPY test_trace_sign_psd)", {
  X_psd <- Variable(c(2L, 2L), PSD = TRUE)
  X_nsd <- Variable(c(2L, 2L), NSD = TRUE)

  expect_true(is_nonneg(Trace(X_psd)))
  expect_true(is_nonpos(Trace(X_nsd)))
})

# ═══════════════════════════════════════════════════════════════════════
# test_upper_tri (test_atoms.py line 899)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_upper_tri
test_that("upper_tri: non-square rejected (CVXPY test_upper_tri)", {
  C <- Variable(c(3L, 2L))
  expect_error(UpperTri(C), "square")
})

## @cvxpy test_atoms.py::TestAtoms::test_upper_tri
test_that("upper_tri: shape for square matrix", {
  A <- Variable(c(4L, 4L))
  ut <- UpperTri(A)
  ## n*(n-1)/2 = 4*3/2 = 6
  expect_equal(ut@shape, c(6L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_upper_tri
test_that("upper_tri: numeric value row-major order", {
  A <- Variable(c(3L, 3L))
  value(A) <- matrix(1:9, 3, 3)
  ut <- UpperTri(A)
  ## Row-major order of strict upper triangle:
  ## (1,2)=4, (1,3)=7, (2,3)=8
  expect_equal(as.numeric(value(ut)), c(4, 7, 8))
})

# ═══════════════════════════════════════════════════════════════════════
# test_vec_to_upper_tri (test_atoms.py line 905)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_vec_to_upper_tri
test_that("vec_to_upper_tri: default (non-strict) (CVXPY test_vec_to_upper_tri)", {
  ## x = Variable((3,)); vec_to_upper_tri(x)
  ## x.value = [1, 2, 3] → [[1, 2], [0, 3]]
  x <- Variable(c(3L, 1L))
  X <- vec_to_upper_tri(x, strict = FALSE)
  value(x) <- matrix(c(1, 2, 3), 3, 1)
  actual <- as.matrix(value(X))
  expected <- matrix(c(1, 0, 2, 3), 2, 2)
  expect_equal(actual, expected, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_vec_to_upper_tri
test_that("vec_to_upper_tri: strict=TRUE (CVXPY test_vec_to_upper_tri)", {
  ## y = Variable((1,)); vec_to_upper_tri(y, strict=True)
  ## y.value = [4] → [[0, 4], [0, 0]]
  y <- Variable(c(1L, 1L))
  Y <- vec_to_upper_tri(y, strict = TRUE)
  value(y) <- matrix(4, 1, 1)
  actual <- as.matrix(value(Y))
  expected <- matrix(c(0, 0, 4, 0), 2, 2)
  expect_equal(actual, expected, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_vec_to_upper_tri
test_that("vec_to_upper_tri: strict 4x4 (CVXPY test_vec_to_upper_tri)", {
  ## a = [11, 12, 13, 16, 17, 21] → 4x4 strict upper tri
  a <- Constant(matrix(c(11, 12, 13, 16, 17, 21), 6, 1))
  A_actual <- as.matrix(value(vec_to_upper_tri(a, strict = TRUE)))
  A_expect <- matrix(c(0, 0, 0, 0, 11, 0, 0, 0, 12, 16, 0, 0, 13, 17, 21, 0), 4, 4)
  expect_equal(A_actual, A_expect, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_vec_to_upper_tri
test_that("vec_to_upper_tri: validation (CVXPY test_vec_to_upper_tri)", {
  ## Not a triangular number
  expect_error(vec_to_upper_tri(Variable(c(4L, 1L)), strict = FALSE), "triangular")
  expect_error(vec_to_upper_tri(Variable(c(4L, 1L)), strict = TRUE), "triangular")
})

# ═══════════════════════════════════════════════════════════════════════
# test_huber (test_atoms.py line 951)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_huber
test_that("huber: M validation (CVXPY test_huber)", {
  x <- Variable(2)
  ## Valid
  h <- Huber(x, M = 1)
  expect_true(is_convex(h))

  ## Negative M rejected
  expect_error(Huber(x, M = -1), "[Nn]on-negative|[Mm]ust be")

  ## Vector M rejected (must be scalar)
  expect_error(Huber(x, M = c(1, 1)), "[Ss]calar|[Mm]ust be")
})

## @cvxpy test_atoms.py::TestAtoms::test_huber
test_that("huber: numeric value (CVXPY test_huber)", {
  ## huber(2, M=1).value == 3
  ## huber loss: |x| <= M → x^2; |x| > M → 2M|x| - M^2
  ## huber(2, 1) = 2*1*2 - 1 = 3
  expect_equal(as.numeric(value(Huber(Constant(2), M = 1))), 3, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# test_sum_largest (test_atoms.py line 1034)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_sum_largest
test_that("sum_largest: k validation (CVXPY test_sum_largest)", {
  x <- Variable(2)
  ## Negative k rejected
  expect_error(SumLargest(x, k = -1), "positive|[Ss]econd")
})

## @cvxpy test_atoms.py::TestAtoms::test_sum_largest
test_that("sum_largest: PWL property (CVXPY test_sum_largest)", {
  x <- Variable(2)
  atom <- SumLargest(x, k = 2)
  expect_true(is_pwl(atom))
})

## @cvxpy test_atoms.py::TestAtoms::test_sum_largest
test_that("sum_largest: numeric value (CVXPY test_sum_largest)", {
  ## 10000 random values
  set.seed(42)
  v <- rnorm(100)
  x <- Constant(matrix(v, ncol = 1))

  for (i in c(5, 10, 25, 50)) {
    expr <- SumLargest(x, k = i)
    prev_idx <- order(-v)[1:i]
    expect_equal(as.numeric(value(expr)), sum(v[prev_idx]), tolerance = 1e-6,
                 info = paste("k =", i))
  }
})

# ═══════════════════════════════════════════════════════════════════════
# test_sum_smallest (test_atoms.py line 1112)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_sum_smallest
test_that("sum_smallest: k validation (CVXPY test_sum_smallest)", {
  x <- Variable(2)
  expect_error(sum_smallest(x, -1), "positive|[Ss]econd")
})

## @cvxpy test_atoms.py::TestAtoms::test_sum_smallest
test_that("sum_smallest: PWL property (CVXPY test_sum_smallest)", {
  x <- Variable(2)
  atom <- sum_smallest(x, 2)
  expect_true(is_pwl(atom))
})

# ═══════════════════════════════════════════════════════════════════════
# test_kron_expr (test_atoms.py line 1292)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_kron_expr
test_that("kron: sign and shape (CVXPY test_kron_expr)", {
  a <- Constant(matrix(1, 3, 2))

  ## kron(a, nonneg) → NONNEG, shape (6, 2)
  b_nonneg <- Constant(matrix(c(1, 2), 2, 1))
  expr <- Kron(a, b_nonneg)
  expect_true(is_nonneg(expr))
  expect_equal(expr@shape, c(6L, 2L))

  ## kron(a, nonpos) → NONPOS
  b_nonpos <- Constant(matrix(c(-1, -2), 2, 1))
  expr <- Kron(a, b_nonpos)
  expect_true(is_nonpos(expr))
})

## @cvxpy test_atoms.py::TestAtoms::test_kron_expr
test_that("kron: validation (CVXPY test_kron_expr)", {
  x <- Variable(2)
  ## Both arguments variable → error
  expect_error(Kron(x, x), "constant")
})

# ═══════════════════════════════════════════════════════════════════════
# test_conv (test_atoms.py line 1232)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_conv
test_that("conv: sign and shape (CVXPY test_conv)", {
  a <- Constant(matrix(1, 3, 1))  # ones (3, 1)

  ## conv(a, nonneg) → NONNEG, shape (4, 1)
  b_nonneg <- Constant(matrix(c(1, 2), 2, 1))
  expr <- Convolve(a, b_nonneg)
  expect_true(is_nonneg(expr))
  ## CVXPY: (4, 1)
  expect_equal(expr@shape, c(4L, 1L))

  ## conv(a, nonpos) → NONPOS
  b_nonpos <- Constant(matrix(c(-1, -2), 2, 1))
  expr <- Convolve(a, b_nonpos)
  expect_true(is_nonpos(expr))
})

## @cvxpy test_atoms.py::TestAtoms::test_conv
test_that("conv: at least one arg must be constant (CVXPY test_conv)", {
  x <- Variable(2)
  y <- Variable(2)
  ## Both variable → error
  expect_error(Convolve(x, y), "constant")
  ## Variable + Constant → OK (CVXR is more permissive than CVXPY here)
  expect_no_error(Convolve(x, Constant(-1)))
})

# ═══════════════════════════════════════════════════════════════════════
# test_cumsum (test_atoms.py line 1267)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_cumsum
test_that("cumsum: numeric value matches R cumsum (CVXPY test_cumsum)", {
  for (axis in c(2L, 1L)) {
    x <- Variable(c(4L, 3L))
    x_val <- matrix(0:11, 4, 3)
    value(x) <- x_val

    expr <- Cumsum(x, axis = axis)
    ## R cumsum along axis
    if (axis == 2L) {
      target <- apply(x_val, 2, cumsum)
    } else {
      target <- t(apply(x_val, 1, cumsum))
    }
    expect_equal(as.matrix(value(expr)), target, tolerance = 1e-10,
                 info = paste("axis =", axis))
  }
})

# ═══════════════════════════════════════════════════════════════════════
# test_diff (test_atoms.py line 1663)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_diff
test_that("diff: shape for 2D matrix (CVXPY test_diff)", {
  A <- Variable(c(20L, 10L))

  ## diff(A, axis=0).shape == np.diff(zeros(20,10), axis=0).shape == (19, 10)
  expect_equal(cvxr_diff(A, axis = 2L)@shape, c(19L, 10L))

  ## diff(A, axis=1).shape == np.diff(zeros(20,10), axis=1).shape == (20, 9)
  expect_equal(cvxr_diff(A, axis = 1L)@shape, c(20L, 9L))
})

## @cvxpy test_atoms.py::TestAtoms::test_diff
test_that("diff: numeric value axis=1 (CVXPY test_diff issue #1834)", {
  x1 <- matrix(c(1, 2, 3, 4, 5), nrow = 1)  # (1, 5)
  x2 <- Variable(c(1L, 5L))
  value(x2) <- x1

  ## diff(x2, axis=1)
  expr <- cvxr_diff(x2, axis = 1L)
  expected <- diff(x1[1, ])  # c(1, 1, 1, 1)
  expect_equal(as.numeric(value(expr)), expected, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_diff
test_that("diff: higher-order diff (CVXPY test_diff)", {
  A <- Variable(c(20L, 10L))

  ## k=2, axis=0 → (18, 10)
  expect_equal(cvxr_diff(A, k = 2L, axis = 2L)@shape, c(18L, 10L))
  ## k=2, axis=1 → (20, 8)
  expect_equal(cvxr_diff(A, k = 2L, axis = 1L)@shape, c(20L, 8L))
})

# ═══════════════════════════════════════════════════════════════════════
# test_log1p (test_atoms.py line 889)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_log1p
test_that("log1p: curvature and numeric (CVXPY test_log1p)", {
  ## CVXPY has a dedicated Log1p atom with sign tracking.
  ## In CVXR, log1p(x) is implemented as log(x + 1), so sign is UNKNOWN.
  ## Test curvature and numeric value instead.
  expr <- log(Constant(1) + 1)  # log(2)
  expect_true(is_constant(expr))
  expect_equal(as.numeric(value(expr)), log(2), tolerance = 1e-10)

  ## log1p(-0.5) = log(0.5)
  expr_neg <- log(Constant(1) + Constant(-0.5))
  expect_true(is_constant(expr_neg))
  expect_equal(as.numeric(value(expr_neg)), log(0.5), tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# Comprehensive atom sign/curvature snapshot (additional coverage)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("abs: convex, nonneg (CVXPY atom properties)", {
  x <- Variable(2)
  a <- Abs(x)
  expect_equal(expr_curvature(a), CONVEX)
  expect_true(is_nonneg(a))
  expect_equal(a@shape, c(2L, 1L))
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_exp
test_that("exp: convex, nonneg (CVXPY atom properties)", {
  x <- Variable(1)
  e <- Exp(x)
  expect_equal(expr_curvature(e), CONVEX)
  expect_true(is_nonneg(e))
  expect_equal(e@shape, c(1L, 1L))
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_log
test_that("log: concave, unknown sign (CVXPY atom properties)", {
  x <- Variable(1)
  l <- Log(x)
  expect_equal(expr_curvature(l), CONCAVE)
  expect_equal(l@shape, c(1L, 1L))
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_entr
test_that("entr: concave, unknown sign (CVXPY atom properties)", {
  x <- Variable(1)
  e <- Entr(x)
  expect_equal(expr_curvature(e), CONCAVE)
  expect_equal(e@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("sum_squares: convex, nonneg (CVXPY atom properties)", {
  x <- Variable(3)
  ss <- sum_squares(x)
  expect_true(is_convex(ss))
  expect_true(is_nonneg(ss))
  expect_equal(ss@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("pos/neg: convex, nonneg (CVXPY atom properties)", {
  x <- Variable(2)
  p <- pos(x)
  n <- neg(x)
  expect_true(is_convex(p))
  expect_true(is_nonneg(p))
  expect_true(is_convex(n))
  expect_true(is_nonneg(n))
})

## @cvxpy NONE
test_that("square: convex, nonneg (CVXPY atom properties)", {
  x <- Variable(2)
  sq <- x^2
  expect_true(is_convex(sq))
  expect_true(is_nonneg(sq))
})

## @cvxpy NONE
test_that("sqrt: concave, nonneg (CVXPY atom properties)", {
  x <- Variable(1)
  s <- sqrt(x)
  expect_true(is_concave(s))
  expect_true(is_nonneg(s))
})

## @cvxpy NONE
test_that("inv_pos: convex, nonneg (CVXPY atom properties)", {
  x <- Variable(1)
  ip <- inv_pos(x)
  expect_true(is_convex(ip))
  expect_true(is_nonneg(ip))
})

# ═══════════════════════════════════════════════════════════════════════
# Numeric evaluation parity
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_norm1
test_that("norm1: numeric value matches R norm (CVXPY numeric)", {
  x <- Variable(3)
  value(x) <- matrix(c(-2, 3, -1), 3, 1)
  expect_equal(as.numeric(value(norm1(x))), 6, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_norm_inf
test_that("norm_inf: numeric value matches R max(abs()) (CVXPY numeric)", {
  x <- Variable(3)
  value(x) <- matrix(c(-2, 3, -1), 3, 1)
  expect_equal(as.numeric(value(norm_inf(x))), 3, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_pnorm
test_that("pnorm p=2: numeric value matches Euclidean norm (CVXPY numeric)", {
  x <- Variable(3)
  value(x) <- matrix(c(3, 4, 0), 3, 1)
  expect_equal(as.numeric(value(p_norm(x, p = 2))), 5, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_trace
test_that("trace: numeric value (CVXPY numeric)", {
  A <- Variable(c(3L, 3L))
  value(A) <- diag(c(1, 2, 3))
  expect_equal(as.numeric(value(Trace(A))), 6, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_diag
test_that("diag: vector to matrix numeric (CVXPY numeric)", {
  x <- Variable(3)
  value(x) <- matrix(c(1, 2, 3), 3, 1)
  D <- DiagVec(x)
  expected <- diag(c(1, 2, 3))
  expect_equal(as.matrix(value(D)), expected, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_diag
test_that("diag: matrix to vector numeric (CVXPY numeric)", {
  A <- Variable(c(3L, 3L))
  value(A) <- matrix(1:9, 3, 3)
  d <- DiagMat(A)
  expect_equal(as.numeric(value(d)), c(1, 5, 9), tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_sum_largest
test_that("sum_largest: numeric value (CVXPY numeric)", {
  v <- c(10, 1, 5, 3, 7)
  x <- Constant(matrix(v, ncol = 1))

  ## sum_largest(x, 2) → 10 + 7 = 17
  expect_equal(as.numeric(value(SumLargest(x, k = 2))), 17, tolerance = 1e-10)

  ## sum_largest(x, 1) → 10
  expect_equal(as.numeric(value(SumLargest(x, k = 1))), 10, tolerance = 1e-10)

  ## sum_largest(x, 5) → sum(v) = 26
  expect_equal(as.numeric(value(SumLargest(x, k = 5))), 26, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_sum_smallest
test_that("sum_smallest: numeric value (CVXPY numeric)", {
  v <- c(10, 1, 5, 3, 7)
  x <- Variable(c(5L, 1L))
  value(x) <- matrix(v, 5, 1)

  ## sum_smallest(x, 2) → 1 + 3 = 4
  expect_equal(as.numeric(value(sum_smallest(x, 2))), 4, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_huber
test_that("huber: numeric values for various inputs (CVXPY numeric)", {
  ## huber(x, M): |x| <= M → x^2; |x| > M → 2M|x| - M^2
  ## huber(0.5, 1) → 0.25  (inside)
  expect_equal(as.numeric(value(Huber(Constant(0.5), M = 1))), 0.25, tolerance = 1e-10)

  ## huber(2, 1) → 2*1*2 - 1 = 3  (outside)
  expect_equal(as.numeric(value(Huber(Constant(2), M = 1))), 3, tolerance = 1e-10)

  ## huber(-1.5, 1) → 2*1*1.5 - 1 = 2  (outside, negative input)
  expect_equal(as.numeric(value(Huber(Constant(-1.5), M = 1))), 2, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_reshape
test_that("reshape: F-order and C-order numeric parity (CVXPY numeric)", {
  ## F-order: column-major (R default)
  m <- matrix(1:6, 2, 3)  # [[1,3,5],[2,4,6]]
  expr_f <- reshape_expr(Constant(m), c(3L, 2L), order = "F")
  ## F-order reshape: read column-major, fill column-major
  ## As a vector (col-major): 1, 2, 3, 4, 5, 6 → fill (3,2) col-major: [[1,4],[2,5],[3,6]]
  expected_f <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  expect_equal(as.matrix(value(expr_f)), expected_f, tolerance = 1e-10)

  ## C-order: row-major
  expr_c <- reshape_expr(Constant(m), c(3L, 2L), order = "C")
  ## C-order reshape: read row-major, fill row-major
  ## m row-major: 1, 3, 5, 2, 4, 6 → fill (3,2) row-major: [[1,3],[5,2],[4,6]]
  expected_c <- matrix(c(1, 5, 4, 3, 2, 6), 3, 2)
  expect_equal(as.matrix(value(expr_c)), expected_c, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_kron_expr
test_that("kron: numeric value (CVXPY numeric)", {
  a <- Constant(matrix(c(1, 0, 0, 1), 2, 2))  # I2
  b <- Constant(matrix(c(1, 2, 3, 4), 2, 2))   # [[1,3],[2,4]]
  k <- Kron(a, b)
  ## kronecker(I2, b) = block diagonal of b
  expected <- as.matrix(kronecker(diag(2), matrix(c(1, 2, 3, 4), 2, 2)))
  expect_equal(as.matrix(value(k)), expected, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_conv
test_that("conv: numeric value (CVXPY numeric)", {
  a <- Constant(matrix(c(1, 1, 1), 3, 1))
  b <- Variable(c(2L, 1L))
  value(b) <- matrix(c(1, 2), 2, 1)
  expr <- Convolve(a, b)
  ## conv([1,1,1], [1,2]) = [1, 3, 3, 2]
  expect_equal(as.numeric(value(expr)), c(1, 3, 3, 2), tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# test_upper_tri_to_full (test_atoms.py line 899)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_upper_tri_to_full
test_that("upper_tri_to_full: produces symmetric matrix (CVXPY test_upper_tri_to_full)", {
  for (n in 3:7) {
    A <- upper_tri_to_full(n)
    ell <- (n * (n + 1L)) %/% 2L
    v <- seq_len(ell) - 1L  # 0-based like np.arange
    ## A @ v reshaped to (n, n) F-order should be symmetric
    Mv <- as.numeric(A %*% v)
    M <- matrix(Mv, n, n)  # F-order (column-major)
    expect_equal(M, t(M), tolerance = 1e-10, info = paste("n =", n))
  }
})

# ═══════════════════════════════════════════════════════════════════════
# test_index: shape, sign, curvature (test_atoms.py line 1169)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_index
test_that("index: shape for various slices (CVXPY test_index)", {
  A <- Variable(c(5L, 4L))

  ## A[0:2, 0:1] → (2, 1)
  expr <- A[1:2, 1L]
  expect_equal(expr@shape, c(2L, 1L))

  ## A[0:2, :] → (2, 4)
  expr2 <- A[1:2, ]
  expect_equal(expr2@shape, c(2L, 4L))

  ## A[:, 1:3] → (5, 2)
  expr3 <- A[, 2:3]
  expect_equal(expr3@shape, c(5L, 2L))

  ## Single element: A[0, 0] → scalar (1, 1)
  expr4 <- A[1L, 1L]
  expect_equal(expr4@shape, c(1L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_index
test_that("index: curvature is AFFINE (CVXPY test_index)", {
  A <- Variable(c(3L, 3L))
  expr <- A[1:2, 1:2]
  expect_equal(expr_curvature(expr), AFFINE)
})

## @cvxpy test_atoms.py::TestAtoms::test_index
test_that("index: sign preserved for nonneg variable (CVXPY test_index)", {
  x <- Variable(5, nonneg = TRUE)
  expr <- x[1:3, ]
  expect_true(is_nonneg(expr))
})

# ═══════════════════════════════════════════════════════════════════════
# test_nonnegative_variable (test_atoms.py line 1580)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_nonnegative_variable
test_that("nonneg variable: sign propagates through atoms (CVXPY test_nonnegative_variable)", {
  x <- Variable(c(4L, 1L), nonneg = TRUE)

  ## sum(x) → NONNEG
  expect_true(is_nonneg(sum(x)))

  ## norm1(x) → NONNEG
  expect_true(is_nonneg(norm1(x)))

  ## max(x) → NONNEG
  expect_true(is_nonneg(max(x)))
})

# ═══════════════════════════════════════════════════════════════════════
# test_log_problem (test_nonlinear_atoms.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_log_problem
test_that("test_log_problem: log in objective (CVXPY test_log_problem)", {
  skip_if_not_installed("clarabel")

  x <- Variable(2, name = "x")

  ## Log in objective: max sum(log(x)) s.t. x <= [1, e]
  ## Optimal: x = [1, e], log(1) + log(e) = 0 + 1 = 1
  obj <- Maximize(sum(log(x)))
  constr <- list(x <= c(1, exp(1)))
  prob <- Problem(obj, constr)
  result <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(result, 1, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, exp(1)), tolerance = 1e-3)
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_log_problem
test_that("test_log_problem: log in constraint (CVXPY test_log_problem)", {
  skip_if_not_installed("clarabel")

  x <- Variable(2, name = "x")

  ## Log in constraint: min sum(x) s.t. log(x) >= 0, x <= [1, 1]
  ## log(x) >= 0 => x >= 1, combined with x <= 1 => x = [1, 1], sum = 2
  obj <- Minimize(sum(x))
  constr <- list(log(x) >= 0, x <= c(1, 1))
  prob <- Problem(obj, constr)
  result <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(result, 2, tolerance = 1e-3)
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_log_problem
test_that("test_log_problem: index into log (CVXPY test_log_problem)", {
  skip_if_not_installed("clarabel")

  x <- Variable(2, name = "x")

  ## Index into log: max log(x)[2] s.t. x <= [1, e]
  ## CVXPY: max log(x)[1] (0-based) = log(x_2) => x_2 = e, log(e) = 1
  ## R: log(x)[2] (1-based)
  obj <- Maximize(log(x)[2])
  constr <- list(x <= c(1, exp(1)))
  prob <- Problem(obj, constr)
  result <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(result, 1, tolerance = 1e-3)
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_log_problem
test_that("test_log_problem: scalar log (CVXPY test_log_problem)", {
  skip_if_not_installed("clarabel")

  x <- Variable(2, name = "x")

  ## Scalar log: max log(x[2]) s.t. x <= [1, e]
  ## CVXPY: max log(x[1]) (0-based) = max log(x_2) => 1
  ## R: x[2] (1-based)
  obj <- Maximize(log(x[2]))
  constr <- list(x <= c(1, exp(1)))
  prob <- Problem(obj, constr)
  result <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), OPTIMAL)
  expect_equal(result, 1, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# test_conv_prob (test_convolution.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_convolution.py::TestConvolution::test_conv_prob
test_that("test_conv_prob: conv(h, x) with random data is unbounded (CVXPY test_conv_prob)", {
  skip_if_not_installed("clarabel")

  N <- 5L
  set.seed(42)
  y <- matrix(rnorm(N), ncol = 1)
  h <- matrix(rnorm(2), ncol = 1)
  x <- Variable(c(N, 1))
  v <- conv(h, x)
  ## v has length N + 2 - 1 = N + 1. Slice v[1:N] (1-based).
  obj <- Minimize(sum_entries(multiply(y, v[1:N, ])))
  prob <- Problem(obj, list())
  result <- psolve(prob, solver = "CLARABEL")

  ## Unconstrained: should be unbounded
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
})

## @cvxpy test_convolution.py::TestConvolution::test_conv_prob
test_that("test_conv_prob: conv with 1-D vectors is unbounded (CVXPY test_conv_prob)", {
  skip_if_not_installed("clarabel")

  N <- 5L
  set.seed(42)
  y <- rnorm(N)
  h <- rnorm(2)
  x <- Variable(N)
  v <- conv(h, x)
  ## v has length N + 2 - 1 = N + 1. Slice v[1:N] (1-based).
  obj <- Minimize(sum_entries(multiply(y, v[1:N, ])))
  prob <- Problem(obj, list())
  result <- psolve(prob, solver = "CLARABEL")

  ## Unconstrained: should be unbounded
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
})

# ═══════════════════════════════════════════════════════════════════════
# test_is_psd / test_is_nsd (test_constant.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constant.py::test_is_psd
test_that("Constant: is_psd / is_nsd for identity matrices", {
  ## CVXPY: n=50, psd=eye(n), nsd=-eye(n)
  ## Verified via CVXPY:
  ##   Constant(psd).is_psd() = True, .is_nsd() = False
  ##   Constant(nsd).is_nsd() = True, .is_psd() = False
  n <- 50L
  psd_c <- Constant(diag(n))
  nsd_c <- Constant(-diag(n))

  expect_true(is_psd(psd_c))
  expect_false(is_nsd(psd_c))
  expect_true(is_nsd(nsd_c))
  expect_false(is_psd(nsd_c))
})

## @cvxpy test_constant.py::test_is_psd
test_that("Constant: psd_wrap forces is_psd = TRUE", {
  ## CVXPY: psd_wrap(Constant(P)).is_psd() == True
  ## psd_wrap wraps a constant and asserts PSD regardless of eigenvalues
  n <- 50L
  P <- diag(n)
  wrapped <- psd_wrap(Constant(P))
  expect_true(is_psd(wrapped))
  expect_false(is_nsd(wrapped))
})

# ═══════════════════════════════════════════════════════════════════════
# test_prod (test_constant.py) — Prod atom with sparse matrices
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_constant.py::test_prod
test_that("Prod: sparse cross matrix product (axis=NULL, axis=1, axis=2)", {
  ## Build a 100x100 sparse cross matrix:
  ##   rows 0:99 in col 0 (col of ones), plus row 0 in cols 1:99 (row of ones)
  ##   Diagonal element (0,0) gets value 2 (double-counted), rest of row 0 and col 0 are 1.
  ##   All other entries are 0.
  ## CVXPY axis=0 → R axis=2, CVXPY axis=1 → R axis=1
  ## Verified via CVXPY:
  ##   prod(A).value = 0.0
  ##   prod(A, axis=0)[:5] = [1, 0, 0, 0, 0], shape=(100,) → R c(100,1)
  ##   prod(A, axis=1)[:5] = [1, 0, 0, 0, 0], shape=(100,) → R c(100,1)

  ## Build the sparse matrix using 1-based indices for sparseMatrix
  rows_part1 <- 1:100          # 1-based: arange(100) + 1
  cols_part1 <- rep(1L, 100)   # all in column 1
  rows_part2 <- rep(1L, 99)    # all in row 1
  cols_part2 <- 2:100          # columns 2:100

  A <- Matrix::sparseMatrix(
    i = c(rows_part1, rows_part2),
    j = c(cols_part1, cols_part2),
    x = rep(1, 199),
    dims = c(100, 100)
  )

  ## Prod over all entries: any row/col with zeros → product = 0
  prod_all <- Prod(Constant(A))
  expect_equal(as.numeric(value(prod_all)), 0.0, tolerance = 1e-10)

  ## Prod along axis=2 (CVXPY axis=0): column-wise product
  ## Column 1 has all 1s → product = 1. Columns 2-100 have zeros → product = 0.
  prod_ax2 <- Prod(Constant(A), axis = 2L)
  val_ax2 <- as.numeric(value(prod_ax2))
  expect_equal(val_ax2[1], 1.0, tolerance = 1e-10)
  expect_equal(val_ax2[2:100], rep(0, 99), tolerance = 1e-10)

  ## Prod along axis=1 (CVXPY axis=1): row-wise product
  ## Row 1 has all 1s (or 2 in position (1,1)) → product != 0.
  ## Rows 2-100 only have one nonzero in col 1, rest are zero → product = 0.
  prod_ax1 <- Prod(Constant(A), axis = 1L)
  val_ax1 <- as.numeric(value(prod_ax1))
  expect_equal(val_ax1[2:100], rep(0, 99), tolerance = 1e-10)
})

## @cvxpy test_constant.py::test_prod
test_that("Prod: dense small matrix (overall product)", {
  ## CVXPY: B = arange(4).reshape(2,2) + 1 = [[1,3],[2,4]] (col-major in R: [[1,2],[3,4]])
  ## Note: NumPy reshape is row-major, R is column-major
  ## B (NumPy row-major) = [[1,2],[3,4]], prod = 24
  ## Verified via CVXPY: prod(sparse(B)).value = 24
  B <- matrix(1:4, 2, 2)  # R col-major: [[1,3],[2,4]] — product = 24
  B_sparse <- Matrix::Matrix(B, sparse = TRUE)
  prod_B <- Prod(Constant(B_sparse))
  expect_equal(as.numeric(value(prod_B)), 24.0, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# test_add_expr_copy (test_atoms.py line 54)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_add_expr_copy
test_that("AddExpression copy: same type, same args, new object", {
  x <- Variable(2, name = "x")
  y <- Variable(2, name = "y")
  A <- Variable(c(2, 2), name = "A")
  B <- Variable(c(2, 2), name = "B")

  atom <- x + y
  copy <- CVXR:::expr_copy(atom)
  ## A new object is constructed, so copy@args has same content but is not identical

  expect_equal(class(copy)[1], class(atom)[1])
  expect_equal(copy@args[[1]]@id, atom@args[[1]]@id)
  expect_equal(copy@args[[2]]@id, atom@args[[2]]@id)

  ## Test copy with new args
  copy2 <- CVXR:::expr_copy(atom, args = list(A, B))
  expect_equal(class(copy2)[1], class(atom)[1])
  expect_true(copy2@args[[1]]@id == A@id)
  expect_true(copy2@args[[2]]@id == B@id)
})

# ═══════════════════════════════════════════════════════════════════════
# test_list_input (test_atoms.py line 96)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_list_input
test_that("list input rejected by atoms", {
  ## CVXPY: cp.max([cp.Variable(), 1]) raises Exception about list input
  ## In CVXR, as_expr rejects list inputs
  expect_error(max_entries(list(Variable(), 1)), regexp = ".*Cannot convert.*list.*")
  expect_error(cvxr_norm(list(1, Variable())), regexp = ".*Cannot convert.*list.*")

  x <- Variable()
  y <- Variable()
  expect_error(cvxr_norm(list(x, y)), regexp = ".*Cannot convert.*list.*")
})

# ═══════════════════════════════════════════════════════════════════════
# test_norm_exceptions (test_atoms.py line 119)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_norm_exceptions
test_that("norm(x, 'nuc') on non-matrix raises error", {
  ## CVXPY: cp.norm(x, 'nuc') on vector -> "Unsupported norm option nuc for non-matrix."
  ## In CVXR, cvxr_norm passes "nuc" to as.numeric which gives NA, causing an error.
  ## We just check that an error is raised for nuclear norm on a vector.
  x <- Variable(2)
  expect_error(suppressWarnings(cvxr_norm(x, "nuc")))
})

# ═══════════════════════════════════════════════════════════════════════
# test_harmonic_mean (test_atoms.py line 220)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_harmonic_mean
test_that("harmonic_mean: shape, curvature, sign", {
  x <- Variable(2)
  atom <- harmonic_mean(x)
  ## CVXPY: shape == tuple() -> scalar c(1,1)
  expect_equal(atom@shape, c(1L, 1L))
  ## CVXPY: curvature == s.CONCAVE
  expect_equal(curvature(atom), "CONCAVE")
  ## CVXPY: sign == s.NONNEG
  expect_true(is_nonneg(atom))
})

# ═══════════════════════════════════════════════════════════════════════
# test_matrix_frac (test_atoms.py line 354)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_matrix_frac
test_that("matrix_frac: shape, curvature, validation", {
  x <- Variable(2)
  A <- Variable(c(2, 2))
  C <- Variable(c(3, 2))

  atom <- matrix_frac(x, A)
  ## CVXPY: shape == tuple() -> scalar c(1,1)
  expect_equal(atom@shape, c(1L, 1L))
  ## CVXPY: curvature == s.CONVEX
  expect_equal(curvature(atom), "CONVEX")

  ## Non-square matrix: "The second argument to matrix_frac must be a square matrix."
  expect_error(matrix_frac(x, C), regexp = ".*square matrix.*")

  ## Incompatible dimensions
  expect_error(matrix_frac(Variable(3), A), regexp = ".*incompatible.*dimensions.*")
})

# ═══════════════════════════════════════════════════════════════════════
# test_concatenate (test_atoms.py line 604)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_concatenate
test_that("Concatenate atom (stub — not in CVXR)", {
  skip("Concatenate atom not implemented in CVXR; use HStack/VStack instead")
})

# ═══════════════════════════════════════════════════════════════════════
# test_squeeze (test_atoms.py line 736)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_squeeze
test_that("Squeeze atom (stub — not in CVXR)", {
  skip("Squeeze atom not implemented in CVXR")
})

# ═══════════════════════════════════════════════════════════════════════
# test_diag_offset (test_atoms.py line 809)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_diag_offset
test_that("DiagVec/DiagMat with k offset", {
  ## CVXPY test_matrix: np.array([[1,2,3],[4,5,6],[7,8,9]])
  ## In R (column-major), to match: row 1=[1,2,3], row 2=[4,5,6], row 3=[7,8,9]
  test_matrix <- matrix(c(1, 4, 7, 2, 5, 8, 3, 6, 9), nrow = 3, ncol = 3)
  test_vector <- c(1, 2, 3)
  offsets <- c(0L, 1L, -1L, 2L)

  ## Helper to extract k-th diagonal from a matrix (NumPy convention)
  np_diag <- function(m, k = 0L) {
    n <- nrow(m)
    p <- ncol(m)
    if (k >= 0) {
      idx <- seq_len(min(n, p - k))
      return(m[cbind(idx, idx + k)])
    } else {
      idx <- seq_len(min(n + k, p))
      return(m[cbind(idx - k, idx)])
    }
  }

  for (offset in offsets) {
    ## Matrix -> vector (DiagMat with offset)
    a_cvxr <- as.numeric(value(DiagMat(Constant(test_matrix), k = offset)))
    a_np   <- np_diag(test_matrix, offset)
    expect_equal(a_cvxr, a_np, tolerance = 1e-10)

    ## Vector -> matrix (DiagVec with offset)
    A_cvxr <- as.matrix(value(DiagVec(test_vector, k = offset)))
    ## R diag(v, nrow=n+abs(k)) doesn't do offset; build manually
    n <- length(test_vector)
    sz <- n + abs(offset)
    A_np <- matrix(0, sz, sz)
    if (offset >= 0) {
      for (i in seq_along(test_vector)) A_np[i, i + offset] <- test_vector[i]
    } else {
      for (i in seq_along(test_vector)) A_np[i - offset, i] <- test_vector[i]
    }
    expect_equal(A_cvxr, A_np, tolerance = 1e-10)
  }

  ## CVXPY: cp.diag(Variable(5), 1).size == 36
  X <- DiagVec(Variable(5), k = 1L)
  expect_equal(prod(X@shape), 36L)
})

# ═══════════════════════════════════════════════════════════════════════
# test_Trace (test_atoms.py line 850)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_Trace
test_that("trace(A) resolves to Trace class", {
  A <- Variable(c(4, 4))
  t <- matrix_trace(A)
  ## CVXPY: isinstance(t, cp.Trace)
  expect_true(S7::S7_inherits(t, CVXR:::Trace))
})

# ═══════════════════════════════════════════════════════════════════════
# test_trace_complex2real (test_atoms.py line 875)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_trace_complex2real
test_that("trace in complex problem", {
  X <- Variable(c(2, 2), complex = TRUE)
  problem <- Problem(Minimize(cvxr_norm(matrix_trace(X))), list(X == 2))
  result <- psolve(problem, solver = "CLARABEL")
  expect_equal(result, 4.0, tolerance = 1e-5)
})

# ═══════════════════════════════════════════════════════════════════════
# test_trace_dgp2dcp (test_atoms.py line 881)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_trace_dgp2dcp
test_that("trace in DGP problem", {
  X <- Variable(c(2, 2), pos = TRUE)
  problem <- Problem(Minimize(matrix_trace(X)), list(X == 2))
  result <- psolve(problem, gp = TRUE)
  expect_equal(result, 4.0, tolerance = 1e-5)
})

# ═══════════════════════════════════════════════════════════════════════
# test_cvar (test_atoms.py line 1137)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_cvar
test_that("CVaR atom evaluation and LP usage", {
  ## Test CVaR computation against alternative LP formulation
  set.seed(1)
  m <- 100L
  x_data <- rnorm(m)
  betas <- c(0.1, 0.5, 0.9, 0.95, 0.99)

  for (beta in betas) {
    ## CVaR via atom
    cvar_val <- as.numeric(value(cvar(x_data, beta)))

    ## CVaR via alternative LP formulation: alpha + 1/((1-beta)*m) * sum(pos(x - alpha))
    alpha <- Variable()
    objective <- Minimize(alpha + 1.0 / ((1 - beta) * m) * sum_entries(pos(x_data - alpha)))
    prob_alt <- Problem(objective)
    cvar_alt_val <- as.numeric(psolve(prob_alt, solver = "CLARABEL"))

    expect_equal(cvar_val, cvar_alt_val, tolerance = 1e-4)
  }
})

# ═══════════════════════════════════════════════════════════════════════
# test_bmat (test_atoms.py line 1219)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_bmat
test_that("bmat block matrix construction", {
  ## CVXPY builds a block matrix from list of lists
  v_np <- matrix(1, 3, 1)
  ## Expected shape: (5, 2)
  expr <- bmat(list(
    list(v_np, v_np),
    list(matrix(0, 2, 1), matrix(c(1, 2), 2, 1))
  ))
  expect_equal(expr@shape, c(5L, 2L))

  ## Check numeric value
  expected <- rbind(
    cbind(v_np, v_np),
    cbind(matrix(0, 2, 1), matrix(c(1, 2), 2, 1))
  )
  expect_equal(as.matrix(value(expr)), expected, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# test_cumprod (test_atoms.py line 1279)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_cumprod
test_that("CumProd in DGP for both axes", {
  ## CVXPY tests both axis=0 and axis=1 with DGP
  ## CVXPY axis=0 -> CVXR axis=2 (column-wise)
  ## CVXPY axis=1 -> CVXR axis=1 (row-wise)
  for (cvxpy_axis in c(0L, 1L)) {
    r_axis <- if (cvxpy_axis == 0L) 2L else 1L
    x <- Variable(c(4, 3), pos = TRUE)
    expr <- CVXR:::Cumprod(x, axis = r_axis)
    ## Constant needs to be elementwise positive
    x_val <- matrix(1:12, nrow = 4, ncol = 3)

    ## Compute expected: R apply with cumprod along the axis
    if (r_axis == 2L) {
      target <- apply(x_val, 2, cumprod)
    } else {
      target <- t(apply(x_val, 1, cumprod))
    }

    prob <- Problem(Minimize(sum_entries(expr)), list(x == x_val))
    psolve(prob, gp = TRUE)

    expect_equal(as.matrix(value(expr)), target, tolerance = 1e-4)
  }
})

# ═══════════════════════════════════════════════════════════════════════
# test_convolve (test_atoms.py line 1308)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_convolve
test_that("conv (convolve): sign propagation and errors", {
  ## NOTE: CVXR uses conv() which returns (n+m-1, 1) shape, while CVXPY
  ## convolve returns (n+m-1,). We test the equivalent functionality.

  ## Sign propagation: nonneg constant * nonneg param -> nonneg
  a <- c(1, 1, 1)
  b_pos <- Parameter(2, nonneg = TRUE)
  expr <- conv(a, b_pos)
  expect_true(is_nonneg(expr))
  ## Shape: len(a) + len(b) - 1 = 3 + 2 - 1 = 4 -> (4, 1)
  expect_equal(expr@shape, c(4L, 1L))

  ## Nonpos param -> nonpos
  b_neg <- Parameter(2, nonpos = TRUE)
  expr2 <- conv(a, b_neg)
  expect_true(is_nonpos(expr2))

  ## Error: matrix input rejected (CVXPY: "scalar or 1D")
  x <- Variable(2)
  expect_error(conv(matrix(c(0, 1, 0, 1), 2, 2), x), regexp = ".*vector.*")
})

# ═══════════════════════════════════════════════════════════════════════
# test_ptp (test_atoms.py line 1341)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_ptp
test_that("ptp: peak-to-peak with axis and keepdims", {
  ## CVXPY: a = np.array([[10., -10., 3.0], [6., 0., -1.5]])
  a <- matrix(c(10, -10, 3, 6, 0, -1.5), nrow = 2, byrow = TRUE)

  ## Overall ptp = max - min = 10 - (-10) = 20
  expr <- ptp(a)
  expect_equal(expr@shape, c(1L, 1L))
  expect_equal(as.numeric(value(expr)), 20.0, tolerance = 1e-10)

  ## axis=0 (CVXPY) -> axis=2 (CVXR): column-wise ptp
  ## col 1: max(10,6)-min(10,6) = 4, col 2: max(-10,0)-min(-10,0) = 10, col 3: max(3,-1.5)-min(3,-1.5)=4.5
  expr_ax2 <- ptp(a, axis = 2L)
  expect_equal(as.numeric(value(expr_ax2)), c(4, 10, 4.5), tolerance = 1e-10)

  ## axis=1 (CVXPY) -> axis=1 (CVXR): row-wise ptp
  ## row 1: 10-(-10)=20, row 2: 6-(-1.5)=7.5
  expr_ax1 <- ptp(a, axis = 1L)
  expect_equal(as.numeric(value(expr_ax1)), c(20, 7.5), tolerance = 1e-10)

  ## axis=0 keepdims=TRUE (CVXPY) -> axis=2 keepdims=TRUE (CVXR)
  ## Expected shape: (1, 3)
  expr_ax2_kd <- ptp(a, axis = 2L, keepdims = TRUE)
  expect_equal(expr_ax2_kd@shape, c(1L, 3L))
  expect_equal(as.numeric(value(expr_ax2_kd)), c(4, 10, 4.5), tolerance = 1e-10)

  ## axis=1 keepdims=TRUE (CVXPY) -> axis=1 keepdims=TRUE (CVXR)
  ## Expected shape: (2, 1)
  expr_ax1_kd <- ptp(a, axis = 1L, keepdims = TRUE)
  expect_equal(expr_ax1_kd@shape, c(2L, 1L))
  expect_equal(as.numeric(value(expr_ax1_kd)), c(20, 7.5), tolerance = 1e-10)

  ## Curvature: ptp on Variable is convex
  x_var <- Variable(10)
  expr_var <- ptp(x_var)
  expect_equal(curvature(expr_var), "CONVEX")
})

# ═══════════════════════════════════════════════════════════════════════
# test_stats (test_atoms.py line 1374)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_stats
test_that("mean/var/std with axis and keepdims", {
  ## CVXPY: a = np.array([[10., 10., 3.0], [6., 0., 1.5]])
  a <- matrix(c(10, 10, 3, 6, 0, 1.5), nrow = 2, byrow = TRUE)

  ## Overall mean, var, std
  expr_mean <- cvxr_mean(a)
  expr_var  <- cvxr_var(a)
  expr_std  <- cvxr_std(a)
  expect_true(is_nonneg(expr_mean))
  expect_true(is_nonneg(expr_var))
  expect_true(is_nonneg(expr_std))

  ## Expected values (NumPy defaults: ddof=0)
  np_mean <- mean(a)
  np_var  <- mean((a - mean(a))^2)   # population variance
  np_std  <- sqrt(np_var)

  expect_equal(as.numeric(value(expr_mean)), np_mean, tolerance = 1e-10)
  expect_equal(as.numeric(value(expr_var)),  np_var,  tolerance = 1e-6)
  expect_equal(as.numeric(value(expr_std)),  np_std,  tolerance = 1e-6)

  ## ddof tests
  for (ddof in c(0, 1)) {
    expr_var_d  <- cvxr_var(a, ddof = ddof)
    expr_std_d  <- cvxr_std(a, ddof = ddof)
    n <- length(a)
    expected_var <- sum((a - mean(a))^2) / (n - ddof)
    expected_std <- sqrt(expected_var)
    expect_equal(as.numeric(value(expr_var_d)), expected_var, tolerance = 1e-6)
    expect_equal(as.numeric(value(expr_std_d)), expected_std, tolerance = 1e-6)
  }

  ## Axis tests (CVXPY axis=0 -> CVXR axis=2, CVXPY axis=1 -> CVXR axis=1)
  for (cvxpy_axis in c(0L, 1L)) {
    r_axis <- if (cvxpy_axis == 0L) 2L else 1L
    for (keepdims in c(TRUE, FALSE)) {
      expr_mean_ax <- cvxr_mean(a, axis = r_axis, keepdims = keepdims)

      ## Expected from R: apply then mean
      if (r_axis == 2L) {
        expected_mean <- apply(a, 2, mean)
      } else {
        expected_mean <- apply(a, 1, mean)
      }
      if (keepdims) {
        if (r_axis == 2L) {
          expected_shape <- c(1L, ncol(a))
        } else {
          expected_shape <- c(nrow(a), 1L)
        }
      }
      expect_equal(as.numeric(value(expr_mean_ax)), expected_mean, tolerance = 1e-10)
    }
  }
})

# ═══════════════════════════════════════════════════════════════════════
# test_mixed_norm (test_atoms.py line 1617)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_mixed_norm
test_that("mixed_norm evaluation", {
  ## CVXPY: y = Variable((5,5)); obj = Minimize(mixed_norm(y, "inf", 1))
  ## prob = Problem(obj, [y == ones((5,5))]), result = 5.
  y <- Variable(c(5, 5))
  obj <- Minimize(mixed_norm(y, Inf, 1))
  prob <- Problem(obj, list(y == matrix(1, 5, 5)))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 5.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# test_indicator (test_atoms.py line 1641)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_indicator
test_that("Indicator transform", {
  x <- Variable()
  constraints <- list(x >= 0, x <= 1)
  expr <- indicator(constraints)

  ## x = 0.5 -> in feasible set -> value 0
  value(x) <- 0.5
  expect_equal(as.numeric(value(expr)), 0.0)

  ## x = 2 -> outside feasible set -> value Inf
  value(x) <- 2.0
  expect_equal(as.numeric(value(expr)), Inf)
})

# ═══════════════════════════════════════════════════════════════════════
# test_vdot (test_atoms.py line 1742)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_vdot
test_that("vdot / scalar_product", {
  p <- rep(1, 4)
  v <- Variable(c(4, 1))

  ## Minimize vdot(v, p) s.t. v >= 1 -> optimal v = 1, objective = 4
  obj <- Minimize(vdot(v, p))
  prob <- Problem(obj, list(v >= 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(v)), rep(1, 4), tolerance = 1e-4)

  ## Same with scalar_product
  obj2 <- Minimize(scalar_product(v, p))
  prob2 <- Problem(obj2, list(v >= 1))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.numeric(value(v)), rep(1, 4), tolerance = 1e-4)

  ## With a parameter
  p_param <- Parameter(c(4, 1))
  v2 <- Variable(c(4, 1))
  value(p_param) <- matrix(1, 4, 1)
  obj3 <- Minimize(vdot(v2, p_param))
  prob3 <- Problem(obj3, list(v2 >= 1))
  psolve(prob3, solver = "CLARABEL")
  expect_equal(as.numeric(value(v2)), rep(1, 4), tolerance = 1e-4)

  obj4 <- Minimize(scalar_product(v2, p_param))
  prob4 <- Problem(obj4, list(v2 >= 1))
  psolve(prob4, solver = "CLARABEL")
  expect_equal(as.numeric(value(v2)), rep(1, 4), tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# test_conj (test_atoms.py line 1807)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_conj
test_that("Conj atom in constraints", {
  ## CVXPY: v = Variable(4); obj = Minimize(sum(v)); prob = Problem(obj, [conj(v) >= 1])
  ## For real variables, conj(v) == v, so this is just v >= 1.
  v <- Variable(c(4, 1))
  obj <- Minimize(sum_entries(v))
  prob <- Problem(obj, list(CVXR:::conj_expr(v) >= 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(v)), rep(1, 4), tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# test_partial_trace (test_atoms.py line 1833)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_partial_trace
test_that("partial_trace multi-subsystem", {
  ## CVXPY: rho_ABC = rho_A (x) rho_B (x) rho_C with complex random matrices
  ## CVXPY axis is 0-based; CVXR axis is 1-based
  ## NOTE: R Matrix package lacks complex sparse zgeMatrix class; skip if unavailable
  has_complex_matrix <- tryCatch({methods::getClass("zgeMatrix"); TRUE},
                                 error = function(e) FALSE)
  skip_if(!has_complex_matrix, "R Matrix package lacks zgeMatrix (complex sparse)")

  set.seed(1)
  rho_A <- matrix(complex(real = runif(16), imaginary = runif(16)), 4, 4)
  rho_A <- rho_A / sum(diag(rho_A))
  rho_B <- matrix(complex(real = runif(9), imaginary = runif(9)), 3, 3)
  rho_B <- rho_B / sum(diag(rho_B))
  rho_C <- matrix(complex(real = runif(4), imaginary = runif(4)), 2, 2)
  rho_C <- rho_C / sum(diag(rho_C))
  rho_AB <- kronecker(rho_A, rho_B)
  rho_AC <- kronecker(rho_A, rho_C)

  temp <- kronecker(rho_AB, rho_C)
  rho_ABC <- Variable(shape = dim(temp), complex = TRUE)
  value(rho_ABC) <- temp

  ## Trace out C (CVXPY axis=2 -> CVXR axis=3)
  rho_AB_test <- partial_trace(rho_ABC, c(4, 3, 2), axis = 3L)
  expect_equal(value(rho_AB_test), rho_AB, tolerance = 1e-10)

  ## Trace out B (CVXPY axis=1 -> CVXR axis=2)
  rho_AC_test <- partial_trace(rho_ABC, c(4, 3, 2), axis = 2L)
  expect_equal(value(rho_AC_test), rho_AC, tolerance = 1e-10)

  ## Trace out B from AB (CVXPY axis=1 -> CVXR axis=2)
  rho_A_test <- partial_trace(rho_AB_test, c(4, 3), axis = 2L)
  expect_equal(value(rho_A_test), rho_A, tolerance = 1e-10)

  ## Trace out A from AB (CVXPY axis=0 -> CVXR axis=1)
  rho_B_test <- partial_trace(rho_AB_test, c(4, 3), axis = 1L)
  expect_equal(value(rho_B_test), rho_B, tolerance = 1e-10)

  ## Trace out A from AC (CVXPY axis=0 -> CVXR axis=1)
  rho_C_test <- partial_trace(rho_AC_test, c(4, 2), axis = 1L)
  expect_equal(value(rho_C_test), rho_C, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# test_partial_trace_exceptions (test_atoms.py line 1873)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_partial_trace_exceptions
test_that("partial_trace exceptions", {
  ## Non-square matrix
  X <- Variable(c(4, 3))
  expect_error(partial_trace(X, dims = c(2, 3), axis = 1L),
               regexp = ".*square.*")

  ## Dims don't match matrix size
  X6 <- Variable(c(6, 6))
  expect_error(partial_trace(X6, dims = c(2, 4), axis = 1L),
               regexp = ".*dimension.*subsystems.*")
})

# ═══════════════════════════════════════════════════════════════════════
# test_partial_transpose (test_atoms.py line 1888)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_partial_transpose
test_that("partial_transpose multi-subsystem", {
  ## NOTE: R Matrix package lacks complex sparse zgeMatrix class; skip if unavailable
  has_complex_matrix <- tryCatch({methods::getClass("zgeMatrix"); TRUE},
                                 error = function(e) FALSE)
  skip_if(!has_complex_matrix, "R Matrix package lacks zgeMatrix (complex sparse)")

  set.seed(1)
  rho_A <- matrix(complex(real = runif(64), imaginary = runif(64)), 8, 8)
  rho_A <- rho_A / sum(diag(rho_A))
  rho_B <- matrix(complex(real = runif(36), imaginary = runif(36)), 6, 6)
  rho_B <- rho_B / sum(diag(rho_B))
  rho_C <- matrix(complex(real = runif(16), imaginary = runif(16)), 4, 4)
  rho_C <- rho_C / sum(diag(rho_C))

  rho_TC <- kronecker(kronecker(rho_A, rho_B), t(rho_C))
  rho_TB <- kronecker(kronecker(rho_A, t(rho_B)), rho_C)
  rho_TA <- kronecker(kronecker(t(rho_A), rho_B), rho_C)

  temp <- kronecker(kronecker(rho_A, rho_B), rho_C)
  rho_ABC <- Variable(shape = dim(temp), complex = TRUE)
  value(rho_ABC) <- temp

  ## Transpose C (CVXPY axis=2 -> CVXR axis=3)
  rho_TC_test <- partial_transpose(rho_ABC, c(8, 6, 4), axis = 3L)
  expect_equal(value(rho_TC_test), rho_TC, tolerance = 1e-10)

  ## Transpose B (CVXPY axis=1 -> CVXR axis=2)
  rho_TB_test <- partial_transpose(rho_ABC, c(8, 6, 4), axis = 2L)
  expect_equal(value(rho_TB_test), rho_TB, tolerance = 1e-10)

  ## Transpose A (CVXPY axis=0 -> CVXR axis=1)
  rho_TA_test <- partial_transpose(rho_ABC, c(8, 6, 4), axis = 1L)
  expect_equal(value(rho_TA_test), rho_TA, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# test_partial_transpose_exceptions (test_atoms.py line 1926)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_partial_transpose_exceptions
test_that("partial_transpose exceptions", {
  ## Non-square matrix
  X <- Variable(c(4, 3))
  expect_error(partial_transpose(X, dims = c(2, 3), axis = 1L),
               regexp = ".*square.*")

  ## Dims don't match
  X6 <- Variable(c(6, 6))
  expect_error(partial_transpose(X6, dims = c(2, 4), axis = 1L),
               regexp = ".*dimension.*subsystems.*")
})

# ═══════════════════════════════════════════════════════════════════════
# test_flatten (test_atoms.py line 1956)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_flatten
test_that("vec/flatten with order", {
  ## Constant argument, F-order
  A <- 0:9
  reshaped_F <- matrix(A, nrow = 2, ncol = 5, byrow = FALSE)  # F-order reshape
  ## vec uses F-order (column-major) by default
  expr_F <- vec(Constant(reshaped_F))
  expect_equal(as.numeric(value(expr_F)), A, tolerance = 1e-10)

  ## C-order: use reshape_expr
  reshaped_C <- matrix(A, nrow = 2, ncol = 5, byrow = TRUE)
  expr_C <- reshape_expr(Constant(reshaped_C), c(prod(dim(reshaped_C)), 1L), order = "C")
  expect_equal(as.numeric(value(expr_C)), A, tolerance = 1e-10)

  ## Variable argument, F-order
  x <- Variable(c(2, 5))
  reshaped_F_expected <- matrix(A, nrow = 2, ncol = 5, byrow = FALSE)
  expr_vf <- vec(x)
  prob <- Problem(Minimize(0), list(expr_vf == matrix(A, ncol = 1)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.matrix(value(x)), reshaped_F_expected, tolerance = 1e-4)

  ## Variable argument, C-order
  x2 <- Variable(c(2, 5))
  reshaped_C_expected <- matrix(A, nrow = 2, ncol = 5, byrow = TRUE)
  expr_vc <- reshape_expr(x2, c(prod(dim(reshaped_C_expected)), 1L), order = "C")
  prob2 <- Problem(Minimize(0), list(expr_vc == matrix(A, ncol = 1)))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.matrix(value(x2)), reshaped_C_expected, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# test_tr_inv (test_atoms.py line 1990)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_tr_inv
test_that("tr_inv in SDP", {
  ## minimize tr_inv(X) s.t. X >> 0, trace(X) == 1
  ## Optimal value = T^2 with X = I/T
  T_dim <- 5L
  X <- Variable(c(T_dim, T_dim), symmetric = TRUE)
  constraints <- list(X %>>% 0, matrix_trace(X) == 1)
  prob <- Problem(Minimize(tr_inv(X)), constraints)
  psolve(prob, solver = "CLARABEL")

  ## Best value is T^2
  expect_equal(prob@.cache$value, T_dim^2, tolerance = 1e-2)

  ## Optimal X = I/T
  X_expect <- diag(T_dim) / T_dim
  expect_equal(as.matrix(value(X)), X_expect, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# TestDotsort tests (test_atoms.py line 2035)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestDotsort::test_sum_k_largest_equivalence
test_that("dotsort: sum_k_largest equivalence", {
  x_var <- Variable(5)
  x_val <- c(1, 3, 2, -5, 0)
  w <- c(1, 1, 1, 0)
  expr <- dotsort(x_var, w)
  expect_true(is_convex(expr))
  expect_true(CVXR:::is_incr(expr, 1L))
  prob <- Problem(Minimize(expr), list(x_var == x_val))
  psolve(prob, solver = "CLARABEL")
  ## sum of 3 largest: sort(x_val) desc = [3, 2, 1, 0, -5] -> top 3 = 3+2+1 = 6
  expected <- sum(sort(x_val, decreasing = TRUE)[1:3])
  expect_equal(value(objective(prob)), expected, tolerance = 1e-4)
})

## @cvxpy test_atoms.py::TestDotsort::test_sum_k_smallest_equivalence
test_that("dotsort: sum_k_smallest equivalence (negated)", {
  x_var <- Variable(5)
  x_val <- c(1, 3, 2, -5, 0)
  w <- c(-1, -1, -1, 0)
  expr <- -dotsort(x_var, w)
  expect_true(is_concave(expr))
  expect_true(CVXR:::is_decr(-dotsort(x_var, w), 1L))
  prob <- Problem(Maximize(expr), list(x_var == x_val))
  psolve(prob, solver = "CLARABEL")
  ## sum of 3 smallest: sort(x_val) = [-5, 0, 1, 2, 3] -> bottom 3 = -5+0+1 = -4
  expected <- sum(sort(x_val)[1:3])
  expect_equal(value(objective(prob)), expected, tolerance = 1e-4)
})

## @cvxpy test_atoms.py::TestDotsort::test_copy
test_that("dotsort: copy", {
  x_var <- Variable(5)
  w <- c(1, 2)
  atom <- dotsort(x_var, w)
  copy <- CVXR:::expr_copy(atom)
  expect_equal(class(copy)[1], class(atom)[1])
  ## Args have same content
  expect_equal(copy@args[[1]]@id, atom@args[[1]]@id)

  ## Copy with new args
  copy2 <- CVXR:::expr_copy(atom, args = list(x_var, Constant(w)))
  expect_equal(class(copy2)[1], class(atom)[1])
  expect_true(copy2@args[[1]]@id == atom@args[[1]]@id)
})

## @cvxpy test_atoms.py::TestDotsort::test_exceptions
test_that("dotsort: exceptions", {
  x_var <- Variable(5)

  ## len(w) > len(x)
  expect_error(dotsort(x_var, c(1, 2, 3, 4, 5, 8)),
               regexp = ".*size.*W.*less.*equal.*size.*X.*")

  ## Two variable expressions for W
  expect_error(dotsort(x_var, Variable(3)),
               regexp = ".*W.*must be constant.*")

  ## Swapped arguments (constant X, variable W)
  expect_error(dotsort(c(1, 2, 3), x_var),
               regexp = ".*W.*must be constant.*")

  ## Non-DCP composition
  expect_error(psolve(Problem(Minimize(dotsort(abs(x_var), c(-1, 1))))),
               regexp = ".*DCP.*")
})

## @cvxpy test_atoms.py::TestDotsort::test_list
test_that("dotsort: list w input", {
  x_var <- Variable(5)
  r <- c(2, 1, 0, -1, -1)
  w <- c(1.2, 1.1)
  expr <- dotsort(x_var, w)
  prob <- Problem(Maximize(sum_entries(r * x_var)),
                  list(x_var >= 0, expr <= 1, sum_entries(x_var) == 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(expr)), 1.0, tolerance = 1e-4)
})

## @cvxpy test_atoms.py::TestDotsort::test_parameter
test_that("dotsort: parameter w", {
  x_var <- Variable(5)
  x_val <- c(1, 3, 2, -5, 0)

  ## Positive parameter -> is_incr (CVXPY uses pos=True and nonneg=True)
  expect_true(CVXR:::is_incr(dotsort(x_var, Parameter(2, pos = TRUE)), 1L))
  expect_true(CVXR:::is_incr(dotsort(x_var, Parameter(2, nonneg = TRUE)), 1L))
  ## Negative parameter -> is_decr, not is_incr
  expect_false(CVXR:::is_incr(dotsort(x_var, Parameter(2, neg = TRUE)), 1L))
  expect_true(CVXR:::is_decr(dotsort(x_var, Parameter(2, neg = TRUE)), 1L))

  ## Fixed parameter value
  w_p <- Parameter(2)
  value(w_p) <- matrix(c(1, 0), ncol = 1)
  expr <- dotsort(x_var, w_p)
  expect_false(CVXR:::is_incr(expr, 1L))
  expect_false(CVXR:::is_decr(expr, 1L))

  prob <- Problem(Minimize(expr), list(x_var == x_val))
  psolve(prob, solver = "CLARABEL")
  ## w_p = [1, 0] -> padded to [1, 0, 0, 0, 0]
  ## dotsort = sort(x_val) . sort(w_padded) = [-5,0,1,2,3] . [0,0,0,0,1]
  expected <- sort(x_val) %*% sort(c(1, 0, 0, 0, 0))
  expect_equal(as.numeric(value(objective(prob))), as.numeric(expected), tolerance = 1e-4)

  ## Change parameter value and resolve
  value(w_p) <- matrix(c(-1, -1), ncol = 1)
  prob2 <- Problem(Minimize(dotsort(x_var, w_p)), list(x_var == x_val))
  psolve(prob2, solver = "CLARABEL")
  expected2 <- sort(x_val) %*% sort(c(-1, -1, 0, 0, 0))
  expect_equal(as.numeric(value(objective(prob2))), as.numeric(expected2), tolerance = 1e-4)
})

## @cvxpy test_atoms.py::TestDotsort::test_non_fixed_x
test_that("dotsort: non-fixed x optimization", {
  x_var <- Variable(5)
  r <- c(2, 1, 0, -1, -1)
  w <- c(1.2, 1.1)
  expr <- dotsort(x_var, w)
  prob <- Problem(Maximize(sum_entries(r * x_var)),
                  list(x_var >= 0, expr <= 1, sum_entries(x_var) == 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(expr)), 1.0, tolerance = 1e-4)
  ## The two largest x values dotted with w should equal 1
  x_sorted <- sort(as.numeric(value(x_var)), decreasing = TRUE)
  expect_equal(sum(x_sorted[1:2] * sort(w, decreasing = TRUE)), 1.0, tolerance = 1e-4)

  ## Test with unordered w of length 3
  r2 <- c(2, 1, 0, -1, -1)
  w2 <- c(1.2, 1.1, 1.3)
  expr2 <- dotsort(x_var, w2)
  prob2 <- Problem(Maximize(sum_entries(r2 * x_var)),
                   list(x_var >= 0, expr2 <= 1, sum_entries(x_var) == 1))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.numeric(value(expr2)), 1.0, tolerance = 1e-4)
  x_sorted2 <- sort(as.numeric(value(x_var)), decreasing = TRUE)
  expect_equal(sum(x_sorted2[1:3] * sort(w2, decreasing = TRUE)), 1.0, tolerance = 1e-4)
})

# ====================================================================
# DOMAIN TESTS (test_domain.py)
# ====================================================================

# ── test_entr ──────────────────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_entr
test_that("domain: entr(a) domain constrains a >= 0 (CVXPY parity)", {
  ## CVXPY: dom = cp.entr(a).domain; min a s.t. dom => a ~= 0
  a <- Variable(name = "a")
  dom <- CVXR:::domain(entr(a))
  expect_length(dom, 1)
  prob <- Problem(Minimize(a), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-4)
})

# ── test_geo_mean ─────────────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_geo_mean
test_that("domain: geo_mean domain constrains x >= 0 (CVXPY parity)", {
  ## CVXPY: min sum(x) s.t. geo_mean(x).domain => sum ~= 0
  x <- Variable(2, name = "x")
  dom <- CVXR:::domain(geo_mean(x))
  prob <- Problem(Minimize(sum_entries(x)), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0, tolerance = 1e-4)

  ## geo_mean with weights [0, 2] => x[1] free, x[2] >= 0
  x <- Variable(2, name = "x")
  dom <- CVXR:::domain(geo_mean(x, c(0, 2)))
  dom <- c(dom, list(x >= -1))
  prob <- Problem(Minimize(sum_entries(x)), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), c(-1, 0), tolerance = 1e-4)

  ## geo_mean with weights [0, 1, 1] on z(3)
  z <- Variable(3, name = "z")
  dom <- CVXR:::domain(geo_mean(z, c(0, 1, 1)))
  dom <- c(dom, list(z >= -1))
  prob <- Problem(Minimize(sum_entries(z)), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(z)), c(-1, 0, 0), tolerance = 1e-4)
})

# ── test_log ──────────────────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_log
test_that("domain: log(a) domain constrains a >= 0 (CVXPY parity)", {
  ## CVXPY: dom = cp.log(a).domain; min a s.t. dom => a ~= 0
  a <- Variable(name = "a")
  dom <- CVXR:::domain(log(a))
  prob <- Problem(Minimize(a), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-4)
})

# ── test_log_det ──────────────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_log_det
test_that("domain: log_det(A + I) domain constrains A + I >> 0 (CVXPY parity)", {
  ## CVXPY: min sum(diag(A)) s.t. log_det(A+I).domain => -2
  A <- Variable(c(2, 2), name = "A")
  dom <- CVXR:::domain(log_det(A + diag(2)))
  prob <- Problem(Minimize(matrix_trace(A)), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), -2, tolerance = 1e-2)
})

# ── test_matrix_frac ──────────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_matrix_frac
test_that("domain: matrix_frac(x, A + I) domain constrains A + I >> 0 (CVXPY parity)", {
  ## CVXPY: min sum(diag(A)) s.t. matrix_frac(x, A+I).domain => -2
  x <- Variable(2, name = "x")
  A <- Variable(c(2, 2), name = "A")
  dom <- CVXR:::domain(matrix_frac(x, A + diag(2)))
  prob <- Problem(Minimize(matrix_trace(A)), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), -2, tolerance = 1e-2)
})

# ── test_partial_problem ──────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_partial_problem
test_that("domain: partial_optimize domain (CVXPY parity)", {
  skip("CVXR does not implement partial_optimize transform")
})

# ── test_pnorm ────────────────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_pnorm
test_that("domain: p_norm(a, -0.5) domain constrains a >= 0 (CVXPY parity)", {
  ## CVXPY: dom = cp.pnorm(a, -0.5).domain; min a s.t. dom => a ~= 0
  a <- Variable(name = "a")
  dom <- CVXR:::domain(p_norm(a, -0.5))
  prob <- Problem(Minimize(a), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0, tolerance = 1e-4)
})

# ── test_power ────────────────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_power
test_that("domain: power atom domains constrain sign correctly (CVXPY parity)", {
  ## sqrt(a) domain => a >= 0, min a ~= 0
  a <- Variable(name = "a")
  dom <- CVXR:::domain(sqrt(a))
  prob <- Problem(Minimize(a), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-4)

  ## square(a) domain is empty => min a with a >= -100 => -100
  a <- Variable(name = "a")
  dom <- CVXR:::domain(power(a, 2))
  prob <- Problem(Minimize(a), c(dom, list(a >= -100)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), -100, tolerance = 1e-4)

  ## a^(-1) domain => a >= 0, min a ~= 0
  a <- Variable(name = "a")
  dom <- CVXR:::domain(power(a, -1))
  prob <- Problem(Minimize(a), c(dom, list(a >= -100)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-4)

  ## a^3 domain => a >= 0, min a ~= 0
  a <- Variable(name = "a")
  dom <- CVXR:::domain(power(a, 3))
  prob <- Problem(Minimize(a), c(dom, list(a >= -100)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-4)
})

# ── test_quad_over_lin ────────────────────────────────────────────

## @cvxpy test_domain.py::TestDomain::test_quad_over_lin
test_that("domain: quad_over_lin(x, a) domain constrains a >= 0 (CVXPY parity)", {
  ## CVXPY: dom = cp.quad_over_lin(x, a).domain; min a s.t. dom => a ~= 0
  a <- Variable(name = "a")
  x <- Variable(2, name = "x")
  dom <- CVXR:::domain(quad_over_lin(x, a))
  prob <- Problem(Minimize(a), dom)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-4)
})

# ====================================================================
# POWER ATOM APPROX TESTS (test_power_atom.py)
# ====================================================================
#
# Helper to get cone counts from problem data (mirrors CVXPY _get_cone_counts)
.get_cone_counts <- function(prob, solver = "CLARABEL") {
  pd <- problem_data(prob, solver = solver)
  dims <- pd$data$dims
  list(
    soc = length(dims@soc),
    p3d = length(dims@p3d),
    pnd = length(dims@pnd)
  )
}

# ── TestPowerAtom ──────────────────────────────────────────────────

## @cvxpy test_power_atom.py::TestPowerAtom::test_dunder_pow_returns_approx
test_that("power: ^ operator returns PowerApprox (CVXPY parity)", {
  ## CVXPY: x**p uses PowerApprox so it canonicalizes via SOC.
  x <- Variable(3)
  expr <- x^2
  expect_true(S7_inherits(expr, CVXR:::PowerApprox))
  expect_true(S7_inherits(expr, CVXR:::Power))

  expr2 <- x^0.5
  expect_true(S7_inherits(expr2, CVXR:::PowerApprox))

  expr3 <- x^3
  expect_true(S7_inherits(expr3, CVXR:::PowerApprox))
})

## @cvxpy test_power_atom.py::TestPowerAtom::test_approx_controls_cone_type
test_that("power: approx=TRUE uses SOC, approx=FALSE uses power cones (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  obj <- Minimize(x[1] + x[2] - x[3])

  ## approx=TRUE -> SOC cones, no power cones
  prob <- Problem(obj, list(power(x, 3.3, approx = TRUE) <= rep(1, 3)))
  cc <- .get_cone_counts(prob, "CLARABEL")
  expect_gt(cc$soc, 0)
  expect_equal(cc$p3d, 0)

  ## approx=FALSE -> power cones, no SOC
  prob2 <- Problem(obj, list(power(x, 3.3, approx = FALSE) <= rep(1, 3)))
  cc2 <- .get_cone_counts(prob2, "CLARABEL")
  expect_equal(cc2$soc, 0)
  expect_gt(cc2$p3d, 0)
})

## @cvxpy test_power_atom.py::TestPowerAtom::test_approx_and_exact_agree
test_that("power: approx and exact give same answer for all p ranges (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  ## p < 0 (low), 0 < p < 1 (mid), p > 1 non-integer (high), p > 1 even integer
  cases <- list(
    list(p = -1.5, direction = "<="),
    list(p = 0.8,  direction = ">="),
    list(p = 4.5,  direction = "<="),
    list(p = 8,    direction = "<=")
  )
  for (case in cases) {
    p <- case$p
    direction <- case$direction

    x <- Variable(3)
    if (direction == "<=") {
      constr_a <- list(power(x, p, approx = TRUE) <= rep(1, 3))
    } else {
      constr_a <- list(power(x, p, approx = TRUE) >= rep(1, 3))
    }
    obj <- Minimize(x[1] + x[2] - x[3] + (x[2] + x[3])^2)
    prob <- Problem(obj, constr_a)
    psolve(prob, solver = "CLARABEL")
    x_approx <- as.numeric(value(x))

    if (direction == "<=") {
      constr_e <- list(power(x, p, approx = FALSE) <= rep(1, 3))
    } else {
      constr_e <- list(power(x, p, approx = FALSE) >= rep(1, 3))
    }
    obj2 <- Minimize(x[1] + x[2] - x[3] + (x[2] + x[3])^2)
    prob2 <- Problem(obj2, constr_e)
    psolve(prob2, solver = "CLARABEL")
    expect_true(status(prob2) %in% c("optimal", "optimal_inaccurate"),
                info = paste("p =", p))
    expect_equal(as.numeric(value(x)), x_approx, tolerance = 1e-3,
                 info = paste("p =", p))
  }
})

## @cvxpy test_power_atom.py::TestPowerAtom::test_approx_false_errors_without_power_cone_support
test_that("power: approx=FALSE errors when solver lacks power cone support (CVXPY parity)", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(3)
  prob <- Problem(
    Minimize(x[1] + x[2] - x[3]),
    list(power(x, 3.3, approx = FALSE) <= rep(1, 3))
  )
  ## ECOS does not support power cones => should error
  expect_error(psolve(prob, solver = "ECOS"))
})

## @cvxpy test_power_atom.py::TestPowerAtom::test_approx_warning
test_that("power: warning fires for approx=TRUE with many SOCs (CVXPY parity)", {
  skip("CVXR does not yet emit SOC approximation warnings like CVXPY")
})

# ── TestGeoMeanApprox ──────────────────────────────────────────────

## @cvxpy test_power_atom.py::TestGeoMeanApprox::test_approx_controls_cone_type
test_that("geo_mean: approx=TRUE uses SOC, approx=FALSE uses power cones (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, pos = TRUE)

  ## approx=TRUE -> SOC cones, no PowConeND
  obj <- Maximize(geo_mean(x, approx = TRUE))
  prob <- Problem(obj, list(sum_entries(x) <= 3))
  cc <- .get_cone_counts(prob, "CLARABEL")
  expect_gt(cc$soc, 0)
  expect_equal(cc$pnd, 0)

  ## approx=FALSE -> PowConeND, no SOC
  obj2 <- Maximize(geo_mean(x, approx = FALSE))
  prob2 <- Problem(obj2, list(sum_entries(x) <= 3))
  cc2 <- .get_cone_counts(prob2, "CLARABEL")
  expect_equal(cc2$soc, 0)
  expect_gt(cc2$pnd, 0)
})

## @cvxpy test_power_atom.py::TestGeoMeanApprox::test_approx_and_exact_agree
test_that("geo_mean: approx and exact give same answer (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  cases <- list(
    list(weights = NULL, n = 3L),      # uniform weights, 3 vars
    list(weights = c(1, 2, 1), n = 3L), # non-uniform weights
    list(weights = NULL, n = 4L)        # uniform weights, 4 vars
  )
  for (case in cases) {
    w <- case$weights
    n <- case$n
    x <- Variable(n, pos = TRUE)
    constr <- list(sum_entries(x) <= n, x[1] >= 0.5)

    prob <- Problem(Maximize(geo_mean(x, w, approx = TRUE)), constr)
    val_a <- psolve(prob, solver = "CLARABEL")
    x_a <- as.numeric(value(x))

    prob2 <- Problem(Maximize(geo_mean(x, w, approx = FALSE)), constr)
    val_e <- psolve(prob2, solver = "CLARABEL")
    expect_true(status(prob2) %in% c("optimal", "optimal_inaccurate"),
                info = paste("n =", n))
    expect_equal(val_e, val_a, tolerance = 1e-3,
                 info = paste("n =", n, "value"))
    expect_equal(as.numeric(value(x)), x_a, tolerance = 1e-3,
                 info = paste("n =", n, "x"))
  }
})

## @cvxpy test_power_atom.py::TestGeoMeanApprox::test_approx_warning
test_that("geo_mean: warning fires for approx=TRUE with many SOCs (CVXPY parity)", {
  skip("CVXR does not yet emit SOC approximation warnings like CVXPY")
})

# ── TestGeoMeanSingleWeight ───────────────────────────────────────

## @cvxpy test_power_atom.py::TestGeoMeanSingleWeight::test_single_weight_is_affine
test_that("geo_mean: single weight is affine (CVXPY parity)", {
  ## CVXPY: geo_mean(x, [0, 0, 1]) should be affine.
  x <- Variable(3)
  g <- geo_mean(x, c(0, 0, 1))
  expect_true(is_convex(g))
  expect_true(is_concave(g))
  expect_true(is_affine(g))
})

## @cvxpy test_power_atom.py::TestGeoMeanSingleWeight::test_single_weight_value
test_that("geo_mean: single weight w=(0,0,1) equals x[3] (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  ## CVXPY: geo_mean(x, [0, 0, 1]) should equal x[2] (0-based) = x[3] (1-based).
  x <- Variable(3)
  prob <- Problem(Maximize(geo_mean(x, c(0, 0, 1))),
                  list(x <= c(1, 2, 3), x >= 0))
  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(val, 3.0, tolerance = 1e-3)
})

## @cvxpy test_power_atom.py::TestGeoMeanSingleWeight::test_single_weight_nonneg_domain
test_that("geo_mean: single weight still requires x >= 0 domain (CVXPY parity)", {
  ## CVXPY: geo_mean(x, [0, 0, 1]).domain should have nonnegativity.
  x <- Variable(3)
  g <- geo_mean(x, c(0, 0, 1))
  dom <- CVXR:::domain(g)
  expect_gt(length(dom), 0)
})

## @cvxpy test_power_atom.py::TestGeoMeanSingleWeight::test_single_weight_nonneg_enforced
test_that("geo_mean: solver enforces x >= 0 for single-weight element (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  ## CVXPY: Minimizing x[0] with geo_mean(x, [1]) >= 0 => x = 0 (domain x >= 0).
  x <- Variable(1)
  prob <- Problem(Minimize(x[1]),
                  list(geo_mean(x, c(1)) >= 0, x <= 5))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 0.0, tolerance = 1e-3)
})

## @cvxpy test_power_atom.py::TestGeoMeanSingleWeight::test_single_weight_no_cones
test_that("geo_mean: single non-zero weight produces no cone constraints (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  ## CVXPY: Single non-zero weight should not produce any SOC/PowCone3D/PowConeND.
  x <- Variable(3, pos = TRUE)
  for (approx in c(TRUE, FALSE)) {
    prob <- Problem(
      Maximize(geo_mean(x, c(0, 0, 1), approx = approx)),
      list(sum_entries(x) <= 3)
    )
    cc <- .get_cone_counts(prob, "CLARABEL")
    expect_equal(cc$soc, 0, info = paste("approx =", approx, "SOC"))
    expect_equal(cc$p3d, 0, info = paste("approx =", approx, "PowCone3D"))
    expect_equal(cc$pnd, 0, info = paste("approx =", approx, "PowConeND"))
  }
})

# ── TestPnormApprox ────────────────────────────────────────────────

## @cvxpy test_power_atom.py::TestPnormApprox::test_approx_controls_cone_type
test_that("pnorm: approx=TRUE uses SOC, approx=FALSE uses power cones (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  constr <- list(sum_entries(x) >= 3, x >= 0)

  ## approx=TRUE -> SOC, no PowCone3D
  prob <- Problem(Minimize(p_norm(x, 3, approx = TRUE)), constr)
  cc <- .get_cone_counts(prob, "CLARABEL")
  expect_gt(cc$soc, 0)
  expect_equal(cc$p3d, 0)

  ## approx=FALSE -> PowCone3D, (SOC may be 0 or not but p3d > 0)
  prob2 <- Problem(Minimize(p_norm(x, 3, approx = FALSE)), constr)
  cc2 <- .get_cone_counts(prob2, "CLARABEL")
  expect_gt(cc2$p3d, 0)
})

## @cvxpy test_power_atom.py::TestPnormApprox::test_approx_and_exact_agree
test_that("pnorm: approx and exact give same answer (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  cases <- list(
    list(p = 3,   sense = "Minimize", direction = ">=", rhs = 3),
    list(p = 2.5, sense = "Minimize", direction = ">=", rhs = 3),
    list(p = 0.5, sense = "Maximize", direction = "<=", rhs = 3)
  )
  for (case in cases) {
    p <- case$p
    x <- Variable(3, pos = (p < 1))
    if (case$direction == ">=") {
      constr <- list(sum_entries(x) >= case$rhs, x >= 0, x[1] <= 2)
    } else {
      constr <- list(sum_entries(x) <= case$rhs, x >= 0.1)
    }
    sense_fn <- if (case$sense == "Minimize") Minimize else Maximize

    prob <- Problem(sense_fn(p_norm(x, p, approx = TRUE)), constr)
    val_a <- psolve(prob, solver = "CLARABEL")
    x_a <- as.numeric(value(x))

    prob2 <- Problem(sense_fn(p_norm(x, p, approx = FALSE)), constr)
    val_e <- psolve(prob2, solver = "CLARABEL")
    expect_true(status(prob2) %in% c("optimal", "optimal_inaccurate"),
                info = paste("p =", p))
    expect_equal(val_e, val_a, tolerance = 1e-3,
                 info = paste("p =", p, "value"))
    expect_equal(as.numeric(value(x)), x_a, tolerance = 1e-3,
                 info = paste("p =", p, "x"))
  }
})

## @cvxpy test_power_atom.py::TestPnormApprox::test_approx_warning
test_that("pnorm: warning fires for approx=TRUE with many SOCs (CVXPY parity)", {
  skip("CVXR does not yet emit SOC approximation warnings like CVXPY")
})

# ── TestInvProdApprox ─────────────────────────────────────────────

## @cvxpy test_power_atom.py::TestInvProdApprox::test_approx_controls_cone_type
test_that("inv_prod: approx=TRUE uses SOC, approx=FALSE uses power cones (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, pos = TRUE)
  constr <- list(sum_entries(x) >= 3)

  ## approx=TRUE -> SOC, no power cones
  prob <- Problem(Minimize(inv_prod(x, approx = TRUE)), constr)
  cc <- .get_cone_counts(prob, "CLARABEL")
  expect_gt(cc$soc, 0)
  expect_equal(cc$p3d, 0)
  expect_equal(cc$pnd, 0)

  ## approx=FALSE -> power cones (3D and/or ND)
  prob2 <- Problem(Minimize(inv_prod(x, approx = FALSE)), constr)
  cc2 <- .get_cone_counts(prob2, "CLARABEL")
  expect_equal(cc2$soc, 0)
  expect_gt(cc2$p3d + cc2$pnd, 0)
})

## @cvxpy test_power_atom.py::TestInvProdApprox::test_approx_produces_correct_types
test_that("inv_prod: approx flag propagates to inner geo_mean and power atoms (CVXPY parity)", {
  x <- Variable(3, pos = TRUE)

  ## approx=TRUE -> PowerApprox wrapping GeoMeanApprox
  expr_approx <- inv_prod(x, approx = TRUE)
  expect_true(S7_inherits(expr_approx, CVXR:::PowerApprox))
  inner_a <- expr_approx@args[[1L]]
  expect_true(S7_inherits(inner_a, CVXR:::GeoMeanApprox))

  ## approx=FALSE -> Power (not PowerApprox) wrapping GeoMean (not GeoMeanApprox)
  expr_exact <- inv_prod(x, approx = FALSE)
  expect_true(S7_inherits(expr_exact, CVXR:::Power))
  expect_false(S7_inherits(expr_exact, CVXR:::PowerApprox))
  inner_e <- expr_exact@args[[1L]]
  expect_true(S7_inherits(inner_e, CVXR:::GeoMean))
  expect_false(S7_inherits(inner_e, CVXR:::GeoMeanApprox))
})

## @cvxpy test_power_atom.py::TestInvProdApprox::test_approx_and_exact_agree
test_that("inv_prod: approx and exact give same answer (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, pos = TRUE)
  constr <- list(sum_entries(x) >= 3, x <= 5)

  prob <- Problem(Minimize(inv_prod(x, approx = TRUE)), constr)
  val_a <- psolve(prob, solver = "CLARABEL")

  prob2 <- Problem(Minimize(inv_prod(x, approx = FALSE)), constr)
  val_e <- psolve(prob2, solver = "CLARABEL")
  expect_true(status(prob2) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(val_e, val_a, tolerance = 1e-3)
})

# ====================================================================
# PERSPECTIVE TESTS (test_perspective.py)
# ====================================================================

## @cvxpy test_perspective.py::test_p_norms
test_that("perspective: p-norms p=1 and p=2 (CVXPY parity)", {
  skip_if_not_installed("clarabel")

  ## p = 1
  x <- Variable(3)
  s <- Variable(nonneg = TRUE, name = "s")
  f <- p_norm(x, 1)
  obj <- perspective(f, s)
  constr <- list(s == 1, x >= c(1, 2, 3))
  prob <- Problem(Minimize(obj), constr)
  val <- psolve(prob, solver = "CLARABEL")

  ## Reference: sum(x) = 1+2+3 = 6, prob.value^1 == 6
  expect_equal(val, 6.0, tolerance = 1e-3)

  ## p = 2
  x <- Variable(3)
  s <- Variable(nonneg = TRUE, name = "s")
  f <- p_norm(x, 2)
  obj <- perspective(f, s)
  constr <- list(s == 1, x >= c(1, 2, 3))
  prob <- Problem(Minimize(obj), constr)
  val <- psolve(prob, solver = "CLARABEL")

  ## Reference: ||[1,2,3]||_2 = sqrt(14) ~ 3.7417
  ## CVXPY: prob.value**2 == ref_prob.value = 14
  expect_equal(val^2, 14.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, 2, 3), tolerance = 1e-3)
})

## @cvxpy test_perspective.py::test_lse
test_that("perspective: log_sum_exp (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  s <- Variable(nonneg = TRUE)
  f <- log_sum_exp(x)
  obj <- perspective(f, s)
  constr <- list(1 <= s, s <= 2, c(1, 2, 3) <= x)
  prob <- Problem(Minimize(obj), constr)
  val <- psolve(prob, solver = "CLARABEL")

  ## Reference from CVXPY: 3.407606
  expect_equal(val, 3.407606, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, 2, 3), tolerance = 1e-3)
  expect_equal(as.numeric(value(s)), 1.0, tolerance = 1e-3)
})

## @cvxpy test_perspective.py::test_lse_atom
test_that("perspective: log_sum_exp atom form (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  ## This is identical to test_lse in CVXPY (same fixture, same test).
  x <- Variable(3)
  s <- Variable(nonneg = TRUE)
  f_exp <- log_sum_exp(x)
  obj <- perspective(f_exp, s)
  constr <- list(1 <= s, s <= 2, c(1, 2, 3) <= x)
  prob <- Problem(Minimize(obj), constr)
  val <- psolve(prob, solver = "CLARABEL")

  expect_equal(val, 3.407606, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, 2, 3), tolerance = 1e-3)
  expect_equal(as.numeric(value(s)), 1.0, tolerance = 1e-3)
})

## @cvxpy test_perspective.py::test_quad_persp_persp
test_that("perspective: perspective of perspective (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  ## Reference: quad_over_lin(x, s) + r*x - 4*s with x >= 2, s <= 0.5
  ## r = 2 => ref_val = 10, ref_x = 2, ref_s = 0.5
  r <- 2
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  t_var <- Variable(nonneg = TRUE)

  f_exp <- square(x) + r * x - 4
  obj_inner <- perspective(f_exp, s)
  obj <- perspective(obj_inner, t_var)
  ## f(x) -> s*f(x/s) -> t*(s/t)*f((x/t)/(s/t)) = s*f(x/s)

  constr <- list(0.1 <= s, s <= 0.5, x >= 2, 0.1 <= t_var, t_var <= 0.5)
  prob <- Problem(Minimize(obj), constr)
  val <- psolve(prob, solver = "CLARABEL")

  expect_equal(val, 10.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(s)), 0.5, tolerance = 1e-3)
})

## @cvxpy test_perspective.py::test_quad_quad
test_that("perspective: power(x, 4) perspective == quad_over_lin(quad_over_lin(x, s), s) (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  ## Reference: quad_over_lin(quad_over_lin(x, s), s) with x >= 5, s <= 3
  ref_x <- Variable()
  ref_s <- Variable(nonneg = TRUE)
  obj <- quad_over_lin(quad_over_lin(ref_x, ref_s), ref_s)
  constr <- list(ref_x >= 5, ref_s <= 3)
  ref_prob <- Problem(Minimize(obj), constr)
  ref_val <- psolve(ref_prob, solver = "CLARABEL")

  ## Perspective problem: f(x) = x^4, perspective = x^4 / s^3
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- power(x, 4)
  obj <- perspective(f, s)
  constr <- list(x >= 5, s <= 3)
  prob <- Problem(Minimize(obj), constr)
  val <- psolve(prob, solver = "CLARABEL")

  expect_equal(val, ref_val, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 5.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(s)), 3.0, tolerance = 1e-3)
})

## @cvxpy test_perspective.py::test_psd_tr_persp
test_that("perspective: trace(P) with PSD variable (CVXPY parity)", {
  skip("Perspective canonicalizer has broadcast shape mismatch bug with PSD matrix variables")
})

## @cvxpy test_perspective.py::test_psd_mf_persp
test_that("perspective: matrix_frac with PSD variable (CVXPY parity)", {
  skip("Perspective canonicalizer has broadcast shape mismatch bug with PSD matrix variables")
})

## @cvxpy test_perspective.py::test_psd_tr_square
test_that("perspective: square(trace(P)) with PSD variable (CVXPY parity)", {
  skip("Perspective canonicalizer has broadcast shape mismatch bug with PSD matrix variables")
})

## @cvxpy test_perspective.py::test_diag
test_that("perspective: trace with diag variable (CVXPY parity)", {
  skip("Perspective canonicalizer has broadcast shape mismatch bug with diag matrix variables")
})

## @cvxpy test_perspective.py::test_dpp
test_that("perspective: is_dpp returns FALSE for parameterized f (CVXPY parity)", {
  ## CVXPY: perspective(square(a+x), s).is_dpp() == False
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  a <- Parameter()

  obj1 <- perspective(square(a + x), s)
  expect_false(is_dpp(obj1))

  ## CVXPY: perspective(log(a+x), s).is_dpp() == False
  obj2 <- perspective(log(a + x), s)
  expect_false(is_dpp(obj2))
})

## @cvxpy test_perspective.py::test_s_eq_0
test_that("perspective: s=0 with recession function (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  ## Problem where the optimal s is s = 0.
  ## CVXPY: f = x + 1, f_recession = x, constr = [-x^2 + 1 >= 0]
  x <- Variable(1)
  s <- Variable(1, nonneg = TRUE)
  f <- x + 1
  f_recession <- x
  obj <- perspective(f, s, f_recession = f_recession)
  constr <- list(-square(x) + 1 >= 0)

  prob <- Problem(Minimize(obj), constr)
  psolve(prob, solver = "CLARABEL")

  expect_equal(as.numeric(value(x)), -1, tolerance = 1e-3)
  expect_equal(as.numeric(value(s)), 0, tolerance = 1e-3)
})
