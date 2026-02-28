## Phase 9b: CVXPY Parity Test Expansion
## Tests ported from CVXPY test files:
##   test_kron_canon.py, test_convolution.py, test_quad_form.py,
##   test_expressions.py
##
## Reference values obtained via `uv run python` against CVXPY 1.8.1

# ═══════════════════════════════════════════════════════════════════
# 1. Kron Solve Tests (test_kron_canon.py)
# ═══════════════════════════════════════════════════════════════════

## ── gen_kronr_const: kron(Constant C, Variable Z) ────────────────
## CVXPY SOURCE: test_kron_canon.py::TestKronRightVar::test_gen_kronr_const

## @cvxpy test_kron_canon.py::TestKronRightVar::test_gen_kronr_const
test_that("kron solve: C(1x1) on left, Variable Z(2x2) on right", {
  ## np.random.seed(0): C=[[0.55]], L=[[0.72,0.60],[0.54,0.42]]
  C <- Constant(matrix(0.55, 1, 1))
  L <- matrix(c(0.72, 0.54, 0.60, 0.42), 2, 2)  # col-major
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Z)),
                  list(kron(C, Z) >= kron(C, Constant(L)), Z >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, sum(L), tolerance = 1e-3)
  expect_equal(value(Z), L, tolerance = 1e-3)
})

## @cvxpy test_kron_canon.py::TestKronRightVar::test_gen_kronr_const
test_that("kron solve: C(2x1) on left, Variable Z(2x2) on right", {
  ## np.random.seed(0): C=[[0.55],[0.72]], L=[[0.60,0.54],[0.42,0.65]]
  C <- Constant(matrix(c(0.55, 0.72), 2, 1))
  L <- matrix(c(0.60, 0.42, 0.54, 0.65), 2, 2)
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Z)),
                  list(kron(C, Z) >= kron(C, Constant(L)), Z >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, sum(L), tolerance = 1e-3)
  expect_equal(value(Z), L, tolerance = 1e-3)
})

## @cvxpy test_kron_canon.py::TestKronRightVar::test_gen_kronr_const
test_that("kron solve: C(1x2) on left, Variable Z(2x2) on right", {
  ## np.random.seed(0): C=[[0.55,0.72]], L=[[0.60,0.54],[0.42,0.65]]
  C <- Constant(matrix(c(0.55, 0.72), 1, 2))
  L <- matrix(c(0.60, 0.42, 0.54, 0.65), 2, 2)
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Z)),
                  list(kron(C, Z) >= kron(C, Constant(L)), Z >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, sum(L), tolerance = 1e-3)
  expect_equal(value(Z), L, tolerance = 1e-3)
})

## @cvxpy test_kron_canon.py::TestKronRightVar::test_gen_kronr_const
test_that("kron solve: C(2x2) on left, Variable Z(2x2) on right", {
  ## np.random.seed(0): C=[[0.55,0.72],[0.60,0.54]], L=[[0.42,0.65],[0.44,0.89]]
  C <- Constant(matrix(c(0.55, 0.60, 0.72, 0.54), 2, 2))
  L <- matrix(c(0.42, 0.44, 0.65, 0.89), 2, 2)
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Z)),
                  list(kron(C, Z) >= kron(C, Constant(L)), Z >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, sum(L), tolerance = 1e-3)
  expect_equal(value(Z), L, tolerance = 1e-3)
})

## ── symvar_kronl_const: kron(symmetric Var X, Constant b) ────────
## CVXPY SOURCE: test_kron_canon.py::TestKronLeftVar::symvar_kronl

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_symvar_kronl_const
test_that("kron solve: symmetric Variable X(2x2) on left, scalar b on right (min)", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  b <- Constant(matrix(1.5, 1, 1))
  L <- matrix(c(0.5, 2, 1, 3), 2, 2)   # col-major: [[0.5,1],[2,3]]
  U <- matrix(c(10, 12, 11, 13), 2, 2)  # col-major: [[10,11],[12,13]]
  kronX <- kron(X, b)
  prob <- Problem(Minimize(sum_entries(X)), list(U >= kronX, kronX >= L))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  ## X_min ≈ [[0.5, 2], [2, 3]] / 1.5 (symmetry forces off-diag equal)
  expect_equal(val, 5.0, tolerance = 1e-3)
  Xval <- value(X)
  expect_equal(Xval[1, 1], 0.5 / 1.5, tolerance = 1e-3)
  expect_equal(Xval[2, 2], 3.0 / 1.5, tolerance = 1e-3)
  expect_equal(Xval[1, 2], Xval[2, 1], tolerance = 1e-6) # symmetric
})

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_symvar_kronl_const
test_that("kron solve: symmetric Variable X(2x2) on left, scalar b on right (max)", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  b <- Constant(matrix(1.5, 1, 1))
  L <- matrix(c(0.5, 2, 1, 3), 2, 2)
  U <- matrix(c(10, 12, 11, 13), 2, 2)
  kronX <- kron(X, b)
  prob <- Problem(Maximize(sum_entries(X)), list(U >= kronX, kronX >= L))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  ## X_max ≈ [[10, 11], [11, 13]] / 1.5 (symmetry forces off-diag to min(11,12)/1.5)
  expect_equal(val, 30.0, tolerance = 1e-3)
  Xval <- value(X)
  expect_equal(Xval[1, 2], Xval[2, 1], tolerance = 1e-6) # symmetric
})

## ── scalar_kronl_const: kron(scalar Var y, Constant A) ───────────
## CVXPY SOURCE: test_kron_canon.py::TestKronLeftVar::scalar_kronl

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_scalar_kronl_const
test_that("kron solve: scalar Variable y on left, matrix A on right (min/max)", {
  y <- Variable(c(1, 1))
  A_val <- matrix(c(1, 3, 2, 4), 2, 2)  # col-major: [[1,2],[3,4]]
  L <- matrix(c(0.5, 2, 1, 3), 2, 2)
  U <- matrix(c(10, 12, 11, 13), 2, 2)
  krony <- kron(y, Constant(A_val))
  ## Min y s.t. U >= kron(y,A) >= L
  prob_min <- Problem(Minimize(y), list(U >= krony, krony >= L))
  val_min <- psolve(prob_min, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_min), OPTIMAL)
  ## y_min = max(L / A_val) = 0.75
  expect_equal(val_min, 0.75, tolerance = 1e-3)
  ## Max y
  prob_max <- Problem(Maximize(y), list(U >= krony, krony >= L))
  val_max <- psolve(prob_max, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_max), OPTIMAL)
  ## y_max = min(U / A_val) = 3.25
  expect_equal(val_max, 3.25, tolerance = 1e-3)
})

## ── gen_kronl_const: kron(Variable Z, Constant C) ────────────────
## CVXPY SOURCE: test_kron_canon.py::TestKronLeftVar::test_gen_kronr_const (var_left=TRUE)

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_gen_kronr_const
test_that("kron solve: Variable Z(2x2) on left, C(1x1) on right", {
  C <- Constant(matrix(0.55, 1, 1))
  L <- matrix(c(0.72, 0.54, 0.60, 0.42), 2, 2)
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Z)),
                  list(kron(Z, C) >= kron(Constant(L), C), Z >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, sum(L), tolerance = 1e-3)
  expect_equal(value(Z), L, tolerance = 1e-3)
})

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_gen_kronr_const
test_that("kron solve: Variable Z(2x2) on left, C(2x2) on right", {
  ## Fixed in v0.12.0: Added KRON_L C++ enum + get_kron_l_mat function.
  ## CVXPY gives Z≈L (sum=2.4); CVXR now matches.
  C <- Constant(matrix(c(0.55, 0.60, 0.72, 0.54), 2, 2))
  L <- matrix(c(0.42, 0.44, 0.65, 0.89), 2, 2)
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Z)),
                  list(kron(Z, C) >= kron(Constant(L), C), Z >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, sum(L), tolerance = 1e-3)
  expect_equal(value(Z), L, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# 2. Convolution Tests (test_convolution.py)
# ═══════════════════════════════════════════════════════════════════

## CVXPY SOURCE: test_convolution.py::test_1D_conv

## @cvxpy test_convolution.py::TestConvolution::test_1D_conv
test_that("conv solve: minimize norm1(conv(f, g)) with g fixed", {
  f <- c(1, 2, 3)
  g <- Variable(3)
  g_val <- c(0, 1, 0.5)
  ## conv(f, g_val) = [0, 1, 2.5, 4, 1.5], sum of abs = 9.0
  prob <- Problem(Minimize(norm1(conv(f, g))), list(g == g_val))
  val <- psolve(prob, solver = SCS_SOLVER)
  expect_equal(val, 9.0, tolerance = 1e-2)
  expect_equal(status(prob), OPTIMAL)
})

## @cvxpy test_convolution.py::TestConvolution::test_1D_conv
test_that("conv: constant conv produces correct value", {
  f <- c(1, 2, 3)
  g <- c(0, 1, 0.5)
  expr <- conv(f, g)
  expect_true(is_constant(expr))
  expect_equal(expr@shape, c(5L, 1L))
  expected <- matrix(c(0, 1, 2.5, 4, 1.5), 5, 1)
  expect_equal(value(expr), expected, tolerance = 1e-10)
})

## CVXPY SOURCE: test_convolution.py::test_convolve / scalar shapes

## @cvxpy test_convolution.py::TestConvolution::test_convolve
test_that("conv: scalar * vector", {
  expr <- conv(Constant(2), Constant(c(0, 1, 0.5)))
  expect_equal(expr@shape, c(3L, 1L))
  expect_equal(as.numeric(value(expr)), c(0, 2, 1), tolerance = 1e-10)
})

## @cvxpy test_convolution.py::TestConvolution::test_convolve
test_that("conv: vector * scalar", {
  expr <- conv(Constant(c(1, 2, 3)), Constant(2))
  expect_equal(expr@shape, c(3L, 1L))
  expect_equal(as.numeric(value(expr)), c(2, 4, 6), tolerance = 1e-10)
})

## @cvxpy test_convolution.py::TestConvolution::test_conv_prob
## CVXPY SOURCE: test_convolution.py::test_conv_prob (unbounded)

## @cvxpy NONE
test_that("conv: unconstrained convolution should return UNBOUNDED", {
  h <- c(1, -1)
  x <- Variable(5)
  v <- conv(h, x)
  expect_equal(v@shape, c(6L, 1L))
  prob <- Problem(Minimize(sum_entries(v[1:5, ])))
  psolve(prob, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(status(prob), UNBOUNDED)
})

## CVXPY SOURCE: test_convolution.py::test_sparse_convolution

## @cvxpy test_convolution.py::TestConvolution::test_sparse_convolution
test_that("conv: sparse convolution constant evaluation", {
  a <- c(0, 1)
  c_val <- c(1, 0, 2)
  r <- conv(Constant(a), Constant(c_val))
  expect_equal(r@shape, c(4L, 1L))
  expect_equal(as.numeric(value(r)), c(0, 1, 0, 2), tolerance = 1e-10)
})

## CVXPY SOURCE: test_convolution.py::test_convolve (2D rejection)

## @cvxpy test_convolution.py::TestConvolution::test_convolve
test_that("conv: 2D input is rejected", {
  A <- Variable(c(3, 2))
  expect_error(conv(c(1, 2), A), "vector")
})

# ═══════════════════════════════════════════════════════════════════
# 3. Quad Form Edge Cases (test_quad_form.py)
# ═══════════════════════════════════════════════════════════════════

## CVXPY SOURCE: test_quad_form.py — non-PSD DCP error

## @cvxpy test_quad_form.py::TestNonOptimal::test_non_psd
test_that("quad_form: indefinite P is not DCP for Minimize", {
  x <- Variable(2)
  P <- matrix(c(1, 0, 0, -1), 2, 2)
  qf <- quad_form(x, P)
  prob <- Problem(Minimize(qf), list(x == c(1, 2)))
  ## CVXPY raises DCPError; in CVXR, psolve should error
  expect_error(psolve(prob, solver = CLARABEL_SOLVER), "DCP")
})

## CVXPY SOURCE: test_quad_form.py — NSD P for Maximize is concave/DCP

## @cvxpy test_quad_form.py::TestNonOptimal::test_non_psd
test_that("quad_form: NSD P is concave (Maximize is DCP)", {
  x <- Variable(2)
  P <- matrix(c(-1, 0, 0, -2), 2, 2)
  qf <- quad_form(x, P)
  expect_true(is_concave(qf))
  prob <- Problem(Maximize(qf), list(x == c(1, 2)))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  ## x'Px = 1*(-1)*1 + 2*(-2)*2 = -1 - 8 = -9
  expect_equal(val, -9.0, tolerance = 1e-3)
})

## CVXPY SOURCE: test_quad_form.py — zero P matrix

## @cvxpy test_quad_form.py::TestNonOptimal::test_zero_matrix
test_that("quad_form: zero P reduces to LP", {
  x <- Variable(3)
  P_zero <- matrix(0, 3, 3)
  c_vec <- c(-1, -1, -1)
  ## 0.5 * quad_form(x, 0) + c'x s.t. 0 <= x <= 1
  prob <- Problem(Minimize(0.5 * quad_form(x, P_zero) + t(c_vec) %*% x),
                  list(x >= 0, x <= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, -3.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-3)
})

## CVXPY SOURCE: test_quad_form.py — zero coefficient on quad_form

## @cvxpy test_quad_form.py::TestNonOptimal::test_zero_term
test_that("quad_form: zero coefficient makes quad_form vanish", {
  x <- Variable(2)
  target <- c(3, 4)
  prob <- Problem(Minimize(0 * quad_form(x, diag(2)) + sum_squares(x - target)))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(3, 4), tolerance = 1e-3)
})

## Objective eval consistency

## @cvxpy NONE
test_that("quad_form: objective value consistency after solve", {
  x <- Variable(3)
  B <- matrix(c(1, 0.5, 0, 0.5, 1, 1), 3, 2)  # 3x2 matrix
  A <- diag(2)
  ## Minimize -sum(B'x) + quad_form(B'x, I)
  obj_expr <- -sum_entries(t(B) %*% x) + quad_form(t(B) %*% x, A)
  prob <- Problem(Minimize(obj_expr), list(x >= 0, x <= 1))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, -0.5, tolerance = 1e-2)
})

# ═══════════════════════════════════════════════════════════════════
# 4. Expression Indexing Tests (test_expressions.py)
# ═══════════════════════════════════════════════════════════════════

## CVXPY SOURCE: test_expressions.py::test_index_expression

## @cvxpy test_expressions.py::TestExpressions::test_index_expression
test_that("indexing: basic shapes for Variable", {
  x <- Variable(5)
  ## Single element
  expect_equal(x[1, ]@shape, c(1L, 1L))
  ## Slice
  expect_equal(x[2:4, ]@shape, c(3L, 1L))

  A <- Variable(c(3, 3))
  expect_equal(A[1, 1]@shape, c(1L, 1L))
  expect_equal(A[1, ]@shape, c(1L, 3L))
  expect_equal(A[, 2]@shape, c(3L, 1L))
})

## @cvxpy test_expressions.py::TestExpressions::test_index_expression
test_that("indexing: curvature preservation", {
  x <- Variable(5)
  ## Variable index is affine
  expect_true(is_affine(x[1, ]))
  ## Constant index is constant
  c_val <- Constant(c(10, 20, 30))
  expect_true(is_constant(c_val[c(1, 3), ]))
})

## @cvxpy test_expressions.py::TestExpressions::test_index_expression
test_that("indexing: sign preservation for nonneg/nonpos", {
  x <- Variable(5, nonneg = TRUE)
  expect_true(is_nonneg(x[1:3, ]))

  y <- Variable(5, nonpos = TRUE)
  expect_true(is_nonpos(y[2:4, ]))
})

## @cvxpy test_expressions.py::TestExpressions::test_index_expression
test_that("indexing: Constant numeric value", {
  c_val <- Constant(c(10, 20, 30))
  res <- value(c_val[c(1, 3), ])
  expect_equal(as.numeric(res), c(10, 30))

  c_mat <- Constant(matrix(1:6, 2, 3))
  ## c_mat[1,2] = 3 (R col-major: [1,1]=1, [2,1]=2, [1,2]=3, ...)
  res2 <- value(c_mat[1, 2])
  expect_equal(as.numeric(res2), 3)
})

## CVXPY SOURCE: test_expressions.py::test_neg_indices

## @cvxpy test_expressions.py::TestExpressions::test_neg_indices
test_that("indexing: negative indices exclude elements", {
  x <- Variable(5)
  ## x[-1, ] excludes first → shape (4,1) in R
  expect_equal(x[-1, ]@shape, c(4L, 1L))
  ## x[-c(1,2), ] excludes first two → shape (3,1)
  expect_equal(x[-c(1, 2), ]@shape, c(3L, 1L))
})

## @cvxpy test_expressions.py::TestExpressions::test_neg_indices
test_that("indexing: negative index on Constant gives correct value", {
  c_val <- Constant(matrix(c(1, 2, 3, 4), 2, 2))
  ## c_val[-1, -1] = last element = 4 (CVXPY: c[-1,-1] with neg indexing)
  ## In R, -1 means "exclude row 1, exclude col 1" → c_val[2,2] = 4
  res <- value(c_val[-1, -1])
  expect_equal(as.numeric(res), 4)
})

## @cvxpy test_expressions.py::TestExpressions::test_logical_indices
## CVXPY SOURCE: test_expressions.py::test_logical_indices

## @cvxpy NONE
test_that("indexing: logical indexing", {
  x <- Variable(4)
  idx <- c(TRUE, FALSE, TRUE, FALSE)
  expect_equal(x[idx, ]@shape, c(2L, 1L))
})

## CVXPY SOURCE: test_expressions.py::test_selector_list_indices

## @cvxpy test_expressions.py::TestExpressions::test_selector_list_indices
test_that("indexing: integer vector indexing", {
  x <- Variable(5)
  expect_equal(x[c(1, 3, 5), ]@shape, c(3L, 1L))
})

# ═══════════════════════════════════════════════════════════════════
# 5. Expression Operator Edge Cases (test_expressions.py)
# ═══════════════════════════════════════════════════════════════════

## CVXPY SOURCE: test_expressions.py::test_add_expression (shape mismatch)

## @cvxpy test_expressions.py::TestExpressions::test_add_expression
test_that("operator: add shape mismatch errors", {
  x <- Variable(2)
  y <- Variable(3)
  expect_error(x + y)
})

## CVXPY SOURCE: test_expressions.py::test_mul_expression (var*var)

## @cvxpy test_expressions.py::TestExpressions::test_mul_expression
test_that("operator: Variable * Variable is quadratic (not affine)", {
  x <- Variable(2)
  y <- Variable(2)
  expr <- x * y
  expect_true(is_quadratic(expr))
  expect_false(is_affine(expr))
})

## CVXPY SOURCE: test_expressions.py::test_div_expression

## @cvxpy test_expressions.py::TestExpressions::test_div_expression
test_that("operator: x/2 is affine", {
  x <- Variable(2)
  expr <- x / 2
  expect_true(is_affine(expr))
  expect_equal(expr@shape, c(2L, 1L))
})

## CVXPY SOURCE: test_expressions.py::test_neg_expression

## @cvxpy test_expressions.py::TestExpressions::test_neg_expression
test_that("operator: negation flips sign", {
  x <- Variable(3, nonneg = TRUE)
  expr <- -x
  expect_true(is_nonpos(expr))
  expect_true(is_affine(expr))
})

## @cvxpy test_expressions.py::TestExpressions::test_neg_expression
test_that("operator: -square(x) is concave", {
  x <- Variable(3)
  expr <- -square(x)
  expect_true(is_concave(expr))
})

## CVXPY SOURCE: test_expressions.py::test_scalar_const_promotion

## @cvxpy test_expressions.py::TestExpressions::test_scalar_const_promotion
test_that("operator: scalar promotion shapes", {
  x <- Variable(3)
  expect_equal((2 + x)@shape, c(3L, 1L))
  expect_equal((x + 2)@shape, c(3L, 1L))

  A <- Variable(c(2, 2))
  expect_equal((A + 1)@shape, c(2L, 2L))
  expect_equal((3 * A)@shape, c(2L, 2L))
})

# ═══════════════════════════════════════════════════════════════════
# 6. is_pwl Query Tests (test_expressions.py)
# ═══════════════════════════════════════════════════════════════════

## CVXPY SOURCE: test_expressions.py::test_is_pwl

## @cvxpy test_expressions.py::TestExpressions::test_is_pwl
test_that("is_pwl: affine expressions are PWL", {
  x <- Variable(3)
  expect_true(is_pwl(x))
  expect_true(is_pwl(Constant(5)))
  A <- matrix(1, 2, 3)
  b <- c(1, 1)
  expr <- A %*% x - b
  expect_true(is_pwl(expr))
})

## @cvxpy test_expressions.py::TestExpressions::test_is_pwl
test_that("is_pwl: PWL atoms", {
  y <- Variable(3)
  expect_true(is_pwl(max_elemwise(1, 3 * y)))
  expect_true(is_pwl(abs(y)))
  expect_true(is_pwl(norm1(y)))
  expect_true(is_pwl(min_entries(y)))
})

## @cvxpy test_expressions.py::TestExpressions::test_is_pwl
test_that("is_pwl: nonlinear atoms are not PWL", {
  y <- Variable(3)
  expect_false(is_pwl(square(y)))
  expect_false(is_pwl(quad_form(y, diag(3))))
})

## @cvxpy test_expressions.py::TestExpressions::test_is_pwl
test_that("is_pwl: sum_largest is PWL", {
  y <- Variable(3)
  expect_true(is_pwl(sum_largest(y, 2)))
})

# ═══════════════════════════════════════════════════════════════════
# 7. Symmetric Variable Tests (test_expressions.py)
# ═══════════════════════════════════════════════════════════════════

## CVXPY SOURCE: test_expressions.py::test_symmetric

## @cvxpy test_expressions.py::TestExpressions::test_symmetric
test_that("symmetric: construction and validation", {
  v <- Variable(c(2, 2), symmetric = TRUE)
  expect_true(is_symmetric(v))
  ## Non-square should error
  expect_error(Variable(c(4, 3), symmetric = TRUE))
})

## @cvxpy test_expressions.py::TestExpressions::test_symmetric
test_that("symmetric: arithmetic preserves symmetry", {
  v <- Variable(c(2, 2), symmetric = TRUE)
  expect_true(is_symmetric(v + v))
  expect_true(is_symmetric(-v))
})

## CVXPY SOURCE: test_expressions.py::test_symmetric (PSD/NSD)

## @cvxpy test_expressions.py::TestExpressions::test_symmetric
test_that("PSD/NSD Variable properties", {
  v_psd <- Variable(c(2, 2), PSD = TRUE)
  expect_true(is_psd(v_psd))
  expect_true(is_symmetric(v_psd))

  v_nsd <- Variable(c(2, 2), NSD = TRUE)
  expect_true(is_nsd(v_nsd))
  expect_true(is_symmetric(v_nsd))
})

# ═══════════════════════════════════════════════════════════════════
# 8. Value Assignment Tests (test_expressions.py)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy test_expressions.py::TestExpressions::test_assign_var_value
test_that("value assignment: Variable stores matrix", {
  x <- Variable(2)
  value(x) <- c(1, 2)
  expect_equal(value(x), matrix(c(1, 2), 2, 1))
})

## @cvxpy test_expressions.py::TestExpressions::test_constants
test_that("value: Constant returns correct matrix", {
  c1 <- Constant(c(1, 2, 3))
  expect_equal(value(c1), matrix(c(1, 2, 3), 3, 1))
  c2 <- Constant(matrix(1:6, 2, 3))
  expect_equal(dim(value(c2)), c(2L, 3L))
})
