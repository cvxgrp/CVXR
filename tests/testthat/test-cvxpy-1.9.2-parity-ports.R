## Ports of CVXPY 1.9.2 tests whose subject was ALREADY correct in CVXR.
##
## Each was probed against the installed package before being written, so none
## of these is a fix — they pin behavior that a future edit could regress. The
## code-level 1.9.2 fixes live in test-cvxpy-1.9.2-atom-fixes.R,
## test-hermitian-psd-dual.R and test-bisection-solver-failure.R.

# ===================================================================
# Complex expressions
# ===================================================================

## @cvxpy test_complex.py::TestComplex::test_lambda_sum_largest_real_arg
test_that("lambda_sum_largest of a real matrix survives Complex2Real unchanged", {
  ## Upstream used to double k for real (symmetric) arguments whenever the
  ## problem contained ANY unrelated complex expression. CVXR's
  ## c2r_lambda_sum_largest_canon already guards on imag_args being non-NULL.
  X <- Variable(c(3, 3), symmetric = TRUE)
  cons <- list(lambda_sum_largest(X, 2) <= 5, X %>>% 0, X %<<% (10 * diag(3)))
  obj <- Maximize(matrix_trace(X))
  ref <- psolve(Problem(obj, cons), solver = "CLARABEL")
  expect_equal(ref, 7.5, tolerance = 1e-4)

  ## The same problem, plus an unrelated complex variable.
  z <- Variable(complex = TRUE)
  X2 <- Variable(c(3, 3), symmetric = TRUE)
  cons2 <- list(lambda_sum_largest(X2, 2) <= 5, X2 %>>% 0,
                X2 %<<% (10 * diag(3)), abs(z) <= 1)
  res <- psolve(Problem(Maximize(matrix_trace(X2)), cons2), solver = "CLARABEL")
  expect_equal(res, ref, tolerance = 1e-4)
  expect_lte(as.numeric(value(lambda_sum_largest(X2, 2))), 5 + 1e-4)
})

## @cvxpy test_complex.py::TestComplex::test_scalar_complex_value
test_that("bare complex scalars are accepted as leaf values", {
  ## Upstream's Leaf.project called .astype on the raw value, which a built-in
  ## Python complex lacks. R has no such distinction, but the contract holds.
  p <- Parameter(complex = TRUE)
  value(p) <- 2 + 3i
  expect_equal(as.complex(value(p)), 2 + 3i)

  z <- Variable(complex = TRUE)
  value(z) <- 2 + 3i
  expect_equal(as.complex(value(z)), 2 + 3i)
})

# ===================================================================
# Cone constraints with matrix arguments: duals must be F-ordered
# ===================================================================

## @cvxpy test_constraints.py::TestConstraints::test_exp_cone_matrix_arg_duals
test_that("ExpCone duals with matrix args match the flattened problem", {
  ## ConeMatrixStuffing flattens matrix args in Fortran (column-major) order.
  ## Recovering them in C order permutes the entries whenever both dimensions
  ## exceed one. R is column-major natively, so as.numeric() IS the F-order
  ## flatten -- that is what makes the two formulations comparable here.
  m <- 2L; k <- 3L
  set.seed(4)
  cc <- matrix(runif(m * k, 0.5, 1), m, k)

  x <- Variable(c(m, k)); y <- Variable(c(m, k)); z <- Variable(c(m, k))
  con <- ExpCone(x, y, z)
  psolve(Problem(Minimize(-sum(cc * x) + sum(y) + sum(z)), list(con)),
         solver = "CLARABEL")

  xf <- Variable(m * k); yf <- Variable(m * k); zf <- Variable(m * k)
  conf <- ExpCone(xf, yf, zf)
  psolve(Problem(Minimize(-sum(as.numeric(cc) * xf) + sum(yf) + sum(zf)),
                 list(conf)), solver = "CLARABEL")

  dv <- dual_value(con); dvf <- dual_value(conf)
  expect_true(is.list(dv))
  expect_length(dv, 3L)
  for (i in seq_along(dv)) {
    expect_equal(as.numeric(dv[[i]]), as.numeric(dvf[[i]]), tolerance = 1e-5)
  }
})

## @cvxpy test_constraints.py::TestConstraints::test_pow_cone_3d_matrix_arg_duals
test_that("PowCone3D duals with matrix args match the flattened problem", {
  m <- 2L; k <- 3L; alpha <- 0.4
  set.seed(7)
  cc <- matrix(runif(m * k, 0.5, 1), m, k)

  x <- Variable(c(m, k)); y <- Variable(c(m, k)); z <- Variable(c(m, k))
  con <- PowCone3D(x, y, z, alpha)
  psolve(Problem(Minimize(sum(x) + sum(y) - sum(cc * z)), list(con)),
         solver = "CLARABEL")

  xf <- Variable(m * k); yf <- Variable(m * k); zf <- Variable(m * k)
  conf <- PowCone3D(xf, yf, zf, alpha)
  psolve(Problem(Minimize(sum(xf) + sum(yf) - sum(as.numeric(cc) * zf)),
                 list(conf)), solver = "CLARABEL")

  dv <- dual_value(con); dvf <- dual_value(conf)
  for (i in seq_along(dv)) {
    expect_equal(as.numeric(dv[[i]]), as.numeric(dvf[[i]]), tolerance = 1e-5)
  }
})

# ===================================================================
# Expressions
# ===================================================================

## @cvxpy test_expressions.py::TestExpressions::test_partial_boolean_sign
test_that("a variable boolean at only some indices has unknown sign", {
  ## Upstream's is_nonneg returned the (truthy) index list itself, so the whole
  ## variable was classified NONNEGATIVE. CVXR's own version of this bug --
  ## isTRUE() on a vector -- was fixed in 030e3c1; this pins the sign contract.
  x <- Variable(2, boolean = c(1))
  expect_false(is_nonneg(x))
  expect_equal(expr_sign(x), UNKNOWN_SIGN)

  ## Fully boolean variables remain nonnegative.
  y <- Variable(2, boolean = TRUE)
  expect_true(is_nonneg(y))
})

## @cvxpy test_expressions.py::TestExpressions::test_sparse_symmetry_checks_positions
test_that("asymmetric sparse matrices whose triangle VALUES match are rejected", {
  ## A[2,1] = 5 and A[1,3] = 5: the two triangles hold the same multiset of
  ## values, so any check that compares values without positions passes it.
  A <- Matrix::sparseMatrix(i = c(2L, 1L), j = c(1L, 3L), x = c(5, 5),
                            dims = c(3L, 3L))
  expect_false(is_symmetric(Constant(A)))
  expect_error(quad_form(Variable(3), A))

  B <- Matrix::sparseMatrix(i = c(2L, 1L), j = c(1L, 3L), x = c(5, -5),
                            dims = c(3L, 3L))
  expect_false(is_skew_symmetric(Constant(B)))

  ## Duplicate entries (summed) and explicit zeros must not break detection:
  ## (2,1) gets 2 + 3 = 5, matching (1,2) = 5, and the explicit 0 at (3,1)
  ## has a matching structural zero at (1,3).
  S <- Matrix::sparseMatrix(i = c(2L, 2L, 1L, 3L), j = c(1L, 1L, 2L, 1L),
                            x = c(2, 3, 5, 0), dims = c(3L, 3L))
  expect_true(is_symmetric(Constant(S)))

  ## Sparse and dense classification must agree on random matrices.
  set.seed(0)
  for (trial in 1:20) {
    M <- Matrix::rsparsematrix(5, 5, density = 0.4)
    if (trial %% 2 == 0) M <- M + Matrix::t(M)
    D <- as.matrix(M)
    expect_equal(is_symmetric(Constant(M)),
                 isTRUE(all.equal(D, t(D), tolerance = 1e-12)))
  }
})

# ===================================================================
# Trigonometric NLP atoms
# ===================================================================

## @cvxpy test_nonlinear_atoms.py::test_trig_atoms_metadata_numeric_and_grad
test_that("sin/cos/tan report the documented metadata, value and gradient", {
  val <- c(0.2, -0.3)
  cases <- list(
    list(cls = Sin, numeric = sin(val),  grad = cos(val),      domain = 0L),
    list(cls = Cos, numeric = cos(val),  grad = -sin(val),     domain = 0L),
    ## tan is the only one with a domain: -pi/2 <= x <= pi/2.
    list(cls = Tan, numeric = tan(val),  grad = 1 / cos(val)^2, domain = 2L)
  )
  for (cs in cases) {
    e <- cs$cls(Variable(2))
    expect_equal(as.numeric(numeric_value(e, list(val))), cs$numeric,
                 tolerance = 1e-12)
    sg <- sign_from_args(e)
    expect_false(sg$is_nonneg)
    expect_false(sg$is_nonpos)
    expect_false(is_atom_convex(e))
    expect_false(is_atom_concave(e))
    expect_true(is_atom_smooth(e))
    expect_false(is_incr(e, 1L))
    expect_false(is_decr(e, 1L))
    expect_length(atom_domain(e), cs$domain)
    expect_equal(diag(as.matrix(.grad(e, list(val))[[1L]])), cs$grad,
                 tolerance = 1e-12)
  }
})

# ===================================================================
# PIQP backend selection
# ===================================================================

## @cvxpy test_qp_solvers.py::TestPiqpInterface::test_piqp_dense_backend
test_that("the PIQP dense backend agrees with the default sparse backend", {
  skip_if_not_installed("piqp")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - 3)), list(x[1] + x[2] <= 4))
  psolve(prob, solver = "PIQP", backend = "dense")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(2, 2), tolerance = 1e-4)
})
