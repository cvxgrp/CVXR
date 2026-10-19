## Ports of CVXPY 1.9.3 tests whose subject needed a CVXR code fix.
##
## Each test here was observed FAILING on the pre-fix build (1.9.2.9102)
## before the fix was written. Plan: notes/cvxpy_1.9.3_plan.md Part 1.

# ===================================================================
# P1 (#3384): quadratic objective over a variable with a dim-reducing
# attribute. CvxAttr2Constr substitutes such a variable with an affine
# expression, so the quad canonicalizers must not hand it to
# SymbolicQuadForm as a bare Variable. Pre-fix CVXR crashed in
# coeff_quad_form ("subscript out of bounds"); the oracle is the analytic
# optimum (X = I, value 2, dual 2), which the conic path already gives.
# ===================================================================

## Python sparsity=[(0, 1), (0, 1)] on a 2x2 is the diagonal (0,0), (1,1).
.p1_dim_reducing_attrs <- list(
  symmetric = list(symmetric = TRUE),
  PSD       = list(PSD = TRUE),
  diag      = list(diag = TRUE),
  sparsity  = list(sparsity = cbind(c(1, 2), c(1, 2)))
)

## @cvxpy test_qp_solvers.py::test_quad_obj_dim_reducing_attr
test_that("P1: sum_squares over a dim-reducing-attribute variable solves (value 2, dual 2)", {
  for (solver in c("CLARABEL", "OSQP")) {
    require_solver(solver)
    for (nm in names(.p1_dim_reducing_attrs)) {
      if (nm == "PSD" && solver == "OSQP") next   # OSQP has no PSD cone
      X <- do.call(Variable, c(list(c(2L, 2L)), .p1_dim_reducing_attrs[[nm]]))
      con <- matrix_trace(X) >= 2
      p <- Problem(Minimize(sum_squares(X)), list(con))
      val <- psolve(p, solver = solver)
      expect_equal(val, 2.0, tolerance = 1e-5, info = paste(solver, nm))
      expect_equal(as.numeric(dual_value(con)), 2.0, tolerance = 1e-4,
                   info = paste(solver, nm))
    }
  }
})

## @cvxpy test_qp_solvers.py::test_quad_obj_nsd_attr
test_that("P1: NSD is dim-reducing too (value 2, dual 2)", {
  X <- Variable(c(2L, 2L), NSD = TRUE)
  con <- matrix_trace(X) <= -2
  p <- Problem(Minimize(sum_squares(X)), list(con))
  val <- psolve(p, solver = "CLARABEL")
  expect_equal(val, 2.0, tolerance = 1e-5)
  expect_equal(as.numeric(dual_value(con)), 2.0, tolerance = 1e-4)
})

## @cvxpy test_qp_solvers.py::test_quad_obj_attr_power_canon
test_that("P1: sum(square(X)) over a symmetric X lowers through power_canon", {
  for (solver in c("CLARABEL", "OSQP")) {
    require_solver(solver)
    X <- Variable(c(2L, 2L), symmetric = TRUE)
    con <- matrix_trace(X) >= 2
    p <- Problem(Minimize(sum(square(X))), list(con))
    val <- psolve(p, solver = solver)
    expect_equal(val, 2.0, tolerance = 1e-5, info = solver)
    expect_equal(as.numeric(dual_value(con)), 2.0, tolerance = 1e-4, info = solver)
  }
})

## @cvxpy test_qp_solvers.py::test_quad_obj_attr_quad_form_canon
test_that("P1: quad_form over a sparsity-constrained vector hits quad_form_canon", {
  for (solver in c("CLARABEL", "OSQP")) {
    require_solver(solver)
    ## Python sparsity=[(0, 2)] on shape (4,) = entries 0 and 2 -> R 1 and 3.
    x <- Variable(4L, sparsity = c(1, 3))
    con <- sum(x) >= 2
    p <- Problem(Minimize(quad_form(x, diag(4))), list(con))
    val <- psolve(p, solver = solver)
    ## Support is {1, 3}: optimum x1 = x3 = 1, value 2, dual 2.
    expect_equal(val, 2.0, tolerance = 1e-5, info = solver)
    expect_equal(as.numeric(dual_value(con)), 2.0, tolerance = 1e-4, info = solver)
  }
})

## @cvxpy test_qp_solvers.py::test_quad_obj_symmetric_mixed_dpp
test_that("P1: DPP quad objective mixing a symmetric and a plain variable", {
  for (solver in c("CLARABEL", "OSQP")) {
    require_solver(solver)
    gamma <- Parameter(nonneg = TRUE)
    X <- Variable(c(2L, 2L), symmetric = TRUE)
    y <- Variable(2L)
    con <- matrix_trace(X) + sum(y) >= 4
    prob <- Problem(Minimize(gamma * sum_squares(X) + sum_squares(y)), list(con))
    expect_true(is_dpp(prob))
    for (gv in c(1.0, 4.0)) {
      value(gamma) <- gv
      val <- psolve(prob, solver = solver)
      ## Analytic optimum: value 8g/(1+g), dual 4g/(1+g).
      expect_equal(val, 8 * gv / (1 + gv), tolerance = 1e-5, info = paste(solver, gv))
      expect_equal(as.numeric(dual_value(con)), 4 * gv / (1 + gv), tolerance = 1e-4,
                   info = paste(solver, gv))
    }
  }
})

# ===================================================================
# P2 (#3492): matrix_trace() rewrites trace(A %*% B) as sum(A * t(B)) only
# for a true matrix product -- Multiply (elementwise) subclasses
# MulExpression and must fall through to Trace -- and wraps Real_() only
# when the rewritten sum is not already known to be real.
# ===================================================================

## @cvxpy test_atoms.py::TestAtoms::test_trace_scalar_real_products
test_that("P2: real scalar-shaped products canonicalize without Real_", {
  X <- Variable(c(1L, 1L))
  k <- matrix(2, 1L, 1L)

  matmul_trace <- matrix_trace(k %*% X)
  expect_true(S7_inherits(matmul_trace, SumEntries))
  expect_no_error(canonical_form(matmul_trace))

  elementwise_trace <- matrix_trace(k * X)
  expect_true(S7_inherits(elementwise_trace, Trace))
  expect_no_error(canonical_form(elementwise_trace))
})

## @cvxpy test_atoms.py::TestAtoms::test_trace_AB_hermitian_complex
test_that("P2: trace of a complex elementwise product with a Hermitian X is real", {
  X <- Variable(c(2L, 2L), hermitian = TRUE)
  ## Python [[2, 1+1j], [1-1j, 3]] (row-major) -> R column-major fill.
  k <- matrix(c(2, 1 - 1i, 1 + 1i, 3), 2L, 2L)
  trace_expr <- matrix_trace(k * X)
  expect_true(S7_inherits(trace_expr, Trace))
  expect_true(is_real(trace_expr))
})

## @cvxpy test_atoms.py::TestAtoms::test_trace_elementwise_multiply
test_that("P2: trace of an elementwise product sums only its diagonal", {
  ## Python [[1, 2], [3, 4]] and [[5, 6], [7, 8]] (row-major).
  A <- matrix(c(1, 3, 2, 4), 2L, 2L)
  B <- matrix(c(5, 7, 6, 8), 2L, 2L)
  trace_expr <- matrix_trace(Constant(A) * Constant(B))
  expect_true(S7_inherits(trace_expr, Trace))
  ## sum(diag(A * B)) = 1*5 + 4*8 = 37; pre-fix CVXR returned sum(A * t(B)) = 69.
  expect_equal(as.numeric(value(trace_expr)), sum(diag(A * B)))
})

# ===================================================================
# P6 (#3342, leaf half): complex/imag leaves reject real-valued sign
# attributes and bounds. Pre-fix CVXR accepted them and the solve silently
# ignored the attribute.
# ===================================================================

.p6_sign_attrs <- list(list(nonneg = TRUE), list(nonpos = TRUE),
                       list(pos = TRUE), list(neg = TRUE))

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_complex_sign_attributes_raise
test_that("P6: complex variables cannot use real-valued sign attributes", {
  for (attrs in .p6_sign_attrs) {
    expect_error(do.call(Variable, c(list(c(2L, 2L), complex = TRUE), attrs)),
                 "Cannot combine .* with complex or imaginary attributes",
                 info = names(attrs))
  }
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_imag_sign_attributes_raise
test_that("P6: imaginary variables cannot use real-valued sign attributes", {
  for (attrs in .p6_sign_attrs) {
    expect_error(do.call(Variable, c(list(c(2L, 2L), imag = TRUE), attrs)),
                 "Cannot combine .* with complex or imaginary attributes",
                 info = names(attrs))
  }
})

## @cvxpy test_complex.py::TestComplex::test_complex_bounds_rejected
test_that("P6: complex variables reject bounds", {
  expect_error(Variable(c(2L, 2L), complex = TRUE, bounds = list(0, 1)),
               "Cannot combine bounds with complex or imaginary attributes")
})

## @cvxpy test_complex.py::TestComplex::test_imag_bounds_rejected
test_that("P6: imaginary variables reject bounds", {
  expect_error(Variable(c(2L, 2L), imag = TRUE, bounds = list(0, 1)),
               "Cannot combine bounds with complex or imaginary attributes")
})

## @cvxpy NONE
test_that("P6: hermitian is not a sign attribute; complex Parameters are checked too", {
  ## CVXPY's check covers only nonneg/nonpos/pos/neg and bounds.
  expect_no_error(Variable(c(2L, 2L), complex = TRUE, hermitian = TRUE))
  expect_error(Parameter(c(2L, 2L), complex = TRUE, nonneg = TRUE),
               "Cannot combine .* with complex or imaginary attributes")
})

# ===================================================================
# P7 (#3342, Complex2Real half): a complex Variable declared PSD/NSD keeps
# its cone. Complex2Real splits it into real and imaginary parts, so the
# reduction must re-impose PSD on [[Re, -Im], [Im, Re]]. Pre-fix CVXR (and
# CVXPY 1.9.2) returned "unbounded" for the bounded problems below.
# ===================================================================

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_complex_psd_attribute_trace_problem
test_that("P7: complex PSD attribute keeps the PSD domain (optimum 1)", {
  X <- Variable(c(2L, 2L), complex = TRUE, PSD = TRUE)
  prob <- Problem(Maximize(Re(X[1, 2] + X[2, 1])), list(Re(matrix_trace(X)) == 1))
  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(val, 1.0, tolerance = 1e-5)
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_complex_psd_attribute_matches_explicit_constraint
test_that("P7: complex PSD attribute matches an explicit X %>>% 0", {
  X_attr <- Variable(c(2L, 2L), complex = TRUE, PSD = TRUE)
  prob_attr <- Problem(Maximize(Re(X_attr[1, 2] + X_attr[2, 1])),
                       list(Re(matrix_trace(X_attr)) == 1))
  X_constr <- Variable(c(2L, 2L), complex = TRUE)
  prob_constr <- Problem(Maximize(Re(X_constr[1, 2] + X_constr[2, 1])),
                         list(Re(matrix_trace(X_constr)) == 1, X_constr %>>% 0))
  val_attr <- psolve(prob_attr, solver = "CLARABEL")
  val_constr <- psolve(prob_constr, solver = "CLARABEL")
  expect_equal(status(prob_attr), "optimal")
  expect_equal(status(prob_constr), "optimal")
  expect_equal(val_attr, val_constr, tolerance = 1e-5)
  expect_equal(val_attr, 1.0, tolerance = 1e-5)
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_complex_nsd_attribute_trace_problem
test_that("P7: complex NSD attribute keeps the NSD domain (optimum -1)", {
  X <- Variable(c(2L, 2L), complex = TRUE, NSD = TRUE)
  prob <- Problem(Minimize(Re(X[1, 2] + X[2, 1])), list(Re(matrix_trace(X)) == -1))
  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(val, -1.0, tolerance = 1e-5)
})

# ===================================================================
# P5 (#3398): cumsum of a parameter-only expression stays parameter-affine
# (no auxiliary Variable), so cumsum(p) %*% x is DPP.
# 0.2 (CVXR defect, decision D2): cumsum(X, axis = NULL) of a matrix is a
# column vector c(size, 1) in COLUMN-MAJOR order -- base R's cumsum() -- and
# the atom, its value, and its canonicalization now agree. CVXPY returns a
# 1-D C-order vector; R deviates deliberately (column-major is the R idiom).
# ===================================================================

## @cvxpy test_dpp.py::TestDcp::test_cumsum_of_parameter
test_that("P5: cumsum(p) %*% x is DPP and re-solves with new parameter values", {
  p <- Parameter(3L)
  value(p) <- c(1, 2, 3)
  x <- Variable(3L)

  ## Upstream form (on CVXR this one already passed: the quadratic term routes
  ## the problem through EvalParams).
  problem <- Problem(Minimize(t(cumsum(p)) %*% x + sum_squares(x)))
  expect_true(is_dpp(problem))
  val <- psolve(problem, solver = "CLARABEL")
  ## Analytic optimum: -||cumsum(p)||^2 / 4 = -(1 + 9 + 36) / 4.
  expect_equal(val, -11.5, tolerance = 1e-6)

  constr_prob <- Problem(Minimize(sum_squares(x)), list(x >= cumsum(p)))
  psolve(constr_prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), cumsum(c(1, 2, 3)), tolerance = 1e-5)

  value(p) <- c(2, -1, 1)
  val <- psolve(problem, solver = "CLARABEL")
  expect_equal(val, -sum(cumsum(c(2, -1, 1))^2) / 4, tolerance = 1e-6)

  ## LINEAR-objective form: no EvalParams, so it exercises the DPP
  ## canonicalization of cumsum(p) itself. Pre-fix: "Product of two
  ## non-constant expressions is not DCP."
  q <- Parameter(3L)
  value(q) <- c(1, 2, 3)
  y <- Variable(3L)
  lin <- Problem(Maximize(t(cumsum(q)) %*% y), list(y >= 0, y <= 1))
  expect_true(is_dpp(lin))
  expect_equal(psolve(lin, solver = "CLARABEL"), 1 + 3 + 6, tolerance = 1e-6)
  value(q) <- c(2, -1, 1)          # cumsum = (2, 1, 2), all positive
  expect_equal(psolve(lin, solver = "CLARABEL"), 2 + 1 + 2, tolerance = 1e-6)
})

## @cvxpy test_dpp.py::TestDcp::test_cumsum_of_parameter_axes
test_that("P5: cumsum of a 2-D parameter along each axis matches base R", {
  ## Python np.arange(6.0).reshape(2, 3) (row-major).
  Pv <- matrix(0:5, 2L, 3L, byrow = TRUE)
  P <- Parameter(c(2L, 3L))
  value(P) <- Pv
  ## R axis 2 = CVXPY axis 0 (down columns); R axis 1 = CVXPY axis 1;
  ## NULL is column-major (R deviation, see the section header).
  expected <- list(
    "2"    = apply(Pv, 2L, cumsum),
    "1"    = t(apply(Pv, 1L, cumsum)),
    "NULL" = matrix(cumsum(as.vector(Pv)), ncol = 1L)
  )
  for (nm in names(expected)) {
    axis <- if (nm == "NULL") NULL else as.integer(nm)
    target <- expected[[nm]]
    y <- Variable(dim(target))
    problem <- Problem(Minimize(0), list(y == cumsum_axis(P, axis = axis)))
    expect_true(is_dpp(problem), info = nm)
    psolve(problem, solver = "CLARABEL")
    expect_equal(value(y), target, tolerance = 1e-6, info = nm)
  }
})

## @cvxpy NONE
test_that("0.2: cumsum_axis(X, NULL) of a matrix variable: atom, value and solver agree", {
  Xv <- matrix(0:5, 2L, 3L, byrow = TRUE)
  X <- Variable(c(2L, 3L))
  cs <- cumsum_axis(X, axis = NULL)
  expect_equal(dim(cs), c(6L, 1L))

  ## The solver-side quantity (through an auxiliary scalar) and the
  ## tree-evaluated value must be the same column-major sum, 40 (C-order
  ## would be 35). Pre-fix the solver saw 35 and value() reported 40.
  z <- Variable()
  pr <- Problem(Minimize(0), list(X == Xv, z == sum(cs)))
  psolve(pr, solver = "CLARABEL")
  expect_equal(as.numeric(value(z)), sum(cumsum(as.vector(Xv))), tolerance = 1e-6)
  expect_equal(as.numeric(value(sum(cs))), sum(cumsum(as.vector(Xv))), tolerance = 1e-6)

  ## An equality against the atom's declared shape (pre-fix: "Cannot
  ## broadcast shapes (2, 3) and (6, 3)").
  Y <- Variable(c(6L, 1L))
  pr2 <- Problem(Minimize(0), list(X == Xv, Y == cumsum_axis(X, axis = NULL)))
  psolve(pr2, solver = "CLARABEL")
  expect_equal(value(Y), matrix(cumsum(as.vector(Xv)), ncol = 1L), tolerance = 1e-6)
})

## cumsum.py:75 `axis: int | None = 0`: the default is down the columns, as
## for cumprod and cummax.
## @cvxpy NONE
test_that("cumsum(X) of a matrix expression is column-wise by default", {
  Xv <- matrix(0:5, 2L, 3L, byrow = TRUE)
  X <- Variable(c(2L, 3L))
  expect_equal(dim(cumsum(X)), c(2L, 3L))
  expect_equal(dim(cumsum_axis(X)), c(2L, 3L))
  expect_equal(dim(cumsum(X)), dim(cumprod(X)))
  expect_equal(value(cumsum(Constant(Xv))), apply(Xv, 2L, cumsum))

  Y <- Variable(c(2L, 3L))
  pr <- Problem(Minimize(0), list(X == Xv, Y == cumsum(X)))
  psolve(pr, solver = "CLARABEL")
  expect_equal(value(Y), apply(Xv, 2L, cumsum), tolerance = 1e-6)

  ## A column vector is unchanged.
  v <- Variable(3)
  expect_equal(dim(cumsum(v)), c(3L, 1L))
  expect_equal(as.numeric(value(cumsum(Constant(c(1, 2, 3))))), c(1, 3, 6))
})

# ===================================================================
# Stage 2 -- pre-existing CVXR defects (not 1.9.3 ports); plan Part 0.
# ===================================================================

## 0.1 (decision D2): CVXR's axis contract is {1, 2, NULL}. A negative axis
## was normalized for the SHAPE only, while numeric_value and
## graph_implementation tested the raw value: sum_entries(A, axis = -1) had
## shape (1, 3), solved as column sums, and value() returned row sums.
## @cvxpy NONE
test_that("0.1: negative axes are rejected by every axis atom", {
  A <- matrix(1:6, 2L, 3L)
  X <- Variable(c(2L, 3L))
  for (ax in c(-1L, -2L)) {
    expect_error(sum_entries(A, axis = ax), "Negative axis", info = ax)
    expect_error(sum_entries(X, axis = ax), "Negative axis", info = ax)
    expect_error(p_norm(X, 2, axis = ax), "Negative axis", info = ax)
    expect_error(max_entries(X, axis = ax), "Negative axis", info = ax)
    expect_error(sum_largest(X, 2, axis = ax), "Negative axis", info = ax)
    expect_error(cumsum_axis(X, axis = ax), "Negative axis", info = ax)
  }
  ## The supported axes are unchanged.
  expect_equal(dim(sum_entries(X, axis = 1L)), c(2L, 1L))
  expect_equal(dim(sum_entries(X, axis = 2L)), c(1L, 3L))
})

## 0.5 (decision D4): Complex2Real forwards diag/sparsity to the split real
## and imaginary variables (CVXPY 1.9.3 variable_canon.py
## _split_complex_attributes). Upstream that function is not on the solve
## path, so complex + diag/sparsity is unenforced in CVXPY 1.9.2 and 1.9.3;
## in CVXR c2r_variable_canon IS the solve path, so CVXR now enforces it.
.c2r_split <- function(X) {
  real2imag <- new.env(parent = emptyenv())
  assign(as.character(X@id), next_expr_id(), envir = real2imag)
  c2r_variable_canon(X, NULL, NULL, real2imag)
}
.attr_of <- function(v, nm) .attributes(v)[[nm]]

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_hermitian_variable_canon_returns_real_and_imag
test_that("0.5: a Hermitian variable splits into a symmetric real part and a full imag part", {
  X <- Variable(c(3L, 3L), hermitian = TRUE)
  parts <- .c2r_split(X)
  expect_equal(dim(parts[[1L]]), c(3L, 3L))
  expect_equal(dim(parts[[2L]]), c(3L, 3L))
  expect_true(isTRUE(.attr_of(parts[[1L]], "symmetric")))
  expect_false(isTRUE(.attr_of(parts[[1L]], "complex")))
  expect_false(isTRUE(.attr_of(parts[[1L]], "hermitian")))
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_complex_diag_attribute_preservation
test_that("0.5: diag is forwarded to both parts of a general complex variable", {
  X <- Variable(c(4L, 4L), complex = TRUE, diag = TRUE)
  parts <- .c2r_split(X)
  expect_true(isTRUE(.attr_of(parts[[1L]], "diag")))
  expect_true(isTRUE(.attr_of(parts[[2L]], "diag")))
  expect_false(isTRUE(.attr_of(parts[[1L]], "complex")))
  expect_false(isTRUE(.attr_of(parts[[2L]], "complex")))
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_complex_psd_real_part_preserves_psd
test_that("0.5: PSD is kept on the real part of a complex PSD variable", {
  X <- Variable(c(3L, 3L), complex = TRUE, PSD = TRUE)
  parts <- .c2r_split(X)
  expect_true(isTRUE(.attr_of(parts[[1L]], "PSD")))
  expect_false(isTRUE(.attr_of(parts[[2L]], "PSD")))
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_imag_diag_attribute_preservation
test_that("0.5: diag is forwarded to the split of an imaginary variable", {
  X <- Variable(c(4L, 4L), imag = TRUE, diag = TRUE)
  parts <- .c2r_split(X)
  expect_null(parts[[1L]])
  expect_true(isTRUE(.attr_of(parts[[2L]], "diag")))
  expect_false(isTRUE(.attr_of(parts[[2L]], "imag")))
})

## Python sparsity=([0, 1], [1, 2]) = entries (0,1), (1,2) -> R (1,2), (2,3).
.p05_pattern <- function() cbind(c(1, 2), c(2, 3))

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_complex_sparsity_attribute_preservation
test_that("0.5: sparsity is forwarded to both parts of a general complex variable", {
  X <- Variable(c(3L, 3L), complex = TRUE, sparsity = .p05_pattern())
  parts <- .c2r_split(X)
  expect_equal(.attr_of(parts[[1L]], "sparsity"), .p05_pattern())
  expect_equal(.attr_of(parts[[2L]], "sparsity"), .p05_pattern())
  expect_false(isTRUE(.attr_of(parts[[1L]], "complex")))
  expect_false(isTRUE(.attr_of(parts[[2L]], "complex")))
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_imag_sparsity_attribute_preservation
test_that("0.5: sparsity is forwarded to the split of an imaginary variable", {
  X <- Variable(c(3L, 3L), imag = TRUE, sparsity = .p05_pattern())
  parts <- .c2r_split(X)
  expect_null(parts[[1L]])
  expect_equal(.attr_of(parts[[2L]], "sparsity"), .p05_pattern())
  expect_false(isTRUE(.attr_of(parts[[2L]], "imag")))
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_hermitian_sparsity_attribute_rejected
test_that("0.5: hermitian and sparsity are mutually exclusive", {
  expect_error(Variable(c(3L, 3L), hermitian = TRUE, sparsity = cbind(1, 2)),
               "cannot have more than one")
})

## @cvxpy test_complex.py::TestComplex2RealAttributes::test_hermitian_diag_attribute_rejected
test_that("0.5: hermitian and diag are mutually exclusive", {
  expect_error(Variable(c(3L, 3L), hermitian = TRUE, diag = TRUE),
               "cannot have more than one")
})

## @cvxpy NONE
test_that("0.5: complex/imag diag and sparsity are enforced in the solve", {
  ## Absolute oracle: the attribute-honoring optimum (the matching real
  ## problem gives the same number). Pre-fix CVXR, like CVXPY, gave -9.
  Xd <- Variable(c(3L, 3L), complex = TRUE, diag = TRUE)
  pd <- Problem(Minimize(Re(sum_entries(Xd))), list(Re(Xd) >= -1))
  expect_equal(psolve(pd, solver = "CLARABEL"), -3, tolerance = 1e-6)
  vd <- value(Xd)
  expect_equal(Mod(vd[row(vd) != col(vd)]), rep(0, 6), tolerance = 1e-8)

  Xs <- Variable(c(3L, 3L), complex = TRUE, sparsity = .p05_pattern())
  ps <- Problem(Minimize(Re(sum_entries(Xs))), list(Re(Xs) >= -1))
  expect_equal(psolve(ps, solver = "CLARABEL"), -2, tolerance = 1e-6)

  Xi <- Variable(c(3L, 3L), imag = TRUE, diag = TRUE)
  pi_ <- Problem(Minimize(Im(sum_entries(Xi))), list(Im(Xi) >= -1))
  expect_equal(psolve(pi_, solver = "CLARABEL"), -3, tolerance = 1e-6)
})

## 0.4 (decision D3): cvxr_norm() of a MATRIX (both dimensions > 1) with
## axis = NULL is a matrix norm, as in CVXPY norm.py and base R norm():
## 1 = max column sum, 2 = largest singular value, Inf = max row sum.
## It returned the ENTRYWISE norm (66 / 22.494 / 11 for the 3x4 below where
## the matrix norms are 21 / 22.409 / 38). Column vectors are unchanged.
## 1.9.3 follow-up (user, 2026-10-14; supersedes the 2026-10-06 choice to
## leave norm() alone): the masking norm() and cvxr_norm() agree with base R
## norm() for every type and shape. A row vector is a 1 x n matrix, as in
## base R and CVXPY 2-D: its 1-norm is the max |entry| and its Inf-norm the
## sum of |entries|.
.p04_C <- function() matrix(c(1, -2, 3, 4,
                              -5, 6, -7, 8,
                              9, -1, 2, -3), 3L, 4L, byrow = TRUE)

## @cvxpy test_atoms.py::TestAtoms::test_matrix_norms
test_that("0.4: cvxr_norm matrix 1-, 2-, Inf-, Frobenius and nuclear norms", {
  set.seed(1)
  for (shp in list(c(2L, 2L), c(3L, 2L))) {
    X <- Variable(shp)
    for (p in list(1, 2, Inf, "fro", "nuc")) {
      atom <- cvxr_norm(X, p)
      expect_equal(dim(atom), c(1L, 1L))
      expect_true(is_convex(atom))
      expect_true(is_nonneg(atom))
      Xv <- matrix(rnorm(prod(shp)), shp[1L], shp[2L])
      value(X) <- Xv
      ref <- switch(as.character(p),
        "1" = norm(Xv, "O"), "2" = norm(Xv, "2"), "Inf" = norm(Xv, "I"),
        "fro" = norm(Xv, "F"), "nuc" = sum(svd(Xv)$d))
      expect_equal(as.numeric(value(atom)), ref, tolerance = 1e-8,
                   info = paste(paste(shp, collapse = "x"), p))
    }
  }
})

## @cvxpy NONE
test_that("0.4: matrix norms solve to the matrix-norm value", {
  C <- .p04_C()
  X <- Variable(c(3L, 4L))
  for (p in list(1, 2, Inf)) {
    ref <- switch(as.character(p), "1" = norm(C, "O"), "2" = norm(C, "2"),
                  "Inf" = norm(C, "I"))
    prob <- Problem(Minimize(cvxr_norm(X, p)), list(X == C))
    expect_equal(psolve(prob, solver = "CLARABEL"), ref, tolerance = 1e-5,
                 info = as.character(p))
  }
  ## A numeric string is the same p (cvxr_docs fast-mixing-mc.qmd writes
  ## cvxr_norm(P - 1/n, "2") for the spectral norm).
  Cx <- Constant(C)
  expect_equal(as.numeric(value(cvxr_norm(Cx, "2"))), norm(C, "2"), tolerance = 1e-8)
  expect_equal(as.numeric(value(cvxr_norm(Cx, "1"))), norm(C, "O"))
  expect_equal(as.numeric(value(cvxr_norm(Cx, "Inf"))), norm(C, "I"))
})

## @cvxpy NONE
test_that("0.4: no regression -- column vectors, axis and other p unchanged", {
  C <- .p04_C()
  X <- Constant(C)
  ## Column vectors: entrywise, as before (equal to the matrix norms).
  v <- matrix(c(3, -4, 12), 3L, 1L)
  V <- Constant(v)
  expect_equal(as.numeric(value(cvxr_norm(V, 1))), 19)
  expect_equal(as.numeric(value(cvxr_norm(V, 2))), 13)
  expect_equal(as.numeric(value(cvxr_norm(V, Inf))), 12)
  ## p outside {1, 2, Inf} on a column vector (CVXPY's 1-D case): the
  ## entrywise p-norm, which is also its induced matrix norm.
  expect_equal(as.numeric(value(cvxr_norm(V, 3))), sum(abs(v)^3)^(1/3))
  expect_equal(as.numeric(value(cvxr_norm(Constant(matrix(-5, 1L, 1L)), 3))), 5)
  ## With an axis: unchanged per-slice norms.
  expect_equal(as.numeric(value(cvxr_norm(X, 1, axis = 2L))), colSums(abs(C)))
  expect_equal(as.numeric(value(cvxr_norm(X, 2, axis = 1L))), sqrt(rowSums(C^2)))
  ## cvxr_std() of a matrix (uses cvxr_norm on a vec) is unchanged.
  expect_equal(as.numeric(value(cvxr_std(X))), sqrt(mean((C - mean(C))^2)))
})

## CVXPY norm.py:66-67: on a 2-D argument with no axis, p outside
## {1, 2, Inf, 'fro', 'nuc'} raises "Unsupported matrix norm." -- the induced
## p-norm has no tractable form, and returning the entrywise p-norm under the
## matrix-norm name was silent (user, 2026-10-14). A row vector is a 1 x n
## matrix (its induced p-norm is the dual q-norm), so it raises too; a column
## vector or scalar is CVXPY's 1-D case and keeps the entrywise norm above.
## The entrywise norm stays available as p_norm(x, p). No upstream test
## asserts this message (checked test_atoms.py v1.9.3), hence NONE.
## @cvxpy NONE
test_that("0.4: cvxr_norm() and norm() raise for an unsupported matrix p", {
  C <- .p04_C()
  r <- matrix(c(3, -4, 12), 1L, 3L)
  for (M in list(C, r)) {
    X <- Constant(M)
    for (p in list(3, 1.5, 0.5, "3")) {
      expect_error(cvxr_norm(X, p), "Unsupported matrix norm",
                   info = paste(dim(M), collapse = "x"))
    }
    expect_error(norm(X, "3"), "Unsupported matrix norm")
    expect_error(cvxr_norm(Variable(dim(M)), 3), "Unsupported matrix norm")
    ## The explicit entrywise atom is unaffected.
    expect_equal(as.numeric(value(p_norm(X, 3))), sum(abs(M)^3)^(1/3),
                 tolerance = 1e-8)
  }
})

## Every base::norm() type on a matrix, a column, a row and a scalar.
.p04_shapes <- function() list(
  matrix = .p04_C(),
  column = matrix(c(3, -4, 12), 3L, 1L),
  row    = matrix(c(3, -4, 12), 1L, 3L),
  scalar = matrix(-5, 1L, 1L))

## @cvxpy NONE
test_that("0.4: masking norm() on an Expression equals base::norm() for every type", {
  for (nm in names(.p04_shapes())) {
    M <- .p04_shapes()[[nm]]
    X <- Constant(M)
    for (ty in c("O", "o", "1", "I", "i", "F", "f", "M", "m", "2", "E", "e")) {
      expect_equal(as.numeric(value(norm(X, ty))), base::norm(M, ty),
                   tolerance = 1e-8, info = paste(nm, ty))
    }
    expect_equal(as.numeric(value(norm(X))), base::norm(M, "2"),
                 tolerance = 1e-8, info = paste(nm, "default"))
  }
})

## @cvxpy test_atoms.py::TestAtoms::test_matrix_norms
test_that("0.4: cvxr_norm() p = 1, 2, Inf, 'fro' equals base::norm() on every shape", {
  ## CVXPY norm.py:54-65 treats every 2-D argument as a matrix; a row
  ## vector's 1-norm is its max |entry| (CVXPY 2-D (1, 3): 12, not 19).
  for (nm in names(.p04_shapes())) {
    M <- .p04_shapes()[[nm]]
    X <- Constant(M)
    for (p in list(1, 2, Inf, "fro")) {
      ty <- switch(as.character(p), "1" = "O", "2" = "2", "Inf" = "I", fro = "F")
      expect_equal(as.numeric(value(cvxr_norm(X, p))), base::norm(M, ty),
                   tolerance = 1e-8, info = paste(nm, p))
    }
  }
})

## @cvxpy NONE
test_that("0.4: norm() type changes the solve where the matrix norm differs", {
  C <- .p04_C()
  X <- Variable(c(3L, 4L))
  for (ty in c("O", "I", "2", "M")) {
    prob <- Problem(Minimize(norm(X, ty)), list(X == C))
    expect_equal(psolve(prob, solver = "CLARABEL"), base::norm(C, ty),
                 tolerance = 1e-5, info = ty)
  }
  r <- matrix(c(3, -4, 12), 1L, 3L)
  y <- Variable(c(1L, 3L))
  for (p in list(1, Inf)) {
    prob <- Problem(Minimize(cvxr_norm(y, p)), list(y == r))
    expect_equal(psolve(prob, solver = "CLARABEL"),
                 base::norm(r, if (p == 1) "O" else "I"), tolerance = 1e-5,
                 info = as.character(p))
  }
})

# ===================================================================
# P8 (#3493): geo_mean over any 2-D shape, with axis and keepdims.
# GeoMean is an AxisAtom; zero-weight entries stay in the argument and
# are dropped through a keep mask in the aligned view. Pre-fix CVXR
# aborted "x must be a row or column vector" on a matrix and had no
# axis argument. R axis 2L = CVXPY axis 0 (per column), 1L = CVXPY 1.
# Solves are checked against ANALYTIC optima, never against cvxpy 1.9.3,
# which gets the U1 problem wrong under its C++ backend (plan §4 U1).
# ===================================================================

## numpy arange(1, 7).reshape(2, 3) is row-major.
.p8_M <- function() matrix(1:6, 2L, 3L, byrow = TRUE)
.p8_gm <- function(v) exp(mean(log(v)))

## @cvxpy test_atoms.py::TestAtoms::test_geo_mean_shapes
test_that("P8: geo_mean of a matrix is the mean over all entries, column-major", {
  M <- .p8_M()
  for (approx in c(TRUE, FALSE)) {
    atom <- geo_mean(Constant(M), approx = approx)
    expect_equal(atom@shape, c(1L, 1L))
    expect_equal(as.numeric(atom@w), rep(1 / 6, 6))
    expect_equal(as.numeric(value(atom)), prod(M)^(1 / 6))
    ## Weights line up with the column-major flattening.
    weighted <- geo_mean(Constant(M), c(2, 0, 0, 0, 0, 1), approx = approx)
    flat <- as.vector(M)
    expect_equal(as.numeric(value(weighted)), (flat[1]^2 * flat[6])^(1 / 3))
  }
  ## The (2, 3, 4) part of the upstream test is N-D: N/A for R.
})

## @cvxpy test_atoms.py::test_geo_mean_axis
test_that("P8: every axis reduces, and the shape follows cvxr_mean()", {
  M <- .p8_M()
  C <- Constant(M)
  for (keepdims in c(FALSE, TRUE)) {
    for (axis in list(NULL, 1L, 2L)) {
      got <- geo_mean(C, axis = axis, keepdims = keepdims)
      expect_equal(got@shape, cvxr_mean(C, axis = axis, keepdims = keepdims)@shape)
      expected <- if (is.null(axis)) .p8_gm(M) else apply(M, axis, .p8_gm)
      expect_equal(as.numeric(value(got)), as.numeric(expected))
      expect_equal(dim(value(got)), got@shape)
    }
  }
  ## The upstream negative axes and axis tuples are N-D / not CVXR's axis
  ## contract (axis is 1L, 2L or NULL; negative axes are refused, 0.1).
})

## @cvxpy test_atoms.py::test_geo_mean_axis_out_of_bounds
test_that("P8: an axis the expression does not have is refused", {
  C <- Constant(.p8_M())
  expect_error(geo_mean(C, axis = 3L), "out of bounds")
  expect_error(geo_mean(C, axis = -1L), "Negative axis")
})

## @cvxpy test_atoms.py::test_geo_mean_axis_with_weights
test_that("P8: weights apply to each slice along the reduced axis", {
  M <- .p8_M()
  got <- geo_mean(Constant(M), c(3, 1), axis = 2L)
  expect_equal(as.numeric(value(got)), (M[1, ]^3 * M[2, ])^(1 / 4))
  ## Row-wise: one weight per column entry.
  got1 <- geo_mean(Constant(M), c(1, 0, 2), axis = 1L)
  expect_equal(as.numeric(value(got1)), (M[, 1] * M[, 3]^2)^(1 / 3))
})

## @cvxpy test_atoms.py::TestAtoms::test_geo_mean_shares_one_set_of_constraints
test_that("P8: every output is covered by one set of constraints, not one each", {
  counts <- vapply(c(2L, 20L, 200L), function(cols) {
    x <- Variable(c(3L, cols), nonneg = TRUE)
    expr <- geo_mean(x, axis = 2L)
    length(CVXR:::geo_mean_approx_canon(expr, expr@args)[[2L]])
  }, integer(1))
  expect_length(unique(counts), 1L)
})

## @cvxpy test_atoms.py::TestAtoms::test_geo_mean_row_vector
test_that("P8: column, row and matrix arguments all reach the solver", {
  require_solver("SCS")
  for (shape in list(c(4L, 1L), c(1L, 4L), c(2L, 2L))) {
    x <- Variable(shape, nonneg = TRUE)
    prob <- Problem(Maximize(geo_mean(x)), list(sum_entries(x) <= 4))
    expect_equal(psolve(prob, solver = "SCS", eps = 1e-6), 1.0,
                 tolerance = 1e-4, info = paste(shape, collapse = "x"))
  }
})

## Plan §4 U1: max sum(geo_mean(X, axis)) s.t. per-slice budgets, X >= 0.
## AM-GM: an even split per slice, so the optimum is sum(budget) / n = 6.
## cvxpy 1.9.3 returns 2.89 here under its C++ backend; CVXR shares that
## backend, so this is the guard against porting the defect. On SCS the
## exact form goes through the PowConeND -> PowCone3D conversion.
## @cvxpy NONE
test_that("P8/U1: per-column geo_mean reaches the analytic optimum 6", {
  b <- c(2, 4, 6)
  for (solver in c("CLARABEL", "SCS")) {
    require_solver(solver)
    for (approx in c(TRUE, FALSE)) {
      info <- paste(solver, if (approx) "approx" else "exact")
      X <- Variable(c(2L, 3L), nonneg = TRUE)
      prob <- Problem(Maximize(sum_entries(geo_mean(X, axis = 2L, approx = approx))),
                      list(sum_entries(X, axis = 2L) <= matrix(b, 1L)))
      expect_equal(psolve(prob, solver = solver), 6, tolerance = 1e-4, info = info)
      expect_equal(value(X), rbind(b / 2, b / 2), tolerance = 1e-3, info = info)
      ## Row-wise on the transposed problem: same optimum.
      Y <- Variable(c(3L, 2L), nonneg = TRUE)
      prob1 <- Problem(Maximize(sum_entries(geo_mean(Y, axis = 1L, approx = approx))),
                       list(sum_entries(Y, axis = 1L) <= matrix(b, ncol = 1L)))
      expect_equal(psolve(prob1, solver = solver), 6, tolerance = 1e-4, info = info)
      expect_equal(value(Y), cbind(b / 2, b / 2), tolerance = 1e-3, info = info)
    }
  }
})

## The aligned view is checked by identity through the C++ lowering
## before any cone is built on it (plan §4 U1 mitigation (ii)).
## @cvxpy NONE
test_that("P8/U1: the aligned view survives the C++ lowering unchanged", {
  require_solver("CLARABEL")
  M <- .p8_M()
  cases <- list(
    list(axis = NULL, p = NULL,          want = matrix(as.vector(M), ncol = 1L)),
    list(axis = 2L,   p = NULL,          want = M),
    list(axis = 1L,   p = NULL,          want = t(M)),
    list(axis = 2L,   p = c(0, 1),       want = M[2, , drop = FALSE]),
    list(axis = 1L,   p = c(1, 0, 2),    want = t(M)[c(1, 3), , drop = FALSE])
  )
  for (cs in cases) {
    X <- Variable(c(2L, 3L))
    gm <- geo_mean(X, p = cs$p, axis = cs$axis)
    A <- CVXR:::.gm_aligned_arg(gm)
    Y <- Variable(A@shape)
    prob <- Problem(Minimize(0), list(X == M, Y == A))
    psolve(prob, solver = "CLARABEL")
    expect_equal(value(Y), cs$want, tolerance = 1e-6,
                 info = paste(c(cs$axis, "|", cs$p), collapse = " "))
  }
})

## R rebuilds atoms through the constructor from get_data() (expr_copy), so
## get_data() must carry the full-length weights even after zero-weight
## entries are dropped from the canonicalization.
## @cvxpy NONE
test_that("P8: copy round-trips axis, keepdims and zero weights", {
  X <- Variable(c(2L, 3L))
  for (approx in c(TRUE, FALSE)) {
    gm <- geo_mean(X, p = c(3, 0), axis = 2L, keepdims = TRUE, approx = approx)
    cp <- CVXR:::expr_copy(gm)
    expect_identical(class(cp), class(gm))
    expect_equal(cp@shape, gm@shape)
    expect_equal(get_data(cp), get_data(gm))
    expect_equal(as.numeric(cp@w), as.numeric(gm@w))
  }
})

## @cvxpy NONE
test_that("P8: the gradient is per slice and skips zero-weight entries", {
  M <- .p8_M()
  X <- Variable(c(2L, 3L)); value(X) <- M
  g <- grad(geo_mean(X, axis = 2L))[[1L]]
  expect_equal(dim(g), c(6L, 3L))
  ## Column j of the Jacobian: (1/2) gm(M[, j]) / M[, j] on rows of column j.
  want <- matrix(0, 6L, 3L)
  for (j in 1:3) want[(2L * j - 1L):(2L * j), j] <- 0.5 * .p8_gm(M[, j]) / M[, j]
  expect_equal(as.matrix(g), want, tolerance = 1e-10)
  ## A zero weight gives a zero derivative, not NaN.
  gw <- grad(geo_mean(X, p = c(1, 0), axis = 2L))[[1L]]
  wantw <- matrix(0, 6L, 3L)
  for (j in 1:3) wantw[2L * j - 1L, j] <- 1
  expect_equal(as.matrix(gw), wantw, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("P8: DGP geo_mean with an axis is the per-slice weighted log average", {
  M <- .p8_M()
  X <- Variable(c(2L, 3L), pos = TRUE)
  prob <- Problem(Minimize(sum_entries(geo_mean(X, axis = 2L))), list(X >= M))
  expect_true(is_dgp(prob))
  expect_equal(psolve(prob, gp = TRUE, solver = "CLARABEL"),
               sum(apply(M, 2L, .p8_gm)), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("P8: DNLP geo_mean with an axis solves per column", {
  b <- c(2, 4, 6)
  X <- Variable(c(2L, 3L), nonneg = TRUE); value(X) <- matrix(1, 2L, 3L)
  prob <- Problem(Maximize(sum_entries(geo_mean(X, axis = 2L))),
                  list(sum_entries(X, axis = 2L) <= matrix(b, 1L)))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(value(X), rbind(b / 2, b / 2), tolerance = 1e-4)
})
