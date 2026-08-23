## Atom-level fixes ported from CVXPY 1.9.2.
##
## Each block is an R port of an upstream 1.9.2 regression test, and each one
## FAILED against CVXR 1.9.1.9014 before the corresponding port.  Measured
## before/after values are recorded beside the assertions, so a future
## regression says which upstream fix was lost.
##
## Sources: atoms/norm1.py, atoms/elementwise/log1p.py, atoms/elementwise/ceil.py,
## atoms/elementwise/hyperbolic.py, atoms/affine/partial_trace.py,
## atoms/affine/partial_transpose.py, atoms/elementwise/elementwise.py,
## atoms/elementwise/entr.py.

# ===================================================================
# norm1 is log-log convex  (atoms/norm1.py:58-66)
# ===================================================================

## @cvxpy test_dgp.py::TestDgp::test_norm1_dgp
test_that("norm1 is log-log convex, so norm(x, 1) is DGP", {
  x <- Variable(2, pos = TRUE)
  ## Before the port: FALSE -- the DGP analyzer inherited Atom's default and
  ## rejected the problem with "Problem is not DGP compliant".
  expect_true(is_log_log_convex(norm1(x)))
  expect_false(is_log_log_concave(norm1(x)))

  ## minimize x1 + x2 s.t. x1 * x2 >= 4 has optimum 4 at x = (2, 2).
  prob <- Problem(Minimize(norm1(x)), list(x[1] * x[2] >= 4))
  expect_true(is_dgp(prob))
  val <- psolve(prob, gp = TRUE, solver = "CLARABEL")
  expect_equal(val, 4.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(2.0, 2.0), tolerance = 1e-4)
})

## @cvxpy test_dgp.py::TestDgp::test_quad_form_dgp
test_that("QuadForm is canonicalized under DGP, not copied verbatim", {
  ## Upstream keyed its dgp2dcp registry by the quad_form FUNCTION rather than
  ## the QuadForm CLASS, so the lookup never matched.  CVXR dispatches on S7
  ## classes and never had the defect; this pins the behavior.
  x <- Variable(2, pos = TRUE)
  P <- matrix(c(2.0, 1.0, 1.0, 3.0), 2, 2)
  prob <- Problem(Minimize(quad_form(x, P)), list(x[1] * x[2] >= 4))
  expect_true(is_dgp(prob))
  val <- psolve(prob, gp = TRUE, solver = "CLARABEL")
  ## With x2 = 4/x1 the objective is 2*x1^2 + 8 + 48/x1^2, minimized at
  ## x1^2 = sqrt(24), optimum 8 + 8*sqrt(6).
  expect_equal(val, 8 + 8 * sqrt(6), tolerance = 1e-4)
  xv <- as.numeric(value(x))
  expect_equal(as.numeric(t(xv) %*% P %*% xv), val, tolerance = 1e-4)
})

# ===================================================================
# log1p bounds are those of log(1 + x)  (elementwise/log1p.py:43-47)
# ===================================================================

## @cvxpy test_bounds.py::TestBoundsPropagation::test_log1p_bounds
test_that("log1p bounds are log(1 + x)'s, not log(x)'s", {
  ## Before the port Log1p inherited Log's bounds_from_args, which forgot the
  ## +1 shift.  Measured before: [1,2] -> (0, log 2); [-1,0] -> (-Inf, +Inf).
  cases <- list(
    list(bounds = c(1, 2),      lb = log(2),  ub = log(3)),
    list(bounds = c(-0.5, 0.5), lb = log(0.5), ub = log(1.5)),
    ## lb on the domain boundary x > -1 gives -Inf; the ub must stay finite.
    list(bounds = c(-1, 0),     lb = -Inf,    ub = 0)
  )
  for (cs in cases) {
    x <- Variable(3, bounds = list(cs$bounds[1], cs$bounds[2]))
    b <- get_bounds(log1p_expr(x))
    expect_equal(as.numeric(b[[1L]]), rep(cs$lb, 3), tolerance = 1e-8,
                 info = paste("lb for bounds", cs$bounds[1], cs$bounds[2]))
    expect_equal(as.numeric(b[[2L]]), rep(cs$ub, 3), tolerance = 1e-8,
                 info = paste("ub for bounds", cs$bounds[1], cs$bounds[2]))
  }
})

## @cvxpy test_bounds.py::TestBoundsPropagation::test_monotonic_atom_bounds
test_that("exp/log/sqrt bounds are unaffected by the log1p fix", {
  b <- get_bounds(exp(Variable(3, bounds = list(0, 1))))
  expect_equal(as.numeric(b[[1L]]), rep(1, 3), tolerance = 1e-8)
  expect_equal(as.numeric(b[[2L]]), rep(exp(1), 3), tolerance = 1e-8)

  b <- get_bounds(log(Variable(3, bounds = list(1, exp(1)))))
  expect_equal(as.numeric(b[[1L]]), rep(0, 3), tolerance = 1e-8)
  expect_equal(as.numeric(b[[2L]]), rep(1, 3), tolerance = 1e-8)

  b <- get_bounds(sqrt(Variable(3, bounds = list(1, 4))))
  expect_equal(as.numeric(b[[1L]]), rep(1, 3), tolerance = 1e-8)
  expect_equal(as.numeric(b[[2L]]), rep(2, 3), tolerance = 1e-8)
})

# ===================================================================
# floor snaps values within solver noise of an integer  (ceil.py:109-112)
# ===================================================================

## @cvxpy test_dqcp.py::TestDqcp::test_ceil_floor_near_integer_tolerance
test_that("floor/ceil treat near-integer values as that integer", {
  x <- Variable()
  ## Before the port floor(3 - 1e-10) was 2: ordinary solver noise moved the
  ## result by a whole integer.  ceil already rounded first.
  for (noisy in c(3 - 1e-10, 3 + 1e-10)) {
    value(x) <- noisy
    expect_equal(as.numeric(value(floor_expr(x))), 3.0)
    expect_equal(as.numeric(value(ceil_expr(x))), 3.0)
  }
  ## Exact integers and clearly-non-integer values are unaffected.
  value(x) <- 4.0
  expect_equal(as.numeric(value(floor_expr(x))), 4.0)
  expect_equal(as.numeric(value(ceil_expr(x))), 4.0)
  value(x) <- 2.5
  expect_equal(as.numeric(value(floor_expr(x))), 2.0)
  expect_equal(as.numeric(value(ceil_expr(x))), 3.0)
})

# ===================================================================
# atanh domain is the closed interval  (elementwise/hyperbolic.py:199)
# ===================================================================

## @cvxpy test_nonlinear_atoms.py::test_hyperbolic_atoms_metadata_numeric_and_unimplemented_grad
test_that("atanh reports a two-constraint domain instead of aborting", {
  x <- Variable(2)
  expr <- Atanh(x)
  ## Through 1.9.1 upstream used STRICT inequalities, which raise, so CVXR
  ## aborted here too.  1.9.2 relaxed both bounds to non-strict.
  dom <- atom_domain(expr)
  expect_length(dom, 2L)
  for (cn in dom) expect_true(S7_inherits(cn, Inequality))

  ## The rest of the upstream metadata block, which never depended on the fix.
  val <- c(0.2, -0.3)
  expect_equal(as.numeric(numeric_value(expr, list(val))), atanh(val),
               tolerance = 1e-12)
  sg <- sign_from_args(expr)
  expect_false(sg$is_nonneg)
  expect_false(sg$is_nonpos)
  expect_false(is_atom_convex(expr))
  expect_false(is_atom_concave(expr))
  expect_true(is_atom_smooth(expr))
  expect_true(is_incr(expr, 1L))
  expect_false(is_decr(expr, 1L))
  expect_error(.grad(expr, list(val)), "Gradient not implemented")
})

# ===================================================================
# partial_trace / partial_transpose propagate DCP attributes
# (affine/partial_trace.py:87-93, affine/partial_transpose.py:88-94)
# ===================================================================

## @cvxpy test_atoms.py::TestAtoms::test_partial_trace_dcp_attributes
test_that("partial_trace propagates PSD/symmetric/Hermitian", {
  ## Before the port EVERY attribute came back FALSE for EVERY input.
  ## PSD input -> PSD output (partial trace preserves PSD).
  X_psd <- Variable(c(4, 4), PSD = TRUE)
  pt_psd <- partial_trace(X_psd, c(2, 2))
  expect_true(is_psd(pt_psd))
  expect_true(is_hermitian(pt_psd))
  expect_true(is_symmetric(pt_psd))

  ## Symmetric input -> symmetric output.
  pt_sym <- partial_trace(Variable(c(4, 4), symmetric = TRUE), c(2, 2))
  expect_true(is_symmetric(pt_sym))
  expect_true(is_hermitian(pt_sym))

  ## Hermitian input -> Hermitian, but not symmetric.
  pt_herm <- partial_trace(Variable(c(4, 4), hermitian = TRUE), c(2, 2))
  expect_true(is_hermitian(pt_herm))
  expect_false(is_symmetric(pt_herm))

  ## Plain input -> no special attributes.
  pt_plain <- partial_trace(Variable(c(4, 4)), c(2, 2))
  expect_false(is_psd(pt_plain))
  expect_false(is_hermitian(pt_plain))

  ## Complex PSD input -> PSD and Hermitian, not symmetric.
  pt_psd_c <- partial_trace(Variable(c(4, 4), PSD = TRUE, complex = TRUE), c(2, 2))
  expect_true(is_psd(pt_psd_c))
  expect_true(is_hermitian(pt_psd_c))
  expect_false(is_symmetric(pt_psd_c))
})

## @cvxpy test_atoms.py::TestAtoms::test_partial_transpose_dcp_attributes
test_that("partial_transpose propagates symmetry but NOT PSD", {
  pp_sym <- partial_transpose(Variable(c(4, 4), symmetric = TRUE), c(2, 2))
  expect_true(is_symmetric(pp_sym))
  expect_true(is_hermitian(pp_sym))

  pp_herm <- partial_transpose(Variable(c(4, 4), hermitian = TRUE), c(2, 2))
  expect_true(is_hermitian(pp_herm))

  ## Partial transpose does NOT preserve PSD -- that is what makes it an
  ## entanglement witness (the PPT criterion).  It must stay Hermitian.
  pp_psd <- partial_transpose(Variable(c(4, 4), PSD = TRUE), c(2, 2))
  expect_false(is_psd(pp_psd))
  expect_true(is_hermitian(pp_psd))

  pp_plain <- partial_transpose(Variable(c(4, 4)), c(2, 2))
  expect_false(is_hermitian(pp_plain))
  expect_false(is_symmetric(pp_plain))

  pp_herm_c <- partial_transpose(
    Variable(c(4, 4), hermitian = TRUE, complex = TRUE), c(2, 2))
  expect_true(is_hermitian(pp_herm_c))
})

# ===================================================================
# Already correct in CVXR -- pinned so a future edit cannot regress them
# ===================================================================

## @cvxpy test_atoms.py::TestAtoms::test_elementwise_is_symmetric
test_that("is_symmetric does not crash for scalar or vector elementwise atoms", {
  ## Upstream raised IndexError on a 1-D shape.  An R shape is always length
  ## 2, so the subscript is always in range; these are the same assertions.
  expect_true(is_symmetric(abs(Variable())))
  expect_false(is_symmetric(abs(Variable(3))))
  expect_true(is_symmetric(abs(Variable(c(2, 2), symmetric = TRUE))))
  expect_false(is_symmetric(abs(Variable(c(2, 2)))))
  expect_false(is_symmetric(abs(Variable(c(2, 3)))))
})

## @cvxpy test_atoms.py::TestAtoms::test_entr
test_that("entr agrees on dense and sparse constants", {
  dense <- matrix(c(0.5, 0, 0, 2), 2, 2)
  expected <- matrix(c(0.5 * log(2), 0, 0, -2 * log(2)), 2, 2)
  expect_equal(as.numeric(value(entr(Constant(dense)))),
               as.numeric(expected), tolerance = 1e-10)
  sparse <- Matrix::Matrix(dense, sparse = TRUE)
  expect_equal(as.numeric(value(entr(Constant(sparse)))),
               as.numeric(expected), tolerance = 1e-10)
})
