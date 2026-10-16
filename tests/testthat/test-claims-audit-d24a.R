## Claims audit decision 24a (notes/audit/claims_audit_item1.md): thirteen small
## exact ports in atoms, constraints and Problem, CVXR 1.9.2.9158. Expected values
## were produced by cvxpy 1.9.3 in a pinned `uv run` process (scratch d24a_cvxpy.py);
## every test here was seen failing on 1.9.2.9157 first.

# -- C0187: power bounds use the ORIGINAL exponent -----------------------------

## @cvxpy NONE
test_that("power bounds: Parameter exponent uses its current value, unset falls back to 1", {
  x <- Variable(2, bounds = list(c(1, 2), c(3, 4)))
  p <- Parameter(nonneg = TRUE)
  value(p) <- 2
  b <- get_bounds(power(x, p))
  expect_equal(as.numeric(b[[1L]]), c(1, 4))
  expect_equal(as.numeric(b[[2L]]), c(9, 16))

  ## power.py:209: an UNSET Parameter exponent defaults to 1.0.
  q <- Parameter(nonneg = TRUE)
  b0 <- get_bounds(power(x, q))
  expect_equal(as.numeric(b0[[1L]]), c(1, 2))
  expect_equal(as.numeric(b0[[2L]]), c(3, 4))
})

## @cvxpy NONE
test_that("power bounds: a fractional constant exponent is the original p, not its rational approximation", {
  x <- Variable(2, bounds = list(c(1, 2), c(3, 4)))
  e <- power(x, pi / 4)
  ## p_used (355/452) drives monotonicity, quadratic/PWL and the gradient; curvature and bounds follow float(self.p.value).
  expect_false(isTRUE(all.equal(as.numeric(e@p_used), pi / 4, tolerance = 1e-12)))
  b <- get_bounds(e)
  expect_equal(as.numeric(b[[1L]]), c(1, 2)^(pi / 4), tolerance = 1e-12)
  expect_equal(as.numeric(b[[2L]]), c(3, 4)^(pi / 4), tolerance = 1e-12)
})

# -- C0216 / C0261: log_det and tr_inv gradients without a conditioning test ---

## @cvxpy test_grad.py::TestSingleVarAtomGradients::test_single_var_atom[log_det]
test_that("log_det gradient: ill-conditioned but positive definite X still yields X^{-T}", {
  X <- Variable(c(2, 2))
  value(X) <- matrix(c(1, 0, 0, 1e-17), 2, 2)
  g <- grad(log_det(X))
  expect_equal(as.numeric(g[[1L]]), c(1, 0, 0, 1e17), tolerance = 1e-9)
})

## @cvxpy test_grad.py::TestSingleVarAtomGradients::test_single_var_atom[tr_inv]
test_that("tr_inv gradient: ill-conditioned but positive definite X still yields -(X^{-T})^2", {
  X <- Variable(c(2, 2))
  value(X) <- matrix(c(1, 0, 0, 1e-17), 2, 2)
  g <- grad(tr_inv(X))
  expect_equal(as.numeric(g[[1L]]), c(-1, 0, 0, -1e34), tolerance = 1e-9)
})

# -- C0216 / C0261: domain() appends the argument domains ----------------------

## @cvxpy test_atoms.py::TestAtoms::test_log_det
test_that("log_det and tr_inv domain(): the atom's PSD constraint plus the argument's domain", {
  y <- Variable(2, nonneg = TRUE)
  Dm <- diag(sqrt(y))
  arg_dom <- domain(Dm)
  expect_gt(length(arg_dom), 0L)
  for (f in list(log_det, tr_inv)) {
    dom <- domain(f(Dm))
    expect_length(dom, 1L + length(arg_dom))
    expect_true(S7_inherits(dom[[1L]], PSD))
  }
})

# -- C0237 / C0244: perspective validation and gradient arity ------------------

## @cvxpy NONE
test_that("perspective rejects complex arguments, as Atom.validate_arguments", {
  z <- Variable(complex = TRUE)
  s <- Variable(nonneg = TRUE)
  expect_error(perspective(Re(z) + Im(z), s), "cannot be complex")
})

## @cvxpy NONE
test_that("perspective gradient: one NULL per argument when f has several variables", {
  a <- Variable(); b <- Variable(); s <- Variable(nonneg = TRUE)
  value(a) <- 1; value(b) <- 2; value(s) <- 1
  p <- perspective(a + b, s)
  expect_length(.args(p), 3L)
  g <- grad(p)  ## indexed out of bounds on 1.9.2.9157
  expect_true(all(vapply(g, is.null, logical(1L))))
})

# -- C0266: value(constraint, tolerance) ---------------------------------------

## @cvxpy NONE
test_that("value(constraint) exposes the absolute tolerance on the violation", {
  v <- Variable()
  value(v) <- 1e-6
  con <- v <= 0
  expect_false(value(con))
  expect_true(value(con, tolerance = 1e-5))
})

# -- C0272 / C0526 / C0963: is_dpp context validation --------------------------

## @cvxpy NONE
test_that("is_dpp aborts on an unsupported context; dcp/dgp are case-insensitive", {
  w <- Variable()
  con <- w >= 0
  prob <- Problem(Minimize(w), list(con))
  t <- Variable()
  inner <- Problem(Minimize(t), list(t >= w))
  pp <- partial_optimize(inner, opt_vars = list(t))

  expect_error(is_dpp(con, "foo"), "Unsupported context")
  expect_error(is_dpp(prob, "foo"), "Unsupported context")
  expect_error(is_dpp(pp, "foo"), "Unsupported context")
  expect_error(.problem_is_dpp(prob, "dcp", quad_form_dpp = "bad"), "Unsupported")

  expect_true(is_dpp(con, "DCP"))
  expect_true(is_dpp(prob, "Dcp"))
  expect_true(is_dpp(pp, "dcp"))
  expect_true(.problem_is_dpp(prob, "dcp", quad_form_dpp = "qp"))
})

# -- C0308: constraint size is that of lhs - rhs -------------------------------

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("Inequality and Equality size follow the constrained expression, whichever side is scalar", {
  u <- Variable(3)
  expect_equal(constr_size(0 <= u), 3)
  expect_equal(constr_size(u >= 0), 3)
  expect_equal(constr_size(0 == u), 3)
  expect_equal(constr_size(u == 0), 3)
  expect_equal(constr_size(Variable() <= 1), 1)
})

# -- C0359 / C0369: PSD and SvecPSD dual_cone argument checks ------------------

## @cvxpy NONE
test_that("dual_cone(PSD) and dual_cone(SvecPSD) check the explicit arguments", {
  M <- Variable(c(2, 2), symmetric = TRUE)
  con <- M %>>% 0
  N3 <- Variable(c(3, 3), symmetric = TRUE)
  expect_true(S7_inherits(dual_cone(con), PSD))
  expect_true(S7_inherits(dual_cone(con, M), PSD))
  expect_error(dual_cone(con, N3), "shapes of the")
  expect_error(dual_cone(con, M, M), "expects 1 argument")

  sv <- SvecPSD(Variable(3), n = 2L)
  expect_true(S7_inherits(dual_cone(sv, Variable(3)), SvecPSD))
  expect_error(dual_cone(sv, Variable(3), Variable(3)), "expects 1 argument")
  expect_error(dual_cone(sv, Variable(6)), "shapes of the")
})
