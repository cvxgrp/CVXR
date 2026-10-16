## Claims audit decision 29 (notes/audit/claims_audit_item1.md): the DPP rules move
## from a recursive `is_dpp` over the arguments into the curvature hooks under the
## DPP scope, as in CVXPY (atom.py:300-308, kron.py:54-69, conv.py:72-86,
## problem.py:400), CVXR 1.9.2.9164. Expected values were produced by cvxpy
## 1.9.3 in a pinned `uv run` process (scratch d29_oracle.py). The first three
## blocks were seen failing on 1.9.2.9162; the pins in the last two blocks assert
## existing behavior and cannot be seen failing.

# -- Indicator: DPP at every level above a direct call -------------------------

## @cvxpy transforms/indicator.py::indicator.is_dpp atoms/atom.py::Atom.is_dpp
test_that("indicator: is_dpp is FALSE on the transform itself, TRUE once it is an argument", {
  s <- Variable()
  ps <- Parameter(value = 2)
  ind <- indicator(list(ps * s <= 1))
  expect_false(is_dpp(ind))               # indicator.py:109-112
  expect_true(is_dpp(ind + s))            # Indicator.is_convex is unconditional
  expect_true(is_dpp(Minimize(ind + s)))
  expect_true(is_dpp(Problem(Minimize(ind - s), list(s >= -5))))
})

## @cvxpy transforms/indicator.py::indicator.is_dpp problems/problem.py::Problem.is_dpp
test_that("indicator: a parametric indicator problem re-solves through the cached parameter program", {
  require_solver("CLARABEL")
  s <- Variable()
  ps <- Parameter(value = 2)
  prob <- Problem(Minimize(indicator(list(ps * s <= 1)) - s), list(s >= -5))
  expect_true(is_dpp(prob))
  expect_no_warning(v1 <- psolve(prob, solver = "CLARABEL"))
  expect_equal(v1, -0.5, tolerance = 1e-6)
  pp <- prob@.cache$param_prog
  expect_false(is.null(pp))
  value(ps) <- 4
  expect_equal(psolve(prob, solver = "CLARABEL"), -0.25, tolerance = 1e-6)
  expect_identical(prob@.cache$param_prog, pp)
})

# -- Variable bounds: checked once, at the problem level ----------------------

## @cvxpy expressions/variable.py::Variable.is_dpp problems/problem.py::Problem.is_dpp
test_that("parametric bounds: a non-affine bound is caught by the Problem, not by every enclosing expression", {
  p1 <- Parameter(value = 1)
  p2 <- Parameter(value = 1)
  bad <- Variable(bounds = list(p1 * p2, NULL))
  expect_false(is_dpp(bad))               # variable.py:96-103, direct call
  expect_true(is_dpp(bad + 1))            # atom.py:300-308: one scoped is_dcp pass
  expect_true(is_dpp(bad <= 10))          # constraint.py:188-194
  expect_true(is_dpp(Minimize(bad)))      # objective.py:173-180
  expect_false(is_dpp(Problem(Minimize(bad), list(bad <= 10))))  # problem.py:400
  good <- Variable(bounds = list(-p1, p1))
  expect_true(is_dpp(Problem(Minimize(good), list(good <= 10))))
})

## Decision 29 side finding 1, CVXR 1.9.2.9170: the gp = TRUE chain asks the same
## question, problem.is_dpp('dgp') (solving_chain.py:193,201). Expected values from
## cvxpy 1.9.3 pinned (notes/audit/repro/30_gp_branch_bounds.py); seen failing on
## 1.9.2.9169, whose chain skipped the bounds.
## @cvxpy expressions/variable.py::Variable.is_dpp problems/problem.py::Problem.is_dpp
test_that("parametric bounds: the gp = TRUE chain checks them, as Problem.is_dpp('dgp') does", {
  require_solver("CLARABEL")
  p <- Parameter(pos = TRUE, value = 2)
  q <- Parameter(pos = TRUE, value = 3)
  bad <- Variable(pos = TRUE, bounds = list(p + q, NULL))  # log-log convex, not affine
  prob <- Problem(Minimize(bad))
  expect_false(is_dpp(prob, "dgp"))                         # variable.py:105-113
  expect_warning(ch <- construct_solving_chain(prob, solver = "CLARABEL", gp = TRUE),
                 class = "DPPWarning")
  expect_true(.s7_is(ch@reductions[[1L]], EvalParams))      # solving_chain.py:207
  expect_error(psolve(prob, gp = TRUE, solver = "CLARABEL", enforce_dpp = TRUE),
               class = "DPPError")                          # solving_chain.py:203-204
  good <- Variable(pos = TRUE, bounds = list(p * q, NULL))  # log-log affine
  gprob <- Problem(Minimize(good))
  expect_no_warning(ch <- construct_solving_chain(gprob, solver = "CLARABEL", gp = TRUE))
  expect_false(.s7_is(ch@reductions[[1L]], EvalParams))
  expect_equal(psolve(gprob, gp = TRUE, solver = "CLARABEL"), 6, tolerance = 1e-6)
  value(p) <- 5
  expect_equal(psolve(gprob, gp = TRUE, solver = "CLARABEL"), 15, tolerance = 1e-6)
})

# -- kron: the rule now lives in is_atom_convex under the DPP scope ------------
## Pins of existing behavior: the D_PORT.5 rule (data operand parameter-free) is
## unchanged, only its placement moved.

## @cvxpy atoms/affine/kron.py::kron.is_atom_convex
test_that("kron: DPP verdicts are unchanged after the rule moved into the curvature hook", {
  C <- Constant(diag(2))
  P <- Parameter(c(2, 2), value = diag(2))
  x <- Variable(c(2, 2))
  p <- Parameter(c(2, 2), value = matrix(1, 2, 2))
  expect_false(is_dpp(kron(P, x)))        # both sides FALSE
  expect_false(is_dpp(kron(P, C)))        # data operand P is parametric
  expect_true(is_dpp(kron(C, x + p)))     # ADR D_PORT.5: upstream says False
  expect_true(is_dpp(kron(x + p, C)))
  expect_false(is_dpp(Problem(Minimize(sum(kron(P, x))), list(x >= 0))))
  expect_true(is_dpp(Problem(Minimize(sum(kron(C, x + p))), list(x >= 0))))
  ## Outside the scope kron is affine whatever the parameters.
  expect_true(is_affine(kron(P, x)))
  expect_true(is_dcp(kron(P, x)))
})

## @cvxpy atoms/affine/conv.py::conv.is_atom_convex
test_that("conv: DPP verdicts are unchanged after the rule moved into the curvature hook", {
  k <- Constant(c(1, 2))
  pk <- Parameter(2, value = c(1, 2))
  y <- Variable(3)
  q <- Parameter(3, value = rep(0, 3))
  expect_true(is_dpp(conv(k, y + q)))     # conv.py:76-79 tests only the kernel
  expect_false(is_dpp(conv(pk, y)))
  expect_false(is_dpp(conv(y, pk)))       # CVXR accepts the kernel in either position
  expect_false(is_dpp(Problem(Minimize(sum(conv(pk, y))), list(y >= 0, y <= 1))))
  expect_true(is_affine(conv(pk, y)))
})
