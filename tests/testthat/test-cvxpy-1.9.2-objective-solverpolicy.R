## Two CVXPY 1.9.2 gaps that needed CODE, not just a test:
##   - arithmetic on Objective objects had no Ops dispatch at all
##   - DEFAULT_TO_COMMERCIAL_SOLVERS (PR #3352) was unported
## plus the two PIQP cases, which are blocked on an upstream piqp-r fix.

## Scoped env-var / option helpers. Deliberately NOT withr: it is not a CVXR
## dependency and this is the only place that would need it.
.with_env <- function(name, value, expr) {
  old <- Sys.getenv(name, unset = NA_character_)
  on.exit(
    if (is.na(old)) Sys.unsetenv(name)
    else do.call(Sys.setenv, stats::setNames(list(old), name)),
    add = TRUE)
  if (is.na(value)) Sys.unsetenv(name)
  else do.call(Sys.setenv, stats::setNames(list(value), name))
  expr
}

.with_opt <- function(name, value, expr) {
  old <- options(stats::setNames(list(value), name))
  on.exit(options(old), add = TRUE)
  expr
}

# ===================================================================
# Objective arithmetic and predicates  (objective.py:61-95, 134-192, 210-268)
# ===================================================================

## @cvxpy test_objectives.py::TestObjectives::test_objective_arithmetic_and_predicate_branches
test_that("objective arithmetic follows the CVXPY rules", {
  ## Before this port every arithmetic line below failed with "non-numeric
  ## argument to binary operator": the helper functions existed and were
  ## correct, but nothing dispatched to them.
  x <- Variable(nonneg = TRUE)
  minimize <- Minimize(square(x))
  maximize <- Maximize(sqrt(x))

  ## Additive identity. NOTE a deliberate deviation from upstream: CVXPY
  ## accepts `0 + obj` (via __radd__) but rejects `obj + 0`, an artifact of
  ## Python's __add__/__radd__ split that exists so sum() can start at 0. R's
  ## Ops handler sees both operands at once, and CVXR's Problem handler already
  ## accepts either side, so both are accepted here.
  expect_identical(0 + minimize, minimize)
  expect_identical(minimize + 0, minimize)
  expect_error(1 + minimize)
  expect_error(minimize - 1)

  ## Negation flips the direction.
  expect_true(S7_inherits(0 - minimize, Maximize))
  expect_true(S7_inherits(-minimize, Maximize))
  expect_true(S7_inherits(-maximize, Minimize))
  expect_error(1 - minimize)

  ## Non-numeric scaling is refused.
  expect_error(minimize * "bad")
  expect_error(minimize / "bad")

  ## Scaling by a negative scalar reverses the direction.
  expect_true(S7_inherits(minimize * 2, Minimize))
  expect_true(S7_inherits(minimize * -2, Maximize))
  expect_true(S7_inherits(maximize * 2, Maximize))
  expect_true(S7_inherits(maximize * -2, Minimize))
  expect_true(S7_inherits(minimize / 2, Minimize))
  ## Scalar on the left is the same operation (upstream's __rmul__ = __mul__).
  expect_true(S7_inherits(2 * minimize, Minimize))
  expect_true(S7_inherits(-2 * minimize, Maximize))

  ## Mixing senses violates DCP in either order.
  expect_error(minimize + maximize)
  expect_error(maximize + minimize)

  ## Predicates.
  expect_null(value(minimize))
  value(x) <- 3
  expect_equal(as.numeric(value(minimize)), 9)
  expect_true(is_quadratic(minimize))
  expect_true(is_qpwa(minimize))
  expect_true(is_dpp(minimize, "dcp"))
  expect_true(is_dpp(maximize, "dcp"))
  expect_false(is_dpp(minimize, "dgp"))
  expect_false(is_dpp(maximize, "dgp"))

  ## An unrecognized context must RAISE, not silently answer the dcp question.
  expect_error(is_dpp(minimize, "bad"), "Unsupported DPP context")
  expect_error(is_dpp(maximize, "bad"), "Unsupported DPP context")

  ## Upstream's final two assertions cover Minimize/Maximize.primal_to_result,
  ## a static that passes the solver primal through or negates it. CVXR has no
  ## such method: the Maximize sign flip happens during canonicalization
  ## (objective.R, `method(canonicalize, Maximize)` -> neg_expr_linop), so the
  ## equivalent observable is that a maximization reports a positive optimum.
  y <- Variable(nonneg = TRUE)
  expect_equal(psolve(Problem(Maximize(y), list(y <= 2.5)), solver = "CLARABEL"),
               2.5, tolerance = 1e-6)
})

# ===================================================================
# DEFAULT_TO_COMMERCIAL_SOLVERS  (settings.py:114-123, problem_form.py:356)
# ===================================================================

## @cvxpy test_problem_form.py::TestPickDefaultSolver::test_default_to_commercial_solvers_env_var
test_that("the commercial-solver env var parses like CVXPY's", {
  ## Mirrors upstream `_env_var_to_bool`: 0/false/no/off are FALSE (case- and
  ## whitespace-insensitive), an unset variable takes the default, and anything
  ## else is TRUE.
  .with_env("CVXR_DEFAULT_TO_COMMERCIAL_SOLVERS", "0",
    expect_false(.env_var_to_bool("CVXR_DEFAULT_TO_COMMERCIAL_SOLVERS", TRUE)))
  .with_env("CVXR_DEFAULT_TO_COMMERCIAL_SOLVERS", "1",
    expect_true(.env_var_to_bool("CVXR_DEFAULT_TO_COMMERCIAL_SOLVERS", FALSE)))

  for (falsey in c("false", "FALSE", " off ", "No", "0"))
    .with_env("CVXR_TEST_FLAG", falsey,
              expect_false(.env_var_to_bool("CVXR_TEST_FLAG", TRUE)))
  for (truthy in c("1", "true", "yes", "anything"))
    .with_env("CVXR_TEST_FLAG", truthy,
              expect_true(.env_var_to_bool("CVXR_TEST_FLAG", FALSE)))

  ## Unset takes the supplied default, either way.
  .with_env("CVXR_TEST_FLAG", NA_character_, {
    expect_true(.env_var_to_bool("CVXR_TEST_FLAG", TRUE))
    expect_false(.env_var_to_bool("CVXR_TEST_FLAG", FALSE))
  })
})

## @cvxpy test_problem_form.py::TestPickDefaultSolver::test_default_to_commercial_solvers_setting
test_that("the commercial-solver preference can be switched off", {
  ## Default is ON: upstream's default, and CVXR's previous unconditional
  ## behavior.
  .with_opt("CVXR.default_to_commercial_solvers", NULL,
    .with_env("CVXR_DEFAULT_TO_COMMERCIAL_SOLVERS", NA_character_,
              expect_true(.default_to_commercial_solvers())))

  ## The R option wins when set.
  .with_opt("CVXR.default_to_commercial_solvers", FALSE,
            expect_false(.default_to_commercial_solvers()))
  .with_opt("CVXR.default_to_commercial_solvers", TRUE,
            expect_true(.default_to_commercial_solvers()))

  ## The environment variable applies when the option is unset...
  .with_opt("CVXR.default_to_commercial_solvers", NULL,
    .with_env("CVXR_DEFAULT_TO_COMMERCIAL_SOLVERS", "0",
              expect_false(.default_to_commercial_solvers())))
  ## ...and the option beats the environment.
  .with_opt("CVXR.default_to_commercial_solvers", TRUE,
    .with_env("CVXR_DEFAULT_TO_COMMERCIAL_SOLVERS", "0",
              expect_true(.default_to_commercial_solvers())))

  ## The flag must actually gate the picker: with it off, a plain LP is never
  ## routed to a commercial solver, whether or not one is installed.
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  .with_opt("CVXR.default_to_commercial_solvers", FALSE, {
    picked <- .pick_default_conic_solver(prob, names(SOLVER_MAP_CONIC),
                                         list(NonNeg))
    expect_false(picked %in% .COMMERCIAL_CONIC_SOLVERS)
    expect_equal(picked, CLARABEL_SOLVER)
  })
})

# ===================================================================
# PIQP -- blocked upstream, not a CVXR defect
# ===================================================================

## @cvxpy test_qp_solvers.py::TestPiqpInterface::test_piqp_invalid_backend
## @v19-pending: piqp-r backend validation [upstream]
test_that("PIQP rejects an unrecognized backend", {
  skip(paste("Blocked on piqp-r. CVXR passes `backend` straight through and the",
             "R piqp package validates it with match.arg, giving \"'arg' should",
             "be one of auto, sparse, dense\" rather than CVXPY's \"backend must",
             "be either dense or sparse\"; it also accepts an `auto` backend that",
             "CVXPY does not. To be fixed in piqp-r, not here (2026-08-12)."))
})

## @cvxpy test_qp_solvers.py::TestPiqpInterface::test_piqp_unknown_setting
## @v19-pending: piqp-r unknown-setting validation [upstream]
test_that("PIQP rejects an unknown solver setting", {
  skip(paste("Blocked on piqp-r. An unrecognized solver setting is SILENTLY",
             "IGNORED rather than raising, so a typo'd option is accepted;",
             "CVXPY raises TypeError('Unrecognized solver setting'). Measured",
             "2026-08-12: psolve(prob, solver = 'PIQP', not_a_real_setting = 1)",
             "returned an optimal solution. To be fixed in piqp-r, not here."))
})
