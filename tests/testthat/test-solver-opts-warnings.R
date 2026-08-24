## The standard solver-parameter layer must not discard silently.
##
## `feastol`/`reltol`/`abstol`/`num_iter` are CVXR's portable spellings,
## translated per solver by .SOLVER_DEFAULT_PARAM. The map is uneven --
## MOSEK/CPLEX/HIGHS carry only num_iter, GLPK/GLPK_MI nothing -- and
## .build_solver_params used to skip an unmapped standard param with no
## signal at all: psolve(prob, solver = "MOSEK", feastol = 1e-9) did
## nothing, silently. Now it warns (class "cvxr_unmapped_solver_param").
##
## These are pure-function tests on .build_solver_params: no solver is
## invoked, so nothing here needs a licence or an installed backend.

## @cvxpy NONE -- the standard-name layer is R-specific (CVXPY passes
## kwargs verbatim and has no portable tolerance names).
test_that("an unmapped standard param warns and is dropped", {
  ## CPLEX: Rcplex's control list exposes no feasibility tolerance at all
  ## (trace/method/itlim/epgap/epagap/tilim/...), so feastol is unmappable.
  opts <- solver_opts(feastol = 1e-9)
  expect_warning(
    params <- .build_solver_params("CPLEX", opts),
    class = "cvxr_unmapped_solver_param"
  )
  expect_false("feastol" %in% names(params))

  ## GLPK has no mapping at all: every standard param set must be named.
  opts_all <- solver_opts(feastol = 1e-8, reltol = 1e-8,
                          abstol = 1e-8, num_iter = 50)
  w <- tryCatch(
    { .build_solver_params("GLPK", opts_all); NULL },
    warning = function(w) w
  )
  expect_s3_class(w, "cvxr_unmapped_solver_param")
  expect_match(conditionMessage(w), "feastol")
  expect_match(conditionMessage(w), "reltol")
  expect_match(conditionMessage(w), "abstol")
  expect_match(conditionMessage(w), "num_iter")
})

## @cvxpy NONE
test_that("mapped standard params translate without warning", {
  expect_no_warning(p <- .build_solver_params("CLARABEL",
                                              solver_opts(feastol = 1e-7)))
  expect_identical(p$tol_feas, 1e-7)

  ## Partially mapped solver: only the unmapped one is named in the warning.
  w <- tryCatch(
    { p2 <- .build_solver_params("SCS",
                                 solver_opts(feastol = 1e-7, reltol = 1e-6))
      NULL },
    warning = function(w) w
  )
  expect_s3_class(w, "cvxr_unmapped_solver_param")
  expect_match(conditionMessage(w), "feastol")
  expect_no_match(conditionMessage(w), "reltol")
})

## @cvxpy NONE
test_that("native names pass through and win over standard names", {
  ## Native passthrough is untouched by the warning path.
  expect_no_warning(
    p <- .build_solver_params("MOSEK",
                              solver_opts(mosek_params = list(a = 1)))
  )
  expect_identical(p$mosek_params, list(a = 1))

  ## Precedence: a directly-set native name beats the standard translation.
  expect_no_warning(
    p2 <- .build_solver_params("CLARABEL",
                               solver_opts(feastol = 1e-7, tol_feas = 1e-3))
  )
  expect_identical(p2$tol_feas, 1e-3)
})

## @cvxpy NONE
test_that("unset standard params never warn", {
  expect_no_warning(.build_solver_params("MOSEK", solver_opts()))
  expect_no_warning(.build_solver_params("GLPK", solver_opts()))
})

## ---- item 2: the filled tolerance rows -------------------------------

## @cvxpy NONE -- CVXR's own portable-name layer; MSK_* names chosen from
## the same families as CVXPY's MOSEK.tolerance_params()
## (mosek_conif.py:713-739).
test_that("MOSEK and HiGHS standard params translate to real native names", {
  p <- .build_solver_params("MOSEK",
                            solver_opts(feastol = 1e-9, reltol = 1e-7,
                                        num_iter = 50))
  expect_identical(p$MSK_DPAR_INTPNT_CO_TOL_PFEAS, 1e-9)
  expect_identical(p$MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1e-7)
  expect_identical(p$MSK_IPAR_INTPNT_MAX_ITERATIONS, 50L)
  ## abstol is deliberately unmapped for MOSEK (MIP-only gap upstream).
  expect_warning(.build_solver_params("MOSEK", solver_opts(abstol = 1e-8)),
                 class = "cvxr_unmapped_solver_param")

  h <- .build_solver_params("HIGHS", solver_opts(feastol = 1e-8))
  expect_identical(h$primal_feasibility_tolerance, 1e-8)
})

## @cvxpy mosek_conif.py::MOSEK::solve_via_data -- parameters must reach
## the solver; through 1.9.2 every one of these calls CRASHED with
## "missing value where TRUE/FALSE needed" because all options were merged
## into Rmosek's interface `opts` instead of problem$dparam/$iparam.
test_that("MOSEK solver parameters reach Rmosek instead of crashing", {
  skip_if_not_installed("Rmosek")
  skip_if_not("MOSEK" %in% installed_solvers(), "MOSEK not available")

  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  ## Rmosek-native dparam list
  v1 <- psolve(prob, solver = "MOSEK",
               dparam = list(INTPNT_CO_TOL_REL_GAP = 1e-3))
  expect_equal(v1, 2, tolerance = 1e-4)

  ## standard portable name, now translated to MSK_DPAR_*
  prob2 <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  v2 <- psolve(prob2, solver = "MOSEK", feastol = 1e-9)
  expect_equal(v2, 2, tolerance = 1e-6)

  ## CVXPY's mosek_params carrier
  prob3 <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  v3 <- psolve(prob3, solver = "MOSEK",
               mosek_params = list(MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-9))
  expect_equal(v3, 2, tolerance = 1e-6)

  ## unknown options raise, as upstream (mosek_conif.py:674-676)
  prob4 <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  expect_error(psolve(prob4, solver = "MOSEK", not_an_option = 1),
               "Invalid option")
  prob5 <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  expect_error(psolve(prob5, solver = "MOSEK",
                      mosek_params = list(BAD_NAME = 1)),
               "Invalid MOSEK parameter")
})

## @cvxpy NONE -- HiGHS is in Imports, so this end-to-end runs everywhere.
test_that("HiGHS receives the translated feasibility tolerance", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  v <- psolve(prob, solver = "HIGHS", feastol = 1e-9)
  expect_equal(v, 2, tolerance = 1e-6)
  expect_true(status(prob) %in% SOLUTION_PRESENT)
})

## ---- item 3: the eps keyword ----------------------------------------

## @cvxpy scs_conif.py::SCS::parse_solver_options -- under SCS >= 3.0,
## `eps` sets both eps_abs and eps_rel (scs_conif.py:206-209). Before the
## port, eps was silently ignored and this solve came back ~1.2e-5 off:
## the tightness assertion below is what fails on the pre-fix build.
test_that("SCS eps keyword tightens both tolerances", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  v <- psolve(prob, solver = "SCS", eps = 1e-10)
  expect_lt(abs(v - 2), 1e-7)
})

## @cvxpy mosek_conif.py::MOSEK::parse_eps_keyword -- eps fans out to the
## 17 tolerance parameters (:690-710); explicit tolerances win. On the
## pre-fix build this call ABORTED as an invalid MOSEK option.
test_that("MOSEK eps keyword fans out and solves", {
  skip_if_not_installed("Rmosek")
  skip_if_not("MOSEK" %in% installed_solvers(), "MOSEK not available")

  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  v <- psolve(prob, solver = "MOSEK", eps = 1e-9)
  expect_equal(v, 2, tolerance = 1e-7)

  ## eps together with an explicit tolerance (the explicit one must win
  ## without error): portable feastol translates to MSK_DPAR_INTPNT_CO_
  ## TOL_PFEAS before the fan-out sees it, so both routes coexist.
  prob2 <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  v2 <- psolve(prob2, solver = "MOSEK", eps = 1e-6, feastol = 1e-9)
  expect_equal(v2, 2, tolerance = 1e-5)
})
