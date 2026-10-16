## Claims audit, decision 25 (CVXR 1.9.2.9161): C0540, enforce_dpp and the
## compile cache key.
##
## CVXPY 1.9.3 `Cache.make_key` (problem.py:118-119) is
## (solver, gp, ignore_dpp, use_quad_obj): the key holds what alters
## compilation (:824-825), and enforce_dpp only decides whether building the
## chain raises (solving_chain.py:203-204). Through 1.9.2.9160 CVXR keyed on
## enforce_dpp too, so every flip of the flag on a DPP problem rebuilt the
## chain and discarded the cached ParamConeProg.
##
## The hybrid the user chose ("Hybrid, go"): the key is upstream's, AND a
## cache hit under enforce_dpp = TRUE still raises DPPError when the cached
## chain was built for a DPP-violating problem. CVXPY 1.9.3 solves silently
## there (measured), against its own contract for the flag (problem.py:624-627
## "a DPPError will be thrown when trying to solve a non-DPP problem").
## Recorded as ADR D_PORT.13.
##
## Seen failing on 1.9.2.9160: the two "cache hit" tests (chain and
## ParamConeProg identity across flips; chain untouched by the re-raise). The
## rest pin behavior that is the same on both builds and cannot be seen
## failing; they guard the hybrid against drifting toward either plain option.

non_dpp_problem <- function() {       # parameter * parameter: not DPP
  x <- Variable()
  p <- Parameter(value = 2)
  Problem(Minimize(square(x - p * p)))
}

dpp_conic_problem <- function() {     # DPP, conic path, ParamConeProg cached
  y <- Variable(3)
  q <- Parameter(3, value = c(1, 2, 3))
  Problem(Minimize(sum_squares(y - q) + p_norm(y, 1)),
          list(y >= 0, sum(y) <= 10))
}

chain_of <- function(prob) prob@.cache$compile_chain

# -- the key is upstream's --------------------------------------------------

## @cvxpy cvxpy/problems/problem.py:118-119
test_that("flipping enforce_dpp on a DPP problem is a cache hit", {
  prob <- dpp_conic_problem()
  expect_true(is_dpp(prob))
  v1 <- psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE)
  ## (solver, gp, use_quad_obj, ignore_dpp): four entries, as make_key.
  expect_length(prob@.cache$compile_key, 4L)
  ch1 <- chain_of(prob)
  pp1 <- prob@.cache$param_prog
  expect_false(is.null(pp1))
  v2 <- psolve(prob, solver = "CLARABEL", enforce_dpp = FALSE)
  expect_identical(chain_of(prob), ch1)
  expect_identical(prob@.cache$param_prog, pp1)
  v3 <- psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE)
  expect_identical(chain_of(prob), ch1)
  expect_identical(prob@.cache$param_prog, pp1)
  expect_equal(v2, v1, tolerance = 1e-6)
  expect_equal(v3, v1, tolerance = 1e-6)
})

## @cvxpy cvxpy/problems/problem.py:118-119
test_that("ignore_dpp stays in the key: it changes the chain that is built", {
  prob <- non_dpp_problem()
  expect_no_warning(psolve(prob, solver = "CLARABEL", ignore_dpp = TRUE))
  ch_ignore <- chain_of(prob)
  ## A different key: rebuilt, and the rebuild refuses as upstream does.
  expect_error(psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE),
               class = "DPPError", regexp = "does not follow DPP rules")
  ## The failed build left nothing behind (problem.py:831-842 ordering).
  expect_null(prob@.cache$compile_key)
  expect_null(chain_of(prob))
  expect_false(identical(chain_of(prob), ch_ignore))
})

# -- the hybrid: the flag keeps its promise on a cache hit --------------------

## @cvxpy cvxpy/problems/problem.py:624-627
test_that("enforce_dpp = TRUE after a plain solve of a non-DPP problem still raises", {
  ## RECORDED DIVERGENCE (ADR D_PORT.13): CVXPY 1.9.3 reuses the cached
  ## non-DPP chain here and solves silently.
  prob <- non_dpp_problem()
  expect_false(is_dpp(prob))
  expect_warning(psolve(prob, solver = "CLARABEL"), class = "DPPWarning")
  ch1 <- chain_of(prob)
  key1 <- prob@.cache$compile_key
  expect_true(isTRUE(ch1@.cache$dpp_violated))
  expect_error(psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE),
               class = "DPPError", regexp = "does not follow DPP rules")
  ## No rebuild: the cache is exactly what the plain solve left.
  expect_identical(chain_of(prob), ch1)
  expect_identical(prob@.cache$compile_key, key1)
  ## And the plain solve still works on that same chain, silently
  ## (the warning belongs to the build, which did not happen again).
  expect_no_warning(psolve(prob, solver = "CLARABEL"))
  expect_identical(chain_of(prob), ch1)
})

## @cvxpy cvxpy/problems/problem.py:624-627
test_that("enforce_dpp = TRUE on a fresh non-DPP problem raises and caches nothing", {
  ## Same on both sides and on both builds: a pin.
  prob <- non_dpp_problem()
  expect_error(psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE),
               class = "DPPError", regexp = "does not follow DPP rules")
  expect_null(prob@.cache$compile_key)
  expect_error(psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE),
               class = "DPPError", regexp = "does not follow DPP rules")
})

## @cvxpy cvxpy/reductions/solvers/solving_chain.py:202-207
test_that("the re-raise keys on the DPP status, not on EvalParams being present", {
  ## A DPP quadratic problem whose chain carries EvalParams for another reason
  ## (the QP-with-parameters route, claims audit C0919): the cached chain must
  ## not be mistaken for a non-DPP one.
  y <- Variable(3)
  q <- Parameter(3, value = c(1, 2, 3))
  prob <- Problem(Minimize(sum_squares(y - q)), list(y >= 0))
  expect_true(is_dpp(prob))
  v1 <- expect_no_warning(psolve(prob, solver = "CLARABEL"))
  ch1 <- chain_of(prob)
  expect_false(isTRUE(ch1@.cache$dpp_violated))
  v2 <- expect_no_error(psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE))
  expect_identical(chain_of(prob), ch1)
  expect_equal(v2, v1, tolerance = 1e-8)
})

## @cvxpy cvxpy/problems/problem.py:817-818
test_that("the contradictory flag pair is still rejected before the key", {
  prob <- dpp_conic_problem()
  psolve(prob, solver = "CLARABEL")
  ch1 <- chain_of(prob)
  expect_error(psolve(prob, solver = "CLARABEL",
                      enforce_dpp = TRUE, ignore_dpp = TRUE),
               class = "DPPError", regexp = "Cannot set")
  expect_identical(chain_of(prob), ch1)
})
