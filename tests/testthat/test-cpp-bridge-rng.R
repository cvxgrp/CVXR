## The C++ canonicalization bridge must not touch R's random number stream.
##
## Through 1.9.2.9514 every C++ LinOp and LinOpVector constructor drew 16
## uniforms from R's RNG (Rcpp::runif) to build a debug-only "almost uuid"
## that nothing read outside #ifdef _R_DEBUG prints, so a cold compile
## advanced the user's stream by 16 draws per expression node and
## `set.seed(); psolve(); rnorm()` differed from `set.seed(); rnorm()`.
## Found by the lifecycle fuzz (notes/lifecycle_fuzz_2026-10-17.md, "Side
## findings"); removed at the user's request.

## @cvxpy NONE
test_that("a cold psolve() leaves .Random.seed unchanged", {
  require_solver("CLARABEL")
  set.seed(20261017)
  before <- .Random.seed
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Minimize(sum(x)), list(x >= c(1, 2, 3)))
  psolve(prob, solver = "CLARABEL")
  expect_identical(.Random.seed, before)
})

## @cvxpy NONE
test_that("a cold problem_data() on a parametric problem leaves .Random.seed unchanged", {
  require_solver("CLARABEL")
  set.seed(20261017)
  before <- .Random.seed
  p <- Parameter(3, value = c(1, 2, 3))
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Minimize(sum(x)), list(x >= p))
  invisible(problem_data(prob, solver = "CLARABEL"))
  expect_identical(.Random.seed, before)
})

## @cvxpy NONE
test_that("draws after a solve equal draws without one", {
  require_solver("CLARABEL")
  set.seed(20261017)
  a <- runif(3)
  set.seed(20261017)
  x <- Variable(2)
  psolve(Problem(Minimize(sum_squares(x - c(1, 2)))), solver = "CLARABEL")
  b <- runif(3)
  expect_identical(a, b)
})
