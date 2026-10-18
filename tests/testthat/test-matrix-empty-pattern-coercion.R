## A problem with no constraints made CVXR build its empty constraint matrix
## with Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = ...) and no
## `x =`, which is a PATTERN matrix (ngCMatrix). Coercing that to dgCMatrix is
## deprecated in Matrix ('as(<ngCMatrix>, "dgCMatrix")' is deprecated), and a
## deprecated coercion eventually becomes defunct. The empty blocks are now
## built numeric (x = numeric(0)), as the other solver interfaces already did.
## Matrix.warnDeprecatedCoerce = 2 turns any deprecated coercion into an error.

## @cvxpy NONE
test_that("unconstrained QP: no deprecated Matrix coercion on any solver", {
  old_opt <- options(Matrix.warnDeprecatedCoerce = 2L); on.exit(options(old_opt), add = TRUE)
  solvers <- intersect(c("OSQP", "HIGHS", "PIQP", "GUROBI", "CPLEX", "XPRESS",
                         "MOSEK", "CLARABEL", "SCS", "ECOS", "CVXOPT", "SCIP"),
                       installed_solvers())
  if (length(solvers) == 0L) skip("no QP-capable solver installed")
  for (solver in solvers) {
    x <- Variable(2)
    prob <- Problem(Minimize(sum_squares(x - c(1, 2))))
    ## The xpress package itself calls as(<dgCMatrix>, "dgTMatrix")
    ## (x_to_sparse_triplet_matrix) on any matrix it is given; not CVXR's.
    options(Matrix.warnDeprecatedCoerce = if (solver == "XPRESS") 0L else 2L)
    expect_equal(psolve(prob, solver = solver), 0, tolerance = 1e-4, info = solver)
    expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-3, info = solver)
  }
})

## Bounds handed to the solver natively leave the LP with no constraint rows.
## @cvxpy NONE
test_that("bounds-only LP: no deprecated Matrix coercion on any solver", {
  old_opt <- options(Matrix.warnDeprecatedCoerce = 2L); on.exit(options(old_opt), add = TRUE)
  solvers <- intersect(c("HIGHS", "OSQP", "GUROBI", "CPLEX", "XPRESS", "MOSEK",
                         "CLARABEL", "SCS", "GLPK", "SCIP", "PIQP"),
                       installed_solvers())
  if (length(solvers) == 0L) skip("no LP solver installed")
  for (solver in solvers) {
    x <- Variable(2, bounds = list(c(1, 2), c(3, 4)))
    prob <- Problem(Minimize(sum_entries(x)))
    expect_equal(psolve(prob, solver = solver), 3, tolerance = 1e-4, info = solver)
  }
})
