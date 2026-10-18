## Claims audit decision 30 (notes/audit/claims_audit_item1.md): the compile cost
## of the DPP tensor path, CVXR 1.9.2.9165. Three changes, none of which changes a
## result: (A) `format_constraints` applies the cone-interleaving row permutation
## to A_tensor as a signed permutation instead of a Kronecker product (ADR
## D_PORT.15, conic_solver.py:275-295); (B) a parameter-free P in
## `coeff_quad_form_tensor` is read as its numeric value (coeff_extractor.py:
## 149-158) instead of being canonicalized; (C) `get_problem_matrix_tensor`
## extracts the whole tensor in one C++ call (cppbackend.py:65-96). Seen failing
## on 9164: the `.restructure_A_tensor` and `getTensor` existence checks and the
## `dsCMatrix` quad_form case ("Need S4 class dgCMatrix"). Every other block pins
## behavior CVXR already had and cannot be seen failing; each holds the new
## mechanism against the reference it replaced (the Kronecker product, the
## per-element getters, or a closed form).

d30_kron_ref <- function(R, At, n_cols) {
  methods::as((Matrix::Diagonal(n_cols) %x% R) %*% At, "CsparseMatrix")
}
d30_same_csc <- function(A, B) {
  A <- methods::as(A, "CsparseMatrix"); B <- methods::as(B, "CsparseMatrix")
  identical(dim(A), dim(B)) && identical(A@p, B@p) && identical(A@i, B@i) &&
    isTRUE(all.equal(A@x, B@x))
}
d30_random_tensor <- function(m, x_len, p, nnz, seed = 1) {
  set.seed(seed)
  Matrix::sparseMatrix(i = sample.int(m * (x_len + 1L), nnz), j = sample.int(p, nnz, replace = TRUE),
                       x = rnorm(nnz), dims = c(m * (x_len + 1L), p))
}
## The cone-interleaving matrix format_constraints builds for a constraint list,
## recovered by formatting an identity A (the numeric A is restruct_mat %*% A).
d30_restruct <- function(cons, m) {
  methods::as(CVXR:::format_constraints(cons, Matrix::Diagonal(m), numeric(m),
                                        exp_cone_order = c(0L, 1L, 2L))$A, "CsparseMatrix")
}

# -- (A) the signed row permutation ---------------------------------------------

## @cvxpy reductions/solvers/conic_solvers/conic_solver.py::ConicSolver.format_constraints
test_that("A: a Zero-then-NonNeg block (signed identity) keeps the CSC pattern and flips the Zero rows", {
  expect_true(is.function(CVXR:::.restructure_A_tensor))      # absent on 9164
  R <- Matrix::bdiag(list(-Matrix::Diagonal(2), Matrix::Diagonal(3)))
  At <- d30_random_tensor(m = 5L, x_len = 2L, p = 3L, nnz = 9L)
  out <- CVXR:::.restructure_A_tensor(R, At, 3L)
  ref <- d30_kron_ref(R, At, 3L)
  expect_true(d30_same_csc(out, ref))
  ## the pattern is the input's, only the values changed
  expect_identical(out@p, methods::as(At, "CsparseMatrix")@p)
  expect_identical(out@i, methods::as(At, "CsparseMatrix")@i)
  ## rows 1-2 of every block are negated, rows 3-5 untouched
  Atc <- methods::as(At, "CsparseMatrix")
  r <- Atc@i %% 5L
  expect_equal(out@x, Atc@x * ifelse(r < 2L, -1, 1))
})

## @cvxpy reductions/solvers/conic_solvers/conic_solver.py::ConicSolver.format_constraints
test_that("A: the ExpCone and PowCone3D interleavings equal the Kronecker product, pattern included", {
  x <- Variable(2L); y <- Variable(2L); z <- Variable(2L)
  for (cons in list(list(ExpCone(x, y, z)), list(PowCone3D(x, y, z, 0.3)),
                    list(NonNeg(x), ExpCone(x, y, z), Zero(z - 1)))) {
    m <- sum(vapply(cons, function(cn) as.integer(prod(CVXR:::.shape(cn))), integer(1L)))
    R <- d30_restruct(cons, m)
    ## a real permutation, not the identity
    expect_false(isTRUE(all.equal(as.matrix(R), diag(m), check.attributes = FALSE)))
    At <- d30_random_tensor(m = m, x_len = 6L, p = 4L, nnz = 2L * m)
    expect_true(d30_same_csc(CVXR:::.restructure_A_tensor(R, At, 7L), d30_kron_ref(R, At, 7L)))
  }
})

## @cvxpy reductions/solvers/conic_solvers/conic_solver.py::ConicSolver.format_constraints
test_that("A: a block that is not a signed permutation falls back to the general product", {
  R <- Matrix::Matrix(c(1, 1, 0, 2), 2, 2, sparse = TRUE)   # column 1 has two nonzeros
  At <- d30_random_tensor(m = 2L, x_len = 3L, p = 2L, nnz = 6L, seed = 2)
  out <- CVXR:::.restructure_A_tensor(R, At, 4L)
  expect_equal(as.matrix(out), as.matrix(d30_kron_ref(R, At, 4L)))
  ## and the empty-constraint edge case goes the same way
  expect_equal(dim(CVXR:::.restructure_A_tensor(Matrix::Matrix(0, 0, 0, sparse = TRUE),
                                                Matrix::Matrix(0, 0, 2, sparse = TRUE), 4L)),
               c(0L, 2L))
})

## @cvxpy reductions/solvers/conic_solvers/conic_solver.py::ConicSolver.format_constraints problems/param_prob.py::ParamConeProg.apply_parameters
test_that("A: an ExpCone DPP problem re-solves through its formatted tensor to the fresh-compile answer", {
  require_solver("CLARABEL")
  n <- 5L
  mk <- function(qv) {
    w <- Variable(n); q <- Parameter(n, nonneg = TRUE); value(q) <- qv
    list(prob = Problem(Minimize(sum(q * exp(w))), list(w >= -1, sum(w) >= 1)), q = q)
  }
  q1 <- seq(1, 2, length.out = n); q2 <- rev(q1)
  a <- mk(q1)
  v1 <- psolve(a$prob, solver = "CLARABEL")
  pp <- a$prob@.cache$param_prog
  expect_true(isTRUE(pp@formatted))
  value(a$q) <- q2
  v2 <- psolve(a$prob, solver = "CLARABEL")
  expect_identical(a$prob@.cache$param_prog, pp)          # reused, not rebuilt
  b <- mk(q2)
  expect_equal(v2, psolve(b$prob, solver = "CLARABEL"), tolerance = 1e-7)
  expect_equal(v1, psolve(mk(q1)$prob, solver = "CLARABEL"), tolerance = 1e-12)
})

# -- (B) a parameter-free P is read as its value ----------------------------------

## @cvxpy utilities/coeff_extractor.py::CoeffExtractor.extract_quadratic_coeffs
test_that("B: dense, dgCMatrix and dsCMatrix P give the same P tensor, whose constant slice is 2P", {
  require_solver("CLARABEL")
  set.seed(1); P <- crossprod(matrix(rnorm(9), 3)) + diag(3)
  Ps <- list(dense = P,
             dgC = methods::as(methods::as(Matrix::Matrix(P, sparse = TRUE), "generalMatrix"), "CsparseMatrix"),
             dsC = Matrix::Matrix(P, sparse = TRUE))              # errored on 9164
  expect_s4_class(Ps$dsC, "dsCMatrix")
  tensors <- list(); objs <- list()
  for (nm in names(Ps)) {
    x <- Variable(3L); q <- Parameter(3L); value(q) <- c(1, 2, 3)
    prob <- Problem(Minimize(quad_form(x, Ps[[nm]]) + t(q) %*% x), list(x >= -2))
    v1 <- psolve(prob, solver = "CLARABEL")
    pp <- prob@.cache$param_prog
    tensors[[nm]] <- methods::as(pp@P_tensor, "CsparseMatrix")
    value(q) <- c(-1, 0, 1)
    objs[[nm]] <- c(v1, psolve(prob, solver = "CLARABEL"))
  }
  Pt <- tensors$dense
  expect_identical(dim(Pt), c(9L, 4L))
  expect_equal(which(Matrix::colSums(Pt != 0) > 0), 4L)          # constant slice only
  expect_equal(matrix(Pt[, 4L], 3L), 2 * P)
  for (nm in c("dgC", "dsC")) {
    expect_true(d30_same_csc(tensors[[nm]], Pt), info = nm)
    expect_equal(objs[[nm]], objs$dense, tolerance = 1e-9, info = nm)
  }
  ## closed form: min x'Px + q'x over x >= -2 (interior): x = -P^{-1} q / 2, obj = -q'P^{-1}q / 4
  cf <- function(qv) -drop(crossprod(qv, solve(P, qv))) / 4
  expect_equal(objs$dense, c(cf(c(1, 2, 3)), cf(c(-1, 0, 1))), tolerance = 1e-6)
})

## @cvxpy utilities/coeff_extractor.py::CoeffExtractor.extract_quadratic_coeffs
test_that("B: a parametric P still goes through the tensor extraction", {
  require_solver("CLARABEL")
  x <- Variable(3L); Pp <- Parameter(c(3L, 3L), PSD = TRUE); value(Pp) <- diag(3) * 2
  prob <- Problem(Minimize(quad_form(x, Pp) - sum(x)), list(x <= 1))
  v1 <- psolve(prob, solver = "CLARABEL")
  pp <- prob@.cache$param_prog
  Pt <- methods::as(pp@P_tensor, "CsparseMatrix")
  ## columns: the parameter's tensor slices, then the constant slice (a symmetric
  ## parameter is carried by its free elements, so not n * n columns)
  const_col <- pp@param_id_to_col[["-1"]] + 1L
  expect_identical(dim(Pt), c(9L, sum(unlist(pp@param_to_size))))
  expect_true(all(Matrix::colSums(Pt != 0)[-const_col] > 0))
  expect_equal(sum(Pt[, const_col] != 0), 0L)
  ## closed form with x interior: x = P^{-1} 1 / 2, obj = -1' P^{-1} 1 / 4
  expect_equal(v1, -3 / 8, tolerance = 1e-6)
  value(Pp) <- diag(3) * 3 + 0.5
  expect_equal(psolve(prob, solver = "CLARABEL"), -1 / 6, tolerance = 1e-6)
  expect_identical(prob@.cache$param_prog, pp)
})

# -- (C) the tensor extraction in one C++ call ------------------------------------

## @cvxpy cvxcore/python/cppbackend.py::build_matrix cvxcore/python/canonInterface.py::get_problem_matrix
test_that("C: get_problem_matrix_tensor lays out an affine parameter tensor as CVXPY does", {
  x <- Variable(3L); q <- Parameter(2L); value(q) <- c(1, 2)
  A <- matrix(1:6, 2)
  lin <- CVXR:::canonical_form(A %*% x + q)[[1L]]
  id_to_col <- c(0L); names(id_to_col) <- as.character(CVXR:::.id(x))
  qid <- as.character(CVXR:::.id(q))
  pts <- list(); pts[[qid]] <- 2L; pts[["-1"]] <- 1L
  pim <- list(); pim[[qid]] <- 0L; pim[["-1"]] <- 2L
  Tn <- CVXR:::get_problem_matrix_tensor(list(lin), id_to_col, var_length = 3L,
                                         param_to_size = pts, param_id_map = pim)
  ## rows: constraint row i + variable column j * m (m = 2), the constant column is
  ## j = 3; columns: the two q elements, then the constant slice.
  expected <- Matrix::sparseMatrix(i = c(7L, 8L, 1:6), j = c(1L, 2L, rep(3L, 6L)),
                                   x = c(1, 1, 1:6), dims = c(8L, 3L))
  expect_true(d30_same_csc(Tn, expected))
})

## @cvxpy cvxcore/python/cppbackend.py::build_matrix
test_that("C: the batched getter returns exactly what the per-element getters return, in order", {
  x <- Variable(3L); q <- Parameter(2L); value(q) <- c(1, 2)
  A <- matrix(1:6, 2)
  lin <- CVXR:::canonical_form(A %*% x + q)[[1L]]
  id_to_col <- c(0L); names(id_to_col) <- as.character(CVXR:::.id(x))
  pts <- c(2L, 1L); names(pts) <- c(as.character(CVXR:::.id(q)), "-1")
  tmp <- CVXR:::make_vec(); lin_vec <- CVXR:::CVXcanon.LinOpVector()
  tree <- CVXR:::build_lin_op_tree(lin, tmp); tmp$push_back(tree); lin_vec$push_back(tree)
  pd <- CVXR:::CVXcanon.build_matrix(lin_vec, id_to_col, 3L, pts)
  expect_true(is.function(pd$getTensor))                        # absent on 9164
  col_offsets <- c(0L, 2L); constr_length <- 2L
  tens <- pd$getTensor(as.integer(names(pts)), unname(pts), col_offsets, constr_length)
  V <- numeric(0); I <- integer(0); J <- integer(0)
  for (k in seq_along(pts)) {
    pd$setParamId(as.integer(names(pts)[k]))
    for (i in seq_len(pts[[k]]) - 1L) {
      pd$setVecIdx(i)
      v <- pd$getV()
      V <- c(V, v)
      I <- c(I, pd$getI() + pd$getJ() * constr_length)
      J <- c(J, rep(i + col_offsets[k], length(v)))
    }
  }
  expect_identical(tens$V, V)
  expect_identical(tens$I, I)
  expect_identical(tens$J, J)
  expect_length(tens$V, 8L)
  ## programming errors stop instead of reading out of bounds
  expect_error(pd$getTensor(as.integer(names(pts)), unname(pts), 0L, constr_length), "same length")
  expect_error(pd$getTensor(as.integer(names(pts)), unname(pts) + 5L, col_offsets, constr_length), "tensor slot")
})
