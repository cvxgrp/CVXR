## PowConeND -> PowCone3D conversion (reductions/cone2cone/exact.R,
## PowNDConversion; CVXPY exact.py:69-305 + cone_tree.py).
##
## CVXPY runs its PowConeND tests on CLARABEL *and* SCS: CLARABEL accepts
## PowConeND natively, SCS only takes 3D power cones, so the SCS leg is what
## exercises the conversion. CVXR's earlier ports kept only the CLARABEL leg
## (test-cvxpy-medium-logic-complex-pownd.R, test-cvxpy-parity-gaps-batch2.R,
## test-cvxpy-complex-pownd-parity.R, test-cvxpy-cone2cone-gaps.R) because,
## before this conversion, PowConeND could not reach SCS or MOSEK at all. This
## file holds the SCS legs, plus MOSEK legs (MOSEK also lacks PowConeND).
## CVXPY axis 0 -> R axis 2L, CVXPY axis 1 -> R axis 1L.

.pnd_solve <- function(prob, solver, ...) {
  psolve(prob, solver = solver, ...)
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"),
              info = solver)
  prob
}

## <primal, dual> over a PowConeND constraint (complementary slackness).
.pnd_complementarity <- function(con) {
  W <- value(con@.W); z <- value(con@.z)
  dW <- value(con@dual_variables[[1L]]); dz <- value(con@dual_variables[[2L]])
  sum(as.numeric(W) * as.numeric(dW)) + sum(as.numeric(z) * as.numeric(dz))
}

## TestPowND.pcp_3 / test_pow_cone_nd_3d: W = [[x0, x2], [x1, 1]].
.pnd_pcp3 <- function(axis, swap = FALSE) {
  x <- Variable(3)
  hypos <- Variable(2)
  W <- if (swap) bmat(list(list(x[1], x[2]), list(x[3], 1.0)))
       else bmat(list(list(x[1], x[3]), list(x[2], 1.0)))
  alpha <- matrix(c(0.2, 0.4, 0.8, 0.6), 2, 2, byrow = TRUE)
  if (axis == 1L) { W <- t(W); alpha <- t(alpha) }
  lin <- if (swap) x[1] + x[3] + 0.5 * x[2] == 2 else x[1] + x[2] + 0.5 * x[3] == 2
  con <- PowConeND(W, hypos, alpha, axis = axis)
  list(prob = Problem(Maximize(sum_entries(hypos) - x[1]), list(lin, con)),
       x = x, con = con)
}

## test_pow_cone_nd: three rows, so the conversion builds a tree.
.pnd_nd5 <- function(axis, swap = FALSE) {
  x <- Variable(5)
  hypos <- Variable(2)
  W <- if (swap) bmat(list(list(x[5], x[4]), list(x[2], x[1]), list(x[3], 1.0)))
       else bmat(list(list(x[1], x[4]), list(x[2], x[5]), list(x[3], 1.0)))
  alpha <- matrix(c(0.2, 0.4, 0.4, 0.3, 0.4, 0.3), 3, 2, byrow = TRUE)
  if (axis == 1L) { W <- t(W); alpha <- t(alpha) }
  lin <- if (swap) x[5] + x[2] + 0.5 * x[3] + 0.5 * x[4] + 0.25 * x[1] == 2
         else x[1] + x[2] + 0.5 * x[3] + 0.5 * x[4] + 0.25 * x[5] == 2
  obj <- if (swap) Maximize(sum_entries(hypos) - x[5]) else Maximize(sum_entries(hypos) - x[1])
  list(prob = Problem(obj, list(lin, PowConeND(W, hypos, alpha, axis = axis))), x = x)
}

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_3d
test_that("PowND conversion: 2-row PowConeND solves on SCS (both axes)", {
  require_solver("SCS")
  for (axis in c(2L, 1L)) {
    p <- .pnd_pcp3(axis)
    .pnd_solve(p$prob, "SCS")
    expect_equal(as.numeric(value(p$x)), c(0.06393515, 0.78320961, 2.30571048),
                 tolerance = 1e-3, info = axis)
  }
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_3d_variable_swap
test_that("PowND conversion: 2-row PowConeND, variables swapped, on SCS", {
  require_solver("SCS")
  for (axis in c(2L, 1L)) {
    p <- .pnd_pcp3(axis, swap = TRUE)
    .pnd_solve(p$prob, "SCS")
    expect_equal(as.numeric(value(p$x)), c(0.06393515, 2.30571048, 0.78320961),
                 tolerance = 1e-3, info = axis)
  }
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd
test_that("PowND conversion: 3-row PowConeND (tree) solves on SCS (both axes)", {
  require_solver("SCS")
  for (axis in c(2L, 1L)) {
    p <- .pnd_nd5(axis)
    .pnd_solve(p$prob, "SCS")
    expect_lt(max(abs(as.numeric(value(p$x)) - c(0, 0, 0, 2.28571379, 3.42857186))),
              1e-3)
  }
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_variable_swap
test_that("PowND conversion: 3-row PowConeND, variables swapped, on SCS", {
  require_solver("SCS")
  for (axis in c(2L, 1L)) {
    p <- .pnd_nd5(axis, swap = TRUE)
    .pnd_solve(p$prob, "SCS")
    expect_lt(max(abs(as.numeric(value(p$x)) - c(3.42857186, 0, 0, 2.28571379, 0))),
              1e-3)
  }
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_single_cone
test_that("PowND conversion: a single 2-row cone solves on SCS", {
  require_solver("SCS")
  for (axis in c(2L, 1L)) {
    x <- Variable(2)
    hypos <- Variable(1)
    W <- bmat(list(list(x[1]), list(x[2])))
    alpha <- matrix(c(0.2, 0.8), 2, 1)
    if (axis == 1L) { W <- t(W); alpha <- t(alpha) }
    prob <- Problem(Maximize(sum_entries(hypos) - x[1]),
                    list(x[1] + x[2] == 2, PowConeND(W, hypos, alpha, axis = axis)))
    .pnd_solve(prob, "SCS")
  }
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_decomposition
test_that("PowND conversion: balanced tree for n = 3 .. 100 on SCS", {
  require_solver("SCS")
  for (n in c(3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 16L, 32L, 100L)) {
    W <- Variable(n, pos = TRUE)
    z <- Variable()
    con <- PowConeND(W, z, matrix(1 / n, n, 1L))
    prob <- Problem(Maximize(z), list(con, W <= 2))
    .pnd_solve(prob, "SCS")
    expect_lt(abs(as.numeric(value(z)) - 2), 1e-3, label = paste("n =", n))
  }
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_nonuniform_alpha
test_that("PowND conversion: balanced tree, non-uniform alpha, on SCS", {
  require_solver("SCS")
  alpha <- c(0.05, 0.1, 0.15, 0.2, 0.2, 0.15, 0.1, 0.05)
  bounds <- c(1, 2, 3, 4, 5, 6, 7, 8)
  W <- Variable(8, pos = TRUE)
  z <- Variable()
  prob <- Problem(Maximize(z), list(PowConeND(W, z, matrix(alpha, ncol = 1L)),
                                    W <= bounds))
  .pnd_solve(prob, "SCS")
  expect_lt(abs(as.numeric(value(z)) - prod(bounds^alpha)), 1e-3)
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_multiple_cones
test_that("PowND conversion: balanced tree, k = 3 cones, on SCS", {
  require_solver("SCS")
  W <- Variable(c(8L, 3L), pos = TRUE)
  z <- Variable(3)
  prob <- Problem(Maximize(sum_entries(z)),
                  list(PowConeND(W, z, matrix(1 / 8, 8L, 3L), axis = 2L), W <= 2))
  .pnd_solve(prob, "SCS")
  expect_lt(max(abs(as.numeric(value(z)) - 2)), 1e-3)
})

.pnd_bounds6 <- function() matrix(c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7), 6, 2, byrow = TRUE)
.pnd_alpha6 <- function() matrix(c(0.1, 0.3, 0.2, 0.2, 0.3, 0.1, 0.15, 0.15,
                                   0.15, 0.1, 0.1, 0.15), 6, 2, byrow = TRUE)

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_balanced_tree_multiple_cones_nonuniform
test_that("PowND conversion: balanced tree, k = 2, non-uniform alpha, on SCS", {
  require_solver("SCS")
  bounds <- .pnd_bounds6(); alpha <- .pnd_alpha6()
  W <- Variable(c(6L, 2L), pos = TRUE)
  z <- Variable(2)
  prob <- Problem(Maximize(sum_entries(z)),
                  list(PowConeND(W, z, alpha, axis = 2L), W <= bounds))
  .pnd_solve(prob, "SCS")
  expect_lt(max(abs(as.numeric(value(z)) - apply(bounds^alpha, 2L, prod))), 1e-3)
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_pow_cone_nd_dual_variables
test_that("PowND conversion: recovered SCS duals match CLARABEL's native duals", {
  require_solver("SCS"); require_solver("CLARABEL")
  cases <- list(c(2, 1e-3), c(3, 1e-2), c(4, 1e-2), c(8, 1e-2), c(16, 5e-2))
  for (cs in cases) {
    n <- as.integer(cs[1L]); rtol <- cs[2L]
    duals <- lapply(c("CLARABEL", "SCS"), function(solver) {
      W <- Variable(c(n, 1L), pos = TRUE)
      z <- Variable(1)
      con <- PowConeND(W, z, matrix(1 / n, n, 1L), axis = 2L)
      .pnd_solve(Problem(Maximize(z[1]), list(con, W <= 2)), solver)
      list(w = as.numeric(value(con@dual_variables[[1L]])),
           z = as.numeric(value(con@dual_variables[[2L]])))
    })
    expect_equal(duals[[2L]]$w, duals[[1L]]$w, tolerance = rtol, info = n)
    expect_equal(duals[[2L]]$z, duals[[1L]]$z, tolerance = rtol, info = n)
  }
})

## @cvxpy test_cone2cone.py::TestPowND::test_pcp_3a test_cone2cone.py::TestPowND::test_pcp_3b
test_that("PowND conversion: TestPowND pcp_3 on SCS, with complementarity", {
  require_solver("SCS")
  for (axis in c(2L, 1L)) {
    p <- .pnd_pcp3(axis)
    .pnd_solve(p$prob, "SCS", eps = 1e-8)
    expect_equal(value(p$prob), 1.8073406786220672, tolerance = 1e-3, info = axis)
    expect_equal(as.numeric(value(p$x)), c(0.06393515, 0.78320961, 2.30571048),
                 tolerance = 1e-3, info = axis)
    expect_lt(abs(.pnd_complementarity(p$con)), 1e-3)
  }
})

## Fisher market (TestPowND.pcp_4): the reference is the log formulation.
.pnd_pcp4 <- function(b, solver) {
  set.seed(0)
  V <- matrix(runif(4L * 6L), 4L, 6L)
  X <- Variable(c(4L, 6L), nonneg = TRUE)
  u <- sum_entries(V * X, axis = 1L)
  log_prob <- Problem(Maximize(sum_entries(b * log(u))),
                      list(sum_entries(X, axis = 2L) <= 1))
  psolve(log_prob, solver = "CLARABEL")
  expect_X <- value(X)
  z <- Variable()
  con <- PowConeND(W = u, z = z, alpha = b)
  pow_prob <- Problem(Maximize(z), list(sum_entries(X, axis = 2L) <= 1, con))
  .pnd_solve(pow_prob, solver, eps = 1e-8)
  expect_equal(value(pow_prob), exp(value(log_prob)), tolerance = 1e-3)
  expect_equal(value(X), expect_X, tolerance = 1e-2)
  expect_lt(abs(.pnd_complementarity(con)), 1e-3)
}

## @cvxpy test_cone2cone.py::TestPowND::test_pcp_4a test_cone2cone.py::TestPowND::test_pcp_4b
test_that("PowND conversion: TestPowND pcp_4 (Fisher market) on SCS", {
  require_solver("SCS"); require_solver("CLARABEL")
  .pnd_pcp4(rep(0.25, 4L), "SCS")
  .pnd_pcp4(c(0.3, 0.15, 0.2, 0.35), "SCS")
})

## The tree the conversion records: n - 1 cones per column, root last, every
## variable index present once.
## @cvxpy NONE
test_that("PowND conversion: the recorded tree covers every entry once", {
  for (n in c(3L, 5L, 8L)) {
    W <- Variable(c(n, 2L)); z <- Variable(2)
    con <- PowConeND(W, z, matrix(1 / n, n, 2L), axis = 2L)
    res <- CVXR:::PowNDConversion$canonicalize(con, con@args)
    tree <- get("pow_tree", envir = res[[1L]]@.cache)
    expect_equal(CVXR:::get_all_cone_ids(tree), seq_len(n - 1L))
    expect_equal(CVXR:::get_root_cone_id(tree), n - 1L)
    expect_equal(.shape(res[[1L]])[2L], 2L * (n - 1L))
    idx <- unlist(lapply(CVXR:::get_leaf_nodes(tree), `[[`, "var_indices"))
    walk <- function(nd) {
      if (inherits(nd, "SingleVarNode")) return(nd$var_index)
      if (inherits(nd, "LeafNode")) return(nd$var_indices)
      c(walk(nd$left), walk(nd$right))
    }
    expect_equal(sort(walk(tree)), seq_len(n))
    expect_true(all(idx %in% seq_len(n)))
  }
})

## Exact geo_mean emits a PowConeND; before the conversion it failed on SCS
## and MOSEK ("does not support the required cone types: PowConeND").
## @cvxpy NONE
test_that("PowND conversion: exact geo_mean solves on SCS (vector and per column)", {
  require_solver("SCS")
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Maximize(geo_mean(x, approx = FALSE)), list(sum_entries(x) <= 3))
  expect_equal(psolve(prob, solver = "SCS"), 1, tolerance = 1e-3)
  ## Plan §4 U1 with the exact canonicalizer: analytic optimum 6.
  b <- c(2, 4, 6)
  X <- Variable(c(2L, 3L), nonneg = TRUE)
  prob2 <- Problem(Maximize(sum_entries(geo_mean(X, axis = 2L, approx = FALSE))),
                   list(sum_entries(X, axis = 2L) <= matrix(b, 1L)))
  expect_equal(psolve(prob2, solver = "SCS"), 6, tolerance = 1e-3)
  expect_equal(value(X), rbind(b / 2, b / 2), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("PowND conversion: PowConeND and exact geo_mean solve on MOSEK", {
  require_solver("MOSEK"); require_solver("CLARABEL")
  b <- c(2, 4, 6)
  X <- Variable(c(2L, 3L), nonneg = TRUE)
  prob <- Problem(Maximize(sum_entries(geo_mean(X, axis = 2L, approx = FALSE))),
                  list(sum_entries(X, axis = 2L) <= matrix(b, 1L)))
  expect_equal(psolve(prob, solver = "MOSEK"), 6, tolerance = 1e-5)
  ## Non-uniform tree, k = 2, with duals against CLARABEL's native ones.
  bounds <- .pnd_bounds6(); alpha <- .pnd_alpha6()
  duals <- lapply(c("CLARABEL", "MOSEK"), function(solver) {
    W <- Variable(c(6L, 2L), pos = TRUE)
    z <- Variable(2)
    con <- PowConeND(W, z, alpha, axis = 2L)
    .pnd_solve(Problem(Maximize(sum_entries(z)), list(con, W <= bounds)), solver)
    expect_lt(max(abs(as.numeric(value(z)) - apply(bounds^alpha, 2L, prod))), 1e-4)
    list(w = value(con@dual_variables[[1L]]), z = value(con@dual_variables[[2L]]))
  })
  expect_equal(duals[[2L]]$w, duals[[1L]]$w, tolerance = 1e-3)
  expect_equal(duals[[2L]]$z, duals[[1L]]$z, tolerance = 1e-3)
})

## save_dual_value(PowConeND) assigned the z dual as a plain vector, so a 1 x k
## z failed with "Invalid dimensions (k, 1)" (plan Stage 3 finding 2).
## CVXPY saves dz into z's shape (power.py:301-314).
## @cvxpy NONE
test_that("PowConeND: a row-vector z gets its dual in z's shape", {
  for (solver in c("CLARABEL", "SCS")) {
    require_solver(solver)
    W <- Variable(c(3L, 2L), pos = TRUE)
    z <- Variable(c(1L, 2L))
    con <- PowConeND(W, z, matrix(1 / 3, 3L, 2L), axis = 2L)
    prob <- Problem(Maximize(sum_entries(z)), list(con, W <= 2))
    expect_equal(psolve(prob, solver = solver), 4, tolerance = 1e-3, info = solver)
    dz <- value(con@dual_variables[[2L]])
    expect_equal(dim(dz), c(1L, 2L), info = solver)
    ## Same problem with a column z: identical duals.
    W2 <- Variable(c(3L, 2L), pos = TRUE)
    z2 <- Variable(c(2L, 1L))
    con2 <- PowConeND(W2, z2, matrix(1 / 3, 3L, 2L), axis = 2L)
    psolve(Problem(Maximize(sum_entries(z2)), list(con2, W2 <= 2)), solver = solver)
    expect_equal(as.numeric(dz), as.numeric(value(con2@dual_variables[[2L]])),
                 tolerance = 1e-3, info = solver)
    expect_equal(value(con@dual_variables[[1L]]), value(con2@dual_variables[[1L]]),
                 tolerance = 1e-3, info = solver)
  }
})

## The constraint shape override read CVXPY's 1-D W as expr_is_vector(W), so a
## 1 x n W with axis = 1 got shape (2, 1) instead of (n + 1, 1) and its duals
## were cut short on a solver that takes PowConeND natively. Checked against
## the same cone written as an n x 1 W with axis = 2.
## @cvxpy NONE
test_that("PowConeND: a 1 x n W with axis = 1 has shape (n + 1, 1) and full duals", {
  require_solver("CLARABEL")
  alpha <- c(0.2, 0.3, 0.5)
  mk <- function(row) {
    x <- Variable(3); h <- Variable(1)
    W <- if (row) t(x) else x
    a <- if (row) matrix(alpha, 1L) else matrix(alpha, ncol = 1L)
    con <- PowConeND(W, h, a, axis = if (row) 1L else 2L)
    lin <- sum_entries(x) == 3
    prob <- Problem(Maximize(h - x[1]), list(lin, con))
    list(prob = prob, con = con, lin = lin)
  }
  r <- mk(TRUE); c <- mk(FALSE)
  expect_equal(r$con@shape, c(4L, 1L))
  expect_equal(psolve(r$prob, solver = "CLARABEL"), psolve(c$prob, solver = "CLARABEL"),
               tolerance = 1e-5)
  expect_equal(as.numeric(value(r$con@dual_variables[[1L]])),
               as.numeric(value(c$con@dual_variables[[1L]])), tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(r$lin)), as.numeric(dual_value(c$lin)),
               tolerance = 1e-4)
})
