## The last six CVXPY 1.9.2 test-port gaps.
##
## Every assertion below was probed against the installed package before being
## written. Five of the six subjects were already correct; porting the sixth
## (the complex2real mapping helpers) exposed a live CVXR bug in
## `project(<Leaf>)` — see "a Hermitian leaf accepts a complex value" below.
##
## Two of these are PARTIAL ports, and both say so at the point of divergence
## rather than in prose here: the ExpCone metadata test cannot cover
## `as_quad_approx` (RelEntrConeQuad is NOT TARGETED in CVXR), and the lin_utils
## test cannot cover the 14 upstream helpers CVXR has no counterpart for.

# ===================================================================
# Canonicalizer registry keys            (test_canon_methods.py)
# ===================================================================

## @cvxpy test_canon_methods.py::test_canon_method_keys_are_classes
test_that("every canonicalizer dispatches on an Expression or Constraint class", {
  ## Upstream keys nine CANON_METHODS dicts by atom class and looks them up via
  ## type(expr), so a key that is a FUNCTION rather than a class can never
  ## match -- the atom silently falls through to default handling. That is a
  ## real bug upstream hit once (dgp2dcp keyed the `quad_form` function instead
  ## of the QuadForm class, producing wrong DGP results).
  ##
  ## CVXR cannot express that particular bug: canonicalizers are S7 method
  ## registrations, and `method(gen, x) <- f` rejects an `x` that is not a
  ## class. The check still has content, though, because S7 accepts base types
  ## and any S7 class -- registering a canonicalizer on, say, `class_numeric`
  ## would compile and then never fire for an atom. So: walk every generic's
  ## dispatch table and require each key to be a CVXR Expression/Constraint
  ## subclass, or the one documented `S7_object` fallback.

  ## Does `cls` have Expression or Constraint among its S7 ancestors?
  descends_from_expr_or_constraint <- function(cls) {
    while (inherits(cls, "S7_class")) {
      if (identical(cls, Expression) || identical(cls, Constraint)) return(TRUE)
      cls <- S7::prop(cls, "parent")
    }
    FALSE
  }

  ## CVXR's counterparts of upstream's nine registries. eliminate_pwl has none
  ## (merged into Dcp2Cone, see eliminate_pwl.R) and discrete2mixedint
  ## registers FiniteSet on dcp_canonicalize (valinvec2mixedint.R:82), so both
  ## are covered by `dcp_canonicalize` here.
  generics <- list(
    dcp_canonicalize    = dcp_canonicalize,
    quad_canonicalize   = quad_canonicalize,
    c2r_canonicalize    = c2r_canonicalize,
    dgp_canonicalize    = dgp_canonicalize,
    smooth_canonicalize = smooth_canonicalize,
    has_dcp_canon       = has_dcp_canon,
    has_quad_canon      = has_quad_canon,
    has_c2r_canon       = has_c2r_canon)

  for (gname in names(generics)) {
    keys <- ls(S7::prop(generics[[gname]], "methods"), all.names = TRUE)
    expect_gt(length(keys), 1L)
    for (key in keys) {
      ## The single permitted non-atom key: the base-class default that stands
      ## in for "absent from the dict" upstream (dcp2cone.R:22, :50).
      if (identical(key, "S7_object")) next
      expect_match(key, "^CVXR::", info = paste(gname, key))
      nm <- sub("^CVXR::", "", key)
      expect_true(exists(nm, envir = asNamespace("CVXR"), inherits = FALSE),
                  info = paste(gname, key, "does not resolve to a CVXR object"))
      cls <- get(nm, envir = asNamespace("CVXR"))
      expect_true(inherits(cls, "S7_class"),
                  info = paste(gname, key, "is not an S7 class"))
      expect_true(descends_from_expr_or_constraint(cls),
                  info = paste(gname, key, "is not an Expression/Constraint subclass"))
    }
  }

  ## The cone2cone conversions are a genuine table rather than S7 dispatch
  ## (S7 class objects cannot be list names -- exact.R:43), so their `source`
  ## and `targets` entries get the same check.
  for (tbl in list(EXACT_CONE_CONVERSIONS, .EXACT_CONVERSIONS)) {
    expect_gt(length(tbl), 0L)
    for (conv in tbl) {
      expect_true(descends_from_expr_or_constraint(conv$source))
      for (tgt in conv$targets)
        expect_true(descends_from_expr_or_constraint(tgt))
    }
  }
})

# ===================================================================
# Complex2Real                                   (test_complex.py)
# ===================================================================

## @cvxpy test_complex.py::TestComplex::test_complex2real_invert_primal_cases
test_that("Complex2Real recombines general and Hermitian primal variables", {
  ## Drives reduction_invert directly with a hand-built Solution, so it covers
  ## the three complex branches (is_imag / hermitian / general complex) without
  ## needing a solver to produce them.
  x <- Variable(2, complex = TRUE)
  H <- Variable(c(2, 2), hermitian = TRUE)
  problem <- Problem(Minimize(sum_squares(Re(x))), list(H %>>% 0))

  reduction <- Complex2Real()
  applied <- reduction_apply(reduction, problem)
  inverse_data <- applied[[2L]]

  ## CVXR's real2imag maps an original leaf id to the id of its IMAGINARY part;
  ## the real part keeps the original id. Upstream mints fresh ids for both
  ## halves (canon_methods._variables); the deliberate divergence is ADR
  ## D_19.5's shared-id design, and the recovered values are identical either
  ## way, which is what this test asserts.
  xid <- as.character(x@id)
  Hid <- as.character(H@id)
  x_imag <- as.character(get(xid, envir = inverse_data$real2imag))
  H_imag <- as.character(get(Hid, envir = inverse_data$real2imag))

  primal <- list()
  primal[[xid]]    <- c(1, 3)
  primal[[x_imag]] <- c(2, 4)
  primal[[Hid]]    <- matrix(c(1, 2, 2, 3), 2, 2)
  primal[[H_imag]] <- 5    # compact strict upper triangle of the imaginary part

  solution <- Solution(status = OPTIMAL, opt_val = 0,
                       primal_vars = primal, dual_vars = list(), attr = list())
  inverted <- reduction_invert(reduction, solution, inverse_data)

  expect_equal(as.vector(inverted@primal_vars[[xid]]), c(1 + 2i, 3 + 4i))
  ## Skew-symmetric imaginary part: [[0, 5], [-5, 0]].
  expect_equal(inverted@primal_vars[[Hid]],
               matrix(c(1 + 0i, 2 - 5i, 2 + 5i, 3 + 0i), 2, 2))
})

test_that("a Hermitian leaf accepts a complex value", {
  ## NOT an upstream test -- an R-only regression found while porting the
  ## mapping-helpers test below, where `hparam.value = [[1, 2+5j], [2-5j, 3]]`
  ## simply would not run.
  ##
  ## leaf.py:318 defines is_complex() as complex || is_imag() || hermitian, and
  ## project() strips the imaginary part only when NOT is_complex(). CVXR's
  ## project() tested the raw `complex`/`imag` attributes and omitted
  ## `hermitian`, so it reached Re(val) first and the hermitian branch below it
  ## never saw a complex value. Assigning any genuinely complex Hermitian value
  ## aborted with "value must be real", and project() silently returned Re(val).
  H <- matrix(c(1, 2 - 5i, 2 + 5i, 3), 2, 2)   # [[1, 2+5i], [2-5i, 3]]

  for (leaf in list(Parameter(c(2, 2), hermitian = TRUE),
                    Variable(c(2, 2), hermitian = TRUE))) {
    expect_true(is_complex(leaf))
    expect_equal(project(leaf, H), H)
    value(leaf) <- H
    expect_equal(value(leaf), H)
  }

  ## A non-Hermitian input is still projected onto (M + Conj(t(M))) / 2, and
  ## the result stays complex.
  p <- Parameter(c(2, 2), hermitian = TRUE)
  M <- matrix(c(1, 0, 2 + 4i, 3), 2, 2)        # [[1, 2+4i], [0, 3]]
  expect_equal(project(p, M),
               matrix(c(1 + 0i, 1 - 2i, 1 + 2i, 3 + 0i), 2, 2))

  ## A real leaf must still have the imaginary part stripped.
  r <- Parameter(c(2, 2))
  expect_false(is_complex(r))
  expect_equal(project(r, H), matrix(c(1, 2, 2, 3), 2, 2))
})

## @cvxpy test_complex.py::TestComplex::test_complex2real_mapping_helpers_for_general_and_hermitian_leaves
test_that("Complex2Real splits and recombines general and Hermitian leaves", {
  ## Upstream reads the real/imag halves out of canon_methods._variables and
  ## ._parameters and checks their values. CVXR has no such dict -- the same
  ## information lives in the reduction's real2imag map -- so this test asserts
  ## the four chain-rule hooks instead, which is where those halves are
  ## actually observable, and checks the same numbers upstream does.
  general   <- Variable(2, complex = TRUE)
  hermitian <- Variable(c(2, 2), hermitian = TRUE)
  param     <- Parameter(2, complex = TRUE)
  hparam    <- Parameter(c(2, 2), hermitian = TRUE)
  value(param)  <- c(1 + 2i, 3 + 4i)
  value(hparam) <- matrix(c(1, 2 - 5i, 2 + 5i, 3), 2, 2)

  gid <- as.character(general@id)
  hid <- as.character(hermitian@id)
  pid <- as.character(param@id)
  qid <- as.character(hparam@id)

  reduction <- Complex2Real()

  ## Before apply() there is no mapping, so every hook is the identity and the
  ## id maps are empty -- upstream asserts exactly this.
  expect_equal(var_id_map(reduction), list())
  expect_equal(param_id_map(reduction), list())
  expect_equal(var_backward(reduction, stats::setNames(list(c(1 + 2i, 3 + 4i)), gid)),
               stats::setNames(list(c(1 + 2i, 3 + 4i)), gid))
  expect_equal(var_forward(reduction, stats::setNames(list(c(1, 2)), gid)),
               stats::setNames(list(c(1, 2)), gid))
  expect_equal(param_backward(reduction, stats::setNames(list(c(1, 2)), pid)),
               stats::setNames(list(c(1, 2)), pid))
  expect_equal(param_forward(reduction, stats::setNames(list(c(1, 2)), pid)),
               stats::setNames(list(c(1, 2)), pid))

  problem <- Problem(Minimize(sum_squares(Re(general - param))),
                     list(hermitian == hparam))
  reduction_apply(reduction, problem)
  update_parameters(reduction, problem)

  ## All four complex leaves are now mapped, plus the complex constraint.
  real2imag <- reduction@.cache$real2imag
  for (id in c(gid, hid, pid, qid))
    expect_true(exists(id, envir = real2imag, inherits = FALSE))
  expect_setequal(ls(reduction@.cache$id_to_var), c(gid, hid))
  expect_setequal(ls(reduction@.cache$id_to_param), c(pid, qid))

  g_imag <- as.character(get(gid, envir = real2imag))
  h_imag <- as.character(get(hid, envir = real2imag))
  p_imag <- as.character(get(pid, envir = real2imag))
  q_imag <- as.character(get(qid, envir = real2imag))

  ## var_id_map() stays empty even after apply(), unlike upstream. That is the
  ## shared-id design of ADR D_19.5 (#3147 part B deferred): the producers keep
  ## leaf ids rather than minting fresh ones, so there is nothing to remap. The
  ## consumer idiom still resolves correctly -- see chain.R:63-80.
  expect_equal(var_id_map(reduction), list())
  expect_equal(param_id_map(reduction), list())

  ## Upstream then reads real_param.value == [1, 3], imag_param.value == [2, 4],
  ## real_hparam.value == [[1, 2], [2, 3]], imag_hparam.value == [5]. In CVXR
  ## those four arrays ARE the output of param_forward, so the numbers below are
  ## upstream's numbers.
  split_params <- param_forward(reduction, stats::setNames(
    list(c(1 + 2i, 3 + 4i), matrix(c(1, 2 - 5i, 2 + 5i, 3), 2, 2)),
    c(pid, qid)))
  expect_equal(as.vector(split_params[[pid]]), c(1, 3))
  expect_equal(as.vector(split_params[[p_imag]]), c(2, 4))
  expect_equal(split_params[[qid]], matrix(c(1, 2, 2, 3), 2, 2))
  expect_equal(as.vector(split_params[[q_imag]]), 5)

  ## var_backward SPLITS variables (outer complex -> inner real/imag)...
  backward <- var_backward(reduction, stats::setNames(
    list(c(1 + 2i, 3 + 4i), matrix(c(1, 2 - 5i, 2 + 5i, 3), 2, 2)),
    c(gid, hid)))
  expect_equal(as.vector(backward[[gid]]), c(1, 3))
  expect_equal(as.vector(backward[[g_imag]]), c(2, 4))
  expect_equal(backward[[hid]], matrix(c(1, 2, 2, 3), 2, 2))
  expect_equal(as.vector(backward[[h_imag]]), 5)

  ## ...and var_forward COMBINES them again.
  forward <- var_forward(reduction, stats::setNames(
    list(c(1, 3), c(2, 4), matrix(c(1, 2, 2, 3), 2, 2), 5),
    c(gid, g_imag, hid, h_imag)))
  expect_equal(as.vector(forward[[gid]]), c(1 + 2i, 3 + 4i))
  expect_equal(forward[[hid]], matrix(c(1 + 0i, 2 - 5i, 2 + 5i, 3 + 0i), 2, 2))
  expect_null(forward[[g_imag]])   # the imaginary id is popped
  expect_null(forward[[h_imag]])

  ## Parameters run the other way round: backward combines, forward splits.
  param_bwd <- param_backward(reduction, stats::setNames(
    list(c(1, 3), c(2, 4), matrix(c(1, 2, 2, 3), 2, 2), 5),
    c(pid, p_imag, qid, q_imag)))
  expect_equal(as.vector(param_bwd[[pid]]), c(1 + 2i, 3 + 4i))
  expect_equal(param_bwd[[qid]], matrix(c(1 + 0i, 2 - 5i, 2 + 5i, 3 + 0i), 2, 2))
})

test_that("a complex Parameter refreshes across re-solves", {
  ## Upstream needs Complex2Real.update_parameters() to push a complex
  ## parameter's value into its real/imag halves before each solve. CVXR keeps
  ## the real part at the ORIGINAL parameter id, so the Reduction default
  ## no-op is correct -- but only behavior can show that, so: solve, change the
  ## value, solve again, and require the answer a freshly built problem gives.
  build <- function(pval) {
    p <- Parameter(2, complex = TRUE)
    value(p) <- pval
    x <- Variable(2, complex = TRUE)
    list(prob = Problem(Minimize(sum_squares(abs(x - p)))), x = x, p = p)
  }
  v1 <- c(1 + 2i, 3 + 4i)
  v2 <- c(-5 + 1i, 2 - 7i)

  fresh <- build(v2)
  psolve(fresh$prob, solver = "CLARABEL")
  reference <- as.vector(value(fresh$x))

  reused <- build(v1)
  psolve(reused$prob, solver = "CLARABEL")
  expect_equal(as.vector(value(reused$x)), v1, tolerance = 1e-6)
  value(reused$p) <- v2
  psolve(reused$prob, solver = "CLARABEL")
  expect_equal(as.vector(value(reused$x)), reference, tolerance = 1e-6)
  expect_equal(as.vector(value(reused$x)), v2, tolerance = 1e-6)
})

# ===================================================================
# ExpCone metadata, duals, validation        (test_constraints.py)
# ===================================================================

## @cvxpy test_constraints.py::TestConstraints::test_exponential_cone_metadata_duals_and_validation
test_that("ExpCone reports its metadata and dual cone correctly", {
  x <- Variable(2); y <- Variable(2); z <- Variable(2)
  con <- ExpCone(x, y, z)

  expect_match(expr_name(con), "^ExpCone\\(")
  expect_match(paste(utils::capture.output(print(con)), collapse = ""), "^ExpCone\\(")
  expect_null(residual(con))
  expect_equal(size(con), 6L)
  expect_equal(num_cones(con), 2L)
  expect_equal(cone_sizes(con), c(3L, 3L))
  expect_equal(con@shape, c(3L, 2L))
  expect_true(is_dcp(con))
  expect_false(is_dgp(con))
  expect_true(is_dqcp(con))
  ## Upstream spells the DPP question `con.is_dcp(dpp=True)`; CVXR has no `dpp`
  ## argument on is_dcp -- DPP is a scope (with_dpp_scope) entered by is_dpp().
  expect_true(is_dpp(con, "dcp"))

  ## PARTIAL, permanently: upstream also checks con.as_quad_approx(2, 3) ->
  ## RelEntrConeQuad. That constraint is NOT TARGETED in CVXR -- see the PARTIAL
  ## PORT sentinel at the top of constraints/exponential.R -- so the slice has no
  ## subject here and never will. A closed decision, not a pending port.

  value(x) <- c(1, 2); value(y) <- c(3, 4); value(z) <- c(5, 6)
  expect_true(S7_inherits(dual_cone(con), ExpCone))
  explicit <- dual_cone(con, x, y, z)
  expect_true(S7_inherits(explicit, ExpCone))
  ## value() keeps dimensions (a 2x1 matrix), so compare as vectors.
  expect_equal(as.vector(value(explicit@args[[1L]])), -as.vector(value(y)))
  expect_equal(as.vector(value(explicit@args[[2L]])), -as.vector(value(x)))
  expect_equal(as.vector(value(explicit@args[[3L]])), exp(1) * as.vector(value(z)))
  ## Upstream raises AssertionError on a shape mismatch; CVXR raises through
  ## the ExpCone constructor's own shape check.
  expect_error(dual_cone(con, Variable(1), y, z), "same shapes")
})

## @cvxpy test_constraints.py::TestConstraints::test_exponential_cone_metadata_duals_and_validation
test_that("ExpCone saves duals in argument order and rejects bad arguments", {
  x <- Variable(2); y <- Variable(2); z <- Variable(2)
  con <- ExpCone(x, y, z)
  save_dual_value(con, 0:5)
  expect_equal(as.vector(value(con@dual_variables[[1L]])), c(0, 3))
  expect_equal(as.vector(value(con@dual_variables[[2L]])), c(1, 4))
  expect_equal(as.vector(value(con@dual_variables[[3L]])), c(2, 5))

  expect_error(ExpCone(square(x), y, z), "affine and real")
  expect_error(ExpCone(Variable(1), y, z), "same shapes")
})

# ===================================================================
# LinOp constructors and helpers                 (test_lin_ops.py)
# ===================================================================

## @cvxpy test_lin_ops.py::test_lin_ops::test_parameter_and_replacement_helpers
test_that("the LinOp constructors build the operators CVXPY's lin_utils does", {
  ## PARTIAL PORT. CVXR's lin_utils.R carries 24 of upstream's 38 helpers; the
  ## other 14 have no counterpart because CVXR assembles and rewrites the
  ## constraint/parameter graph at the Expression level (eval_params.R,
  ## canonicalization.R) or in the C++ backend, never on R-level LinOps:
  ##   sub_expr, promote_lin_ops_for_mul, broadcast_to, concatenate,
  ##   get_constr_expr, create_eq, create_leq, create_geq, get_expr_vars,
  ##   get_expr_params, copy_constr, replace_new_vars, check_param_val,
  ##   replace_params_with_consts (an Expression-level `.replace_params_with_consts`
  ##   in reductions/eval_params.R plays the last one's role).
  ## Upstream's replacement half of this test therefore has no subject here.
  ## CVXR also suffixes the builders `_linop` to keep them clear of the atom
  ## namespace (`trace`, `index`, `conv`, ... are all taken).
  x <- create_var(c(2L, 2L), var_id = 10L)
  y <- create_var(c(2L, 2L), var_id = 11L)
  param <- Parameter(c(2, 2))
  value(param) <- matrix(0:3, 2, 2)

  expect_equal(x$type, LINOP_VARIABLE)
  expect_equal(x$shape, c(2L, 2L))
  expect_equal(x$data, 10L)
  expect_true(is_scalar_linop(create_var(c(1L, 1L), 12L)))
  expect_true(is_const_linop(create_const(1, c(1L, 1L))))

  ## create_param takes either the Parameter or just its id.
  p_obj <- create_param(c(2L, 2L), param)
  expect_equal(p_obj$type, LINOP_PARAM)
  expect_true(S7_inherits(p_obj$data, Parameter))
  p_op <- create_param(c(2L, 2L), param@id)
  expect_equal(p_op$type, LINOP_PARAM)
  expect_equal(p_op$data, param@id)

  expect_equal(sum_expr_linop(list(x, y))$type, "sum")
  expect_equal(neg_expr_linop(x)$type, "neg")
  expect_equal(mul_expr_linop(p_op, x, c(2L, 2L))$type, "mul_expr")
  expect_equal(rmul_expr_linop(x, p_op, c(2L, 2L))$type, "rmul_expr")
  expect_equal(multiply_linop(x, y)$type, "mul_elem")
  expect_equal(kron_r_linop(p_op, x, c(4L, 4L))$type, "kron")
  expect_equal(kron_l_linop(x, p_op, c(4L, 4L))$type, "kron_l")
  expect_equal(div_expr_linop(x, create_const(2, c(1L, 1L)))$type, "div")
  expect_equal(promote_linop(create_const(1, c(1L, 1L)), c(2L, 2L))$type, "promote")
  expect_equal(sum_entries_linop(x, c(1L, 1L))$type, "sum_entries")
  expect_equal(trace_linop(x)$type, "trace")
  expect_equal(index_linop(x, c(1L, 2L), list(1L, NULL))$type, "index")
  expect_equal(conv_linop(p_op, x, c(3L, 1L))$type, "conv")
  expect_equal(reshape_linop(x, c(4L, 1L))$type, "reshape_expr")
  expect_equal(hstack_linop(list(x, y), c(2L, 4L))$type, "hstack")
  expect_equal(vstack_linop(list(x, y), c(4L, 2L))$type, "vstack")

  ## Shapes, where upstream checks a shape rather than a type. CVXR has no 1-D
  ## shapes -- upstream's transpose(create_var((3,))) -> (3,) is the degenerate
  ## case of (3, 1) -> (1, 3) here.
  expect_equal(transpose_linop(create_var(c(3L, 1L), 13L))$shape, c(1L, 3L))
  expect_equal(transpose_linop(x)$shape, c(2L, 2L))
  expect_equal(diag_vec_linop(create_var(c(2L, 1L), 14L))$shape, c(2L, 2L))
  expect_equal(diag_mat_linop(x)$shape, c(2L, 1L))
  expect_equal(upper_tri_linop(x)$shape, c(1L, 1L))
})

# ===================================================================
# power_tools edge branches                   (test_power_tools.py)
# ===================================================================

## @cvxpy test_power_tools.py::test_power_tools_edge_branches
test_that("the power_tools helpers handle their edge branches", {
  ## gmp is never attached -- it masks %*%, matrix and apply.
  Q <- function(n, d = 1L) gmp::as.bigq(n, d)

  ## A single unit weight collapses to t == x.
  t_var <- Variable(); xv <- Variable()
  constraints <- gm_constrs(t_var, list(xv), Q(1))
  expect_length(constraints, 1L)
  expect_true(S7_inherits(constraints[[1L]], Equality))
  expect_identical(constraints[[1L]]@args[[1L]], t_var)
  expect_identical(constraints[[1L]]@args[[2L]], xv)

  ## Power-splitting helpers.
  high <- pow_high(2, approx = FALSE)
  expect_equal(as.numeric(high[[1L]]), 2)
  expect_equal(as.numeric(high[[2L]]), c(0.5, 0.5))
  expect_equal(as.character(pow_high(Q(3, 2))[[1L]]), "3/2")
  neg <- pow_neg(-2, approx = FALSE)
  expect_equal(as.numeric(neg[[1L]]), -2)
  expect_equal(as.numeric(neg[[2L]]), c(2 / 3, 1 / 3), tolerance = 1e-12)
  expect_equal(as.numeric(pow_neg(-2)[[1L]]), -2)
  mid <- pow_mid(Q(1, 2), approx = FALSE)
  expect_equal(as.character(mid[[1L]]), "1/2")

  ## Predicates.
  expect_true(is_dyad(Q(2)))
  expect_false(is_dyad(Q(1, 3)))
  expect_true(is_weight(Q(c(0, 0, 1))))
  expect_true(is_power2(8))
  expect_false(is_power2(6))
  expect_equal(next_pow2(0), 1L)
  expect_equal(get_max_denom(Q(c(1, 2), c(3, 3))), 3L)

  ## fracify rejects its three bad inputs.
  expect_error(fracify(c(-1, 2)), "nonnegative")
  expect_error(fracify(c(1, 2), max_denom = 0L), "max_denom")
  expect_error(fracify(Q(c(1, 2), c(3, 3)), max_denom = 2L), "reliably represent")

  ## ...and returns a weight plus its dyadic completion.
  fr <- fracify(c(0.1, 0.9), max_denom = 16L)
  expect_equal(as.character(fr[[1L]]), c("1/10", "9/10"))
  expect_equal(as.character(fr[[2L]]), c("1/16", "9/16", "3/8"))
  expect_true(is_weight(fr[[1L]]))
  expect_true(is_dyad_weight(fr[[2L]]))

  expect_equal(as.character(make_frac(c(0.2, 0.8), 4L)), c("1/4", "3/4"))
  expect_equal(as.character(dyad_completion(Q(c(1, 2), c(3, 3)))),
               c("1/4", "1/2", "1/4"))
  expect_equal(approx_error(c(1, 1), Q(c(1, 1), c(2, 2))), 0)

  half <- Q(c(1, 1), c(2, 2))
  expect_true(check_dyad(half, half))
  expect_false(check_dyad(Q(c(2, 1), c(3, 1)), Q(c(2, 1), c(3, 1))))
  expect_length(.split_dyad(Q(c(1, 0))), 0L)

  ## Decomposition and its bounds / pretty-printing.
  tree <- decompose(half)
  expect_type(tree, "list")
  expect_equal(over_bound(half, tree), 0)
  expect_equal(lower_bound(half), 1)
  expect_equal(prettytuple(half), "(1/2, 1/2)")
  expect_match(prettydict(tree), "(1/2, 1/2)", fixed = TRUE)
})
