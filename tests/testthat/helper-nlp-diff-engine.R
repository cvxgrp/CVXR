## R-SPECIFIC: finite-difference derivative checker for the NLP diff engine.
##
## Mirrors nlp_tests/derivative_checker.DerivativeChecker.run_and_assert: build
## the diff-engine C_problem for a (DNLP) Problem, then at the problem's initial
## point compare the engine's analytic objective gradient, constraint Jacobian,
## and Lagrangian Hessian against central finite differences of objective_forward
## / constraint_forward. CVXPY's checker applies Dnlp2Smooth first and evaluates
## at the Bounds x0 (the variable values, kept in-domain by each test); we do the
## same. Random Lagrange multipliers are used for the Hessian, exactly as in
## DerivativeChecker.check_hessian (`duals = np.random.rand(len(cl))`).

## Some CVXPY NLP tests sweep many sizes/factors (e.g. abs lasso x 20 factors,
## ML_Gaussian x 20 sizes). By default we run a small representative subset to
## keep R CMD check fast; set CVXR_NLP_FULL=1 to run CVXPY's full loop counts
## (e.g. an overnight `CVXR_NLP_FULL=1 R CMD check`). Only the number of points
## changes -- the assertions are identical.
.nlp_full <- function() nzchar(Sys.getenv("CVXR_NLP_FULL"))

## Is an R NLP backend available to actually solve a DNLP?
.de_have_backend <- function() {
  requireNamespace("sparsediff", quietly = TRUE) &&
    (requireNamespace("ipopt", quietly = TRUE) ||
       requireNamespace("Uno", quietly = TRUE))
}

## Solve a DNLP, skipping if no backend. `print_level = 0` quiets IPOPT's
## per-iteration dump (ignored by other backends).
.de_nlp_solve <- function(prob, ...) {
  if (!.de_have_backend()) skip("No R NLP backend (sparsediff + ipopt/Uno)")
  invisible(capture.output(psolve(prob, nlp = TRUE, print_level = 0L, ...),
                           type = "output"))
  prob
}

## Solve a DNLP specifically with the IPOPT backend (the robust interior-point
## reference; several CVXPY NLP examples are skipped for UNO upstream). Skips if
## the ipopt R package is absent. Quiets IPOPT's per-iteration dump.
.de_ipopt_solve <- function(prob, ...) {
  if (!requireNamespace("sparsediff", quietly = TRUE) ||
        !requireNamespace("ipopt", quietly = TRUE)) {
    skip("IPOPT NLP backend (sparsediff + ipopt) not installed")
  }
  invisible(capture.output(
    psolve(prob, nlp = TRUE, solver = "IPOPT", print_level = 0L, ...),
    type = "output"))
  prob
}

## Build x0 (column-major flatten of variable values) in the diff engine's
## variable order, via the same Bounds machinery the NLP solver uses.
.nlp_x0 <- function(problem) {
  Bounds <- get("Bounds", envir = asNamespace("CVXR"))
  Bounds(problem)@x0
}

## Full derivative check at one or more in-domain points. By default the single
## point is the problem's x0 (mirrors DerivativeChecker, which checks at x0).
nlp_derivative_check <- function(problem, x = NULL, seed = 0L,
                                 grad_tol = 1e-4, jac_tol = 1e-4,
                                 hess_tol = 1e-4, eps = 1e-6,
                                 canonicalize = TRUE) {
  skip_if_not_installed("sparsediff")
  ns <- asNamespace("CVXR")
  de_C_problem <- get(".de_C_problem", envir = ns)
  if (canonicalize) {
    Dnlp2Smooth <- get("Dnlp2Smooth", envir = ns)
    reduction_apply <- get("reduction_apply", envir = ns)
    canon <- reduction_apply(Dnlp2Smooth(), problem)
    problem <- canon[[1L]]
  }
  if (is.null(x)) x <- .nlp_x0(problem)

  cp <- de_C_problem(problem)
  cp$init_jacobian_coo()
  cp$init_hessian_coo_lower_tri()
  nv <- cp$n_vars
  js <- cp$jacobian_sparsity()      # rows/cols 0-based, nrow = #scalar constraints
  hs <- cp$hessian_sparsity()       # rows/cols 0-based, lower triangle
  ncons <- js$nrow

  set.seed(seed)
  duals <- if (ncons > 0L) stats::rnorm(ncons) else numeric(0)

  lagrangian_gradient <- function(xe) {
    cp$objective_forward(xe)
    gf <- as.numeric(cp$gradient())
    if (ncons > 0L) {
      cp$constraint_forward(xe)
      jv <- as.numeric(cp$jacobian_values())
      J <- matrix(0, ncons, nv)
      J[cbind(js$rows + 1L, js$cols + 1L)] <- jv
      gf + as.numeric(t(J) %*% duals)
    } else {
      gf
    }
  }

  ## --- objective gradient ---
  cp$objective_forward(x)
  g <- as.numeric(cp$gradient())
  gfd <- numeric(nv)
  for (j in seq_len(nv)) {
    xp <- x; xp[j] <- xp[j] + eps
    xm <- x; xm[j] <- xm[j] - eps
    gfd[j] <- (cp$objective_forward(xp) - cp$objective_forward(xm)) / (2 * eps)
  }
  expect_equal(g, gfd, tolerance = grad_tol)

  ## --- constraint Jacobian ---
  if (ncons > 0L) {
    cp$constraint_forward(x)
    jv <- as.numeric(cp$jacobian_values())
    Jc <- matrix(0, ncons, nv)
    Jc[cbind(js$rows + 1L, js$cols + 1L)] <- jv
    Jn <- matrix(0, ncons, nv)
    for (j in seq_len(nv)) {
      xp <- x; xp[j] <- xp[j] + eps
      xm <- x; xm[j] <- xm[j] - eps
      Jn[, j] <- (cp$constraint_forward(xp) - cp$constraint_forward(xm)) / (2 * eps)
    }
    expect_equal(Jc, Jn, tolerance = jac_tol)
  }

  ## --- Lagrangian Hessian (obj_factor = 1, random duals) ---
  cp$objective_forward(x)
  cp$gradient()
  if (ncons > 0L) cp$constraint_forward(x)
  hv <- as.numeric(cp$hessian_values(1, duals))
  Hc <- matrix(0, nv, nv)
  Hc[cbind(hs$rows + 1L, hs$cols + 1L)] <- hv
  off <- hs$rows != hs$cols
  if (any(off)) Hc[cbind(hs$cols[off] + 1L, hs$rows[off] + 1L)] <- hv[off]
  Hn <- matrix(0, nv, nv)
  for (j in seq_len(nv)) {
    xp <- x; xp[j] <- xp[j] + eps
    xm <- x; xm[j] <- xm[j] - eps
    Hn[, j] <- (lagrangian_gradient(xp) - lagrangian_gradient(xm)) / (2 * eps)
  }
  Hn <- (Hn + t(Hn)) / 2
  expect_equal(Hc, Hn, tolerance = hess_tol)

  invisible(NULL)
}
