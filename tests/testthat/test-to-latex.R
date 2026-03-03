## Tests for to_latex() — LaTeX rendering of CVXR objects

# ═══════════════════════════════════════════════════════════════════
# Name conversion helpers
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that(".name_to_latex: single letters pass through", {
  expect_equal(.name_to_latex("x"), "x")
  expect_equal(.name_to_latex("A"), "A")
})

## @cvxpy NONE
test_that(".name_to_latex: Greek letters", {
  expect_equal(.name_to_latex("gamma"), "\\gamma")
  expect_equal(.name_to_latex("Sigma"), "\\Sigma")
  expect_equal(.name_to_latex("alpha"), "\\alpha")
  expect_equal(.name_to_latex("lambda"), "\\lambda")
})

## @cvxpy NONE
test_that(".name_to_latex: Greek with trailing digits", {
  expect_equal(.name_to_latex("beta1"), "\\beta_{1}")
  expect_equal(.name_to_latex("gamma42"), "\\gamma_{42}")
})

## @cvxpy NONE
test_that(".name_to_latex: suffix patterns", {
  expect_equal(.name_to_latex("x_hat"), "\\hat{x}")
  expect_equal(.name_to_latex("x_bar"), "\\bar{x}")
  expect_equal(.name_to_latex("x_tilde"), "\\tilde{x}")
})

## @cvxpy NONE
test_that(".name_to_latex: multi-char names", {
  expect_equal(.name_to_latex("myvar"), "\\mathit{myvar}")
  expect_equal(.name_to_latex("var42"), "\\mathit{var42}")
})


# ═══════════════════════════════════════════════════════════════════
# Leaf rendering
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: named variable", {
  x <- Variable(3, name = "x")
  expect_equal(to_latex(x), "x")
})

## @cvxpy NONE
test_that("to_latex: variable with latex_name", {
  x <- Variable(3, name = "x", latex_name = "\\mathbf{x}")
  expect_equal(to_latex(x), "\\mathbf{x}")
})

## @cvxpy NONE
test_that("to_latex: Greek-named variable", {
  g <- Variable(name = "gamma")
  expect_equal(to_latex(g), "\\gamma")
})

## @cvxpy NONE
test_that("to_latex: auto-named variable uses mathit", {
  v <- Variable(2)
  ltx <- to_latex(v)
  expect_match(ltx, "\\\\mathit\\{var")
})

## @cvxpy NONE
test_that("to_latex: named parameter", {
  p <- Parameter(name = "alpha")
  expect_equal(to_latex(p), "\\alpha")
})

## @cvxpy NONE
test_that("to_latex: parameter with latex_name", {
  p <- Parameter(name = "p", latex_name = "\\hat{p}")
  expect_equal(to_latex(p), "\\hat{p}")
})

## @cvxpy NONE
test_that("to_latex: scalar constant", {
  c1 <- Constant(5)
  expect_equal(to_latex(c1), "5")
  c2 <- Constant(-3.14)
  expect_match(to_latex(c2), "-3")
})

## @cvxpy NONE
test_that("to_latex: matrix constant (abbreviated)", {
  c_mat <- Constant(matrix(1:6, 2, 3))
  expect_match(to_latex(c_mat), "2x3")
})


# ═══════════════════════════════════════════════════════════════════
# Affine operators
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: addition", {
  x <- Variable(name = "x")
  y <- Variable(name = "y")
  expr <- x + y
  expect_equal(to_latex(expr), "x + y")
})

## @cvxpy NONE
test_that("to_latex: subtraction renders as minus", {
  x <- Variable(name = "x")
  y <- Variable(name = "y")
  expr <- x - y
  expect_equal(to_latex(expr), "x - y")
})

## @cvxpy NONE
test_that("to_latex: negation", {
  x <- Variable(name = "x")
  expr <- -x
  expect_equal(to_latex(expr), "-x")
})

## @cvxpy NONE
test_that("to_latex: scalar multiply uses juxtaposition, not circ", {
  x <- Variable(3, name = "x")
  expr <- 2 * x
  ltx <- to_latex(expr)
  expect_match(ltx, "2 x")
  expect_no_match(ltx, "\\\\circ")
})

## @cvxpy NONE
test_that("to_latex: Hadamard product uses circ for non-scalar operands", {
  x <- Variable(3, name = "x")
  y <- Variable(3, name = "y")
  expr <- Multiply(x, y)
  ltx <- to_latex(expr)
  expect_match(ltx, "\\\\circ")
})

## @cvxpy NONE
test_that("to_latex: matrix multiply", {
  x <- Variable(3, name = "x")
  A <- Parameter(c(5, 3), name = "A")
  expr <- A %*% x
  ltx <- to_latex(expr)
  expect_match(ltx, "A.*x")
})

## @cvxpy NONE
test_that("to_latex: division", {
  x <- Variable(name = "x")
  expr <- x / 2
  expect_match(to_latex(expr), "\\\\frac")
})

## @cvxpy NONE
test_that("to_latex: transpose", {
  x <- Variable(3, name = "x")
  expr <- t(x)
  expect_match(to_latex(expr), "x\\\\T")
})

## @cvxpy NONE
test_that("to_latex: indexing", {
  x <- Variable(5, name = "x")
  expr <- x[2]
  ltx <- to_latex(expr)
  expect_match(ltx, "x_\\{2")
})

## @cvxpy NONE
test_that("to_latex: sum_entries of column vector", {
  x <- Variable(3, name = "x")
  expr <- sum_entries(x)
  ltx <- to_latex(expr)
  ## Column vector: 1^T x only, no trailing 1
  expect_equal(ltx, "\\ones\\T x")
})

## @cvxpy NONE
test_that("to_latex: sum_entries of matrix", {
  X <- Variable(c(3, 4), name = "X")
  expr <- sum_entries(X)
  ltx <- to_latex(expr)
  ## Matrix: 1^T X 1
  expect_equal(ltx, "\\ones\\T X \\ones")
})

## @cvxpy NONE
test_that("to_latex: trace", {
  X <- Variable(c(3, 3), name = "X")
  expr <- matrix_trace(X)
  ltx <- to_latex(expr)
  expect_match(ltx, "\\\\tr")
})

## @cvxpy NONE
test_that("to_latex: kron", {
  A <- Constant(matrix(1:4, 2, 2))
  B <- Variable(c(2, 2), name = "B")
  expr <- kron(A, B)
  expect_match(to_latex(expr), "\\\\otimes")
})

## @cvxpy NONE
test_that("to_latex: hstack", {
  x <- Variable(c(2, 1), name = "x")
  y <- Variable(c(2, 1), name = "y")
  expr <- hstack(x, y)
  expect_match(to_latex(expr), "bmatrix")
})

## @cvxpy NONE
test_that("to_latex: vstack", {
  x <- Variable(c(1, 2), name = "x")
  y <- Variable(c(1, 2), name = "y")
  expr <- vstack(x, y)
  expect_match(to_latex(expr), "bmatrix")
})

## @cvxpy NONE
test_that("to_latex: diag (vector to matrix)", {
  x <- Variable(3, name = "x")
  expr <- DiagVec(x)
  expect_match(to_latex(expr), "\\\\diag")
})

## @cvxpy NONE
test_that("to_latex: reshape", {
  x <- Variable(6, name = "x")
  expr <- reshape_expr(x, c(2, 3))
  expect_match(to_latex(expr), "reshape")
})


# ═══════════════════════════════════════════════════════════════════
# Elementwise atoms
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: abs", {
  x <- Variable(name = "x")
  expr <- abs(x)
  expect_match(to_latex(expr), "\\\\cvxabs")
})

## @cvxpy NONE
test_that("to_latex: power (square)", {
  x <- Variable(name = "x")
  expr <- x^2
  expect_match(to_latex(expr), "\\^\\{2\\}")
})

## @cvxpy NONE
test_that("to_latex: power (sqrt)", {
  x <- Variable(name = "x", nonneg = TRUE)
  expr <- sqrt(x)
  expect_match(to_latex(expr), "\\\\sqrt")
})

## @cvxpy NONE
test_that("to_latex: power (inverse)", {
  x <- Variable(name = "x", nonneg = TRUE)
  expr <- Power(x, -1)
  expect_match(to_latex(expr), "\\\\frac\\{1\\}")
})

## @cvxpy NONE
test_that("to_latex: exp", {
  x <- Variable(name = "x")
  expr <- exp(x)
  expect_match(to_latex(expr), "e\\^")
})

## @cvxpy NONE
test_that("to_latex: log", {
  x <- Variable(name = "x", nonneg = TRUE)
  expr <- log(x)
  expect_match(to_latex(expr), "\\\\log")
})

## @cvxpy NONE
test_that("to_latex: log1p", {
  x <- Variable(name = "x")
  expr <- log1p_expr(x)
  expect_match(to_latex(expr), "\\\\log.*1 \\+")
})

## @cvxpy NONE
test_that("to_latex: maximum", {
  x <- Variable(name = "x")
  y <- Variable(name = "y")
  expr <- max_elemwise(x, y)
  expect_match(to_latex(expr), "\\\\max")
})

## @cvxpy NONE
test_that("to_latex: minimum", {
  x <- Variable(name = "x")
  y <- Variable(name = "y")
  expr <- min_elemwise(x, y)
  expect_match(to_latex(expr), "\\\\min")
})

## @cvxpy NONE
test_that("to_latex: entropy", {
  x <- Variable(name = "x", nonneg = TRUE)
  expr <- Entr(x)
  ltx <- to_latex(expr)
  expect_match(ltx, "-")
  expect_match(ltx, "\\\\log")
})

## @cvxpy NONE
test_that("to_latex: rel_entr", {
  x <- Variable(name = "x", nonneg = TRUE)
  y <- Variable(name = "y", nonneg = TRUE)
  expr <- RelEntr(x, y)
  ltx <- to_latex(expr)
  expect_match(ltx, "\\\\log")
  expect_match(ltx, "\\\\frac")
})

## @cvxpy NONE
test_that("to_latex: kl_div", {
  x <- Variable(name = "x", nonneg = TRUE)
  y <- Variable(name = "y", nonneg = TRUE)
  expr <- KlDiv(x, y)
  ltx <- to_latex(expr)
  expect_match(ltx, "\\\\log")
})

## @cvxpy NONE
test_that("to_latex: logistic", {
  x <- Variable(name = "x")
  expr <- Logistic(x)
  ltx <- to_latex(expr)
  expect_match(ltx, "\\\\log.*1 \\+ e")
})

## @cvxpy NONE
test_that("to_latex: huber", {
  x <- Variable(name = "x")
  expr <- Huber(x, M = 2)
  ltx <- to_latex(expr)
  expect_match(ltx, "\\\\huber")
})

## @cvxpy NONE
test_that("to_latex: ceil and floor", {
  x <- Variable(name = "x")
  expect_match(to_latex(Ceil(x)), "\\\\cvxceil")
  expect_match(to_latex(Floor(x)), "\\\\cvxfloor")
})


# ═══════════════════════════════════════════════════════════════════
# Norm atoms
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: p_norm (L2)", {
  x <- Variable(3, name = "x")
  expr <- p_norm(x, 2)
  expect_match(to_latex(expr), "\\\\cvxnorm\\{x\\}_\\{2\\}")
})

## @cvxpy NONE
test_that("to_latex: norm1", {
  x <- Variable(3, name = "x")
  expr <- norm1(x)
  expect_match(to_latex(expr), "\\\\cvxnorm\\{x\\}_1")
})

## @cvxpy NONE
test_that("to_latex: norm_inf", {
  x <- Variable(3, name = "x")
  expr <- norm_inf(x)
  expect_match(to_latex(expr), "\\\\cvxnorm\\{x\\}_\\\\infty")
})

## @cvxpy NONE
test_that("to_latex: norm_nuc", {
  X <- Variable(c(3, 3), name = "X")
  expr <- norm_nuc(X)
  expect_match(to_latex(expr), "\\\\cvxnorm.*\\*")
})

## @cvxpy NONE
test_that("to_latex: sigma_max", {
  X <- Variable(c(3, 3), name = "X")
  expr <- sigma_max(X)
  expect_match(to_latex(expr), "\\\\sigmamax")
})


# ═══════════════════════════════════════════════════════════════════
# Reduction atoms
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: max_entries", {
  x <- Variable(5, name = "x")
  expr <- max_entries(x)
  expect_match(to_latex(expr), "\\\\max")
})

## @cvxpy NONE
test_that("to_latex: min_entries", {
  x <- Variable(5, name = "x")
  expr <- min_entries(x)
  expect_match(to_latex(expr), "\\\\min")
})

## @cvxpy NONE
test_that("to_latex: log_sum_exp", {
  x <- Variable(5, name = "x")
  expr <- LogSumExp(x)
  expect_match(to_latex(expr), "\\\\logsumexp")
})


# ═══════════════════════════════════════════════════════════════════
# Matrix / spectral atoms
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: quad_form", {
  x <- Variable(3, name = "x")
  P <- Parameter(c(3, 3), name = "P", PSD = TRUE)
  expr <- quad_form(x, P)
  ltx <- to_latex(expr)
  expect_match(ltx, "x\\\\T")
  expect_match(ltx, "P")
})

## @cvxpy NONE
test_that("to_latex: quad_over_lin", {
  x <- Variable(3, name = "x")
  y <- Variable(name = "y", nonneg = TRUE)
  expr <- quad_over_lin(x, y)
  ltx <- to_latex(expr)
  expect_match(ltx, "\\\\frac")
  expect_match(ltx, "\\\\cvxnorm")
})

## @cvxpy NONE
test_that("to_latex: log_det", {
  X <- Variable(c(3, 3), name = "X", PSD = TRUE)
  expr <- log_det(X)
  expect_match(to_latex(expr), "\\\\logdet")
})

## @cvxpy NONE
test_that("to_latex: lambda_max", {
  X <- Variable(c(3, 3), name = "X", symmetric = TRUE)
  expr <- lambda_max(X)
  expect_match(to_latex(expr), "\\\\lambdamax")
})

## @cvxpy NONE
test_that("to_latex: lambda_sum_largest", {
  X <- Variable(c(3, 3), name = "X", symmetric = TRUE)
  expr <- lambda_sum_largest(X, 2)
  ltx <- to_latex(expr)
  expect_match(ltx, "\\\\lambda")
  expect_match(ltx, "2")
})


# ═══════════════════════════════════════════════════════════════════
# Constraints
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: inequality (x >= 0 flips to x >= 0)", {
  x <- Variable(3, name = "x")
  constr <- x >= 0  # R dispatches as Inequality(0, x); we flip to x >= 0
  ltx <- to_latex(constr)
  expect_match(ltx, "x \\\\geq 0")
})

## @cvxpy NONE
test_that("to_latex: inequality (x <= 0 stays x <= 0)", {
  x <- Variable(3, name = "x")
  constr <- x <= 0
  ltx <- to_latex(constr)
  expect_match(ltx, "x \\\\leq 0")
})

## @cvxpy NONE
test_that("to_latex: Inequality (lhs <= rhs)", {
  x <- Variable(3, name = "x")
  b <- Parameter(3, name = "b")
  constr <- (x <= b)
  ltx <- to_latex(constr)
  expect_match(ltx, "\\\\leq")
})

## @cvxpy NONE
test_that("to_latex: equality", {
  x <- Variable(3, name = "x")
  constr <- (sum_entries(x) == 1)
  ltx <- to_latex(constr)
  expect_match(ltx, "=")
})

## @cvxpy NONE
test_that("to_latex: PSD constraint (X >> 0) decomposes cleanly", {
  X <- Variable(c(3, 3), name = "X", symmetric = TRUE)
  constr <- X %>>% 0
  ltx <- to_latex(constr)
  ## Should render as "X \psd 0" without subtraction artifacts
  expect_equal(ltx, "X \\psd 0")
})

## @cvxpy NONE
test_that("to_latex: PSD constraint (X >> Y) decomposes to X psd Y", {
  X <- Variable(c(3, 3), name = "X", symmetric = TRUE)
  Y <- Variable(c(3, 3), name = "Y", symmetric = TRUE)
  constr <- X %>>% Y
  ltx <- to_latex(constr)
  expect_equal(ltx, "X \\psd Y")
})


# ═══════════════════════════════════════════════════════════════════
# Problem-level rendering
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: simple minimization problem", {
  x <- Variable(3, name = "x")
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  ltx <- to_latex(prob)
  expect_match(ltx, "\\\\begin\\{mini\\*\\}")
  expect_match(ltx, "\\\\end\\{mini\\*\\}")
  expect_match(ltx, "\\\\addConstraint")
})

## @cvxpy NONE
test_that("to_latex: maximization problem", {
  x <- Variable(3, name = "x")
  prob <- Problem(Maximize(sum_entries(x)), list(x <= 1))
  ltx <- to_latex(prob)
  expect_match(ltx, "\\\\begin\\{maxi\\*\\}")
})

## @cvxpy NONE
test_that("to_latex: unconstrained problem", {
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x^2))
  ltx <- to_latex(prob)
  expect_match(ltx, "\\\\begin\\{mini\\*\\}")
  expect_no_match(ltx, "addConstraint")
})

## @cvxpy NONE
test_that("to_latex: problem with multiple constraints", {
  x <- Variable(3, name = "x")
  A <- Parameter(c(5, 3), name = "A")
  b <- Parameter(5, name = "b")
  gamma <- Parameter(name = "gamma", nonneg = TRUE)
  prob <- Problem(
    Minimize(p_norm(x, 2) + gamma * sum_entries(x)),
    list(A %*% x <= b, x >= 0)
  )
  ltx <- to_latex(prob)
  expect_match(ltx, "mini\\*")
  expect_match(ltx, "\\\\cvxnorm")
  expect_match(ltx, "\\\\gamma")
  expect_match(ltx, "\\\\addConstraint")
})


# ═══════════════════════════════════════════════════════════════════
# Name collision handling
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: name collision warning", {
  x1 <- Variable(2, name = "x")
  x2 <- Variable(3, name = "x")
  prob <- Problem(Minimize(sum_entries(x1) + sum_entries(x2)), list(x1 >= 0, x2 >= 0))
  expect_warning(to_latex(prob), "disambiguating")
})


# ═══════════════════════════════════════════════════════════════════
# Atom fallback
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: fallback for unknown atom uses operatorname", {
  ## ConditionNumber is a class with no specific pretty LaTeX;
  ## but we DID register a method for it, so test with Dotsort instead
  x <- Variable(3, name = "x")
  w <- Constant(c(3, 2, 1))
  expr <- Dotsort(x, w)
  ltx <- to_latex(expr)
  expect_match(ltx, "dotsort")
})


# ═══════════════════════════════════════════════════════════════════
# Logic atoms
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: logic atoms", {
  x <- Variable(name = "x", boolean = TRUE)
  y <- Variable(name = "y", boolean = TRUE)

  expect_match(to_latex(Not(x)), "\\\\neg")
  expect_match(to_latex(And(x, y)), "\\\\land")
  expect_match(to_latex(Or(x, y)), "\\\\lor")
  expect_match(to_latex(Xor(x, y)), "\\\\oplus")
})


# ═══════════════════════════════════════════════════════════════════
# Complex atoms
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("to_latex: conj, real, imag", {
  z <- Variable(name = "z", complex = TRUE)
  expect_match(to_latex(Conj_(z)), "\\\\overline")
  expect_match(to_latex(Real_(z)), "\\\\operatorname\\{Re\\}")
  expect_match(to_latex(Imag_(z)), "\\\\operatorname\\{Im\\}")
})


# ═══════════════════════════════════════════════════════════════════
# dcp.sty file exists
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("dcp.sty is accessible via system.file", {
  sty <- system.file("sty", "dcp.sty", package = "CVXR")
  expect_true(nzchar(sty))
  expect_true(file.exists(sty))
  contents <- readLines(sty)
  expect_true(any(grepl("optidef", contents)))
  expect_true(any(grepl("DeclarePairedDelimiter", contents)))
})
