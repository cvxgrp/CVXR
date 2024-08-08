context("test_perspective")

quad_example <- function(r_val) {
  # Reference expected output.
  x <- Variable()
  s <- Variable()
  r <- r_val
  
  obj <- quad_over_lin(x, s) + r*x - 4*s
  constraints <- list(x >= 2, s <= 0.5)
  prob_ref <- Problem(Minimize(obj), constraints)
  result <- solve(prob_ref, solver = "ECOS")
  
  return(list(result$getValue(prob_ref), result$getValue(s), result$getValue(x), r))
}

test_p_norms <- function(p) {
  x <- Variable(3)
  s <- Variable(nonneg = TRUE, name = "s")
  f <- p_norm(x, p)
  obj <- perspective(f, s)
  constraints <- list(1 == s, x >= c(1,2,3))
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "ECOS")
  
  # Reference problem.
  ref_x <- Variable(3, pos = TRUE)
  ref_s <- Variable(pos = TRUE)
  
  obj <- sum(power(ref_x, p) / power(ref_s, p-1))
  
  ref_constraints <- list(ref_x >= c(1,2,3), ref_s == 1)
  ref_prob <- Problem(Minimize(obj), ref_constraints)
  result_ref <- solve(ref_prob, gp = TRUE)
  
  expect_true(np.isclose(result$value^p, result_ref$value))
  expect_true(is.allclose(result$getValue(x), result_ref$getValue(ref_x)))
  if(p != 1)   # s is not used when denominator is s^0.
    expect_true(np.isclose(result$getValue(s), result_ref$getValue(ref_s)))
}

test_rel_entr <- function(cvx) {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- log(x)*ifelse(cvx, -1, 1)
  obj <- perspective(f, s)
  constraints <- list(1 <= s, s <= 2, 1 <= x, x <= 2)
  prob <- Problem(ifelse(cvx, Minimize(obj), Maximize(obj)))
  result <- solve(prob, solver = "ECOS")
  
  # Reference problem.
  ref_x <- Variable()
  ref_s <- Variable()
  obj <- rel_entr(ref_s, ref_x)*ifelse(cvx, 1, -1)
  
  ref_constraints <- list(1 <= ref_x, ref_x <= 2, 1 <= ref_s, ref_s <= 2)
  ref_prob <- Problem(ifelse(cvx, Minimize(obj), Maximize(obj)), ref_constraints)
  result_ref <- solve(ref_prob, solver = "ECOS")
  
  expect_true(np.isclose(result$value, result_ref$value))
  expect_true(is.allclose(result$getValue(x), result_ref$getValue(ref_x)))
  expect_true(np.isclose(result$getValue(s), result_ref$getValue(ref_s)))
}

lse_example <- function() {
  # Reference problem.
  ref_x <- Variable(3)
  ref_s <- Variable()
  ref_z <- Variable(3)
  ref_t <- Variable()
  
  ref_constraints <- list(ref_s >= sum(ref_z), 
                          c(1,2,3) <= ref_x, 1 <= ref_s, ref_s <= 2)
  ref_constraints <- c(ref_constraints, lapply(1:3, function(i) { ExpCone(ref_x[i] - ref_t, ref_s, ref_z[i]) }))
  ref_prob <- Problem(Minimize(ref_t), ref_constraints)
  result <- solve(ref_prob, solver = "ECOS")
  return(list(result$value, result$getValue(ref_x), result$getValue(ref_s)))
}

test_evaluate_persp <- function(x_val, s_val) {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f_exp <- square(x) + 3*x - 5
  obj <- perspective(f_exp, s)
  
  val_array <- matrix(c(s_val, x_val))
  value(x) <- x_val
  value(s) <- s_val   # Currently assumes variables have values before querying.
  val <- as.numeric(val_array)
  
  # True val.
  ref_val <- x_val^2/s_val + 3*x_val - 5*s_val
  expect_true(np.isclose(val, ref_val))
}

test_power <- function(n) {
  # Reference problem.
  ref_x <- Variable(pos = TRUE)
  ref_s <- Variable(pos = TRUE)
  
  # f(x) = x^n -> persp(f)(x,s) = x^n / s^(n-1)
  obj <- power(ref_x, n) / power(ref_s, n-1)
  constraints <- list(ref_x >= 1, ref_s <= 0.5)
  ref_prob <- Problem(Minimize(obj), constraints)
  result_ref <- solve(ref_prob, gp = TRUE)
  
  # Perspective problem.
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- power(x, n)
  obj <- perspective(f, s)
  constraints <- list(x >= 1, s <= 0.5)
  
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "ECOS")
  
  expect_true(np.isclose(result$value, result_ref$value))
  expect_true(np.isclose(result$getValue(x), result_ref$getValue(ref_x)))
  expect_true(np.isclose(result$getValue(s), result_ref$getValue(ref_s)))
}

test_psd_mf_persp <- function(n) {
  # Reference problem.
  ref_x <- Variable(n)
  ref_P <- Variable(n, n, PSD = TRUE)
  
  # matrix_frac is homogenous of degree 1 so its perspective is itself.
  obj <- matrix_frac(ref_x, ref_P)
  constraints <- list(ref_x == 5, ref_P == diag(n))
  ref_prob <- Problem(Minimize(obj), constraints)
  result_ref <- solve(ref_prob, solver = "SCS")
  
  # Perspective problem.
  x <- Variable(n)
  P <- Variable(n, n, PSD = TRUE)
  s <- Variable(nonneg = TRUE)
  
  f <- matrix_frac(x, P)
  obj <- perspective(f, s)
  constraints <- list(x == 5, P == diag(n), s == 1)
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "SCS")
  
  expect_equal(result$status, OPTIMAL)
  expect_true(np.isclose(result$value, result_ref$value, atol = 1e-2))
  expect_true(is.allclose(result$getValue(x), result_ref$getValue(ref_x), atol = 1e-2))
}

test_psd_tr_square <- function(n) {
  # Reference problem.
  ref_s <- Variable(nonneg = TRUE)
  ref_P <- Variable(n, n, PSD = TRUE)
  
  # Tr(X)^2 perspective is quad over line of Tr(X).
  obj <- quad_over_lin(matrix_trace(ref_P), ref_s)
  constraints <- list(ref_s <= 5, ref_P %>>% diag(n))
  ref_prob <- Problem(Minimize(obj), constraints)
  result_ref <- solve(ref_prob, solver = "SCS")
  
  # Perspective problem.
  P <- Variable(n, n, PSD = TRUE)
  s <- Variable(nonneg = TRUE)
  
  f <- perspective(square(matrix_trace(P)), s)
  obj <- perspective(f, s)
  constraints <- list(s <= 5, P %>>% diag(n))
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "SCS")
  
  expect_equal(result$status, OPTIMAL)
  expect_true(np.isclose(result$value, result_ref$value, atol = 1e-3))
  expect_true(is.allclose(result$getValue(P), result_ref$getValue(ref_P), atol = 1e-4))
}

test_that("test monotonicity", {
  x <- Variable(nonneg = TRUE)
  f <- exp(x)
  s <- Variable(nonneg = TRUE)
  p <- perspective(f, s)
  
  expect_true(is_nonneg(p))
  expect_false(is_nonpos(p))
  
  expect_false(is_incr(p, 1))
  expect_false(is_incr(p, 2))
  
  expect_false(is_decr(p, 1))
  expect_false(is_decr(p, 2))
})

test_that("test p-norms", {
  for(p in c(1,2))
    test_p_norms(p)
})

test_that("test rel_entr", {
  for(cvx in c(TRUE, FALSE))
    test_rel_entr(cvx)
})

test_that("test exp", {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- exp(x)
  obj <- perspective(f, s)
  constraints <- list(s >= 1, 1 <= x)
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "ECOS")
  
  # Reference problem.
  ref_x <- Variable()
  ref_s <- Variable()
  ref_z <- Variable()
  
  obj <- ref_z
  ref_constraints <- list(ExpCone(ref_x, ref_s, ref_z), ref_x >= 1, ref_s >= 1)
  ref_prob <- Problem(Minimize(obj), ref_constraints)
  result_ref <- solve(ref_prob, solver = "ECOS")
  
  expect_true(np.isclose(result$value, result_ref$value))
  expect_true(np.isclose(result$getValue(x), result_ref$getValue(ref_x)))
  expect_true(np.isclose(result$getValue(s), result_ref$getValue(ref_s)))
})

test_that("test lse", {
  x <- Variable(3)
  s <- Variable(nonneg = TRUE)
  f <- log_sum_exp(x)
  
  obj <- perspective(f, s)
  constraints <- list(1 <= s, s <= 2, c(1,2,3) <= x)
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "ECOS")
  
  ref <- lse_example()
  ref_prob <- ref[[1]]
  ref_x <- ref[[2]]
  ref_s <- ref[[3]]
  
  expect_true(np.isclose(result$value, ref_prob))
  expect_true(is.allclose(result$getValue(x), ref_x))
  expect_true(np.isclose(result$getValue(s), ref_s))
})

test_that("test evaluate persp", {
  xs_vals <- list(c(1,2), c(5,0.25, c(0.5,7)))
  for(xs_val in xs_vals) {
    x_val <- xs_val[1]
    s_val <- xs_val[2]
    test_evaluate_persp(x_val, s_val)
  }
})

test_that("test quad atom", {
  for(r_val in c(2, 3, 4, -2, 0)) {
    ref <- quad_example(r_val)
    ref_val <- ref[[1]]
    ref_s <- ref[[2]]
    ref_x <- ref[[3]]
    r <- ref[[4]]
    
    # Form objective, introduce original variable.
    x <- Variable()
    s <- Variable(nonneg = TRUE)
    f_exp <- square(x) + r*x - 4
    obj <- perspective(f_exp, s)
    
    constraints <- list(s <= 0.5, x >= 2)
    prob <- Problem(Minimize(obj), constraints)
    result <- solve(prob, verbose = TRUE)
    
    expect_true(np.isclose(result$value, ref_val))
    expect_true(np.isclose(result$getValue(x), ref_x))   # Assuming the solutions are unique...
    expect_true(np.isclose(result$getValue(s), ref_s))
  }
})

test_that("test quad persp persp", {
  for(r_val in c(2, 3, 4, -2, 0)) {
    ref <- quad_example(r_val)
    ref_val <- ref[[1]]
    ref_s <- ref[[2]]
    ref_x <- ref[[3]]
    r <- ref[[4]]
    
    # Form objective, introduce original variable.
    x <- Variable()
    s <- Variable(nonneg = TRUE)
    t <- Variable(nonneg = TRUE)
    
    # f(x) -> sf(x/s) -> t(s/t)*f(xt/ts) - > sf(x/s)
    f_exp <- square(x) + r*x - 4
    obj_inner <- perspective(f_exp, s)
    obj <- perspective(obj_inner, t)
    
    constraints <- list(0.1 <= s, s <= 0.5, x >= 2, 0.1 <= t, t <= 0.5)
    prob <- Problem(Minimize(obj), constraints)
    result <- solve(prob, verbose = TRUE)
    
    expect_true(np.isclose(result$value, ref_val))
    expect_true(np.isclose(result$getValue(x), ref_x))
    expect_true(np.isclose(result$getValue(s), ref_s))
  }
})

test_that("test quad quad", {
  # Reference problem.
  ref_x <- Variable()
  ref_s <- Variable(nonneg =  TRUE)
  
  # f(x) = x^4 -> persp(f)(x,s) = x^4 / s^3 = (x^2/s)^2 / s.
  obj <- quad_over_lin(quad_over_lin(ref_x, ref_s), ref_s)
  constraints <- list(ref_x >= 5, ref_s <= 3)
  ref_prob <- Problem(Minimize(obj), constraints)
  result_ref <- solve(ref_prob, solver = "ECOS")
  
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "ECOS")
  
  expect_true(np.isclose(result$value, result_ref$value))
  expect_true(np.isclose(result$getValue(x), result_ref$getValue(ref_x)))
  expect_true(np.isclose(result$getValue(s), result_ref$getValue(ref_s)))
})

test_that("test power", {
  for(n in c(4, 5, 7, 11))
    test_power(n)
})

test_that("test psd tr persp", {
  # Reference problem.
  ref_P <- Variable(2, 2, PSD = TRUE)
  
  obj <- matrix_trace(ref_P)
  constraints <- list(ref_P == diag(2))
  ref_prob <- Problem(Minimize(obj), constraints)
  result_ref <- solve(ref_prob, solver = "SCS")
  
  # Perspective problem.
  P <- Variable(2, 2, PSD = TRUE)
  s <- Variable(nonneg = TRUE)
  
  f <- matrix_trace(P)
  
  obj <- perspective(f, s)
  constraints <- list(P == diag(2), s == 1)
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "SCS")
  
  expect_equal(result$status, OPTIMAL)
  expect_true(np.isclose(result$value, result_ref$value))
})

test_that("test psd mf persp", {
  for(n in c(2, 3, 11))
    test_psd_mf_persp(n)
})

test_that("test psd tr square", {
  for(n in c(2, 3, 11))
    test_psd_tr_square(n)
})

test_that("test diag", {
  X_ref <- Variable(2, 2, diag = TRUE)
  obj <- matrix_trace(X_ref)
  constraints <- list(diag(X_ref) >= c(1,2))
  ref_prob <- Problem(Minimize(obj), constraints)
  result_ref <- solve(ref_prob)
  
  X <- Variable(2, 2, diag = TRUE)
  f <- matrix_trace(X)
  s <- Variable(nonneg = TRUE)
  obj <- perspective(f, s)
  constraints <- list(diag(X) >= c(1,2), s == 1)
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob)
  
  expect_equal(result$status, OPTIMAL)
  expect_true(np.isclose(result$value, result_ref$value, atol = 1e-3))
  expect_true(is.allclose(result$getValue(X), result_ref$getValue(X_ref), atol = 1e-4))
})

test_that("test scalar x", {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  obj <- perspective(x-1, s)
  
  prob <- Problem(Minimize(obj), list(x >= 3.14, s <= 1))
  result <- solve(prob)
  expect_true(np.isclose(result$value, 3.14 - 1))
})

test_that("test assert s nonzero", {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  obj <- perspective(x+1, s)
  
  prob <- Problem(Minimize(obj), list(x >= 3.14))
  expect_error(result <- solve(prob), "There are valid cases")
})

test_that("test parameter", {
  p <- Parameter(nonneg = TRUE)
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- p*square(x)
  
  obj <- perspective(f, s)
  prob <- Problem(Minimize(obj), list(s <= 1, x >= 2))
  value(p) <- 99
  
  result <- solve(prob)
  
  expect_true(np.isclose(result$value, 4*value(p)))
})

test_that("test affine s", {
  # Test requiring affine s to be nonneg.
  x <- Variable()
  s <- Variable(2)
  expect_error(perspective(square(x), sum(s)), "s must be a variable")
})

test_that("test dpp", {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  a <- Parameter()
  
  obj <- perspective(square(a + x), s)
  expect_false(is_dpp(obj))
  
  obj <- perspective(log(a + x), s)
  expect_false(is_dpp(obj))
})
