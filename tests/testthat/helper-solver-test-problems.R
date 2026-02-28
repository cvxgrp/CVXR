## Standard solver test problems ported from CVXPY's solver_test_helpers.py
## Each function returns a list with:
##   prob       — a CVXR Problem
##   expect_obj — expected optimal value
##   expect_x   — expected primal (NULL if not checked)
##   con_duals  — list of expected dual values (NULL entries = not checked)
##   cone_type  — string label ("LP", "QP", "SOCP", "SDP", "ExpCone", "PowCone3D", etc.)
##
## Expected values were verified via `uv run python` against CVXPY + Clarabel.

# ── LP helpers ────────────────────────────────────────────────────

sth_lp_0 <- function() {
  x <- Variable(2)
  objective <- Minimize(norm1(x) + 1.0)
  constraints <- list(x == 0)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = 1,
    expect_x = c(0, 0),
    con_duals = list(NULL),
    cone_type = "LP"
  )
}

sth_lp_1 <- function() {
  x <- Variable(2)
  objective <- Minimize(-4 * x[1] + -5 * x[2])
  constraints <- list(
    2 * x[1] + x[2] <= 3,
    x[1] + 2 * x[2] <= 3,
    x[1] >= 0,
    x[2] >= 0
  )
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -9,
    expect_x = c(1, 1),
    con_duals = list(1, 2, 0, 0),
    cone_type = "LP"
  )
}

sth_lp_2 <- function() {
  x <- Variable(2)
  objective <- Minimize(x[1] + 0.5 * x[2])
  constraints <- list(
    x[1] >= -100,
    x[1] <= -10,
    x[2] == 1
  )
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -99.5,
    expect_x = c(-100, 1),
    con_duals = list(1, 0, -0.5),
    cone_type = "LP"
  )
}

sth_lp_3 <- function() {
  x <- Variable(5)
  objective <- Minimize(sum_entries(x))
  constraints <- list(x <= 1)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -Inf,
    expect_x = NULL,
    con_duals = list(NULL),
    cone_type = "LP"
  )
}

sth_lp_4 <- function() {
  x <- Variable(5)
  objective <- Minimize(sum_entries(x))
  constraints <- list(x <= 0, x >= 1)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = Inf,
    expect_x = NULL,
    con_duals = list(NULL, NULL),
    cone_type = "LP"
  )
}

sth_lp_5 <- function() {
  ## Deterministic A, b, c from numpy seed=0 (cannot replicate in R)
  A <- matrix(c(
     1.76405235,  0.40015721,  0.97873798,  2.2408932 ,  1.86755799,
    -0.97727788,  0.95008842, -0.15135721, -0.10321885,  0.4105985 ,
     0.14404357,  1.45427351,  0.76103773,  0.12167502,  0.44386323,
     0.33367433,  1.49407907, -0.20515826,  0.3130677 , -0.85409574,
    -2.55298982,  0.6536186 ,  0.8644362 , -0.74216502,  2.26975462,
    -1.45436567,  0.04575852, -0.18718385,  1.53277921,  1.46935877,
     0.15494743,  0.37816252, -0.88778575, -1.98079647, -0.34791215,
     0.15634897,  1.23029068,  1.20237985, -0.38732682, -0.30230275,
    -1.85237813,  1.68045187,  1.41161128, -0.37105888,  3.01286545,
    -1.52280044,  1.45380953, -0.09858135,  1.58000735,  1.03164353,
     0.4692577 ,  0.91108765,  0.63477278,  0.47314843,  1.04942518,
    -0.33228632,  1.20398975,  0.0298212 ,  0.19370282, -0.12409246
  ), nrow = 6, ncol = 10, byrow = TRUE)
  b <- c(2.47092245, -2.09373443, 11.89809328, -1.38621545,
         11.05577432, 0.97091495)
  cvec <- c(-4.92644682,  4.67328509,  1.69191946, -5.67035615,
             5.63813844, -2.43137261,  5.20204355,  1.61746403,
             3.34401775,  1.35955943)

  x <- Variable(10)
  objective <- Minimize(t(cvec) %*% x)
  constraints <- list(x >= 0, A %*% x == b)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = 21.275425028808698,
    expect_x = NULL,  # only feasibility checked, not exact primal
    con_duals = list(NULL, NULL),
    cone_type = "LP"
  )
}

# ── QP helpers ────────────────────────────────────────────────────

sth_qp_0 <- function() {
  x <- Variable(1)
  objective <- Minimize(power(x, 2))
  constraints <- list(x[1] >= 1)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = 1,
    expect_x = 1,
    con_duals = list(2),
    cone_type = "QP"
  )
}

# ── SOCP helpers ──────────────────────────────────────────────────

sth_socp_0 <- function() {
  x <- Variable(2)
  objective <- Minimize(p_norm(x, 2) + 1)
  constraints <- list(x == 0)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = 1,
    expect_x = c(0, 0),
    con_duals = list(NULL),
    cone_type = "SOCP"
  )
}

sth_socp_1 <- function() {
  x <- Variable(3)
  y <- Variable(1)
  soc <- SOC(y, x)
  constraints <- list(
    soc,
    x[1] + x[2] + 3 * x[3] >= 1.0,
    y <= 5
  )
  objective <- Minimize(3 * x[1] + 2 * x[2] + x[3])
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -13.548638904065102,
    expect_x = c(-3.87462, -2.12979, 2.33480),
    expect_y = 5,
    con_duals = list(
      list(c(2.86560), c(2.22063, 1.22063, -1.33812)),  # SOC: [t_dual, x_dual]
      0.779397,
      2.865603
    ),
    cone_type = "SOCP"
  )
}

sth_socp_2 <- function() {
  x <- Variable(2)
  objective <- Minimize(-4 * x[1] + -5 * x[2])
  expr <- reshape_expr(x[1] + 2 * x[2], c(1L, 1L))
  constraints <- list(
    2 * x[1] + x[2] <= 3,
    SOC(Constant(c(3)), expr),
    x[1] >= 0,
    x[2] >= 0
  )
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -9,
    expect_x = c(1, 1),
    con_duals = list(1, list(c(2.0), matrix(-2.0, 1, 1)), 0, 0),
    cone_type = "SOCP"
  )
}

sth_socp_3_ax0 <- function() {
  x <- Variable(2)
  cvec <- c(-1, 2)
  root2 <- sqrt(2)
  u <- matrix(c(1/root2, 1/root2, -1/root2, 1/root2), 2, 2)
  mat1 <- diag(c(root2, 1/root2)) %*% t(u)
  mat2 <- diag(c(1, 1))
  mat3 <- diag(c(0.2, 1.8))
  ## In CVXPY, Variable(2) is 1D so mat@x is 1D, vstack gives (3,2).
  ## In R, Variable(2) is (2,1) so mat%*%x is (2,1), need t() to get (1,2) rows.
  X <- vstack(t(mat1 %*% x), t(mat2 %*% x), t(mat3 %*% x))
  tvec <- Constant(rep(1, 3))
  objective <- Minimize(t(cvec) %*% x)
  con <- SOC(tvec, t(X), axis = 2L)
  prob <- Problem(objective, list(con))
  list(
    prob = prob,
    expect_obj = -1.932105,
    expect_x = c(0.83666003, -0.54772256),
    con_duals = list(
      list(
        c(0, 1.16452, 0.76756),
        matrix(c(0, -0.97431, -0.12845,
                 0,  0.63783,  0.75674), nrow = 2, byrow = TRUE)
      )
    ),
    cone_type = "SOCP"
  )
}

sth_socp_3_ax1 <- function() {
  x <- Variable(2)
  cvec <- c(-1, 2)
  root2 <- sqrt(2)
  u <- matrix(c(1/root2, 1/root2, -1/root2, 1/root2), 2, 2)
  mat1 <- diag(c(root2, 1/root2)) %*% t(u)
  mat2 <- diag(c(1, 1))
  mat3 <- diag(c(0.2, 1.8))
  ## Same transpose trick as socp_3_ax0: each row is a (1,2) vector
  X <- vstack(t(mat1 %*% x), t(mat2 %*% x), t(mat3 %*% x))
  tvec <- Constant(rep(1, 3))
  objective <- Minimize(t(cvec) %*% x)
  con <- SOC(tvec, X, axis = 1L)
  prob <- Problem(objective, list(con))
  list(
    prob = prob,
    expect_obj = -1.932105,
    expect_x = c(0.83666003, -0.54772256),
    con_duals = list(
      list(
        c(0, 1.16452, 0.76756),
        matrix(c(0, 0, -0.97431, 0.63783, -0.12845, 0.75674),
               nrow = 3, byrow = TRUE)
      )
    ),
    cone_type = "SOCP"
  )
}

# ── SDP helpers ───────────────────────────────────────────────────

sth_sdp_1_min <- function() {
  rho <- Variable(c(4, 4), symmetric = TRUE)
  constraints <- list(
    0.6 <= rho[1, 2], rho[1, 2] <= 0.9,
    0.8 <= rho[1, 3], rho[1, 3] <= 0.9,
    0.5 <= rho[2, 4], rho[2, 4] <= 0.7,
    -0.8 <= rho[3, 4], rho[3, 4] <= -0.4,
    rho[1, 1] == 1, rho[2, 2] == 1, rho[3, 3] == 1, rho[4, 4] == 1,
    PSD(rho)
  )
  objective <- Minimize(rho[1, 4])
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -0.39,
    expect_x = NULL,
    con_duals = replicate(length(constraints), NULL, simplify = FALSE),
    cone_type = "SDP"
  )
}

sth_sdp_1_max <- function() {
  rho <- Variable(c(4, 4), symmetric = TRUE)
  constraints <- list(
    0.6 <= rho[1, 2], rho[1, 2] <= 0.9,
    0.8 <= rho[1, 3], rho[1, 3] <= 0.9,
    0.5 <= rho[2, 4], rho[2, 4] <= 0.7,
    -0.8 <= rho[3, 4], rho[3, 4] <= -0.4,
    rho[1, 1] == 1, rho[2, 2] == 1, rho[3, 3] == 1, rho[4, 4] == 1,
    PSD(rho)
  )
  objective <- Maximize(rho[1, 4])
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = 0.23,
    expect_x = NULL,
    con_duals = replicate(length(constraints), NULL, simplify = FALSE),
    cone_type = "SDP"
  )
}

sth_sdp_2 <- function() {
  X1 <- Variable(c(2, 2), symmetric = TRUE)
  X2 <- Variable(c(4, 4), symmetric = TRUE)
  C1 <- matrix(c(1, 0, 0, 6), 2, 2)
  A1 <- matrix(c(1, 1, 1, 2), 2, 2)
  C2 <- matrix(c(1, -3, 0, 0, -3, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 4, 4)
  A2 <- matrix(c(0, 1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3), 4, 4)

  constraints <- list(
    matrix_trace(A1 %*% X1) + matrix_trace(A2 %*% X2) == 23,
    X2[1, 2] <= -3,
    PSD(X1),
    PSD(X2)
  )
  objective <- Minimize(matrix_trace(C1 %*% X1) + matrix_trace(C2 %*% X2))
  prob <- Problem(objective, constraints)

  expect_X1 <- matrix(c(21.04761, 4.07692, 4.07692, 0.78970), 2, 2)
  expect_X2 <- matrix(c(5.05374, -3, 0, 0,
                         -3, 1.78086, 0, 0,
                         0, 0, 0, 0,
                         0, 0, 0, 0), 4, 4, byrow = TRUE)

  list(
    prob = prob,
    expect_obj = 52.4013,
    expect_X1 = expect_X1,
    expect_X2 = expect_X2,
    con_duals = list(
      -0.83772,
      11.04455,
      matrix(c(0.16228, -0.83772, -0.83772, 4.32456), 2, 2),
      matrix(c(1, 1.68455, 0, 0,
               1.68455, 2.83772, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 2.51317), 4, 4, byrow = TRUE)
    ),
    cone_type = "SDP"
  )
}

# ── ExpCone helpers ───────────────────────────────────────────────

sth_expcone_1 <- function() {
  x <- Variable(c(3, 1))
  ## CVXPY: ExpCone(x[2], x[1], x[0]) → z >= y*exp(x/y)
  ## R 1-indexed: ExpCone(x[3], x[2], x[1])
  cone_con <- ExpCone(x[3], x[2], x[1])
  constraints <- list(
    sum_entries(x) <= 1.0,
    sum_entries(x) >= 0.1,
    x >= 0,
    cone_con
  )
  objective <- Minimize(3 * x[1] + 2 * x[2] + x[3])
  prob <- Problem(objective, constraints)

  expect_x <- matrix(c(0.05463, 0.02609, 0.01928), 3, 1)

  list(
    prob = prob,
    expect_obj = 0.23535,
    expect_x = expect_x,
    con_duals = list(0, 2.35348, matrix(0, 3, 1),
      list(matrix(-1.35348, 1, 1), matrix(-0.35348, 1, 1), matrix(0.64652, 1, 1))
    ),
    cone_type = "ExpCone"
  )
}

sth_expcone_socp_1 <- function() {
  sigma <- matrix(c(1.83, 1.79, 3.22,
                     1.79, 2.18, 3.18,
                     3.22, 3.18, 8.69), 3, 3)
  L <- t(chol(sigma))  # R's chol gives upper-tri; transpose for lower
  cc <- 0.75
  tt <- Variable(1)
  x <- Variable(3)
  s <- Variable(3)
  e <- Constant(rep(1, 3))
  objective <- Minimize(tt - cc * t(e) %*% s)
  con1 <- p_norm(t(L) %*% x, 2) <= tt
  con2 <- ExpCone(s, e, x)
  prob <- Problem(objective, list(con1, con2))
  list(
    prob = prob,
    expect_obj = 4.07512,
    expect_x = c(0.57608, 0.54315, 0.28037),
    expect_s = c(-0.55150, -0.61036, -1.27161),
    con_duals = list(
      1.0,
      list(c(-0.75, -0.75, -0.75),
           c(-1.16363, -1.20777, -1.70371),
           c(1.30190, 1.38082, 2.67496))
    ),
    cone_type = "Mixed_ExpCone_SOC"
  )
}

# ── PowCone helpers ──────────────────────────────────────────────

sth_pcp_1 <- function() {
  x <- Variable(3)
  y_square <- Variable(1)
  epis <- Variable(3)
  constraints <- list(
    PowCone3D(Constant(rep(1, 3)), epis, x, Constant(c(0.5, 0.5, 0.5))),
    sum_entries(epis) <= y_square,
    x[1] + x[2] + 3 * x[3] >= 1.0,
    y_square <= 25
  )
  objective <- Minimize(3 * x[1] + 2 * x[2] + x[3])
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -13.5486,
    expect_x = c(-3.87462, -2.12979, 2.33480),
    con_duals = list(
      list(c(4.30209, 1.29985, 1.56212),
           c(0.28656, 0.28656, 0.28656),
           c(2.22063, 1.22063, -1.33811)),
      0.28656,
      0.77940,
      0.28656
    ),
    cone_type = "PowCone3D"
  )
}

sth_pcp_2 <- function() {
  x <- Variable(3)
  ## CVXPY Variable(shape=(2,)) is 1D; hstack gives (2,).
  ## In R, hstack(x[1], x[3]) gives (1,2), so hypos must be (1,2) too.
  hypos <- Variable(c(1, 2))
  objective <- Minimize(-sum_entries(hypos) + x[1])
  arg1 <- hstack(x[1], x[3])
  arg2 <- hstack(x[2], Constant(1.0))
  pc_con <- PowCone3D(arg1, arg2, hypos, matrix(c(0.2, 0.4), 1, 2))
  constraints <- list(
    x[1] + x[2] + 0.5 * x[3] == 2,
    pc_con
  )
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -1.80734,
    expect_x = c(0.06395, 0.78321, 2.30571),
    con_duals = list(
      0.48466,
      list(c(1.48466, 0.24233), c(0.48466, 0.83801), c(-1, -1))
    ),
    cone_type = "PowCone3D"
  )
}

sth_pcp_3 <- function() {
  D <- matrix(c(
    -1.0856306,   0.99734545,
     0.2829785,  -1.50629471,
    -0.57860025,  1.65143654,
    -2.42667924, -0.42891263,
     1.26593626, -0.8667404,
    -0.67888615, -0.09470897,
     1.49138963, -0.638902
  ), nrow = 7, ncol = 2, byrow = TRUE)

  p <- 1 / 0.4  # 2.5
  TT <- nrow(D)
  w <- Variable(c(2, 1))
  tt <- Variable(1)
  d <- Variable(c(TT, 1))
  ones_T <- matrix(1, TT, 1)

  powcone <- PowCone3D(d, tt * Constant(ones_T), D %*% w, 1/p)
  constraints <- list(
    sum_entries(w) == 1,
    w >= 0,
    powcone,
    sum_entries(d) == tt
  )
  objective <- Minimize(tt)
  prob <- Problem(objective, constraints)

  w_opt <- matrix(c(0.39336, 0.60664), 2, 1)

  list(
    prob = prob,
    expect_obj = 1.51430,
    expect_w = w_opt,
    con_duals = list(-1.51430, matrix(0, 2, 1), NULL, 0.40001),
    cone_type = "PowCone3D"
  )
}

# ── MIP LP helpers ──────────────────────────────────────────────

sth_mi_lp_0 <- function() {
  x <- Variable(2)
  bool_var <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(norm1(x) + 1.0),
                  list(x == bool_var, bool_var == 0))
  list(
    prob = prob,
    expect_obj = 1,
    expect_x = c(0, 0),
    cone_type = "MIP_LP"
  )
}

sth_mi_lp_1 <- function() {
  x <- Variable(2)
  boolvar <- Variable(boolean = TRUE)
  intvar <- Variable(integer = TRUE)
  prob <- Problem(
    Minimize(-4 * x[1] - 5 * x[2]),
    list(
      2 * x[1] + x[2] <= intvar,
      x[1] + 2 * x[2] <= 3 * boolvar,
      x >= 0,
      intvar == 3 * boolvar,
      intvar == 3
    )
  )
  list(
    prob = prob,
    expect_obj = -9,
    expect_x = c(1, 1),
    expect_boolvar = 1,
    expect_intvar = 3,
    cone_type = "MIP_LP"
  )
}

sth_mi_lp_2 <- function() {
  ## Knapsack instance "knapPI_1_50_1000_1"
  ## from http://www.diku.dk/~pisinger/genhard.c
  n <- 50L
  c_cap <- 995
  z_opt <- 8373
  profits <- c(94, 506, 416, 992, 649, 237, 457, 815, 446, 422,
               791, 359, 667, 598, 7, 544, 334, 766, 994, 893,
               633, 131, 428, 700, 617, 874, 720, 419, 794, 196,
               997, 116, 908, 539, 707, 569, 537, 931, 726, 487,
               772, 513, 81, 943, 58, 303, 764, 536, 724, 789)
  weights <- c(485, 326, 248, 421, 322, 795, 43, 845, 955, 252,
               9, 901, 122, 94, 738, 574, 715, 882, 367, 984,
               299, 433, 682, 72, 874, 138, 856, 145, 995, 529,
               199, 277, 97, 719, 242, 107, 122, 70, 98, 600,
               645, 267, 972, 895, 213, 748, 487, 923, 29, 674)
  X <- Variable(n, boolean = TRUE)
  prob <- Problem(
    Maximize(sum_entries(multiply(profits, X))),
    list(sum_entries(multiply(weights, X)) <= c_cap)
  )
  list(
    prob = prob,
    expect_obj = z_opt,
    expect_x = NULL,
    cone_type = "MIP_LP"
  )
}

# ── MIP SOCP helpers ────────────────────────────────────────────

sth_mi_socp_1 <- function() {
  x <- Variable(3)
  y <- Variable(2, integer = TRUE)
  prob <- Problem(
    Minimize(3 * x[1] + 2 * x[2] + x[3] + y[1] + 2 * y[2]),
    list(
      cvxr_norm(x, 2) <= y[1],
      cvxr_norm(x, 2) <= y[2],
      x[1] + x[2] + 3 * x[3] >= 0.1,
      y <= 5
    )
  )
  list(
    prob = prob,
    expect_obj = 0.21364,
    expect_x = c(-0.78510, -0.43565, 0.44025),
    expect_y = c(1, 1),
    cone_type = "MIP_SOCP"
  )
}

sth_mi_socp_2 <- function() {
  x <- Variable(2)
  bool_var <- Variable(boolean = TRUE)
  int_var <- Variable(integer = TRUE)
  prob <- Problem(
    Minimize(-4 * x[1] - 5 * x[2]),
    list(
      2 * x[1] + x[2] <= int_var,
      square(x[1] + 2 * x[2]) <= 9 * bool_var,
      x >= 0,
      int_var == 3 * bool_var,
      int_var == 3
    )
  )
  list(
    prob = prob,
    expect_obj = -9,
    expect_x = c(1, 1),
    expect_boolvar = 1,
    expect_intvar = 3,
    cone_type = "MIP_SOCP"
  )
}

# ── MIP LP helpers (continued) ───────────────────────────────────

sth_mi_lp_3 <- function() {
  ## Infeasible boolean MIP: 4 boolean vars, sum==2, pairwise constraints
  ## CVXPY SOURCE: solver_test_helpers.py mi_lp_3()
  x <- Variable(4, boolean = TRUE)
  prob <- Problem(
    Maximize(Constant(1)),
    list(
      x[1] + x[2] + x[3] + x[4] <= 2,
      x[1] + x[2] + x[3] + x[4] >= 2,
      x[1] + x[2] <= 1,
      x[1] + x[3] <= 1,
      x[1] + x[4] <= 1,
      x[3] + x[4] <= 1,
      x[2] + x[4] <= 1,
      x[2] + x[3] <= 1
    )
  )
  list(
    prob = prob,
    expect_obj = Inf,  # infeasible
    expect_x = NULL,
    cone_type = "MIP_LP"
  )
}

sth_mi_lp_4 <- function() {
  ## Boolean abs value MIP
  ## CVXPY SOURCE: solver_test_helpers.py mi_lp_4()
  x <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(abs(x)), list())
  list(
    prob = prob,
    expect_obj = 0,
    expect_x = 0,
    cone_type = "MIP_LP"
  )
}

sth_mi_lp_5 <- function() {
  ## Multiple integer vars with simple constraints
  ## CVXPY SOURCE: solver_test_helpers.py mi_lp_5()
  x <- Variable(2, integer = TRUE)
  prob <- Problem(
    Minimize(sum_entries(x)),
    list(x[1] + x[2] >= 3, x >= 0)
  )
  list(
    prob = prob,
    expect_obj = 3,
    expect_x = NULL,  # multiple optima (e.g., (0,3), (1,2), (2,1), (3,0))
    cone_type = "MIP_LP"
  )
}

sth_lp_6 <- function() {
  ## LP with dual variables checked
  ## CVXPY SOURCE: solver_test_helpers.py lp_6()
  x <- Variable(2)
  prob <- Problem(
    Minimize(x[1] + 2 * x[2]),
    list(x[1] + x[2] >= 1, x >= 0)
  )
  list(
    prob = prob,
    expect_obj = 1,
    expect_x = c(1, 0),
    con_duals = list(1, c(0, 2)),
    cone_type = "LP"
  )
}

# ── PowCone helpers (continued) ─────────────────────────────────

sth_mi_pcp_0 <- function() {
  ## Mixed-integer power cone problem
  ## CVXPY SOURCE: solver_test_helpers.py mi_pcp_0()
  x <- Variable(3)
  y <- Variable(integer = TRUE)
  constraints <- list(
    PowCone3D(x[1], x[2], x[3], 0.25),
    x[1] >= 2,
    x[2] >= 1,
    y >= x[3],
    y <= -1
  )
  objective <- Minimize(y)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = -1,
    expect_x = NULL,
    expect_y = -1,
    cone_type = "MIP_PCP"
  )
}

# ── SDP+PCP mixed helper ────────────────────────────────────────

sth_sdp_pcp_1 <- function() {
  ## Mixed SDP + power cone: minimize PSD matrix trace with pow cone constraint
  ## CVXPY SOURCE: solver_test_helpers.py sdp_pcp_1()
  X <- Variable(c(2, 2), symmetric = TRUE)
  t_var <- Variable(1)
  constraints <- list(
    PSD(X),
    X[1, 1] >= 1,
    X[2, 2] >= 1,
    PowCone3D(X[1, 1], X[2, 2], t_var, 0.5)
  )
  objective <- Minimize(matrix_trace(X) + t_var)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = 1,  # trace(I) + (-1) = 2 + (-1) = 1
    expect_x = NULL,
    cone_type = "Mixed_SDP_PCP"
  )
}

# ── QP with linear objective ────────────────────────────────────

sth_qp_0_linear_obj <- function() {
  ## QP problem with purely linear objective (quadratic only in constraints)
  ## CVXPY SOURCE: test_conic_solvers.py clarabel_qp_0_linear_obj
  x <- Variable(1)
  objective <- Minimize(x[1])
  constraints <- list(x[1] >= 1)
  prob <- Problem(objective, constraints)
  list(
    prob = prob,
    expect_obj = 1,
    expect_x = 1,
    cone_type = "LP"
  )
}
