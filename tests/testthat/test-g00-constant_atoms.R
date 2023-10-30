context("test-constant_atoms")

ROBUST_CVXOPT <- "robust_cvxopt"
SOLVER_TO_TOL <- list(SCS = 1e-2, ECOS = 1e-7, OSQP = 1e-1)
# SOLVER_TO_TOL <- list(SCS = 1e-2, ECOS = 5e-7, OSQP = 1e-1)   # ECOS = 1e-7 fails a few tests.
SOLVERS_TO_TRY <- c("ECOS", "SCS", "OSQP")

# Test CVXOPT if installed.
if("CVXOPT" %in% installed_solvers()) {
    SOLVERS_TO_TRYS <- c(SOLVERS_TO_TRY, "CVXOPT", ROBUST_CVXOPT)
    SOLVERS_TO_TOL$CVXOPT <- 1e-7
    SOLVERS_TO_TOL[[ROBUST_CVXOPT]] <- 1e-7
}

# Test MOSEK if installed.
if("MOSEK" %in% installed_solvers()) {
    SOLVERS_TO_TRY <- c(SOLVERS_TO_TRY, "MOSEK")
    SOLVER_TO_TOL$MOSEK <- 1e-6
}

v_np <- matrix(c(-1, 2, -2), nrow = 1, ncol = 3)
log_sum_exp_axis_1 <- function(x) { log_sum_exp(x, axis = 1) }
log_sum_exp_axis_2 <- function(x) { log_sum_exp(x, axis = 2, keepdims = TRUE) }

# Map from solver name to a list of strings for atoms that fail.
KNOWN_SOLVER_ERRORS <- list()
KNOWN_SOLVER_ERRORS$MOSEK <- c("xexp")

atoms_minimize <- list(
            list(abs, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(5,2), c(3,1)))),
            list(function(x) { cumsum_axis(x, axis=1) }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(-5,2), c(-8,3)))),
            list(function(x) { cumsum_axis(x, axis=2) }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(-5,-3), c(-3,-2)))),
            list(function(x) { cummax_axis(x, axis=1) }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(-5,2), c(-3,2)))),
            list(function(x) { cummax_axis(x, axis=2) }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(-5,2), c(-3,1)))),
            list(diag, c(2,1), list(cbind(c(-5,2), c(-3,1))), Constant(c(-5,1))),
            list(diag, c(2,2), list(matrix(c(-5,1))), Constant(cbind(c(-5,0), c(0,1)))),
            list(exp, c(2,2), list(cbind(c(1,0), c(2,-1))), Constant(cbind(c(exp(1),1), c(exp(2), exp(-1))))),
            list(huber, c(2,2), list(cbind(c(0.5,-1.5), c(4,0))), Constant(cbind(c(0.25,2), c(7,0)))),
            list(function(x) { huber(x,2.5) }, c(2,2), list(cbind(c(0.5,-1.5), c(4,0))), Constant(cbind(c(0.25,2.25), c(13.75,0)))),
            list(inv_pos, c(2,2), list(cbind(c(1,2), c(3,4))), Constant(cbind(c(1,1.0/2), c(1.0/3,1.0/4)))),
            list(function(x) { (x + Constant(0))^-1 }, c(2,2), list(cbind(c(1,2), c(3,4))), Constant(cbind(c(1,1.0/2), c(1.0/3,1.0/4)))),
            list(kl_div, c(1,1), list(exp(1), 1), Constant(1)),
            list(kl_div, c(1,1), list(exp(1), exp(1)), Constant(0)),
            list(kl_div, c(2,1), list(matrix(c(exp(1), 1)), 1), Constant(c(1,0))),
            list(rel_entr, c(1,1), list(exp(1), 1), Constant(exp(1))),
            list(rel_entr, c(1,1), list(exp(1), exp(1)), Constant(0)),
            list(rel_entr, c(2,1), list(matrix(c(exp(1),1)), 1), Constant(c(exp(1),0))),
            
            # kronecker with variable in the right operand.
            list(function(x) { kronecker(cbind(c(1,2), c(3,4)), x) }, c(4,4), list(cbind(c(5,6), c(7,8))),
                 Constant(kronecker(cbind(c(1,2), c(3,4)), cbind(c(5,6), c(7,8))))),
            list(function(x) { kronecker(cbind(c(1,2), c(3,4), c(5,6)), x) }, c(6,4), list(cbind(c(5,6), c(7,8))),
                 Constant(kronecker(cbind(c(1,2), c(3,4), c(5,6)), cbind(c(5,6), c(7,8))))),
            list(function(x) { kronecker(cbind(c(1,2), c(3,4)), x) }, c(6,4), list(cbind(c(5,6), c(7,8), c(9,10))),
                 Constant(kronecker(cbind(c(1,2), c(3,4)), cbind(c(5,6), c(7,8), c(9,10))))),
            
            # kronecker with variable in the left operand.
            list(function(x) { kronecker(x, cbind(c(1,2), c(3,4))) }, c(4,4), list(cbind(c(5,6), c(7,8))),
                 Constant(kronecker(cbind(c(5,6), c(7,8)), cbind(c(1,2), c(3,4))))),
            list(function(x) { kronecker(x, cbind(c(1,2), c(3,4), c(5,6))) }, c(6,4), list(cbind(c(5,6), c(7,8))),
                 Constant(kronecker(cbind(c(5,6), c(7,8)), cbind(c(1,2), c(3,4), c(5,6))))),
            list(function(x) { kronecker(x, cbind(c(1,2), c(3,4))) }, c(6,4), list(cbind(c(5,6), c(7,8), c(9,10))),
                 Constant(kronecker(cbind(c(5,6), c(7,8), c(9,10)), cbind(c(1,2), c(3,4))))),
            
            list(lambda_max, c(1,1), list(cbind(c(2,0), c(0,1))), Constant(2)),
            list(lambda_max, c(1,1), list(cbind(c(2,0,0), c(0,3,0), c(0,0,1))), Constant(3)),
            list(lambda_max, c(1,1), list(cbind(c(5,7), c(7,-3))), Constant(9.06225775)),
            list(function(x) { lambda_sum_largest(x,2) }, c(1,1), list(cbind(c(1,2,3), c(2,4,5), c(3,5,6))), Constant(11.51572947)),
            list(log_sum_exp, c(1,1), list(cbind(c(5,7), c(0,-3))), Constant(7.1277708268)),
            list(log_sum_exp_axis_1, c(3,1), list(cbind(c(5,7,1), c(0,-3,6))), Constant(c(5.00671535, 7.0000454, 6.0067153))),
            list(log_sum_exp_axis_2, c(1,2), list(cbind(c(5,7,1), c(0,-3,6))), t(Constant(c(7.12910890, 6.00259878)))),
            list(logistic, c(2,2), list(cbind(c(log(5), log(7)), c(0, log(0.3)))), Constant(cbind(c(log(6),log(8)), c(log(2),log(1.3))))),
            list(matrix_frac, c(1,1), list(matrix(1:3), diag(3)), Constant(14)),
            list(matrix_frac, c(1,1), list(matrix(1:3), cbind(c(67,78,90), c(78,94,108), c(90,108,127))), Constant(0.46557377049180271)),
            list(matrix_frac, c(1,1), list(cbind(1:3, 4:6), cbind(c(67,78,90), c(78,94,108), c(90,108,127))), Constant(0.768852459016)),
            list(max_elemwise, c(2,1), list(matrix(c(-5,2)), matrix(c(-3,1)), 0, matrix(c(-1,2))), Constant(c(0,2))),
            list(max_elemwise, c(2,2), list(cbind(c(-5,2), c(-3,1)), 0, cbind(c(5,4), c(-1,2))), Constant(cbind(c(5,4), c(0,2)))),
            list(max_entries, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(2)),
            list(max_entries, c(1,1), list(matrix(c(-5,-10))), Constant(-5)),
            list(function(x) { max_entries(x, axis=1) }, c(2,1), list(cbind(c(-5,2), c(-3,1))), Constant(c(-3,2))),
            list(function(x) { max_entries(x, axis=2, keepdims=TRUE) }, c(1,2), list(cbind(c(-5,2), c(-3,1))), t(Constant(c(2,1)))),
            list(function(x) { norm(x,"2") }, c(1,1), list(v_np), Constant(3)),
            list(function(x) { norm(x,"F") }, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(5.47722557)),
            list(function(x) { norm1(x) }, c(1,1), list(v_np), Constant(5)),
            list(function(x) { norm1(x) }, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(7)),
            list(function(x) { norm_inf(x) }, c(1,1), list(v_np), Constant(2)),
            list(function(x) { norm_inf(x) }, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(6)),
            list(function(x) { norm_nuc(x) }, c(1,1), list(cbind(c(2,0), c(0,1))), Constant(3)),
            list(function(x) { norm_nuc(x) }, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(23.173260452512931)),
            list(function(x) { norm_nuc(x) }, c(1,1), list(cbind(3:5, 6:8)), Constant(14.618376738088918)),
            list(function(x) { sum_largest(abs(x),3) }, c(1,1), list(matrix(c(1,2,3,-4,-5))), Constant(5+4+3)),
            list(function(x) { dotsort(abs(x), matrix(c(2,1,0.1))) }, c(1,1), matrix(c(1,2,3,-4,-5)), Constant(10+4+0.3)),
            list(function(x) { mixed_norm(x,1,1) }, c(1,1), list(cbind(c(1,2), c(3,4), c(5,6))), Constant(21)),
            list(function(x) { mixed_norm(x,1,1) }, c(1,1), list(cbind(1:3, 4:6)), Constant(21)),
            # list(function(x) { mixed_norm(x,2,1) }, c(1,1), list(cbind(c(3,1), c(4,sqrt(3)))), Constant(7)),
            list(function(x) { mixed_norm(x,1,Inf) }, c(1,1), list(cbind(c(1,4), c(5,6))), Constant(10)),

            list(p_norm, c(1,1), list(matrix(1:3)), Constant(3.7416573867739413)),
            list(function(x) { p_norm(x,1) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(6.1)),
            list(function(x) { p_norm(x,2) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(3.7696153649941531)),
            list(function(x) { p_norm(x,2,axis=2) }, c(2,1), list(cbind(c(1,2), c(3,4))), Constant(c(sqrt(5), 5))),
            list(function(x) { p_norm(x,2,axis=1) }, c(2,1), list(cbind(c(1,2), c(4,5))), Constant(c(sqrt(17), sqrt(29)))),
            list(function(x) { p_norm(x,Inf) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(3)),
            list(function(x) { p_norm(x,3) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(3.3120161866074733)),
            list(function(x) { p_norm(x,5.6) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(3.0548953718931089)),
            list(function(x) { p_norm(x,1.2) }, c(1,1), list(cbind(1:3, 4:6)), Constant(15.971021676279573)),

            list(pos, c(1,1), list(8), Constant(8)),
            list(pos, c(2,1), list(matrix(c(-3,2))), Constant(c(0,2))),
            list(neg, c(2,1), list(matrix(c(-3,3))), Constant(c(3,0))),

            list(function(x) { power(x,1) }, c(1,1), list(7.45), Constant(7.45)),
            list(function(x) { power(x,2) }, c(1,1), list(7.45), Constant(55.502500000000005)),
            list(function(x) { power(x,-1) }, c(1,1), list(7.45), Constant(0.1342281879194631)),
            list(function(x) { power(x,-0.7) }, c(1,1), list(7.45), Constant(0.24518314363015764)),
            list(function(x) { power(x,-1.34) }, c(1,1), list(7.45), Constant(0.06781263100321579)),
            list(function(x) { power(x,1.34) }, c(1,1), list(7.45), Constant(14.746515290825071)),

            list(quad_over_lin, c(1,1), list(cbind(c(-1,2,-2), c(-1,2,-2)), 2), Constant(2*4.5)),
            list(quad_over_lin, c(1,1), list(v_np,2), Constant(4.5)),
            list(function(x) { norm(x,"2") }, c(1,1), list(cbind(c(2,0), c(0,1))),  Constant(2)),
            list(function(x) { norm(x,"2") }, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(22.368559552680377)),
            list(function(x) { scalene(x,2,3) }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(15,4), c(9,2)))),
            list(function(x) { power(x,2) }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(25,4), c(9,1)))),
            list(sum_entries, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(-5)),
            list(function(x) { sum_entries(x,axis=1) }, c(2,1), list(cbind(c(-5,2), c(-3,1))), Constant(c(-8,3))),
            list(function(x) { sum_entries(x,axis=2) }, c(2,1), list(cbind(c(-5,2), c(-3,1))), Constant(c(-3,-2))),
            list(function(x) { (x + Constant(0))^2 }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(25,4), c(9,1)))),
            list(function(x) { sum_largest(x,3) }, c(1,1), list(matrix(1:5)), Constant(5+4+3)),
            list(function(x) { sum_largest(x,3) }, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(9+10+11)),
            list(sum_squares, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(30)),
            list(matrix_trace, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(3+7+11)),
            list(matrix_trace, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(-5+1)),
            list(tv, c(1,1), list(matrix(c(1,-1,2))), Constant(5)),
            list(tv, c(1,1), list(t(matrix(c(1,-1,2)))), Constant(5)),
            list(tv, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(sqrt(53))),
            list(tv, c(1,1), list(cbind(c(-5,2), c(-3,1)), cbind(c(6,5), c(-4,3)), cbind(c(8,0), c(15,9))),
                 Constant(base::norm(c(7,-1,-8,2,-10,7), "2"))),
            list(tv, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(4*sqrt(10))),
            list(upper_tri, c(3,1), list(cbind(3:5, 6:8, 9:11)), Constant(c(6,9,10))),

            ## Advanced indexing
            list(function(x) { x[cbind(c(2,3), c(1,3))] }, c(2,1), list(cbind(3:5, 6:8, 9:11)), Constant(c(4,11))),
            list(function(x) { x[c(2,3),] }, c(2,2), list(cbind(3:5, 6:8)), Constant(cbind(c(4,5), c(7,8)))),
            list(function(x) { x[cbind(3:5, 6:8) %% 2 == 0] }, c(2,1), list(cbind(3:5, 6:8)), Constant(c(6,4,8))),
            list(function(x) { x[seq(3,2,-1)] }, c(2,1), list(matrix(3:5)), Constant(c(5,4))),
            list(function(x) { x[seq(3,1,-1)] }, c(3,1), list(matrix(3:5)), Constant(c(5,4,3))),
            list(function(x) { x[seq(4,2,-1)] }, c(2,1), list(matrix(3:5)), Constant(c(5,4))),
            list(function(x) { x[seq(4,1,-1)] }, c(3,1), list(matrix(3:5)), Constant(c(5,4,3)))
          )

atoms_maximize <- list(
            list(entr, c(2,2), list(cbind(c(1,exp(1)), c(exp(2), exp(-1)))), Constant(cbind(c(0,-exp(1)), c(-2*exp(2),exp(-1))))),
            list(log_det, c(1,1), list(cbind(c(20, 8, 5, 2),
                                             c(8, 16, 2, 4),
                                             c(5, 2, 5, 2),
                                             c(2, 4, 2, 4))), Constant(7.7424020218157814)),
            list(geo_mean, c(1,1), list(matrix(c(4,1))), Constant(2)),
            list(geo_mean, c(1,1), list(matrix(c(0.01,7))), Constant(0.2645751311064591)),
            list(geo_mean, c(1,1), list(matrix(c(63,7))), Constant(21)),
            list(geo_mean, c(1,1), list(matrix(c(1,10))), Constant(sqrt(10))),
            list(function(x) { geo_mean(x, c(1,1)) }, c(1,1), list(matrix(c(1,10))), Constant(sqrt(10))),
            list(function(x) { geo_mean(x, c(0.4,0.8,4.9)) }, c(1,1), list(matrix(c(0.5,1.8,17))), Constant(10.04921378316062)),
            list(harmonic_mean, c(1,1), list(matrix(c(1,2,3))), Constant(1.6363636363636365)),
            list(harmonic_mean, c(1,1), list(matrix(c(2.5,2.5,2.5,2.5))), Constant(2.5)),
            list(harmonic_mean, c(1,1), list(matrix(c(1e-8,1,2))), Constant(0)),

            list(function(x) { diff(x, differences = 0) }, c(3,1), list(matrix(c(1,2,3))), Constant(c(1,2,3))),
            list(diff, c(2,1), list(matrix(c(1,2,3))), Constant(c(1,1))),
            list(diff, c(1,1), list(matrix(c(1.1,2.3))), Constant(1.2)),
            list(function(x) { diff(x, differences = 2) }, c(1,1), list(matrix(c(1,2,3))), Constant(0)),
            list(diff, c(3,1), list(matrix(c(2.1,1,4.5,-0.1))), Constant(c(-1.1,3.5,-4.6))),
            list(function(x) { diff(x, differences = 2) }, c(2,1), list(matrix(c(2.1,1,4.5,-0.1))), Constant(c(4.6,-8.1))),
            list(function(x) { diff(x, differences = 1, axis = 1) }, c(2,1), list(cbind(c(-5,-3), c(2,1))), Constant(c(7,4))),
            list(function(x) { diff(x, differences = 1, axis = 2) }, c(1,2), list(cbind(c(-5,-3), c(2,1))), Constant(matrix(c(2,-1), nrow = 1))),

            list(function(x) { p_norm(x,0.5) }, c(1,1), list(matrix(c(1.1,2,0.1))), Constant(7.724231543909264)),
            list(function(x) { p_norm(x,-0.4) }, c(1,1), list(matrix(c(1.1,2,0.1))), Constant(0.02713620334)),
            list(function(x) { p_norm(x,-1) }, c(1,1), list(matrix(c(1.1,2,0.1))),  Constant(0.0876494023904)),
            list(function(x) { p_norm(x,-2.3) }, c(1,1), list(matrix(c(1.1,2,0.1))), Constant(0.099781528576)),

            list(lambda_min, c(1,1), list(cbind(c(2,0), c(0,1))), Constant(1)),
            list(lambda_min, c(1,1), list(cbind(c(5,7), c(7,-3))), Constant(-7.06225775)),
            list(function(x) { lambda_sum_smallest(x,2) }, c(1,1), list(cbind(c(1,2,3), c(2,4,5), c(3,5,6))), Constant(-0.34481428)),
            list(log, c(2,2), list(cbind(c(1,exp(1)), c(exp(2), exp(-1)))), Constant(cbind(c(0,1), c(2,-1)))),
            list(log1p, c(2,2), list(cbind(c(0,exp(1)-1), c(exp(2)-1,exp(-1)-1))), Constant(cbind(c(0,1), c(2,-1)))),
            list(min_elemwise, c(2,1), list(matrix(c(-5,2)), matrix(c(-3,1)), 0, matrix(c(1,2))), Constant(c(-5,0))),
            list(min_elemwise, c(2,2), list(cbind(c(-5,2), c(-3,-1)), 0, cbind(c(5,4), c(-1,2))), Constant(cbind(c(-5,0), c(-3,-1)))),
            list(min_entries, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(-5)),
            list(min_entries, c(1,1), list(matrix(c(-5,-10))), Constant(-10)),
            list(function(x) { x^0.25 }, c(1,1), list(7.45), Constant(7.45^0.25)),
            list(function(x) { x^0.32 }, c(2,1), list(matrix(c(7.45,3.9))), Constant(matrix(c(7.45,3.9), nrow = 2, ncol = 1)^0.32)),
            list(function(x) { x^0.9 }, c(2,2), list(cbind(c(7.45,2.2), c(4,7))), Constant(cbind(c(7.45,2.2), c(4,7))^0.9)),
            list(sqrt, c(2,2), list(cbind(c(2,4), c(16,1))), Constant(cbind(c(1.414213562373095,2), c(4,1)))),
            list(function(x) { sum_smallest(x,3) }, c(1,1), list(matrix(c(-1,2,3,4,5))), Constant(-1+2+3)),
            list(function(x) { sum_smallest(x,4) }, c(1,1), list(cbind(c(-3,-4,5), c(6,7,8), c(9,10,11))), Constant(-3-4+5+6)),
            list(function(x) { (x + Constant(0))^0.5 }, c(2,2), list(cbind(c(2,4), c(16,1))), Constant(cbind(c(1.414213562373095,2), c(4,1))))
          )


# Can the solver solver the problem?
check_solver <- function(prob, solver_name) {
  atom_str <- as.character(prob@objective@args[[1]])
  if(solver_name %in% names(KNOW_SOLVER_ERRORS)) {
    for(bad_atom_name in KNOWN_SOLVER_ERRORS[[solver_name]]) {
      if(grepl(bad_atom_name, atom_str, fixed = TRUE))
        return(FALSE)
    }
  }
  
  tryCatch({
    if(solver_name == ROBUST_CVXOPT)
      solver_name <- "CVXOPT"
    chains <- CVXR:::.construct_chains(prob, solver = solver_name)
    return(TRUE)
  }, error = function(e) {
    # warning(e)
    return(FALSE)
  })
}

# Test numeric version of atoms.
run_atom <- function(atom, problem, obj_val, solver, verbose = FALSE) {
  expect_true(is_dcp(problem))
  print(problem)
  
  if(verbose) {
    print(problem@objective)
    print(problem@constraints)
    print(paste("solver", solver))
  }
  
  if(check_solver(problem, solver)) {
    tolerance <- SOLVERS_TO_TOL[[solver]]
    
    tryCatch({
      if(solver == ROBUST_CVXOPT)
        result <- solve(problem, solver = "CVXOPT", verbose = verbose, kktsolver = ROBUST_KKTSOLVER)
      else
        result <- solve(problem, solver = solver, verbose = verbose)
    }, error = function(e) {
      if(solver %in% names(KNOWN_SOLVER_ERRORS) && atom %in% KNOWN_SOLVER_ERRORS[[solver]])
        return()
      stop(e)
    })
    
    if(result$status %in% SOLUTION_PRESENT) {
      if(verbose) {
        print(result$value)
        print(obj_val)
      }
      expect_true(-tolerance <= (result$value - obj_val) / (1 + abs(obj_val)) <= tolerance)
    } else
      stop("No solution present. Problem status: ", result$status)
  }
}

get_indices <- function(size) {
  # Get indices for dimension.
  if(length(size) == 0)
    return(c(1))
  else if(length(size) == 1)
    return(seq(size[1]))
  else
    expand.grid(seq(size[1]), seq(size[2]))
}

atoms_minimize <- lapply(atoms_minimize, function(a) { list(a, Minimize) })
atoms_maximize <- lapply(atoms_maximize, function(a) { list(a, Maximize) })

test_constant_atoms <- function(atom_info, objective_type) {
  atom <- atom_info[[1]]
  size <- atom_info[[2]]
  args <- atom_info[[3]]
  obj_val <- atom_info[[4]]
  
  for(indexer in get_indices(size)) {
    for(solver in SOLVERS_TO_TRY) {
      # Atoms with Constant arguments.
      prob_val <- value(obj_val[indexer])
      const_args <- lapply(args, Constant)
      if(length(size) != 0)
        objective <- objective_type(do.call(atom, const_args)[indexer])
      else
        objective <- objective_type(do.call(atom, const_args))
      problem <- Problem(objective)
      run_atom(atom, problem, prob_val, solver)
      
      # Atoms with Variable arguments.
      variables <- list()
      constraints <- list()
      for(idx in seq_along(args)) {
        expr <- args[[idx]]
        variables <- c(variables, new("Variable", dim = intf_dim(expr)))
        constraints <- c(constraints, variables[[length(variables)]] == expr)
      }
      if(length(size) != 0)
        objective <- objective_type(do.call(atom, variables)[indexer])
      else
        objective <- objective_type(do.call(atom, variables))
      problem <- Problem(objective, constraints)
      run_atom(atom, problem, prob_val, solver)
      
      # Atoms with Parameter arguments.
      parameters <- list()
      for(expr in args) {
        parameters <- c(parameters, new("Parameter", dim = intf_dim(expr)))
        value(parameters[[length(parameters)]]) <- intf_const_to_matrix(expr)
      }
      if(length(size) != 0)
        objective <- objective_type(do.call(atom, parameters)[indexer])
      else
        objective <- objective_type(do.call(atom, parameters))
      run_atom(atom, Problem(objective), prob_val, solver)
    }
  }
}

test_that("Test all constant atoms", {
  skip_on_cran()
  atoms_list <- c(atoms_minimize, atoms_maximize)
  for(a in atoms_list) {
    atom_info <- a[[1]]
    objective_type <- a[[2]]
    test_constant_atom(atom_info, objective_type)
  }
})
