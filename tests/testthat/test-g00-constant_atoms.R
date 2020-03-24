context("test-constant_atoms")

ROBUST_CVXOPT <- "robust_cvxopt"
SOLVERS_TO_TRY <- c("ECOS", "SCS", "OSQP")
SOLVERS_TO_TOL <- list(ECOS = 5e-7, SCS = 1e-2, OSQP = 1e-1) ## ECOS = 1e-7 fails one or two tests!

# Test CVXOPT if installed.
if("CVXOPT" %in% installed_solvers()) {
    SOLVERS_TO_TRYS <- c(SOLVERS_TO_TRY, "CVXOPT", ROBUST_CVXOPT)
    SOLVERS_TO_TOL$CVXOPT <- 1e-7
    SOLVERS_TO_TOL[[ROBUST_CVXOPT]] <- 1e-7
}

# Test MOSEK if installed.
if("MOSEK" %in% installed_solvers()) {
    SOLVERS_TO_TRY <- c(SOLVERS_TO_TRY, "MOSEK")
    SOLVERS_TO_TOL$MOSEK <- 1e-6
}

v_np <- matrix(c(-1, 2, -2), nrow = 1, ncol = 3)
log_sum_exp_axis_1 <- function(x) { log_sum_exp(x, axis = 1) }
log_sum_exp_axis_2 <- function(x) { log_sum_exp(x, axis = 2, keepdims = TRUE) }

atoms <- list(
    list(
        list(
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
            list(function(x) { kronecker(cbind(c(1,2), c(3,4)), x) }, c(4,4), list(cbind(c(5,6), c(7,8))),
                 Constant(kronecker(cbind(c(1,2), c(3,4)), cbind(c(5,6), c(7,8))))),
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
            list(function(x) { norm1(x) }, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(10)),
            list(function(x) { norm_inf(x) }, c(1,1), list(v_np), Constant(2)),
            list(function(x) { norm_inf(x) }, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(4)),
            list(function(x) { norm_nuc(x) }, c(1,1), list(cbind(c(2,0), c(0,1))), Constant(3)),
            list(function(x) { norm_nuc(x) }, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(23.173260452512931)),
            list(function(x) { norm_nuc(x) }, c(1,1), list(cbind(3:5, 6:8)), Constant(14.618376738088918)),
            list(function(x) { sum_largest(abs(x),3) }, c(1,1), list(matrix(c(1,2,3,-4,-5))), Constant(5+4+3)),
            list(function(x) { mixed_norm(x,1,1) }, c(1,1), list(cbind(c(1,2), c(3,4), c(5,6))), Constant(21)),
            list(function(x) { mixed_norm(x,1,1) }, c(1,1), list(cbind(1:3, 4:6)), Constant(21)),
            list(function(x) { mixed_norm(x,2,1) }, c(1,1), list(cbind(c(3,1), c(4,sqrt(3)))), Constant(7)),
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
            list(function(x) { x[c(2,3),] }, c(2,1), list(cbind(3:5, 6:8)), Constant(cbind(c(4,5), c(7,8)))),
            list(function(x) { x[cbind(3:5, 6:8) %% 2 == 0] }, c(3,1), list(cbind(3:5, 6:8)), Constant(c(4,6,8))),
            list(function(x) { x[seq(3,2,-1)] }, c(2,1), list(matrix(3:5)), Constant(c(5,4))),
            list(function(x) { x[seq(3,1,-1)] }, c(3,1), list(matrix(3:5)), Constant(c(5,4,3)))
        ),
        Minimize),
    list(
        list(
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
            list(harmonic_mean, c(1,1), list(matrix(c(0,1,2))), Constant(0)),

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
        ),
        Maximize
    )
)

check_solver <- function(prob, solver_name) {
    tryCatch({
        if(solver_name == ROBUST_CVXOPT)
            solver_name <- "CVXOPT"

        chains <- CVXR:::.construct_chains(prob, solver = solver_name)
        return(TRUE)
    }, error = function(e) {
        return(FALSE)
    })
}

ecnt  <- 0
run_atom <- function(atom, problem, obj_val, solver, verbose = FALSE) {
    expect_true(is_dcp(problem))
    ##print(problem)
    if(verbose) {
        print(problem@objective)
        print(problem@constraints)
        print(paste("solver", solver))
    }

    if(check_solver(problem, solver)) {
        tolerance <- SOLVERS_TO_TOL[[solver]]
        result <- solve(problem, solver = solver, verbose = verbose)

        if(result$status %in% c("optimal", "optimal_inaccurate")) {
            if(verbose) {
                print(result$value)
                print(obj_val)
            }

            obj_diff <- (result$value - obj_val)/(1+abs(obj_val))
            expect_true(abs(obj_diff) <= tolerance)

            if(abs(obj_diff) > tolerance) {
                ecnt  <- ecnt + 1
                sink("test_constant_atoms_out.txt", append = TRUE)
                cat(sprintf("Solver: %s\n", solver))
                print(atom)
                cat(sprintf("Result: %.10f \t Expected: %.10f \t Obj_diff: %.10f\n", result$value, obj_val, obj_diff))
                sink()
                saveRDS(list(problem  = problem, solver = "solver"), file = sprintf("error-%02d.RDS", ecnt))
            }
        } else
            stop("Problem status is sub-optimal: ", result$status)
    }
}

test_that("Test all constant atoms", {
    skip_on_cran()
    ## if(file.exists("test_constant_atoms_out.txt"))
    ##  file.remove("test_constant_atoms_out.txt")
    ##counter  <- 0 ## list item counter
    for(a in atoms) {
        atom_list <- a[[1]]
        objective_type <- a[[2]]
        for(al in atom_list) {
            ##counter  <- counter + 1 ## list item counter
            ##cat(sprintf("Item No: %d\n", counter))
            atom <- al[[1]]
            dims <- al[[2]]
            args <- al[[3]]
            obj_val <- al[[4]]
            for(row in 1:dims[1]) {
                for(col in 1:dims[2]) {
                    for(solver in SOLVERS_TO_TRY) {
                        ## Atoms with Constant arguments
                        const_args <- lapply(args, Constant)
                        run_atom(atom, Problem(objective_type(do.call(atom, const_args)[row, col])),
                                 value(obj_val[row, col]), solver)

                        ## Atoms with Variable arguments
                        variables <- list()
                        constraints <- list()
                        for(expr in args) {
                            expr_dim <- CVXR:::intf_dim(expr)
                            variables <- c(variables, Variable(expr_dim[1], expr_dim[2]))
                            constraints <- c(constraints, variables[[length(variables)]] == expr)
                        }
                        objective <- objective_type(do.call(atom, variables)[row, col])
                        ## print(atom)
                        ## print(value(obj_val[row, col]))
                        run_atom(atom, Problem(objective, constraints), value(obj_val[row, col]), solver)
                        ## cat("Index is", ind, "\n")
                        ## print("ATOM is")
                        ## print(atom)
                        ## print("Args is")
                        ## print(args)
                        ## ##Atoms with Parameter arguments
                        ## parameters <- list()
                        ## for(expr in args) {
                        ##     expr_dim <- dim(expr)
                        ##     parameters <- c(parameters, Parameter(expr_dim[1], expr_dim[2]))
                        ##     value(parameters[[length(parameters)]]) <- as.matrix(expr)
                        ## }
                        ## objective <- objective_type(do.call(atom, parameters)[row, col])
                        ## run_atom(atom, Problem(objective), value(obj_val[row, col]), solver)
                    }
                }
            }
        }
    }
})
